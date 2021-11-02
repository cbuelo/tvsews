if(getRversion() >= "2.15.1"){
  utils::globalVariables(c("grouped_data"))
}

#' Calculate Moran's I from a data frame
#'
#' Computes spatial autocorrelation for a single sample and single variable (sssv) from a data frame with three columns (longitude, latitude, and measurement). Uses inverse of distance between points in Moran's I calculation.
#'
#' @param df_sssv data frame with spatial measurements of a single sample and single variable
#' @param lon_col character string, column name for longitude column
#' @param lat_col character string, column name for latitude column
#' @param meas_col character string, column name for measurement column
#'
#' @return numeric, the Moran's I (spatial autocorrelation)
#' @export
#'
#' @examples
#' example_data <- data.frame(
#'   lon = runif(n = 10, min = -90, max = -85),
#'   lat = runif(n = 10, min = 45, max = 46),
#'   meas = rnorm(10)
#' )
#' calc_moransI(example_data, lon_col = "lon", lat_col = "lat", meas_col = "meas")
#'
#' peter_lake_data <- flame_data %>%
#'   dplyr::filter(Lake == "R" & Date == as.Date("2019-07-22")) %>%
#'   dplyr::select(latitude, longitude, BGApc_ugL_tau)
#' calc_moransI(peter_lake_data, meas_col = "BGApc_ugL_tau")
calc_moransI <- function(df_sssv, lon_col = "longitude", lat_col = "latitude", meas_col = "Measurement") {
  # get rid of data with missing coordinates
  df_sssv_0 <- df_sssv %>%
    dplyr::filter(!is.na(.data[[lon_col]]) & !is.na(.data[[lat_col]]))
  if (sum(!is.na(df_sssv_0[[meas_col]])) < 2) {
    return(NA)
  }
  # calc distances
  distances <- geosphere::distm(x = as.matrix(df_sssv_0[, c(lon_col, lat_col)]))
  distances_inv <- 1 / distances # exp(-distances) #
  diag(distances_inv) <- 0
  distances_inv[is.infinite(distances_inv)] <- 0
  # calc Moran's I
  hold_MI <- ape::Moran.I(df_sssv_0[[meas_col]], distances_inv, na.rm = T)
  return(hold_MI$observed)
}


#' Calculate spatial statistics
#'
#' Calculates spatial early warning statistics (spatial standard deviation and spatial autocorrelation = Moran's I) from a data frame of (optionally multiple) sampling events and variables. By default uses multiple cores to speed up spatial autocorrelation calculation, which can take quite a bit of time depending on how many data points you have.
#'
#' Note that time estimates printed out are very rough estimates and are based on using the default data and variables, and run on a AMD Ryzen 1800X with 8 cores/16 threads.
#'
#' @param spatial_data data frame containing spatial locations and measurements of one or more variables; see \code{\link{flame_data}} for default and formatting
#' @param lon_col character string, column name for longitude column
#' @param lat_col character string, column name for latitude column
#' @param var_cols character vector, column names that hold measurements of different variables
#' @param id_cols character vector, columns that identify unique sampling events to separate data into before calculating stats independently on, defaults are "Lake", "Year", and "DOY"
#' @param statistics character vector, statistics to run; options are "SD", "moransI", "skew", and "kurt"
#' @param multiple_cores TRUE (default) or FALSE, should spatial autocorrelation calculations be run on multiple cores to speed up calculations?
#'
#' @return data frame with calculated spatial stats (SD and Moran's I) for each sample and variable
#' @import dplyr
#' @export
#'
#' @examples
#' spat_stats <- calc_spatial_stats(flame_data)
#' library(ggplot2)
#' ggplot(spat_stats %>% dplyr::filter(Stat == "SD"), aes(x = DOY, y = Value, color = Lake)) +
#'   geom_line() +
#'   facet_grid(rows = vars(Variable), cols = vars(Year), scales = "free_y") +
#'   theme_bw()
calc_spatial_stats <- function(spatial_data = flame_data, lat_col = "latitude", lon_col = "longitude", var_cols = c("BGApc_ugL_tau", "ODO_percent_tau", "pH_tau"), id_cols = c("Lake", "Year", "DOY"), statistics = c("SD", "moransI"), multiple_cores = TRUE) {

  results_list = list()
  # calculate spatial standard deviation
  if("SD" %in% statistics){
    spat_SD <- spatial_data %>%
      group_by(across(all_of(id_cols))) %>%
      select(all_of(var_cols)) %>%
      summarise_all(stats::sd, na.rm = TRUE) %>%
      tidyr::pivot_longer(cols = all_of(var_cols), names_to = "Variable", values_to = "Value") %>%
      mutate(Stat = "SD") %>%
      ungroup()

    results_list[["SD"]] = spat_SD
  }

  if("skew" %in% statistics){
    spat_skew <- spatial_data %>%
      group_by(across(all_of(id_cols))) %>%
      select(all_of(var_cols)) %>%
      summarise_all(moments::skewness, na.rm = TRUE) %>%
      tidyr::pivot_longer(cols = all_of(var_cols), names_to = "Variable", values_to = "Value") %>%
      mutate(Stat = "skew") %>%
      ungroup()

    results_list[["skew"]] = spat_skew
  }

  if("kurt" %in% statistics){
    spat_kurt <- spatial_data %>%
      group_by(across(all_of(id_cols))) %>%
      select(all_of(var_cols)) %>%
      summarise_all(moments::kurtosis, na.rm = TRUE) %>%
      tidyr::pivot_longer(cols = all_of(var_cols), names_to = "Variable", values_to = "Value") %>%
      mutate(Stat = "kurt") %>%
      ungroup()

    results_list[["kurt"]] = spat_kurt
  }

  if("moransI" %in% statistics){
    # calculate moran's I
    # warn this will take a long time if not using multiple cores
    if (multiple_cores == FALSE) {
      message("multiple_cores set to FALSE, this will probably take awhile (~10 minutes for default data)")
    } else {
      cores <- parallel::detectCores()
      use_cores <- max(cores - 2, 1)
      # give a rough estimate of time it will take
      message(paste("Using", use_cores, "cores; this will probably take roughly", round(15 / use_cores, 1), "minutes (assuming using default data)"))
      cluster <- multidplyr::new_cluster(use_cores)
    }
    # re-format data: wide -> long, group, and nest
    flame_data_partitioned <- spatial_data %>%
      tidyr::pivot_longer(cols = all_of(var_cols), names_to = "Variable", values_to = "Measurement") %>%
      select(all_of(c(id_cols, "Variable", lat_col, lon_col, "Measurement"))) %>%
      group_by(across(all_of(c(id_cols, "Variable")))) %>%
      tidyr::nest(grouped_data = - group_cols())

    # partition data and set up cluster if using multiple cores
    if (multiple_cores == TRUE) {
      # partition data
      flame_data_partitioned <- flame_data_partitioned %>%
        multidplyr::partition(cluster)
      # copy info to cluster
      multidplyr::cluster_copy(cluster, "calc_moransI")
      multidplyr::cluster_library(cluster, c("ape", "dplyr", "magrittr"))
    }

    # do the calculation
    morans_I0 <- flame_data_partitioned %>%
      mutate(Value = purrr::map(grouped_data, calc_moransI))

    if (multiple_cores == TRUE) {
      morans_I0 <- morans_I0 %>%
        collect()
    }

    morans_I <- morans_I0 %>%
      tidyr::unnest(.data$Value) %>%
      select(-grouped_data) %>%
      mutate(Stat = "Moran's I")

    results_list[["moransI"]] = morans_I
  }

  spat_stats_comb <- bind_rows(results_list)

  return(spat_stats_comb)
}
