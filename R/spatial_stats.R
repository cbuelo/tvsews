library(ape)
library(geosphere)
library(parallel)
library(multidplyr)
library(magrittr)
library(dplyr)
library(devtools)
library(tidyr)

calc_moransI <- function(df_sssv, lon_col = "longitude", lat_col="latitude", meas_col="Measurement"){
  # get rid of data with missing coordinates
  df_sssv_0 = df_sssv %>%
    dplyr::filter(!is.na(.data[[lon_col]]) & !is.na(.data[[lat_col]]))
  if(sum(!is.na(df_sssv_0[[meas_col]])) < 2){
    return(NA)
  }
  # calc distances
  distances = geosphere::distm(x=as.matrix(df_sssv_0[, c(lon_col, lat_col)]))
  distances_inv =  1/distances # exp(-distances) #
  diag(distances_inv) = 0
  distances_inv[is.infinite(distances_inv)] = 0
  # calc Moran's I
  hold_MI = ape::Moran.I(df_sssv_0[[meas_col]], distances_inv, na.rm = T)
  return(hold_MI$observed)
}


calc_spatial_stats <- function(spatial_data = flame_data, lat_col = "latitude", lon_col = "longitude", var_cols=c("BGApc_ugL_tau", "ODO_percent_tau", "pH_tau"), id_cols = c("Lake", "Year", "DOY"), multiple_cores = TRUE){
  # calculate spatial standard deviation
  spat_SD = flame_data %>%
    group_by(across(all_of(id_cols))) %>%
    select(all_of(var_cols)) %>%
    summarise_all(sd, na.rm=TRUE) %>%
    tidyr::pivot_longer(cols=all_of(var_cols), names_to = "Variable", values_to = "Value") %>%
    mutate(Stat = "SD") %>%
    ungroup()

  # calculate moran's I
  # warn this will take a long time if not using multiple cores
  if(multiple_cores == FALSE){
    message("multiple_cores set to FALSE, this will probably take awhile (~15 minutes for default data)")
  }else{
    cores=parallel::detectCores()
    use_cores = max(cores - 2, 1)
    # give a rough estimate of time it will take
    message(paste("Using", use_cores, "cores; this will probably take roughly", round(15/use_cores, 1), "minutes (assuming using default data)"))
    cluster = new_cluster(use_cores)
  }
  # re-format data: wide -> long, group, and nest
  flame_data_partitioned = flame_data %>%
    tidyr::pivot_longer(cols=all_of(var_cols), names_to = "Variable", values_to = "Measurement") %>%
    select(all_of(c(id_cols, "Variable", lat_col, lon_col, "Measurement"))) %>%
    group_by(across(all_of(c(id_cols, "Variable")))) %>%
    nest()

  # partition data and set up cluster if using multiple cores
  if(multiple_cores == TRUE){
    # partition data
    flame_data_partitioned = flame_data_partitioned %>%
      partition(cluster)
    # copy info to cluster
    cluster_copy(cluster, "calc_moransI")
    cluster_library(cluster, c("ape", "dplyr", "magrittr"))
  }

  # do the calculation
  morans_I0 = flame_data_partitioned %>%
    mutate(Value = purrr::map(data, calc_moransI))

  if(multiple_cores == TRUE){
    morans_I0 = morans_I0 %>%
      collect()
    }

  morans_I = morans_I0 %>%
    unnest(Value) %>%
    select(-data) %>%
    mutate(Stat = "Moran's I")

  spat_stats_comb = bind_rows(spat_SD, morans_I)

  return(spat_stats_comb)
}

tStart = Sys.time()
spst = calc_spatial_stats()
tEnd = Sys.time()

tEnd - tStart

tStart = Sys.time()
spst = calc_spatial_stats(multiple_cores = FALSE)
tEnd = Sys.time()

tEnd - tStart
