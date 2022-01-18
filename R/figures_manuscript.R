if(getRversion() >= "2.15.1"){
  flame_vars = c("Lake", "Year", "DOY", "Date", "date_time", "latitude", "longitude", "BGApc_ugL_tau", "ODO_percent_tau", "pH_tau" )
  ts_vars = c("Lake", "Year", "DOY", "BGA_HYLB", "Manual_Chl", "DO_Sat", "pH" )
  fert_vars = c("Year", "Lake", "fertStartDOY", "fertEndDOY", "bloomStartDOY")
  bounds_vars = c("AREA", "PERIMETER", "LAKE_TYPE", "NAME", "UNIQUE_ID", "COUNTY", "NOTE24", "NEW_KEY", "HECTARES", "ACRES_GIS", "FMU", "Lake", "FID_1", "Lake_Name", "geometry")
  df_vars = c("bloom_fert_dates", "flame_data", "peter_paul_bounds", "ts_data")
  plot_vars = c("Start Of", "Value", "DOY_lab", "DOYtrunc", "Variable", "AlarmType", "Rate", "PositiveType")
  utils::globalVariables(unique(c(flame_vars, ts_vars, fert_vars, bounds_vars, df_vars, plot_vars)))
}

#' Plot state variables
#'
#' @param time_series data frame containing time series of state variables, , see \code{\link{ts_data}} for default and formatting
#' @param bloom_fert_df data frame of fertilization start and end and bloom start dates, see \code{\link{bloom_fert_dates}} for default and formatting
#' @param var_rename_vec named vector, used to change variable names displayed in plot, format is c("new name" = "original name")
#' @param lakes_rename_vec named vector, used to change lake names displayed in plot, format is c("new name" = "original name")
#' @param lake_colors named vector, used to specify colors for each lake, format is c("Lake" = "color")
#' @param legend_location numerical vector of length 2, has c(x, y) position for placing legend on plot
#'
#' @return ggplot object plotting the time series
#' @import ggplot2
#' @export
#'
#' @examples
#' # all data:
#' plot_state_vars()
#'
#' # just latest experiment
#' plot_state_vars(
#'   time_series = ts_data %>% dplyr::filter(Year >= 2018),
#'   bloom_fert_df = bloom_fert_dates %>% dplyr::filter(Lake == "R" & Year >= 2018)
#' )
plot_state_vars <- function(time_series = ts_data,
                      bloom_fert_df = bloom_fert_dates,
                      var_rename_vec = c(`Chl-a (ug/L)` = "Manual_Chl", `BGA (cells/mL)` = "BGA_HYLB", `D.O. sat. (%)` = "DO_Sat", pH = "pH"),
                      lakes_rename_vec = c("Reference" = "L", "Experimental" = "R"),
                      lake_colors = c("Experimental" = "firebrick3", "Reference" = "royalblue3"),
                      legend_location = c(0.13, 0.9)) {
  # rename variables and pivot time series longer
  ts_formatted <- time_series %>%
    dplyr::rename(!!!var_rename_vec) %>%
    tidyr::pivot_longer(cols = names(var_rename_vec), names_to = "Variable", values_to = "Value") %>%
    dplyr::mutate(Variable = factor(.data$Variable, levels = names(var_rename_vec), ordered = TRUE))

  # fix lake names
  if (all(unique(ts_formatted$Lake) %in% lakes_rename_vec)) {
    lakes_rename_vec_rev <- names(lakes_rename_vec)
    names(lakes_rename_vec_rev) <- unname(lakes_rename_vec)
    ts_formatted <- ts_formatted %>%
      dplyr::mutate(Lake = unname(lakes_rename_vec_rev[.data$Lake]))
  }

  # format bloom dates
  dates_formatted <- bloom_fert_df %>%
    dplyr::select(Year, Lake, fertStartDOY, bloomStartDOY) %>%
    dplyr::rename(Nutrients = fertStartDOY, Bloom = bloomStartDOY) %>%
    tidyr::pivot_longer(cols = c("Nutrients", "Bloom"), names_to = "Start Of", values_to = "DOY") %>%
    dplyr::mutate(`Start Of` = factor(`Start Of`, levels = c("Nutrients", "Bloom"), ordered = TRUE))

  # fix lake names
  if (all(unique(dates_formatted$Lake) %in% lakes_rename_vec)) {
    lakes_rename_vec_rev <- names(lakes_rename_vec)
    names(lakes_rename_vec_rev) <- unname(lakes_rename_vec)
    dates_formatted <- dates_formatted %>%
      dplyr::mutate(Lake = unname(lakes_rename_vec_rev[.data$Lake]))
  }

  # make the plot
  p1 <- ts_formatted %>%
    ggplot(aes(x = DOY, y = Value, color = Lake)) +
    geom_line(size = 2) +
    facet_grid(rows = vars(Variable), cols = vars(Year), scales = "free_y") +
    theme_bw() +
    labs(x = "Day of Year", y = "") +
    # add legend
    theme(legend.position = legend_location, legend.background = element_blank(), legend.key = element_blank())

  # do some conditional modifications to the plot
  if (length(unique(ts_formatted$Lake)) == length(lake_colors) & all(names(lake_colors) %in% unique(ts_formatted$Lake))) {
    p1 <- p1 +
      scale_color_manual(values = lake_colors)
  }
  # add bloom lines
  colorMap <- ifelse(length(unique(dates_formatted$Lake)) > 1, "Lake", "'black'")
  if (length(unique(dates_formatted$Lake)) > 1) {
    p1 <- p1 +
      geom_vline(data = dates_formatted, aes(xintercept = DOY, linetype = `Start Of`, color = Lake)) +
      scale_linetype_manual(
        values = c("Bloom" = "solid", "Nutrients" = "dashed"),
        guide = guide_legend(title.position = "top", direction = "horizontal")
      )
  } else {
    p1 <- p1 +
      geom_vline(data = dates_formatted, aes(xintercept = DOY, linetype = `Start Of`)) +
      scale_linetype_manual(
        values = c("Bloom" = "solid", "Nutrients" = "dashed"),
        guide = guide_legend(title.position = "top", direction = "horizontal")
      )
  }


  return(p1)
}


#' Plot example maps of spatial data
#'
#' @param spatial_data data frame containing spatial locations and measurements of one or more variables; see \code{\link{flame_data}} for default and formatting
#' @param lake_boundaries sf data frame containing boundaries for lakes; see \code{\link{peter_paul_bounds}} and str(peter_paul_bounds) for default and formatting
#' @param samples_to_plot a data frame specifying which info from \code{spatial_data} should be used for plotting
#' @param var_cols character vector, column names that hold measurements of different variables
#'
#' @return a ggplot object with the specified maps of spatial data
#' @import ggplot2
#' @export
#'
#' @examples
#' plot_spatial_data(
#'   spatial_data = flame_data,
#'   samples_to_plot = data.frame(Year = 2019, DOY = c(165, 210, 228), stringsAsFactors = FALSE),
#'   var_cols = c("BGApc_ugL_tau")
#' )
plot_spatial_data <- function(
                      spatial_data = flame_data,
                      lake_boundaries = peter_paul_bounds,
                      samples_to_plot = data.frame(Year = 2019, DOY = c(165, 210, 228), stringsAsFactors = FALSE),
                      var_cols = c("BGApc_ugL_tau", "ODO_percent_tau", "pH_tau")) {
  flame_to_plot <- samples_to_plot %>%
    left_join(spatial_data) %>%
    select(c(Year, Lake, DOY, "latitude", "longitude", all_of(var_cols))) %>%
    mutate(DOY_lab = paste("DOY", DOY)) # %>%
  # TODO: make this work with a function argument
  var_labels <- data.frame(
    Var = c("BGApc_ugL_tau", "chlor_ugL_tau", "ODO_percent_tau", "pH_tau", "specCond_tau", "turb_FNU_tau", "fDOM_QSU_tau"),
    Var_form = c("BGA (ug/L)    ", "Chl-a", "D.O. sat. (%)", "pH            ", "Sp. Cond", "Turb.", "fDOM")
  )
  # convert data to sf
  flame_to_plot_sf <- sf::st_as_sf(flame_to_plot, coords = c("longitude", "latitude"), crs = 4326, agr = "constant")

  # make plot using apply
  names_plot <- names(flame_to_plot)[!names(flame_to_plot) %in% c("Year", "Lake", "DOY", "DOY_lab", "latitude", "longitude")]
  plot_list <- lapply(names_plot, function(Var, plotDF, varLabs) {
    ggplot() +
      facet_wrap(vars(DOY_lab), nrow = 1) +
      geom_sf(data = lake_boundaries, inherit.aes = FALSE, fill = "white", size = 0.8) +
      geom_sf(data = plotDF, aes_string(color = Var), size = 0.5) +
      theme_bw() +
      viridis::scale_color_viridis(option = "viridis") +
      scale_x_continuous(breaks = c(-89.503, -89.504, -89.505), labels = c("-89.503", "-89.504", "-89.505")) +
      scale_y_continuous(breaks = seq(46.251, 46.2535, 0.0005), labels = as.character(seq(46.251, 46.2535, 0.0005))) +
      labs(x = "", y = "", color = varLabs[varLabs$Var == Var, "Var_form"])
  },
  plotDF = flame_to_plot_sf,
  varLabs = var_labels
  )
  Margin <- theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
  out_plot <- cowplot::plot_grid(plotlist = lapply(plot_list, "+", Margin), align = "vh", nrow = length(var_cols))
  return(out_plot)
}


#' Plot rolling window stats and quickest detection alarms for all variables
#'
#' @param rolling_window_stats data frame, output from call to \code{\link{calc_rolling_stats}}
#' @param qd_alarms data frame, output from call to \code{\link{format_qd}}
#' @param bloom_fert_df data frame of fertilization start and end and bloom start dates, see \code{\link{bloom_fert_dates}} for default and formatting
#' @param var_rename_vec named vector, used to change variable names displayed in plot, format is c("new name" = "original name")
#' @param legend_location numerical vector of length 2, has c(x, y) position for placing legend on plot
#'
#' @return ggplot object showing rolling window statistics and quickest detection alarms
#' @import ggplot2
#' @import dplyr
#' @export
#'
#' @examples
#' use_vars <- c("Manual_Chl", "BGA_HYLB", "DO_Sat", "pH")
#' rw_stats <- calc_rolling_stats(ts_data, var_cols = use_vars)
#' qd_s <- run_qd(rw_stats, var_cols = use_vars)
#' qd_a <- format_qd(qd_s, bloom_fert_df = bloom_fert_dates)
#' plot_temporal_EWS(rolling_window_stats = rw_stats,
#'   qd_alarms = qd_a, bloom_fert_df = bloom_fert_dates)
plot_temporal_EWS <- function(rolling_window_stats, qd_alarms, bloom_fert_df = bloom_fert_dates, var_rename_vec = c("Chl-a" = "Manual_Chl", "Phyco" = "BGA_HYLB", "D.O. sat." = "DO_Sat", "pH" = "pH"), legend_location = c(0.13, 0.9)) {
  # rename variables
  rolling_window_stats <- rolling_window_stats %>% rename(!!!var_rename_vec)
  # pivot long for plotting
  rw_stats_long <- rolling_window_stats %>%
    filter(Stat %in% c("SD", "Ar1")) %>%
    tidyr::pivot_longer(cols = names(var_rename_vec), names_to = "Variable", values_to = "Value") %>%
    mutate(Variable = factor(.data$Variable, levels = names(var_rename_vec), ordered = TRUE))
  # fromate bloom dates for plotting
  bloom_dates <- bloom_fert_df %>%
    tidyr::pivot_longer(cols = c("fertStartDOY", "fertEndDOY", "bloomStartDOY"), names_to = "Start Of", values_to = "DOYtrunc") %>%
    mutate(
      `Start Of` = ifelse(`Start Of` == "fertStartDOY", "nutrients", `Start Of`),
      `Start Of` = ifelse(`Start Of` == "bloomStartDOY", "bloom", `Start Of`)
    ) %>%
    filter(`Start Of` %in% c("nutrients", "bloom"))
  # format alarms and get just the alarms DOYS
  var_rename_vec_rev <- names(var_rename_vec)
  names(var_rename_vec_rev) <- unname(var_rename_vec)
  qd_alarms_true <- qd_alarms %>%
    filter(!is.na(.data$Alarm) & .data$Alarm == 1) %>%
    mutate(Variable = unname(var_rename_vec_rev[.data$Var])) %>%
    mutate(Variable = factor(.data$Variable, levels = names(var_rename_vec), ordered = TRUE))

  ## make SD plot
  # plot rw stats
  sd_plot <- rw_stats_long %>%
    filter(Stat == "SD") %>%
    ggplot(aes(x = DOYtrunc, y = Value, color = Lake)) +
    geom_line(size = 1.25) +
    theme_bw() +
    facet_grid(rows = vars(Variable), cols = vars(Year), scales = "free_y")

  if (all(unique(rw_stats_long$Lake) %in% c("R", "L"))) {
    sd_plot <- sd_plot +
      scale_color_manual(values = c("R" = "firebrick3", "L" = "royalblue3"), guide = "none")
  }

  sd_plot <- sd_plot +
    labs(x = "", y = "") +
    # add alarms
    geom_point(
      data = qd_alarms_true %>% filter(Stat == "SD"),
      mapping = aes(fill = AlarmType),
      shape = 21,
      size = 3,
      color = "black"
    ) +
    scale_fill_manual(values = c("#FFE77AFF", "#5F9260FF", "darkgrey")) +
    ggtitle("SD") +
    theme(
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.position = legend_location,
      plot.title = element_text(hjust = 0.5, size = 18),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 16),
      strip.text = element_text(size = 16)
    ) +
    # add fert/bloom lines
    geom_vline(data = bloom_dates, aes(xintercept = DOYtrunc, linetype = `Start Of`)) +
    scale_linetype_manual(breaks = c("nutrients", "bloom"), values = c("dashed", "solid"), guide = "none")

  ## make AR1 plot
  # plot rw stats
  ar1_plot <- rw_stats_long %>%
    filter(Stat == "Ar1") %>%
    ggplot(aes(x = DOYtrunc, y = Value, color = Lake)) +
    geom_line(size = 1.25) +
    theme_bw() +
    facet_grid(rows = vars(Variable), cols = vars(Year))

  if (all(unique(rw_stats_long$Lake) %in% c("R", "L"))) {
    ar1_plot <- ar1_plot +
      scale_color_manual(values = c("R" = "firebrick3", "L" = "royalblue3"), guide = "none")
  }
  ar1_plot <- ar1_plot +
    labs(x = "", y = "") +
    # add alarms
    geom_point(
      data = qd_alarms_true %>% filter(Stat == "Ar1"),
      mapping = aes(fill = AlarmType),
      shape = 21,
      size = 3,
      color = "black"
    ) +
    scale_fill_manual(values = c("False" = "#FFE77AFF", "True" = "#5F9260FF", "Late" = "darkgrey"), guide = "none") +
    ggtitle("AR(1)") +
    theme(
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.position = legend_location,
      plot.title = element_text(hjust = 0.5, size = 18),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 16),
      strip.text = element_text(size = 16)
    ) +
    # add fert/bloom lines
    geom_vline(data = bloom_dates, aes(xintercept = DOYtrunc, linetype = `Start Of`)) +
    scale_linetype_manual(breaks = c("nutrients", "bloom"), values = c("dashed", "solid"), guide = "none")


  ## combine
  out_fig <- gridExtra::grid.arrange(sd_plot, ar1_plot, bottom=grid::textGrob("Day of Year", gp=grid::gpar(fontsize=18), vjust=-.25), ncol=2, padding=unit(1, "line"))
  return(out_fig)
}

#' Plot spatial statistics time series
#'
#' @param spatial_stats data frame, output from call to \code{\link{calc_spatial_stats}}
#' @param bloom_fert_df data frame of fertilization start and end and bloom start dates, see \code{\link{bloom_fert_dates}} for default and formatting
#' @param var_rename_vec named vector, used to change variable names displayed in plot, format is c("new name" = "original name")
#'
#' @return a ggplot object showing the change in spatial statistics over the experiment
#' @import ggplot2
#' @import dplyr
#' @export
#'
#' @examples
#' spat_stats <- calc_spatial_stats(
#'   spatial_data = flame_data
#' )
#' plot_spatial_EWS(spat_stats)
plot_spatial_EWS <- function(spatial_stats, bloom_fert_df = bloom_fert_dates, var_rename_vec = c("BGA (ug/L)" = "BGApc_ugL_tau", "DO sat." = "ODO_percent_tau", "pH" = "pH_tau")) {

  # re-name the variables
  var_rename_vec_rev <- names(var_rename_vec)
  names(var_rename_vec_rev) <- unname(var_rename_vec)

  spatial_stats <- spatial_stats %>%
    mutate(Variable = unname(var_rename_vec_rev[.data$Variable])) %>%
    mutate(Variable = factor(.data$Variable, levels = names(var_rename_vec), ordered = TRUE))

  # format bloom dates for plotting
  bloom_dates <- bloom_fert_df %>%
    tidyr::pivot_longer(cols = c("fertStartDOY", "fertEndDOY", "bloomStartDOY"), names_to = "Start Of", values_to = "DOYtrunc") %>%
    mutate(
      `Start Of` = ifelse(`Start Of` == "fertStartDOY", "nutrients", `Start Of`),
      `Start Of` = ifelse(`Start Of` == "bloomStartDOY", "bloom", `Start Of`)
    ) %>%
    filter(`Start Of` %in% c("nutrients", "bloom") & Lake == "R" & Year >= 2018) # TODO don't hard code lake and Year here

  # make the plots
  ar1_spatial_plot <- spatial_stats %>%
    filter(Stat == "Moran's I") %>% #
    ggplot(aes(x = DOY, y = Value, color = Lake)) +
    geom_line(size = 1.25) +
    theme_bw() +
    facet_grid(rows = vars(Variable), cols = vars(Year)) + # , scales="free_y"
    scale_color_manual(values = c("R" = "firebrick3", "L" = "royalblue3"), guide = guide_legend(title.position = "top")) +
    labs(x = "", y = "") +
    # theme(legend.position=c(0.25, 0.65), legend.background = element_blank(), legend.key = element_blank(), legend.direction ="horizontal") +
    geom_vline(data = bloom_dates, aes(xintercept = DOYtrunc, linetype = `Start Of`)) +
    scale_linetype_manual(breaks = c("nutrients", "bloom"), values = c("dashed", "solid"), guide = guide_legend(title.position = "top")) +
    ggtitle("Moran's I") +
    theme(
      legend.title=element_blank(),
      legend.position = "none",
      plot.title=element_text(hjust=0.5, size=18),
      strip.text = element_text(size=16),
      legend.text = element_text(size=12),
      axis.text=element_text(size=12),
      panel.spacing.x = unit(1, "lines")
    )

  sd_spatial_plot <- spatial_stats %>%
    filter(Stat == "SD") %>% #
    ggplot(aes(x = DOY, y = Value, color = Lake)) +
    geom_line(size = 1.25) +
    theme_bw() +
    facet_grid(rows = vars(Variable), cols = vars(Year), scales = "free_y") +
    scale_color_manual(values = c("R" = "firebrick3", "L" = "royalblue3"), guide = guide_legend(title.position = "top")) +
    labs(x = "", y = "") +
    # theme(legend.position=c(0.25, 0.65), legend.background = element_blank(), legend.key = element_blank(), legend.direction ="horizontal") +
    geom_vline(data = bloom_dates, aes(xintercept = DOYtrunc, linetype = `Start Of`)) +
    scale_linetype_manual(breaks = c("nutrients", "bloom"), values = c("dashed", "solid"), guide = guide_legend(title.position = "top")) +
    ggtitle("SD")  +
    theme(
      legend.title=element_blank(),
      legend.position = "none",
      plot.title=element_text(hjust=0.5, size=18),
      strip.text = element_text(size=16),
      legend.text = element_text(size=12),
      axis.text=element_text(size=12),
      panel.spacing.x = unit(1, "lines")
    )

  out_plot <- gridExtra::grid.arrange(sd_spatial_plot, ar1_spatial_plot, bottom = grid::textGrob("Day of Year", gp=grid::gpar(fontsize=18), vjust=-.25), ncol = 2, padding = unit(1, "line"))

  return(out_plot)
}


#' Plot true and false positive quickest detection alarm rates
#'
#' @param qd_alarm_rates data frame, output from call to \code{\link{calc_alarm_rates}}
#' @param var_rename_vec named vector, used to change variable names displayed in plot, format is c("new name" = "original name")
#' @param y_lim numeric vector of length 2, optionally specify y limits manually
#' @param plot_title character string, optional main title for plot
#' @param legend_title character string, optional title for legend
#'
#' @return ggplot2 object, bar plot giving the true and false positive rates
#' @export
#' @import dplyr
#' @import ggplot2
#'
#' @examples
#' use_vars <- c("Manual_Chl", "BGA_HYLB", "DO_Sat", "pH")
#' rw_data <- ts_data %>% dplyr::filter(Year >= 2011)
#' rw_stats <- calc_rolling_stats(rw_data, var_cols = use_vars)
#' qd_stats <- run_qd(rw_stats, var_cols = use_vars)
#' qd_alarms <- format_qd(qd_stats, bloom_fert_df = bloom_fert_dates)
#' qd_rates <- calc_alarm_rates(qd_alarms)
#' plot_alarm_rates(qd_alarm_rates = qd_rates)
plot_alarm_rates <- function(qd_alarm_rates, var_rename_vec = c("Chl-a" = "Manual_Chl", "Phyco" = "BGA_HYLB", "D.O. sat." = "DO_Sat", "pH" = "pH"), y_lim = NULL, plot_title = NULL, legend_title = "Positive Alarm Rate") {
  # format data
  var_rename_vec_rev <- names(var_rename_vec)
  names(var_rename_vec_rev) <- unname(var_rename_vec)
  qd_alarm_rates <- qd_alarm_rates %>%
    mutate(Variable = unname(var_rename_vec_rev[.data$Var])) %>%
    mutate(Variable = factor(.data$Variable, levels = names(var_rename_vec), ordered = TRUE)) %>%
    mutate(PositiveType = factor(.data$PositiveType, levels = c("TRUE", "FALSE"), ordered = TRUE)) %>%
    mutate(Stat = ifelse(Stat == "Ar1", "AR(1)", "SD")) %>%
    mutate(Stat = factor(.data$Stat, levels = c("SD", "AR(1)"), ordered = TRUE))

  out_plot <- qd_alarm_rates %>%
    ggplot(aes(x = Variable, y = Rate, fill = PositiveType)) +
    geom_bar(stat = "identity", position = "dodge", color = "black") +
    facet_grid(cols = vars(Stat)) +
    theme_bw() +
    scale_fill_manual(legend_title, values = c("FALSE" = "#FFE77AFF", "TRUE" = "#5F9260FF")) +
    labs(x = "Variable", y = "Rate (alarms / day)") +
    theme(
      legend.position = c(0.82, 0.76),
      legend.background = element_blank(),
      legend.key = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 18),
      axis.text.y = element_text(size = 12),
      axis.text.x = element_text(size = 16),
      axis.title = element_text(size = 16),
      strip.text = element_text(size = 16)
    )
  if (!is.null(y_lim)) {
    out_plot <- out_plot + lims(y = y_lim)
  }
  if (!is.null(plot_title)) {
    out_plot <- out_plot +
      ggtitle(plot_title)
  }
  return(out_plot)
}
