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
#' plot_fig1()
#'
#' # just latest experiment
#' plot_fig1(time_series = ts_data %>% dplyr::filter(Year >= 2018),
#'    bloom_fert_df = bloom_fert_dates %>% dplyr::filter(Lake == "R" & Year >= 2018))
plot_fig1 <- function(time_series = ts_data,
                      bloom_fert_df = bloom_fert_dates,
                      var_rename_vec = c(`Chl-a (ug/L)` = "Manual_Chl", `BGA (cells/mL)` = "BGA_HYLB", `D.O. sat. (%)` = "DO_Sat", pH = "pH"),
                      lakes_rename_vec = c("Reference" = "L", "Experimental" = "R"),
                      lake_colors = c("Experimental" = "firebrick3", "Reference"= "royalblue3"),
                      legend_location = c(0.13, 0.9)){
  # rename variables and pivot time series longer
  ts_formatted = time_series %>%
    dplyr::rename(!!!var_rename_vec) %>%
    tidyr::pivot_longer(cols = names(var_rename_vec), names_to = "Variable", values_to = "Value") %>%
    dplyr::mutate(Variable = factor(Variable, levels = names(var_rename_vec), ordered = TRUE))

  # fix lake names
  if(all(unique(ts_formatted$Lake) %in% lakes_rename_vec)){
    lakes_rename_vec_rev = names(lakes_rename_vec)
    names(lakes_rename_vec_rev) = unname(lakes_rename_vec)
    ts_formatted = ts_formatted %>%
      dplyr::mutate(Lake = unname(lakes_rename_vec_rev[.data$Lake]))
  }

  # format bloom dates
  dates_formatted = bloom_fert_df %>%
    dplyr::select(Year, Lake, fertStartDOY, bloomStartDOY) %>%
    dplyr::rename(Nutrients = fertStartDOY, Bloom = bloomStartDOY) %>%
    tidyr::pivot_longer(cols = c("Nutrients", "Bloom"), names_to="Start Of", values_to = "DOY") %>%
    dplyr::mutate(`Start Of` = factor(`Start Of`, levels = c("Nutrients", "Bloom"), ordered = TRUE))

  # fix lake names
  if(all(unique(dates_formatted$Lake) %in% lakes_rename_vec)){
    lakes_rename_vec_rev = names(lakes_rename_vec)
    names(lakes_rename_vec_rev) = unname(lakes_rename_vec)
    dates_formatted = dates_formatted %>%
      dplyr::mutate(Lake = unname(lakes_rename_vec_rev[.data$Lake]))
  }

  # make the plot
  p1 = ts_formatted %>%
    ggplot(aes(x=DOY, y=Value, color=Lake)) +
    geom_line(size=2) +
    facet_grid(rows=vars(Variable), cols=vars(Year), scales="free_y") +
    theme_bw() +
    labs(x = "Day of Year", y="") +
    # add legend
    theme(legend.position=legend_location, legend.background = element_blank(), legend.key = element_blank())

    # do some conditional modifications to the plot
    if(length(unique(ts_formatted$Lake)) == length(lake_colors) & all(names(lake_colors) %in% unique(ts_formatted$Lake))){
      p1 = p1 +
        scale_color_manual(values = lake_colors)
    }
    # add bloom lines
    colorMap = ifelse(length(unique(dates_formatted$Lake)) > 1, "Lake", "'black'")
    if(length(unique(dates_formatted$Lake)) > 1){
      p1 = p1 +
        geom_vline(data = dates_formatted, aes(xintercept=DOY, linetype=`Start Of`, color=Lake)) +
        scale_linetype_manual(values=c("Bloom" = "solid", "Nutrients" = "dashed"),
                              guide = guide_legend(title.position = "top", direction="horizontal")
                              )
    }else{
      p1 = p1 +
        geom_vline(data = dates_formatted, aes(xintercept=DOY, linetype=`Start Of`)) +
        scale_linetype_manual(values=c("Bloom" = "solid", "Nutrients" = "dashed"),
                              guide = guide_legend(title.position = "top", direction="horizontal")
        )
    }


    return(p1)
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
#' use_vars = c("Manual_Chl", "BGA_HYLB", "DO_Sat", "pH")
#' rw_stats = calc_rolling_stats(ts_data, var_cols=use_vars)
#' qd_s = run_qd(rw_stats, var_cols=use_vars)
#' qd_a = format_qd(qd_s, bloom_fert_df = bloom_fert_dates)
#' plot_fig3(rolling_window_stats = rw_stats, qd_alarms = qd_a, bloom_fert_df = bloom_fert_dates)
plot_fig3 <- function(rolling_window_stats, qd_alarms, bloom_fert_df, var_rename_vec = c("Chl-a" = "Manual_Chl", "BGA" = "BGA_HYLB", "D.O. sat." = "DO_Sat", "pH" = "pH"), legend_location = c(0.13, 0.9)){
  # rename variables
  rolling_window_stats = rolling_window_stats %>% rename(!!!var_rename_vec)
  # pivot long for plotting
  rw_stats_long = rolling_window_stats %>%
    filter(Stat %in% c("SD", "Ar1")) %>%
    tidyr::pivot_longer(cols = names(var_rename_vec), names_to = "Variable", values_to = "Value") %>%
    mutate(Variable = factor(.data$Variable, levels = c("Chl-a", "BGA", "D.O. sat.", "pH"), ordered = TRUE))
  # fromate bloom dates for plotting
  bloom_dates = bloom_fert_df %>%
    tidyr::pivot_longer(cols = c("fertStartDOY", "fertEndDOY", "bloomStartDOY"), names_to = "Start Of", values_to = "DOYtrunc") %>%
    mutate(
      `Start Of` = ifelse(`Start Of` == "fertStartDOY", "nutrients", `Start Of`),
      `Start Of` = ifelse(`Start Of` == "bloomStartDOY", "bloom", `Start Of`)
      ) %>%
    filter(`Start Of` %in% c("nutrients", "bloom"))
  # format alarms and get just the alarms DOYS
  var_rename_vec_rev = names(var_rename_vec)
  names(var_rename_vec_rev) = unname(var_rename_vec)
  qd_alarms_true = qd_alarms %>%
    filter(!is.na(Alarm) & Alarm == 1) %>%
    mutate(Variable = unname(var_rename_vec_rev[.data$Var]))  %>%
    mutate(Variable = factor(.data$Variable, levels = c("Chl-a", "BGA", "D.O. sat.", "pH"), ordered = TRUE))

  ## make SD plot
  # plot rw stats
  sd_plot = rw_stats_long %>%
    filter(Stat == "SD") %>%
    ggplot(aes(x=DOYtrunc, y=Value, color=Lake)) +
    geom_line(size=1.25) +
    theme_bw() +
    facet_grid(rows=vars(Variable), cols=vars(Year), scales="free_y") +
    scale_color_manual(values=c("R" = "firebrick3", "L" = "royalblue3"), guide = FALSE) +
    labs(x="Day of Year", y="") +
    # add alarms
    geom_point(
      data=qd_alarms_true %>% filter(Stat == "SD"),
      mapping=aes(fill=AlarmType),
      shape=21,
      size=3,
      color="black") +
    scale_fill_manual(values = c("#FFE77AFF", "#5F9260FF", "darkgrey")) +
    ggtitle("SD") +
    theme(legend.background = element_blank(), legend.key = element_blank(), legend.position = legend_location, plot.title=element_text(hjust=0.5)) +
    # add fert/bloom lines
    geom_vline(data = bloom_dates, aes(xintercept=DOYtrunc, linetype=`Start Of`)) +
    scale_linetype_manual(breaks=c("nutrients", "bloom"), values = c("dashed", "solid"), guide = FALSE)

  ## make AR1 plot
  # plot rw stats
  ar1_plot = rw_stats_long %>% 
    filter(Stat == "Ar1") %>%
    ggplot(aes(x=DOYtrunc, y=Value, color=Lake)) +
    geom_line(size=1.25) +
    theme_bw() +
    facet_grid(rows=vars(Variable), cols=vars(Year)) +
    scale_color_manual(values=c("R" = "firebrick3", "L" = "royalblue3"), guide = FALSE) +
    labs(x="Day of Year", y="") +
    # add alarms
    geom_point(
      data=qd_alarms_true %>% filter(Stat == "Ar1"),
      mapping=aes(fill=AlarmType),
      shape=21,
      size=3,
      color="black") +
    scale_fill_manual(values = c("False" = "#FFE77AFF", "True" = "#5F9260FF", "Late" = "darkgrey"), guide=FALSE) +
    ggtitle("Ar1") +
    theme(legend.background = element_blank(), legend.key = element_blank(), legend.position = legend_location, plot.title=element_text(hjust=0.5)) +
  # add fert/bloom lines
  geom_vline(data = bloom_dates, aes(xintercept=DOYtrunc, linetype=`Start Of`)) +
    scale_linetype_manual(breaks=c("nutrients", "bloom"), values = c("dashed", "solid"), guide = FALSE)


  ## combine
  out_fig = gridExtra::grid.arrange(sd_plot, ar1_plot, ncol=2)
  return(out_fig)
}

#' Plot true and false positive quickest detection alarm rates
#'
#' @param qd_alarm_rates data frame, output from call to \code{\link{calc_alarm_rates}}
#' @param var_rename_vec named vector, used to change variable names displayed in plot, format is c("new name" = "original name")
#' @param y_lim numeric vector of length 2, optionally specify y limits manually
#' @param title character string, optional main title for plot
#'
#' @return ggplot2 object, bar plot giving the true and false positive rates
#' @export
#' @import dplyr
#' @import ggplot2
#'
#' @examples
#' use_vars = c("Manual_Chl", "BGA_HYLB", "DO_Sat", "pH")
#' rw_data = ts_data %>% dplyr::filter(Year >= 2011)
#' rw_stats = calc_rolling_stats(rw_data, var_cols=use_vars)
#' qd_stats = run_qd(rw_stats, var_cols=use_vars)
#' qd_alarms = format_qd(qd_stats, bloom_fert_df = bloom_fert_dates)
#' qd_rates = calc_alarm_rates(qd_alarms)
#' plot_fig5(qd_alarm_rates = qd_rates)
plot_fig5 <- function(qd_alarm_rates, var_rename_vec = c("Chl-a" = "Manual_Chl", "BGA" = "BGA_HYLB", "D.O. sat." = "DO_Sat", "pH" = "pH"), y_lim=NULL, title="Positive Alarm Rat"){
  # format data
  var_rename_vec_rev = names(var_rename_vec)
  names(var_rename_vec_rev) = unname(var_rename_vec)
  qd_alarm_rates = qd_alarm_rates %>%
    mutate(Variable = unname(var_rename_vec_rev[.data$Var]))  %>%
    mutate(Variable = factor(.data$Variable, levels = c("Chl-a", "BGA", "D.O. sat.", "pH"), ordered = TRUE)) %>%
    mutate(PositiveType = factor(.data$PositiveType, levels = c("TRUE", "FALSE"), ordered = TRUE))

  out_plot = qd_alarm_rates %>%
    ggplot(aes(x=Variable, y=Rate, fill=PositiveType)) +
    geom_bar(stat="identity", position="dodge", color="black") +
    facet_grid(cols=vars(Stat)) +
    theme_bw() +
    scale_fill_manual(title, values=c("FALSE" = "#FFE77AFF", "TRUE" = "#5F9260FF")) +
    labs(x="Variable", y="Rate (alarms / day)") +
    theme(legend.position=c(0.82, 0.76),
          legend.background = element_blank(),
          legend.key = element_blank(),
          strip.text = element_text(size=18))
  if(!is.null(y_lim)){
    out_plot = out_plot + lims(y = y_lim)
  }
  if(!is.null(title)){
    out_plot = out_plot +
      ggtitle(title)
  }
  return(out_plot)
}
