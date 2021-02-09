#' Plot rolling window stats and quickest detection alarms for all variables
#'
#' @param rolling_window_stats data frame, output from call to
#' @param qd_alarms data frame, output from call to
#' @param bloom_dates data frame, formatted as
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
#' plot_fig3(rolling_window_stats = rw_stats, qd_alarms = qd_a, bloom_dates = bloom_fert_dates)
plot_fig3 <- function(rolling_window_stats, qd_alarms, bloom_dates, rename_vec = c("Chl-a" = "Manual_Chl", "BGA" = "BGA_HYLB", "D.O. sat." = "DO_Sat", "pH" = "pH")){
  # rename variables
  rolling_window_stats = rolling_window_stats %>% rename(!!!rename_vec)
  # pivot long for plotting
  rw_stats_long = rolling_window_stats %>%
    filter(Stat %in% c("SD", "Ar1")) %>%
    tidyr::pivot_longer(cols = names(rename_vec), names_to = "Variable", values_to = "Value") %>%
    mutate(Variable = factor(.data$Variable, levels = c("Chl-a", "BGA", "D.O. sat.", "pH"), ordered = TRUE))
  # just get the Peter bloom dates for the experiment
  bloom_dates_peter = bloom_dates %>%
    filter(Lake == "R" & Year >= 2018) %>%
    tidyr::pivot_longer(cols = c("fertStartDOY", "fertEndDOY", "bloomStartDOY"), names_to = "Start Of", values_to = "DOYtrunc") %>%
    mutate(
      `Start Of` = ifelse(`Start Of` == "fertStartDOY", "nutrients", `Start Of`),
      `Start Of` = ifelse(`Start Of` == "bloomStartDOY", "bloom", `Start Of`)
      ) %>%
    filter(`Start Of` %in% c("nutrients", "bloom"))
  # format alarms and get just the alarms DOYS
  rename_vec_rev = names(rename_vec)
  names(rename_vec_rev) = unname(rename_vec)
  qd_alarms_true = qd_alarms %>%
    filter(!is.na(Alarm) & Alarm == 1) %>%
    mutate(Variable = unname(rename_vec_rev[.data$Var]))  %>%
    mutate(Variable = factor(.data$Variable, levels = c("Chl-a", "BGA", "D.O. sat.", "pH"), ordered = TRUE))

  ## make SD plot
  # plot rw stats
  sd_plot = rw_stats_long %>%
    filter(Year >= 2018 & Stat == "SD") %>%
    ggplot(aes(x=DOYtrunc, y=Value, color=Lake)) +
    geom_line(size=1.25) +
    theme_bw() +
    facet_grid(rows=vars(Variable), cols=vars(Year), scales="free_y") +
    scale_color_manual(values=c("R" = "firebrick3", "L" = "royalblue3"), guide = FALSE) +
    labs(x="Day of Year", y="") +
    # add alarms
    geom_point(
      data=qd_alarms_true %>% filter(Year >= 2018 & Stat == "SD"),
      mapping=aes(fill=AlarmType),
      shape=21,
      size=3,
      color="black") +
    scale_fill_manual(values = c("#FFE77AFF", "#5F9260FF", "darkgrey")) +
    ggtitle("SD") +
    theme(legend.background = element_blank(), legend.key = element_blank(), legend.position = c(0.12, 0.92), plot.title=element_text(hjust=0.5)) +
    # add fert/bloom lines
    geom_vline(data = bloom_dates_peter, aes(xintercept=DOYtrunc, linetype=`Start Of`)) +
    scale_linetype_manual(breaks=c("nutrients", "bloom"), values = c("dashed", "solid"), guide = FALSE)

  ## make AR1 plot
  # plot rw stats
  ar1_plot = rw_stats_long %>%
    filter(Year >= 2018 & Stat == "Ar1") %>%
    ggplot(aes(x=DOYtrunc, y=Value, color=Lake)) +
    geom_line(size=1.25) +
    theme_bw() +
    facet_grid(rows=vars(Variable), cols=vars(Year)) +
    scale_color_manual(values=c("R" = "firebrick3", "L" = "royalblue3"), guide = FALSE) +
    labs(x="Day of Year", y="") +
    # add alarms
    geom_point(
      data=qd_alarms_true %>% filter(Year >= 2018 & Stat == "Ar1"),
      mapping=aes(fill=AlarmType),
      shape=21,
      size=3,
      color="black") +
    scale_fill_manual(values = c("False" = "#FFE77AFF", "True" = "#5F9260FF", "Late" = "darkgrey"), guide=FALSE) +
    ggtitle("Ar1") +
    theme(legend.background = element_blank(), legend.key = element_blank(), legend.position = c(0.12, 0.92), plot.title=element_text(hjust=0.5)) +
  # add fert/bloom lines
  geom_vline(data = bloom_dates_peter, aes(xintercept=DOYtrunc, linetype=`Start Of`)) +
    scale_linetype_manual(breaks=c("nutrients", "bloom"), values = c("dashed", "solid"), guide = FALSE)


  ## combine
  out_fig = gridExtra::grid.arrange(sd_plot, ar1_plot, ncol=2)
  return(out_fig)
}

#' Plot true and false positive quickest detection alarm rates
#'
#' @param qd_alarm_rates data frame, output from call to
#'
#' @return
#' @export
#'
#' @examples
#' # need to add
plot_fig5 <- function(qd_alarm_rates){

}
