#' Exact density for Pearson's correlation coefficient
#'
#' @param r numeric, sample correlation
#' @param rho numeric, "true" correlation coefficient
#' @param N numeric, sample size
#'
#' @return numeric, the probability density at *r*
#' @export
#'
#' @examples
#' f_corr_coef(r = 0.5, rho = 0.7, N = 10)
f_corr_coef <- function(r, rho, N){
  numerator_part = (N-2) * gamma(N-1) * (1 - rho^2)^((N - 1)/2) * (1 - r^2)^((N - 4)/2)
  denominator_part = (2*pi)^0.5 * gamma(N - (1/2)) * (1-rho*r)^(N - (3/2))
  hypergeo_part = hypergeo::hypergeo(1/2, 1/2, (2*N - 1)/2, (rho*r + 1)/2)
  out = (numerator_part/denominator_part) * hypergeo_part
  return(Re(out))
}



#' Quickest Detection calculation
#'
#' @param rolling_window_stats data frame of rolling window statistics that is output from \code{\link{calc_rolling_stats}}
#' @param var_cols character vector, variable names that are columns of *data* containing time series to calculate quickest detection stats on
#' @param widths numeric vector, rolling window width(s) to be used, default is c(21)
#' @param stats_to_qd character vector, rolling window statistics to run quickest detection method on, defaults to all available: c("sd", "Ar1")
#' @param A_adj numeric, the threshold in the S-R statistic at which an alarm is triggered
#' @param exp_lakes character vector, values in *Lake* column of *rolling_window_stats* that correspond to experimental lakes
#' @param ref_lake character vector, values in *Lake* column of *rolling_window_stats* that correspond to the reference lake
#' @param ar1_alarm_rho numeric value between 0 and 1 (exclusive), 'true' value of lag-1 autocorrelation in alarm state, default is 0.95
#'
#' @return a data frame of values from the quickest detection method; TODO add more details on columns
#' @export
#'
#' @examples
#' rw_stats = calc_rolling_stats(ts_data, var_cols=c("DO_Sat"))
#' qd_stats = run_qd(rw_stats, var_cols=c("DO_Sat"))
run_qd <- function(rolling_window_stats, var_cols, widths = c(21), stats_to_qd = c("SD", "Ar1"), A_adj = 1E7, exp_lakes = c("R", "T"), ref_lake="L", ar1_alarm_rho = 0.95){

  # QD data frame setup
  stats_to_qd_completed = c(stats_to_qd, paste0(stats_to_qd, "_sd"))
  rolling_window_stats_Exp = rolling_window_stats[rolling_window_stats$Lake %in% c(exp_lakes) & rolling_window_stats$Stat %in% stats_to_qd_completed & rolling_window_stats$Width %in% widths, ]
  rolling_window_stats_Ref = rolling_window_stats[rolling_window_stats$Lake %in% c(ref_lake) & rolling_window_stats$Stat %in% stats_to_qd_completed & rolling_window_stats$Width %in% widths, ]
  rolling_window_stats_Exp$Exp_or_Ref = "Exp"
  rolling_window_stats_Ref$Exp_or_Ref = "Ref"

  rolling_window_stats_Exp_pivot = tidyr::pivot_wider(rolling_window_stats_Exp,
                                                      id_cols = c("Year", "DOYtrunc", "Width", "Lake"),
                                                      values_from=dplyr::all_of(var_cols),
                                                      names_from=c("Stat", "Exp_or_Ref"),
                                                      names_glue="{.value}_{Stat}_{Exp_or_Ref}")
  rolling_window_stats_Ref_pivot = tidyr::pivot_wider(rolling_window_stats_Ref,
                                                      id_cols = c("Year", "DOYtrunc", "Width"),
                                                      values_from=dplyr::all_of(var_cols),
                                                      names_from=c("Stat", "Exp_or_Ref"),
                                                      names_glue="{.value}_{Stat}_{Exp_or_Ref}")

  rolling_window_stats_forCalcs0 = merge(rolling_window_stats_Exp_pivot, rolling_window_stats_Ref_pivot)
  # get rid of rows where all variable columns are NA
  varsAllColnames = colnames(rolling_window_stats_forCalcs0)[!(colnames(rolling_window_stats_forCalcs0) %in% c("Lake", "Year", "DOYtrunc", "Width"))]
  inds_allNA = apply(rolling_window_stats_forCalcs0[, varsAllColnames], MARGIN=1, FUN=function(x) all(is.na(x)))

  rolling_window_stats_forCalcs = rolling_window_stats_forCalcs0[!inds_allNA, ]

  # create data frame to hold stats
  QD_statCols = paste(stats_to_qd, rep(c("Lambda", "R_vec", "R_squeals",  "Pooled_sd", "Alarm"), each=length(stats_to_qd)), sep="_")
  newCols = paste(rep(var_cols, each=length(QD_statCols)), QD_statCols, sep="_")
  rolling_window_stats_forCalcs[, newCols] = NA

  # iterate through each variable, do set-up calculations for QD calculation
  for(i in 1:length(var_cols)){
    curVar = var_cols[i]
    #Calculate pooled sd for SD
    if("SD" %in% stats_to_qd){
      rolling_window_stats_forCalcs[, paste0(curVar, "_SD_Pooled_sd")] = sqrt(rolling_window_stats_forCalcs[, paste0(curVar, "_SD_sd_Exp")]^2 + rolling_window_stats_forCalcs[, paste0(curVar, "_SD_sd_Ref")]^2)
    }

    #Calculate lambdas for AC
      f.L_AC_hold = rep(NA, nrow(rolling_window_stats_forCalcs))
      g.L_AC_hold = rep(NA, nrow(rolling_window_stats_forCalcs))
      inds_allInfoForHGM = !is.na(rolling_window_stats_forCalcs[, paste0(curVar, "_Ar1_Exp")]) &
        !is.na(rolling_window_stats_forCalcs[, paste0(curVar, "_Ar1_Ref")])
      f.L_AC_hold[inds_allInfoForHGM] = log(
        f_corr_coef(
          r = rolling_window_stats_forCalcs[inds_allInfoForHGM, paste0(curVar, "_Ar1_Exp")],
          rho = rolling_window_stats_forCalcs[inds_allInfoForHGM, paste0(curVar, "_Ar1_Ref")],
          N = rolling_window_stats_forCalcs[inds_allInfoForHGM, "Width"]
        )
      )
      g.L_AC_hold[inds_allInfoForHGM] = log(
        f_corr_coef(
          r = rolling_window_stats_forCalcs[inds_allInfoForHGM, paste0(curVar, "_Ar1_Exp")],
          rho = ar1_alarm_rho,
          N = rolling_window_stats_forCalcs[inds_allInfoForHGM, "Width"]
        )
      )

    rolling_window_stats_forCalcs[, paste0(curVar, "_Ar1_Lambda")] = exp(g.L_AC_hold - f.L_AC_hold)

    #calculate lambdas for SD
    f.L_SD_hold = stats::dnorm(rolling_window_stats_forCalcs[, paste0(curVar, "_SD_Exp")], mean=rolling_window_stats_forCalcs[, paste0(curVar, "_SD_Ref")], sd=rolling_window_stats_forCalcs[, paste0(curVar, "_SD_Pooled_sd")], log=TRUE)
    g.L_SD_hold = stats::dnorm(rolling_window_stats_forCalcs[, paste0(curVar, "_SD_Exp")], mean=rolling_window_stats_forCalcs[, paste0(curVar, "_SD_Ref")]+2* rolling_window_stats_forCalcs[, paste0(curVar, "_SD_Pooled_sd")], sd=rolling_window_stats_forCalcs[, paste0(curVar, "_SD_Pooled_sd")], log=TRUE)
    rolling_window_stats_forCalcs[, paste0(curVar, "_SD_Lambda")] = exp(g.L_SD_hold - f.L_SD_hold)
  }
  # order data by Lake, Width, Year, and then DOY
  rolling_window_stats_forCalcs = rolling_window_stats_forCalcs[order(rolling_window_stats_forCalcs$Lake, rolling_window_stats_forCalcs$Width, rolling_window_stats_forCalcs$Year, rolling_window_stats_forCalcs$DOYtrunc), ]
  #Compute the S-R statistic for each lake, each day that has a rolling window stat calculated for it
  widths = unique(rolling_window_stats_forCalcs$Width)
  for(w in 1:length(widths)){
    curWidth = widths[w]
    for(v in 1:length(var_cols)){
      curVar = var_cols[v]
      print(paste("Current variable QD-ing:", curVar))
      indices_Stat_notNA = !is.na(rolling_window_stats_forCalcs[, paste0(curVar, "_Ar1_Lambda")])
      years = unique(rolling_window_stats_forCalcs$Year)
      for(y in 1:length(years)){
        Lakes = exp_lakes
        for(l in 1:length(Lakes)){
          doys = rolling_window_stats_forCalcs$DOYtrunc[rolling_window_stats_forCalcs$Year == years[y] & rolling_window_stats_forCalcs$Lake == Lakes[l] & rolling_window_stats_forCalcs$Width == widths[w] & indices_Stat_notNA]
          if(length(doys)>0){
            for(d in 1:length(doys)){
              curInd = rolling_window_stats_forCalcs$Year == years[y] & rolling_window_stats_forCalcs$Lake == Lakes[l] & rolling_window_stats_forCalcs$DOYtrunc == doys[d] & rolling_window_stats_forCalcs$Width == widths[w]
              if(d == 1){
                rolling_window_stats_forCalcs[curInd, paste0(curVar, "_Ar1_R_vec")] = 0
                rolling_window_stats_forCalcs[curInd, paste0(curVar, "_Ar1_R_squeals")] = 0
                rolling_window_stats_forCalcs[curInd, paste0(curVar, "_Ar1_Alarm")] = 0
                rolling_window_stats_forCalcs[curInd, paste0(curVar, "_SD_R_vec")] = 0
                rolling_window_stats_forCalcs[curInd, paste0(curVar, "_SD_R_squeals")] = 0
                rolling_window_stats_forCalcs[curInd, paste0(curVar, "_SD_Alarm")] = 0
              }else{
                prevInd = rolling_window_stats_forCalcs$Year == years[y] & rolling_window_stats_forCalcs$Lake == Lakes[l] & rolling_window_stats_forCalcs$DOYtrunc == doys[d-1] & rolling_window_stats_forCalcs$Width == widths[w]
                #calculate S-R stat and see if alarm for AC
                cur_R_AC = (1+rolling_window_stats_forCalcs[prevInd, paste0(curVar, "_Ar1_R_vec")])*rolling_window_stats_forCalcs[curInd, paste0(curVar, "_Ar1_Lambda")]
                rolling_window_stats_forCalcs[curInd, paste0(curVar, "_Ar1_R_squeals")] = cur_R_AC
                if(!is.na(cur_R_AC) & cur_R_AC <=A_adj){
                  rolling_window_stats_forCalcs[curInd, paste0(curVar, "_Ar1_R_vec")] = cur_R_AC
                  rolling_window_stats_forCalcs[curInd, paste0(curVar, "_Ar1_Alarm")] = 0
                }else if(!is.na(cur_R_AC) & cur_R_AC  > A_adj){
                  rolling_window_stats_forCalcs[curInd, paste0(curVar, "_Ar1_R_vec")] = 0
                  rolling_window_stats_forCalcs[curInd, paste0(curVar, "_Ar1_Alarm")] = 1
                }else{
                  rolling_window_stats_forCalcs[curInd, paste0(curVar, "_Ar1_R_vec")] = rolling_window_stats_forCalcs[prevInd, paste0(curVar, "_Ar1_R_vec")]
                  rolling_window_stats_forCalcs[curInd, paste0(curVar, "_Ar1_Alarm")] = 0
                }
                #calculate S-R stat and see if alarm for SD
                cur_R_SD = (1+rolling_window_stats_forCalcs[prevInd, paste0(curVar, "_SD_R_vec")])*rolling_window_stats_forCalcs[curInd, paste0(curVar, "_SD_Lambda")]
                rolling_window_stats_forCalcs[curInd, paste0(curVar, "_SD_R_squeals")] = cur_R_SD
                if(!is.na(cur_R_SD) & cur_R_SD <=A_adj){
                  rolling_window_stats_forCalcs[curInd, paste0(curVar, "_SD_R_vec")] = cur_R_SD
                  rolling_window_stats_forCalcs[curInd, paste0(curVar, "_SD_Alarm")] = 0
                }else if(!is.na(cur_R_SD) & cur_R_SD > A_adj){
                  rolling_window_stats_forCalcs[curInd, paste0(curVar, "_SD_R_vec")] = 0
                  rolling_window_stats_forCalcs[curInd, paste0(curVar, "_SD_Alarm")] = 1
                }else{
                  rolling_window_stats_forCalcs[curInd, paste0(curVar, "_SD_R_vec")] = rolling_window_stats_forCalcs[prevInd, paste0(curVar, "_SD_R_vec")]
                  rolling_window_stats_forCalcs[curInd, paste0(curVar, "_SD_Alarm")] = 0
                }
              }
            }
          }
        }
      }
    }
  }
  return(rolling_window_stats_forCalcs)
}


#' Format quickest detection results for plotting
#'
#' @param qd_stats data frame, output from \code{\link{run_qd}}
#' @param var_cols character vector, variable names that are columns of *data* containing time series to calculate rolling window stats on
#' @param stats character vector, statistics to format alarms for, default = c("SD", "Ar1)
#' @param bloom_fert_df data frame of fertilization start and end dates, see \code{\link{bloom_fert_dates}} for default and formatting
#' @param exp_lakes character vector, values in *Lake* column of qd_stats that correspond to experimental lakes
#' @param just_alarm_DOYs TRUE or FALSE (default), should results be filtered to just DOYs that have an alarm (in any variable)?
#'
#' @return data frame with quickest detection alarm info
#' @export
#'
#' @examples
#' rw_stats = calc_rolling_stats(ts_data, var_cols=c("DO_Sat"))
#' qd_stats = run_qd(rw_stats, var_cols=c("DO_Sat"))
#' qd_alarms = format_qd(qd_stats, var_cols=c("DO_Sat"))
format_qd <- function(qd_stats, var_cols = c("Manual_Chl", "BGA_HYLB", "DO_Sat", "pH"), stats=c("SD", "Ar1"), bloom_fert_df = bloom_fert_dates, exp_lakes = c("T", "R"), just_alarm_DOYs = FALSE){
  # transform QD df to make it easy to add points
  alarmCols = paste0(rep(var_cols, times=length(stats)), rep(paste0("_", stats, "_Alarm"), each = length(var_cols)))
  if(just_alarm_DOYs == TRUE){
    inds_anyAlarm = apply(qd_stats[, alarmCols], MARGIN = 1, FUN=function(x) !all(is.na(x)) & any(x==1, na.rm=TRUE))
    QD_long0 = as.data.frame(
      tidyr::pivot_longer(qd_stats[inds_anyAlarm, ], cols = dplyr::starts_with(var_cols), names_to="Var_Stat", values_to = "Value")
    )
  }else{
    QD_long0 = as.data.frame(
      tidyr::pivot_longer(qd_stats, cols = dplyr::starts_with(var_cols), names_to="Var_Stat", values_to = "Value")
    )
  }

  # separate variables and stats into separate columns
  vars_to_remove = paste(paste0(var_cols, "_"), collapse="|")
  QD_long0 = QD_long0 %>%
    dplyr::mutate(Stat = gsub(pattern=vars_to_remove, replacement="", x = .data$Var_Stat))
  stats_to_remove = paste(paste0("_", unique(QD_long0$Stat)), collapse="|")
  QD_long0 = QD_long0 %>%
    dplyr::mutate(Var = gsub(pattern=stats_to_remove, replacement="", x = .data$Var_Stat))

  # pull out just the stats we want for plotting
  QD_long = QD_long0[QD_long0[, "Stat"] %in% c(paste0(stats, "_Exp"), paste0(stats, "_Alarm")), c("Year", "Lake", "DOYtrunc", "Width", "Var", "Stat", "Value")]
  Alarms = QD_long[grep(pattern="Alarm", QD_long$Stat) & !is.na(QD_long[, "Value"]) & QD_long[, "Value"] %in% c(0,1), ]
  Alarms[, "Stat"] = gsub(pattern="_Alarm", replacement="", x=Alarms[, "Stat"])
  colnames(Alarms)[colnames(Alarms) == "Value"] = "Alarm"
  statsAll = QD_long[QD_long[, "Stat"] %in% paste0(stats, "_Exp"), ]
  statsAll[, "Stat"] = gsub(pattern="_Exp", replacement="", x=statsAll[, "Stat"])

  # combine stats and alarms
  statsAlarms = dplyr::full_join(statsAll, Alarms)
  statsAlarms$Stat = factor(statsAlarms$Stat, levels=stats)

  # classify alarms
  statsAlarms = dplyr::left_join(statsAlarms, bloom_fert_df)
  statsAlarms = statsAlarms %>%
    dplyr::mutate(fertStartDOY = ifelse(is.na(fertStartDOY), 367, fertStartDOY)) %>%
    dplyr::mutate(bloomStartDOY = ifelse(is.na(bloomStartDOY), 368, bloomStartDOY))

  statsAlarms$AlarmType = "NA"
  statsAlarms = statsAlarms %>%
    dplyr:: mutate(AlarmType = ifelse(DOYtrunc < fertStartDOY & !is.na(Alarm) & Alarm == 1, "False", AlarmType)) %>%
    dplyr::mutate(AlarmType = ifelse(DOYtrunc >= fertStartDOY & DOYtrunc < bloomStartDOY & !is.na(Alarm) & Alarm == 1, "True", AlarmType)) %>%
    dplyr::mutate(AlarmType = ifelse(DOYtrunc >= bloomStartDOY & !is.na(Alarm) & Alarm == 1, "Late", AlarmType))

  statsAlarms$AlarmType = factor(statsAlarms$AlarmType, levels = c("False", "True", "Late"), ordered=TRUE)

  statsAlarms_out = statsAlarms %>%
    dplyr::filter(!is.na(Value) & !is.na(Alarm)) %>%
    dplyr::select(-c(fertStartDOY, fertEndDOY, bloomStartDOY))

  return(statsAlarms_out)
}



#' Plot rolling window statistics and qd alarms
#'
#' @param rolling_window_stats data frame of rolling window statistics that is output from \code{\link{calc_rolling_stats}}
#' @param qd_alarms data frame of quickest detection alarms, output from \code{\link{format_qd}}
#' @param var_col character vector, variable names that are columns of *data* containing time series to calculate rolling window stats on
#' @param stats statistics to plot, default and options are c("SD", "Ar1")
#' @param title title to put on the plot
#' @param widths numeric vector, rolling window width(s) to be used, default is c(21)
#' @param plot_bloom_lines TRUE or FALSE (default), should lines for the bloom dates be plotted?
#' @param bloom_dates if plot_bloom_lines = TRUE, this data frame is used to determine where
#'
#' @return a ggplot object
#' @import ggplot2
#' @export

#'
#' @examples
#' use_vars = c("Manual_Chl", "BGA_HYLB", "DO_Sat", "pH")
#' rw_stats = calc_rolling_stats(ts_data, var_cols=use_vars)
#' qd_stats = run_qd(rw_stats, var_cols=use_vars)
#' qd_alarms = format_qd(qd_stats, bloom_fert_df = bloom_fert_dates)
#' plot_qd(rw_stats, qd_alarms, use_vars[1], stats = c("SD", "Ar1"),
#'    plot_bloom_lines = TRUE, bloom_dates = bloom_fert_dates)
plot_qd <- function(rolling_window_stats, qd_alarms, var_col, stats=c("SD", "Ar1"), title="R.W. width = 21", widths=c(21), plot_bloom_lines = FALSE, bloom_dates){

  # get just the alarms
  statsAlarms = qd_alarms %>%
    dplyr::filter(!is.na(Alarm) & Alarm == 1 & Var == var_col)


  # add the points
  p = ggplot(rolling_window_stats[rolling_window_stats$Stat %in% stats & rolling_window_stats$Width %in% widths, ], aes(x=DOYtrunc, y=get(var_col), color=Lake)) +
    facet_grid(cols=vars(Year), rows=vars(Stat), scales="free_y") +
    geom_line(size=1.25) +
    geom_point(data=statsAlarms, mapping=aes(x=DOYtrunc, y=Value, color=Lake, fill=AlarmType), shape=21, size=3, color="black") + # , fill="red", color="black"
    scale_fill_manual(values = c("#FFE77AFF", "#5F9260FF", "darkgrey")) +
    # geom_point(data=statsAlarms, aes(x=DOYtrunc, y=Value, color=Lake), size = 2) +
    theme_bw(base_size=16) + labs(y=var_col, x="Day of Year")
  if(length(unique(rolling_window_stats$Lake)) == 2){
    p = p + scale_color_manual(values=c("dodgerblue", "orangered3"))
  }else if(length(unique(rolling_window_stats$Lake)) == 3){
    p = p + scale_color_manual(values=c("dodgerblue", "orangered3", "grey10"))
  }
  if(plot_bloom_lines){
    p = p + geom_vline(data = bloom_dates, aes(xintercept=bloomStartDOY, color=Lake), linetype="solid") +
      geom_vline(data = bloom_dates, aes(xintercept=fertStartDOY), color="black", linetype="dashed")
  }
  if(!is.na(title)){
    p = p + ggtitle(title)
  }
  return(p)
}


#' Calculate rates of true and false alarms
#'
#' @param qd_alarms data frame of quickest detection alarms, output from \code{\link{format_qd}}
#' @param bloom_fert_df data frame of fertilization start and end dates, see \code{\link{bloom_fert_dates}} for default and formatting
#'
#' @return data frame with true and false positive rates for combinations of variable and statistic
#' @export
#'
#' @examples
#' use_vars = c("Manual_Chl", "BGA_HYLB", "DO_Sat", "pH")
#' rw_stats = calc_rolling_stats(ts_data, var_cols=use_vars)
#' qd_stats = run_qd(rw_stats, var_cols=use_vars)
#' qd_alarms = format_qd(qd_stats, bloom_fert_df = bloom_fert_dates)
#' calc_alarm_rates(qd_alarms)
#'
calc_alarm_rates <- function(qd_alarms, bloom_fert_df = bloom_fert_dates){
  # add DOYs of fert and bloom start
  qd_alarms_wDates = dplyr::left_join(qd_alarms, bloom_fert_df) %>%
    dplyr::mutate(fertStartDOY = ifelse(is.na(fertStartDOY), 367, fertStartDOY)) %>%
    dplyr::mutate(bloomStartDOY = ifelse(is.na(bloomStartDOY), 368, bloomStartDOY))

  # assign DOYs to time periods
  qd_alarms_wDates$TimePeriod = NA
  qd_alarms_wDates = qd_alarms_wDates %>%
    dplyr::mutate(TimePeriod = ifelse(DOYtrunc < fertStartDOY, "preFert", TimePeriod)) %>%
    dplyr::mutate(TimePeriod = ifelse(DOYtrunc >= fertStartDOY & DOYtrunc < bloomStartDOY, "alarmPeriod", TimePeriod)) %>%
    dplyr::mutate(TimePeriod = ifelse(DOYtrunc >= bloomStartDOY, "postBloom", TimePeriod))

  # calculate rates of positive alarms
  qd_alarm_rates = qd_alarms_wDates %>%
    dplyr::filter(TimePeriod != "postBloom") %>%
    dplyr::group_by(Var, Stat, TimePeriod) %>%
    dplyr::summarise(nAlarms = sum(Alarm == 1), nDays = sum(Alarm %in% c(0,1))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Rate = nAlarms / nDays, PositiveType = ifelse(TimePeriod == "preFert", FALSE, TRUE))
}
