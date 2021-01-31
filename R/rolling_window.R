#' Mean with detrending
#'
#' Calculate the mean of values with optional (linear) detrending, and perform check that there's enough data
#'
#' @param x numeric
#' @param detrend TRUE or FALSE, default = FALSE
#' @param prop_req numeric value between 0 and 1, default = 0.99
#'
#' @return numeric
#' @export
#'
#' @examples
#' y = rnorm(10)
#' Mean(y)
#' Mean(y, detrend = TRUE)
#' y[1] = NA
#' Mean(y)
#' Mean(y, prop_req = 0.5)
Mean <- function(x, detrend = FALSE, prop_req = 0.99){
  N = sum(!is.na(x))
  if((N / length(x)) < prop_req){
    warning("Warning: more than threshold of rolling window is NA, returning NA")
    Mean = NA
  }else{
    if(detrend == TRUE){
      x_t = 1:length(x)
      lm_x = stats::lm(x~x_t, na.action="na.exclude")
      xd = stats::residuals(lm_x)
    }else{
      xd = x
    }
    Mean = mean(xd, na.rm=TRUE)
  }
  return(Mean)
}



#' Lag-1 autocorrelation with detrending
#'
#' @param x numeric
#' @param detrend TRUE or FALSE, default = FALSE
#' @param prop_req numeric value between 0 and 1, default = 0.99
#'
#' @return numeric
#' @export
#' @importFrom stats na.pass
#'
#' @examples
#' y = rnorm(10)
#' Ar1(y)
#' Ar1(y, detrend = TRUE)
#' y[1] = NA
#' Ar1(y)
#' Ar1(y, prop_req = 0.5)
Ar1 <-function(x, detrend=FALSE, prop_req = 0.99){
  N = sum(!is.na(x))
  if((N / length(x)) < prop_req){
    warning("Warning: more than threshold of rolling window is NA, returning NA")
    ac = NA
    ac_sd = NA
  }else{
    if(detrend == TRUE){
      x_t = 1:length(x)
      lm_x = stats::lm(x~x_t, na.action="na.exclude")
      xd = stats::residuals(lm_x)
    }else{
      xd = x
    }
    ac = stats::acf(xd, lag=1, na.action=na.pass, plot=FALSE)$acf[2]
  }
  return(ac)
}

#' SD of lag-1 autocorrelation with detrending
#'
#' @param x numeric
#' @param detrend TRUE or FALSE, default = FALSE
#' @param prop_req numeric value between 0 and 1, default = 0.99
#'
#' @return numeric
#' @export
#' @importFrom stats na.pass
#'
#' @examples
#' y = rnorm(10)
#' Ar1_sd(y)
#' Ar1_sd(y, detrend = TRUE)
#' y[1] = NA
#' Ar1_sd(y)
#' Ar1_sd(y, prop_req = 0.5)
Ar1_sd <-function(x, detrend=FALSE, prop_req = 0.99){
  N = sum(!is.na(x))
  if((N / length(x)) < prop_req){
    warning("Warning: more than threshold of rolling window is NA, returning NA")
    ac = NA
    ac_sd = NA
  }else{
    if(detrend == TRUE){
      x_t = 1:length(x)
      lm_x = stats::lm(x~x_t, na.action="na.exclude")
      xd = stats::residuals(lm_x)
    }else{
      xd = x
    }
    ac = stats::acf(xd, lag=1, na.action=na.pass, plot=FALSE)$acf[2]
    var_ac = ((1 - ac^2)^2)/(sum(!is.na(x)) - 3)
    ac_sd = sqrt(var_ac)
  }
  return(ac_sd)
}

#' SD with detrending
#'
#' Calculate the standard deviation of values with optional (linear) detrending, and perform check that there's enough data
#'
#' @param x numeric
#' @param detrend TRUE or FALSE, default = FALSE
#' @param prop_req numeric value between 0 and 1, default = 0.99
#'
#' @return numeric
#' @export
#'
#' @examples
#' y = rnorm(10)
#' SD(y)
#' SD(y, detrend = TRUE)
#' y[1] = NA
#' SD(y)
#' SD(y, prop_req = 0.5)
SD <- function(x, detrend=FALSE, prop_req = 0.99){
  N = sum(!is.na(x))
  if((N / length(x)) < prop_req){
    warning("Warning: more than threshold of rolling window is NA, returning NA")
    sd.out = NA
    sd.sd = NA
  }else{
    if(detrend == TRUE){
      x_t = 1:length(x)
      lm_x = stats::lm(x~x_t, na.action="na.exclude")
      xd = stats::residuals(lm_x)
    }else{
      xd = x
    }
    n.x = sum(!is.na(xd))
    sd.out = stats::sd(xd, na.rm = TRUE)
  }
  return(sd.out)
}

#' SD of SD autocorrelation with detrending
#'
#' @param x numeric
#' @param detrend TRUE or FALSE, default = FALSE
#' @param prop_req numeric value between 0 and 1, default = 0.99
#'
#' @return numeric
#' @export
#'
#' @examples
#' y = rnorm(10)
#' SD_sd(y)
#' SD_sd(y, detrend = TRUE)
#' y[1] = NA
#' SD_sd(y)
#' SD_sd(y, prop_req = 0.5)
SD_sd <- function(x, detrend=FALSE, prop_req = 0.99){
  N = sum(!is.na(x))
  if((N / length(x)) < prop_req){
    warning("Warning: more than threshold of rolling window is NA, returning NA")
    sd.out = NA
    sd.sd = NA
  }else{
    if(detrend == TRUE){
      x_t = 1:length(x)
      lm_x = stats::lm(x~x_t, na.action="na.exclude")
      xd = stats::residuals(lm_x)
    }else{
      xd = x
    }
    n.x = sum(!is.na(xd))
    sd.out = stats::sd(xd, na.rm = TRUE)
    kurt.x = moments::kurtosis(xd, na.rm=TRUE)
    var.s2 = (sd.out^4)*((2/(n.x-1)) + (kurt.x / n.x))
    var.sd = ( 1/(4*sd.out*sd.out))*var.s2
    sd.sd = sqrt(var.sd)
  }
  return(sd.sd)
}


#' Calculate rolling window statistics
#'
#' Calculate rolling window statistics for multiple combinations of statistics, variables, rolling window widths, and detrending options.
#'
#' @param data data frame containing time series to calculate rolling window stats on
#' @param var_cols character vector, variable names that are columns of *data* containing time series to calculate rolling window stats on
#' @param id_cols character vector, two columns specifying "groups" to separate data into before calculating stats independently on, defaults are "Lake" and "Year"#'
#' @param time_col character, column name that defines time step of data within the *id_cols* groupings
#' @param statistics character vector, rolling window statistics to calculate, defaults to all available: c("Mean", "SD", "SD_sd", "Ar1", "Ar1_sd")
#' @param vars_to_log character vector, optional variables names to log10 transform and create new columns before calculating rolling window stats
#' @param log_offset positive numeric value, if any variables are being log transformed that contain negative values, all values will be shifted so that the minimum (before transformation) equals this value
#' @param widths numeric vector, rolling window width(s) to be used, default is c(21)
#' @param detrend TRUE or FALSE, should data within rolling windows be linearly detrended before calculating *statistics*
#' @param min_prop numeric value between 0 and 1, proportion of data within a rolling window to calculate rolling window statistic, default 0.99

#'
#' @return a data frame with the results for the specified stats, widths, variables, etc.
#' @export
#'
#' @examples
#' rw_stats = calc_rolling_stats(data = ts_data, var_cols = c("DO_Sat"))
calc_rolling_stats <- function(data, var_cols, id_cols=c("Lake", "Year"), time_col="DOY", statistics=c("Mean", "SD", "SD_sd", "Ar1", "Ar1_sd"), vars_to_log=NULL, log_offset=.1, widths=c(21), detrend=FALSE, min_prop=0.99){
  if(!is.null(vars_to_log)){
    data[, paste0(vars_to_log, "_log10")] = NA
    for(v in 1:length(vars_to_log)){
      min_val = min(data[, vars_to_log[v]], na.rm=TRUE)
      if(min_val < 0){
        add_offset = abs(min_val) + log_offset
      }else{
        add_offset = 0
      }
      data[, paste0(vars_to_log[v], "_log10")] = log10(data[, vars_to_log[v]] + add_offset)
    }
    var_cols = c(var_cols, paste0(vars_to_log, "_log10"))
  }

  lake_years = unique(data[, id_cols])
  stats_combinations = expand.grid(Stat=statistics, Width=widths, Detrend=detrend, stringsAsFactors = FALSE)
  all_lakeyear_all_stats = list()
  for(i in 1:nrow(lake_years)){
    ly_data = data %>%
      dplyr::filter(!!as.symbol(id_cols[1]) == dplyr::pull(lake_years[i, id_cols[1]]) & !!as.symbol(id_cols[2]) == dplyr::pull(lake_years[i, id_cols[2]])) %>%
      dplyr::select(dplyr::all_of(c(time_col, var_cols))) %>%
      dplyr::arrange(!!as.symbol(time_col))

      # data[, id_cols[1]] == lake_years[i, id_cols[1]] & data[, id_cols[2]] == lake_years[i, id_cols[2]] # TODO: generalize this to work with any number of columns that need to match; include check that all columns are there
    ly_stats_list = list()
    for(s in 1:nrow(stats_combinations)){
      hold_results_lys = as.data.frame(
        zoo::rollapplyr(
          data=ly_data %>% dplyr::select(dplyr::all_of(var_cols)) , #%>% dplyr::pull()
          width=stats_combinations[s, "Width"],
          FUN=get(stats_combinations[s, "Stat"]),
          detrend=stats_combinations[s, "Detrend"],
          prop_req=min_prop,
          fill=NA
        )
      )
      hold_results_lys$DOYtrunc = zoo::rollapplyr(
        data=ly_data %>% dplyr::select(dplyr::all_of(time_col)) %>% dplyr::pull(),
        width=stats_combinations[s, "Width"],
        FUN=max,
        partial=TRUE
        )
      hold_results_lys$Stat = stats_combinations[s, "Stat"]
      hold_results_lys$Width = stats_combinations[s, "Width"]
      ly_stats_list[[s]] = hold_results_lys
    }
    ly_all_stats = dplyr::bind_rows(ly_stats_list)
    ly_all_stats$Lake = dplyr::pull(lake_years[i, "Lake"])
    ly_all_stats$Year = dplyr::pull(lake_years[i, "Year"])
    all_lakeyear_all_stats[[i]] = ly_all_stats
  }
  out_stats = dplyr::bind_rows(all_lakeyear_all_stats)
  return(out_stats)
}
