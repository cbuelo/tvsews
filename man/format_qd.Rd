% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quickest_detection.R
\name{format_qd}
\alias{format_qd}
\title{Format quickest detection results for plotting}
\usage{
format_qd(
  qd_stats,
  var_cols = c("Manual_Chl", "BGA_HYLB", "DO_Sat", "pH"),
  stats = c("SD", "Ar1"),
  bloom_fert_df = bloom_fert_dates,
  exp_lakes = c("T", "R"),
  just_alarm_DOYs = FALSE
)
}
\arguments{
\item{qd_stats}{data frame, output from \code{\link{run_qd}}}

\item{var_cols}{character vector, variable names that are columns of *data* containing time series to calculate rolling window stats on}

\item{stats}{character vector, statistics to format alarms for, default = c("SD", "Ar1)}

\item{bloom_fert_df}{data frame of fertilization start and end dates, see \code{\link{bloom_fert_dates}} for default and formatting}

\item{exp_lakes}{character vector, values in *Lake* column of qd_stats that correspond to experimental lakes}

\item{just_alarm_DOYs}{TRUE or FALSE (default), should results be filtered to just DOYs that have an alarm (in any variable)?}
}
\value{
data frame with quickest detection alarm info
}
\description{
Format quickest detection results for plotting
}
\examples{
rw_stats <- calc_rolling_stats(ts_data, var_cols = c("DO_Sat"))
qd_stats <- run_qd(rw_stats, var_cols = c("DO_Sat"))
qd_alarms <- format_qd(qd_stats, var_cols = c("DO_Sat"))
}
