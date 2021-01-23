#' Daily time series data of lake bloom variables
#'
#' A dataset containing daily data of bloom variables
#'
#' @format A data frame with 1957 rows and 7 variables:
#' \describe{
#'   \item{Lake}{character indicating lake; "L" is Paul Lake, "R" is Peter Lake, "T" is Tuesday Lake}
#'   \item{Year}{year}
#'   \item{DOY}{"Day Of Year" or Julian day of year, integer from 1 to 366}
#'   \item{BGA_HYLB}{blue green algae fluorescence measured by Hydrolab DS5X sonde, in units cells/mL}
#'   \item{Manual_Chl}{chlorophyll concentration measured by manual extraction and lab fluorometry, units ug/L}
#'   \item{DO_Sat}{dissolved oxygen, as percent of saturation in equilibrium with atmosphere}
#'   \item{pH}{pH measured by sonde}
#' }
#' @source \url{Cascade research group; TODO add URL}
"ts_data"

#' Spatial data of lake bloom variables
#'
#' A dataset containing spatial data of bloom variables collected using FLAMe spatial sampling system (TODO: add FLAMe url). Measurements are 'tau' corrected, corrected for the hydrology residence time of the system and each sensor's response time.
#'
#' @format A data frame with 176291 rows and 10 variables:
#' \describe{
#'   \item{Lake}{character indicating lake; "L" is Paul Lake, "R" is Peter Lake, "T" is Tuesday Lake}
#'   \item{Year}{year}
#'   \item{DOY}{"Day Of Year" or Julian day of year that sample was collected on, integer from 1 to 366}
#'   \item{Date}{Date that sample was collected on, Date (YYYY-MM-DD) format}
#'   \item{date_time}{precise data-time that sample was collected at, POSIXct format in UTC time zone}
#'   \item{latitude}{latitude of sample}
#'   \item{longitude}{longitude of sample}
#'   \item{BGApc_ugL_tau}{blue green algae fluorescence measured by YSI EXO3 sonde, in units ug/L of phycocyanin}
#'   \item{ODO_percent_tau}{dissolved oxygen, as percent of saturation in equilibrium with atmosphere}
#'   \item{pH_tau}{pH measured by sonde}
#' }
#' @source \url{Cascade research group; TODO add URL}
"flame_data"
