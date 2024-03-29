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

#' Dates that experimental fertilization started and ended on, and that blooms started
#'
#' A dataset containing dates (Year and Day Of Year = DOY) that algal blooms started (crossed lake-specific thresholds) and that fertilizations started and ended.
#'
#' @format A data frame with 9 rows and 5 variables:
#' \describe{
#'   \item{Lake}{character indicating lake; "L" is Paul Lake, "R" is Peter Lake, "T" is Tuesday Lake}
#'   \item{Year}{year}
#'   \item{fertStartDOY}{"Day Of Year" or Julian day of year that nutrient additions to a lake started on}
#'   \item{fertEndDOY}{"Day Of Year" or Julian day of year that nutrient additions to a lake ended on}
#'   \item{bloomStartDOY}{"Day Of Year" or Julian day of year that bloom threshold was first crossed}
#' }
#' @source \url{Cascade research group; TODO add URL}
"bloom_fert_dates"


#' Boundaries for study lakes
#'
#' An sf data.frame with boundaries for study lakes, Peter and Paul lakes.
#'
#' @format An sf data frame with 2 rows and 15 variables:
#' \describe{
#'   \item{AREA}{numeric, lake area in square meters}
#'   \item{PERIMETER}{numeric, lake perimenter in meters meters}
#'   \item{LAKE_TYPE}{TODO: document}
#'   \item{NAME}{TODO: document}
#'   \item{UNIQUE_ID}{TODO: document}
#'   \item{COUNTY}{TODO: document}
#'   \item{NOTE24}{TODO: document}
#'   \item{NEW_KEY}{TODO: document}
#'   \item{HECTARES}{TODO: document}
#'   \item{ACRES_GIS}{TODO: document}
#'   \item{FMU}{TODO: document}
#'   \item{Lake}{TODO: document}
#'   \item{FID_1}{TODO: document}
#'   \item{Lake_Name}{TODO: document}
#'   \item{geometry}{TODO: document}
#'
#' }
#' @source \url{TODO add where this came from}
"peter_paul_bounds"
