
<!-- README.md is generated from README.Rmd. Please edit that file -->
# tvsews

<!-- badges: start -->
<!-- badges: end -->
tvsews is to provides code and data for a manuscript (in preparation) associated with a whole-lake fertilization experiment that compared the performance of early warning statistics (EWS) in temporal and spatial data.

## Installation

You can install the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("cbuelo/tvsews")
```

## Example

Results from most recent experiment:

### Set up and visualize data

``` r
library(tvsews)
library(dplyr)
ex_data <- ts_data %>% 
  filter(Year >= 2018) %>% 
  select(Lake, Year, DOY, Manual_Chl)

ex_bloom_fert_dates <- bloom_fert_dates %>% 
  filter(Year >= 2018 & Lake == "R")

plot_fig1(time_series = ex_data, 
          bloom_fert_df = ex_bloom_fert_dates,
          variables_rename_vec = c(`Chl-a (ug/L)` = "Manual_Chl")
          )
```

<img src="man/figures/README-example-1.png" width="100%" />

### Calculate rolling window statistics and quickest detection alarms

``` r
ex_rw_stats <- calc_rolling_stats(data = ex_data,
                               var_cols = "Manual_Chl",
                               widths = c(21)
                               ) # warnings are missing data, okay
ex_qd_all <- run_qd(rolling_window_stats = ex_rw_stats,
                    var_cols = "Manual_Chl",
                    widths = c(21),
                    stats_to_qd = c("SD", "Ar1"),
                    exp_lakes = c("R"),
                    ref_lake = "L",
                    ar1_alarm_rho = 0.95
                    )
#> [1] "Current variable QD-ing: Manual_Chl"
ex_qd_formatted <- format_qd(qd_stats = ex_qd_all,
                             var_cols = c("Manual_Chl")
                             )
```

``` r
# plot the rolling window statistics and alarms
time_series_plot = plot_fig3(rolling_window_stats = ex_rw_stats,
          qd_alarms = ex_qd_formatted,
          bloom_fert_df = ex_bloom_fert_dates,
          rename_vec = c(`Chl-a` = "Manual_Chl")
          )
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/master/examples>.
