#' Simulated Biweekly Carabid Pitfall Trapping Dataset
#'
#' A lightweight synthetic dataset of biweekly pitfall trap captures across
#' five sites along a latitudinal gradient (`51.0` to `57.0` degrees North)
#' over three sampling years (`2018`--`2020`), exhibiting bimodal spring and
#' autumn activity modulated by latitude and cumulative year-to-date thermal
#' anomalies (`GDD5_YTD`). Included for quick demonstrations and automated testing.
#'
#' @format A `data.frame` with 390 rows and 7 variables:
#' \describe{
#'   \item{Count}{Integer trap capture count over a 14-day interval.}
#'   \item{DOY}{Integer midpoint Day of Year (`6` to `359`).}
#'   \item{Latitude}{Numeric site latitude in decimal degrees (`51.0` to `57.0`).}
#'   \item{GDD5_YTD}{Numeric standardized year-to-date growing degree-day anomaly.}
#'   \item{LogEffort}{Numeric natural logarithm of trapping effort in trap-days (`log(14)`).}
#'   \item{SubSite}{Factor identifying the 5 sampling sites (`S01`--`S05`).}
#'   \item{fYear}{Factor identifying the 3 sampling years (`2018`--`2020`).}
#' }
#' @usage data(sim_carabids)
#' @examples
#' data(sim_carabids)
#' head(sim_carabids)
"sim_carabids"
