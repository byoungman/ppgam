#' Locations of windstorm peaks and tracks over the North Atlantic
#'
#' A dataset in windstorm peaks between 1st January 1979 and 31st December 2014
#' occuring in [-50, 33] longitude and [36, 77] latitude.
#'
#' @format A data frame with 3133 rows and 4 variables
#' 
#' The variables are as follows:
#'
#' \describe{
#'   \item{date}{date of peak, as class "Date"}
#'   \item{lon}{longitude, in degrees}
#'   \item{lat}{latitude, in degrees}
#'   \item{NAO}{North Atlantic Oscillation index}
#' }
#'
#' @docType data
#' @keywords datasets
#' @name windstorm
#' @usage data(windstorm)
#'
#' @references
#'
#' Youngman, B. D., & Economou, T. (2017). Generalised additive point process models for natural 
#' hazard occurrence. Environmetrics, 28(4), e2444.
#'
#' @examples
#' 
#' data(windstorm)
#' plot(windstorm[,c("lon", "lat")])
#'
NULL

#' Times of landfalling US hurricanes
#'
#' A data frame:
#'
#' @format A data frame with 61129 rows and 2 variables
#' 
#' The variables are as follows:
#'
#' \describe{
#'   \item{date}{date of landfall, as class "Date"}
#'   \item{landfall}{an integer: did a hurricane make landfall on this day?}
#' }
#'
#' @docType data
#' @keywords datasets
#' @name USlandfall
#' @usage data(USlandfall)
#'
#' @references
#'
#' https://www.nhc.noaa.gov/data/
#'
#' @examples
#' 
#' data(USlandfall)
#' plot(USlandfall, type="h")
#'
NULL
