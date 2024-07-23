#' Predictions from a fitted \code{ppgam} object
#'
#' @param object a fitted \code{ppgam} object
#' @param newdata a data frame
#' @param type a character string giving the type of prediction sought; see Details. Defaults to \code{"link"}
#' @param se.fit a logical: should estimated standard errors be returned? Defaults to \code{FALSE}
#' @param ... passed to \code{mgcv::predict()}
#'
#' @details
#'
#' This calls \link[mgcv]{predict.gam} and gives predictions of the intensity function of 
#' the Poisson process on the original scale if \code{type = "response"}, on log 
#' scale if \code{type = "link"} (default), and of the design matrix if 
#' \code{type = "lpmatrix"}.
#'
#' @references 
#'
#' Youngman, B. D., & Economou, T. (2017). Generalised additive point process models for natural 
#' hazard occurrence. Environmetrics, 28(4), e2444.
#' 
#' @seealso \link[mgcv]{predict.gam}
#'
#' @return A data frame or list of predictions
#'
#' @examples
#' 
#' # Times of landfalling US hurricanes
#' data(USlandfall)
#' 
#' # convert dates to years, as a continuous variable
#' year <- as.integer(format(USlandfall$date, "%Y"))
#' day <- as.integer(format(USlandfall$date, "%j"))
#' USlandfall$year <- year + pmin(day / 365, 1)
#' hits <- subset(USlandfall, landfall == 1)
#' 
#' # this creates nodes in the default way
#' m1 <- ppgam( ~ s(year), hits)
#' predict(m1)
#' predict(m1, type = "response")
#' predict(m1, type = "lpmatrix")
#' predict(m1, newdata = data.frame(year = c(2000, 2001)))
#' predict(m1, se.fit = TRUE)
#' 
#' @export
#'
predict.ppgam <- function(object, newdata, type = "link", se.fit = FALSE, ...) {
  mgcv::predict.gam(object, newdata, type, se.fit, ...)
}

#' Extract Model Fitted Values
#'
#' @param object a fitted \code{ppgam} object
#' @param ... not used
#'
#' @examples
#'
#' # Times of landfalling US hurricanes
#' data(USlandfall)
#' 
#' # convert dates to years, as a continuous variable
#' year <- as.integer(format(USlandfall$date, "%Y"))
#' day <- as.integer(format(USlandfall$date, "%j"))
#' USlandfall$year <- year + pmin(day / 365, 1)
#' hits <- subset(USlandfall, landfall == 1)
#' 
#' m1 <- ppgam( ~ s(year), hits)
#' fitted(m1)
#' 
#' @return Fitted values extracted from the object `object'.
#' 
#' @export
#' 
fitted.ppgam <- function(object, ...) {
  predict(object)
}

#' Extract Model Coefficients
#'
#' @param object a fitted \code{ppgam} object
#' @param ... not used
#'
#' @examples
#'
#' # Times of landfalling US hurricanes
#' data(USlandfall)
#' 
#' # convert dates to years, as a continuous variable
#' year <- as.integer(format(USlandfall$date, "%Y"))
#' day <- as.integer(format(USlandfall$date, "%j"))
#' USlandfall$year <- year + pmin(day / 365, 1)
#' hits <- subset(USlandfall, landfall == 1)
#' 
#' m1 <- ppgam( ~ s(year), hits)
#' coef(m1)
#' 
#' @return Model coefficients extracted from the object `object'.
#' 
#' @export
#' 
coef.ppgam <- function(object, ...) {
  object$coefficients
}
#' Plots smooths of a fitted \code{ppgam} object
#'
#' @param x a fitted \code{ppgam} object
#' @param ... arguments to be passed to \code{mgcv::plot.gam}
#'
#' @return Simulations of parameters
#'
#' @seealso \link[mgcv]{plot.gam}
#'
#' @examples
#'
#' # Times of landfalling US hurricanes
#' data(USlandfall)
#' 
#' # convert dates to years, as a continuous variable
#' year <- as.integer(format(USlandfall$date, "%Y"))
#' day <- as.integer(format(USlandfall$date, "%j"))
#' USlandfall$year <- year + pmin(day / 365, 1)
#' hits <- subset(USlandfall, landfall == 1)
#' 
#' # this creates nodes in the default way
#' m1 <- ppgam( ~ s(year), hits)
#' plot(m1)
#'
#' @export
#' 
plot.ppgam <- function(x, ...) {
  mgcv::plot.gam(x, ...)
}

#' Print a fitted \code{ppgam} object
#'
#' @param x a fitted \code{ppgam} object
#' @param ... other arguments passed to \link[mgcv]{print.gam}
#' 
#' @details
#' 
#' Calls \link[mgcv]{print.gam}.
#'
#' @return Prints a details of a fitted ppgam object
#'
#' @seealso \link[mgcv]{print.gam}
#'
#' @examples
#'
#' # Times of landfalling US hurricanes
#' data(USlandfall)
#' 
#' # convert dates to years, as a continuous variable
#' year <- as.integer(format(USlandfall$date, "%Y"))
#' day <- as.integer(format(USlandfall$date, "%j"))
#' USlandfall$year <- year + pmin(day / 365, 1)
#' hits <- subset(USlandfall, landfall == 1)
#' 
#' # this creates nodes in the default way
#' m1 <- ppgam( ~ s(year), hits)
#' print(m1)
#'
#' @export
#' 
print.ppgam <- function(x, ...) {
  mgcv::print.gam(x, ...)
}

#' Simulations from a fitted \code{ppgam} object
#'
#' @param object a fitted \code{ppgam} object
#' @param nsim an integer giving the number of simulations
#' @param seed an integer giving the seed for simulations
#' @param newdata a data frame
#' @param type a character string, as in \code{predict.ppgam}
#' @param ... arguments to be passed to \code{predict.ppgam}
#'
#' @return Simulations of parameters
#'
#' @seealso \link{predict.ppgam}
#'
#' @examples
#'
#' # Times of landfalling US hurricanes
#' data(USlandfall)
#' 
#' # convert dates to years, as a continuous variable
#' year <- as.integer(format(USlandfall$date, "%Y"))
#' day <- as.integer(format(USlandfall$date, "%j"))
#' USlandfall$year <- year + pmin(day / 365, 1)
#' hits <- subset(USlandfall, landfall == 1)
#' 
#' # this creates nodes in the default way
#' m1 <- ppgam( ~ s(year), hits)
#' simulate(m1)
#' simulate(m1, type = "response")
#' simulate(m1, newdata = data.frame(year = c(2000, 2001)))
#'
#' @export
#' 
simulate.ppgam <- function(object, nsim = 1, seed = NULL, newdata,
                           type = "link", ...) {
  if(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1) # initialize the RNG if necessary
  if(is.null(seed)) {
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  } else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  family <- object$family
  V.type <- "Vp"
  B <- .pivchol_rmvn(nsim, object$coefficients, object[[V.type]])
  X <- predict.ppgam(object, newdata, type = "lpmatrix")
  X <- X %*% B
  if (type == "response")
    X <- exp(X)
  return(X)
}
