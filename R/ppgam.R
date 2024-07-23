#' Fit a generalised additive point process model
#'
#' @param formula a formula for a Poisson process log intensity function (compatible with \link[mgcv]{gam})
#' @param data a data frame
#' @param nodes a list or data frame; see `Details'
#' @param weights a scalar, list or vector; see `Details'
#' @param nquad a scalar giving the number of quadrature nodes for each variable
#' @param approx a length 2 character string; see `Details'
#' @param knots spline knots to pass to \link[mgcv]{gam}
#' @param use.data should splines should be constructed from \code{data} (otherwise uses \code{nodes})?
#' @param trace integers controlling what's reported. Defaults to 0
#' @param weight.non.numeric should nodes for non-numeric variables be weighted according how often they occur? Defaults to \code{FALSE}
#' 
#' @details
#'
#' \code{ppgam} fits a Poisson process with intensity function \eqn{\lambda({\bf x})} for covariate
#' \eqn{{\bf x} = (x_1, \ldots, x_d)}. The likelihood for this model with events occurring at \eqn{{\bf x}_i}, 
#' for \eqn{i=1, \ldots, n}, is approximated by quadrature with 
#'
#' \deqn{\exp\bigg[-\sum_{j=1}^m w_j \lambda({\bf x}_j^*)\bigg] \prod_{i=1}^n \lambda({\bf x}_i)}
#'
#' where \eqn{{\bf x}_j^*} and \eqn{w_j} are quadrature nodes and weights, for \eqn{j=1, \ldots, m},
#' defined with \code{nodes} and \code{weights}.
#'
#' \code{formula} gives the formula for the log intensity function of a Poisson process.
#' It is passed to \link[mgcv]{gam}. If \code{formula} has no response, i.e. \code{ ~ s(...)},
#' then \code{data} is assumed to give the times at which events occur. Then \code{nodes}
#' is used to control integration of the intensity function. If \code{formula} has a response,
#' e.g. \code{y ~ s(...)}, then \code{y} is assumed binary, comprising only zeros and ones.
#' Then \code{data} is assumed to give the state space of the Poisson process,
#' (e.g. daily time steps if occurrences of events are measured in days)
#' and ones in \code{y} identify when events occur. Note that
#' if \code{formula} has no response, \code{data} will have \eqn{n} rows, and \eqn{m} rows otherwise.
#'
#' \code{nodes} is used to supply nodes for integrating the Poisson
#' process intensity function by quadrature. It is supplied as a list or data 
#' frame. 
#' 
#' If \code{nodes} is a list, its names must correspond to variables on 
#' the r.h.s. of formula. Elements of the list, \code{x}, say, can be a vector 
#' or 2-column matrix, where \code{length(x) > 1} or \code{nrow(x) > 1}. If a 
#' matrix, its first and second columns are taken as integration nodes and
#' weights, respectively. If a vector of length 2, it is assumed to give the 
#' range of the \code{nquad} midpoints used as integration nodes. If a longer vector,
#' it is assumed to be the integration nodes, and \code{nquad} is ignored.
#' 
#' If \code{nodes} is a data frame, it is assumed to give the integration nodes.
#' 
#' \code{nquad} specifies the number of integration nodes per variable, unless
#' nodes are specified in \code{nodes}. If a single integer and 
#' \code{is.null(names(nquad))} it is used for all variables. Otherwise, 
#' names are matched to variables. An error is returned if any variables do not 
#' have values specified.
#'
#' \code{weights} controls the quadrature weights. If \code{nodes} is a list,
#' a scalar multiplies any weights calculated alongside nodes, i.e. node separations.
#' If \code{nodes} is a data frame, weights can be a scalar that is repeated
#' \code{nrow(nodes)}, or a vector of length \code{nrow(nodes)} that gives the weights
#' for each row of \code{nodes}.
#'  
#' \code{approx} controls quadrature details. Its first term controls the 
#' integration method, which uses either midpoint (\code{"midpoint"}, default), 
#' Simpson's (\code{"Simpson"}) or Gauss-Legendre (\code{"Gauss"}) rules. The second
#' term of \code{approx} controls the integration range, which is either the
#' range of the variable (\code{"exact"}), or by calling \code{pretty()} (\code{"pretty"}).
#'
#' \code{trace} controls what is reported. Details of convergence are printed 
#' with \code{trace = 1}, of nodes with \code{trace = 2}, and \code{trace = 3}
#' prints both.
#' 
#' \code{weight.non.numeric} applied to any non-numeric variables, and gives non-equal 
#' quadrature weights to different nodes if \code{TRUE}. So nodes get weights according
#' to their frequency of occurrence in \code{data}. If \code{nquad} is invoked, only a
#' subset of the unique values of the non-numeric variable are used, which are the \code{nquad}
#' with largest weights.
#'
#' @references
#' 
#' Wood, S. N., Pya, N., & Safken, B. (2016). Smoothing parameter and model selection for general 
#' smooth models. Journal of the American Statistical Association, 111(516), 1548-1563.
#'
#' Youngman, B. D., & Economou, T. (2017). Generalised additive point process models for natural 
#' hazard occurrence. Environmetrics, 28(4), e2444.
#' 
#' @return 
#' 
#' An object of class \code{gam}, as returned by \code{mgcv::gam}, with parameters,
#' covariance matrices and a few other things  swapped
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
#' 
#' # some examples of providing nodes
#' nodes.year <- list(year=pretty(USlandfall$year, 20))
#' # as 2 is in trace, nodes and weights are printed
#' m2 <- ppgam( ~ s(year), hits, nodes = nodes.year, trace = 2)
#' 
#' # alternatively, we might just want to specify how many nodes to use
#' m3 <- ppgam( ~ s(year), hits, nquad = 30)
#'
#' \donttest{
#'
#' data(windstorm)
#' m4 <- ppgam(~ s(lon, lat, k=20), windstorm)
#'
#' ## Storm peak locations, given the North Atlantic Oscillation (NAO) index
#' # NAO values from https://crudata.uea.ac.uk/cru/data/nao/nao.dat
#' # NAO midpoints and weights based on `hist'
#'
#' NAO.mids <- c(-2.75, -2.25, -1.75, -1.25, -0.75, -0.25, 0.25, 0.75, 1.25, 1.75, 2.25)
#' NAO.wts <- c(0.002, 0.014, 0.057, 0.145, 0.302, 0.427, 0.463, 0.364, 0.171, 0.047, 0.007)
#'
#' m5 <- ppgam(~ te(lat, lon, NAO, d = 2:1, k = c(40, 8), bs = c("ts", "cr")), windstorm, 
#'   nodes = list(NAO = cbind(NAO.mids, NAO.wts)))
#'
#' }
#'
#' @export
ppgam <- function(formula, data, nodes = NULL, weights = 1, nquad, 
approx = c("midpoint", "exact"), knots = NULL, use.data = TRUE, trace = 0,
weight.non.numeric = FALSE) {

# convert and objects of class "Date" to integer
data.class <- sapply(data, class)
if (any(data.class == "Date")) for (i in which(data.class == "Date")) data[,i] <- as.integer(data[,i])

# identify whether response is present
gam.info <- mgcv::interpret.gam(formula)
rsp <- gam.info$response
with.response <- !is.null(rsp)
if (with.response) {
# check whether response okay
  if(!all(data[,rsp] %in% 0:1))
    stop(paste("Some data$", rsp, " not in {0, 1}", sep=""))
  nodes <- data
  data <- subset(data, data[,rsp] == 1)
} else {
# add in fake response for mgcv::gam
  formula0 <- formula
  rsp <- "fake_response"
  data$fake_response <- sample(0:1, nrow(data), replace=TRUE)
  formula <- as.character(formula)
  if (length(formula) > 2) {
    formula <- paste(rsp, "~", paste(formula[-1], collpase=" + "))
  } else {
    formula <- paste(rsp, "~", formula[-1])
  }
  formula <- as.formula(formula)
  gam.info <- mgcv::interpret.gam(formula)
}

# identify covariates
vars <- gam.info$fake.names
nvars <- length(vars)

nquad <- .make_nquad(nquad, vars, data[,vars], trace)

## set up quadrature nodes and weights
# identify node type
# node.type == 0 is data.frame, so nodes assumed
# node.type == 1 is list, so nodes produced
if (is.data.frame(nodes)) {
  node.type <- 0
  # nothing else needs to be done to nodes
} else {
  node.type <- 1
}

w0 <- weights

if (!is.data.frame(nodes)) {
  approx <- .check.approx(approx)
  nodes <- lapply(vars, function(x) nodes[[x]])
  nodes <- lapply(seq_len(nvars), function(i) .make_nodes(nodes[[i]], data[,vars[i]], nquad[i], approx, trace, weight.non.numeric))
  names(nodes) <- vars
  if (trace %in% 2:3) {
    cat("\nQuadrature nodes and weights:\n")
    lapply(seq_along(nodes), function(i) .print.nodes(nodes[[i]], names(nodes)[i]))
  }
  weights <- lapply(nodes, function(x) x[,2])
  weights <- Reduce("*", expand.grid(weights))
  if (length(w0) > 1)
    stop("'weights' of length > 1 not compatible with `nodes' supplied as list. See Details for `ppgam'.")
  weights <- w0 * weights
  nodes <- lapply(nodes, function(x) x[,1])
  nodes <- expand.grid(nodes)
  n.nodes <- nrow(nodes)
} else {
  n.nodes <- nrow(nodes)
  if (length(w0) > 1 & length(w0) < nrow(nodes)) {
    stop("length `weights' must be 1 or nrow(nodes)")
  } else {
    if (length(w0) == 1) {
      weights <- rep(w0, n.nodes)
    }
  } 
}

if (use.data) {
  G <- mgcv::gam(formula, data=data, knots=knots, fit=FALSE, method="REML", family=poisson)
} else {
  G <- mgcv::gam(formula, data=nodes, knots=knots, fit=FALSE, method="REML", family=poisson)
}

class(G) <- c("gam", "glm", "lm")
G$nb <- ncol(G$X)
G$coefficients <- numeric(G$nb)
G$wts <- weights

# now sort design matrices
if (use.data) {
  G$XT <- predict(G, newdata=nodes, type="lpmatrix")
} else {
  G$XT <- G$X
  G$X <- predict(G, newdata=data, type="lpmatrix")
}

G$control <- .control.ppgam()

beta <- .give_beta0(G)
rho0 <- numeric(length(G$sp))
G$Sd <- evgam:::.joinSmooth(list(G))
G$S <- evgam:::.makeS(G$Sd, rho0)
G$null.deviance <- NA#.f0(beta, G)
attr(rho0, "beta") <- beta

fit.reml <- evgam:::.BFGS(rho0, .reml0, .reml1, dat=G, control=G$control$outer, trace=trace %in% c(1, 3))

G$coefficients <- attr(fit.reml$objective, "beta")
H <- attr(fit.reml$objective, "Hessian")
G$Vp <- MASS::ginv(H)
G$S <- evgam:::.makeS(G$Sd, exp(fit.reml$par))
G$deviance <- .f0(G$coefficients, G)
H0 <- H - G$S
G$edf <- colSums(G$Vp %*% H0)
G$logLik <-  -.f(G$coefficients, G)
G$scale.estimated <- 0
G$df.residual <- ncol(G$X) - sum(G$edf)
G$aic <- 2 * (G$df.residual - G$logLik)
G$nobs <- nrow(data)
old <- c("mf", "pP", "cl")
new <- c("model", "paraPen", "call")
is.in <- !is.na(match(old, names(G)))
if (any(is.in)) names(G)[match(old[is.in], names(G))] <- new[is.in]
G$na.action <- attr(G$mf, "na.action")
G$sig2 <- G$scale.estimated <- 1
G$method <- "REML"
if (!with.response) G$formula <- formula0
G$family$no.r.sq <- TRUE
G$gcv.ubre <- as.vector(fit.reml$objective)
n.samp <- max(1e3, 2 * G$nb)
if (n.samp < nrow(G$X)) {
  G$R <- G$X[sample(nrow(G$X), n.samp, replace=FALSE),]
} else {
  G$R <- G$X
}
G$family$family <- 'Poisson process'

class(G) <- c('ppgam', class(G))
return(G)

}
