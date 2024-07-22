## functions for quadrature

.make_nquad <- function(x, vars, data, trace) {
# x the `nquad' argument supplied to `ppgam'
# v the names of variables in `formula' supplied to `ppgam' (except response)
nvars <- length(vars)
if (missing(x)) {
  nquad <- rep(15, nvars)
  names(nquad) <- vars
  fc <- sapply(data, inherits, c("factor", "character"))
  if (any(fc)) {
    nquad[fc] <- sapply(data[, fc, drop=FALSE], function(x) length(unique(x)))
  }
} else {
  nquad <- rep(NA, nvars)
  names(nquad) <- vars
  if (is.null(names(x))) {
    if (length(x) == 1) {
      nquad[] <- x
    } else {
      stop("Can't provide NULL-named `nquad' of length > 1. See Details for `ppgam'.")
    }
  } else {
    nquad[names(x)] <- x
  }
}
set2default <- is.na(nquad)
if (any(set2default)) {
  nquad[set2default] <- 15
  if (trace %in% 2:3)
    message(paste("'nquad' set to 15 for: ", paste(vars[set2default], collapse=", "), ".", sep=""))
}
nquad
}

.mids <- function(x, n) {
# function to gives midpoints for vector
# put into n bins
rx <- range(x)
h <- diff(rx) / (n + 1)
rx <- rx + .5 * c(1, -1) * h
cbind(seq(rx[1], rx[2], l=n), h)
}

.simpson <- function(x, n) {
# function to gives midpoints for vector
# put into n bins
if (n %% 2 == 0) n <- n + 1
rx <- range(x)
dx <- diff(rx) 
h <- dx / (n - 1)
x <- seq(rx[1], rx[2], by=h)
wts <- c(1, rep(c(4, 2), (n - 3) / 2), 4, 1)
wts <- dx * wts / sum(wts)
cbind(x, wts)
}

.gaussian <- function(x, n) {
# modified from pracma::gaussLegendre
# original scale stuff
rx <- range(x)
h <- diff(rx)
# unit stuff
i <- seq_len(n - 1)
d <- i / sqrt(4*i^2 - 1)
D <- matrix(0, n, n)
off <- cbind(i, i + 1)
D[off] <- D[off[,2:1]] <- d
E <- eigen(D, symmetric=TRUE)
L <- E$values
V <- E$vectors
inds <- order(L)
x <- L[inds]
V <- t(V[, inds])
wts <- 2 * V[, 1]^2
x <- h * x + sum(rx)
wts <- 2 * h * wts / sum(wts)
.5 * cbind(x, wts)
}

.check.approx <- function(x) {
x <- tolower(x)
if (!all(x %in% c("midpoint", "simpson", "gauss", "exact", "pretty")))
  stop("Unrecognized entries in `approx'.")
if (length(x) == 1) {
  if (x %in% c("exact", "pretty")) {
    x <- c("midpoint", x)
  } else {
    x <- c(x, "exact")
  }
}
x
}

.non.numeric <- function(x, n, wnn) {
ux <- sort(unique(x))
nx <- length(ux)
if (wnn) {
  wts <- as.vector(table(x))
} else {
  wts <- rep(1/n, n)
}
if (nx > n) {
  if (wnn) {
    samp <- sort(order(wts, decreasing=TRUE)[seq_len(n)])
    wts <- wts[samp]
    wts <- wts / sum(wts)
  } else {
    samp <- round(nx * (1:n - .5) / n)
  }
  ux <- ux[samp]
}
data.frame(ux, wts)
}

.approx <- function(x, n, type, wnn) {
if (inherits(x, c("factor", "character"))) {
  out <- .non.numeric(x, n, wnn)
} else {
  if (type == "midpoint") {
    out <- .mids(x, n)
  } else {
    if (type == "simpson") {
      out <- .simpson(x, n)
    } else {
      out <- .gaussian(x, n)
    }
  }
}
out
}

.make_nodes <- function(x, y, n, approx, trace, wnn) {
# function to give quadrature nodes and weights
# x can be ...
# -- NULL, then y and n are used
# -- an integer, then it overrides n
# -- a 2-vector, then it's the range for the nodes
# -- a n-vector, n > 2, then it's the nodes themselves
# -- a n x 2 matrix, then column 1 is as above, and column 2 is the weights
# n controls the number of nodes
# y is used to give x if x is NULL
# approx[1] is the approximation from c("midpoint", "Simpson", "Gaussian") to use
# approx[2] is from c("exact", "hist", "pretty") is the function to create node range
# return a 2-column matrix of nodes and weights
if (!is.null(x)) {
  out <- as.data.frame(x)
  if (ncol(out) == 1) {
    if (trace %in% 2:3) 
      cat("Nodes taken as breaks from which midpoints derived.\n")
    if (inherits(out[,1], c("factor", "character"))) {
      if (ncol(out) == 1) {
        browser()
        if (!all(out[,1] %in% x))
          stop("Some nodes gives for non-numeric variable not in `data'.")
        out[,2] <- 1 / length(out[,1])
      }
    } else {
      h <- diff(out[,1])
      h <- h / sum(h)
      out <- data.frame(x[-1] - .5 * h, h)
    }
  }
  if (trace %in% 2:3) {
    if (sum(out[,2]) != 1)
      cat("Note: weights not scaled to sum to one.\n")
  }
} else {
  if (approx[2] == "pretty") {
    x <- pretty(y, n)
    n <- length(x) - 1
    out <- .approx(x, n, approx[1], wnn)  
  } else {
    out <- .approx(y, n, approx[1], wnn) 
  }
  out[,2] <- out[,2] / sum(out[,2])
}
return(out)
}

.print.nodes <- function(x, nm) {
cat("\n")
cat(paste("**", nm, "**\n"))
cat("- Nodes:\n")
print(x[,1])
cat("- Weights:\n")
print(x[,2])
}

## function for initial basis function coefficients

.give_beta0 <- function(G) {
p <- ncol(G$X)
here <- which.min(colSums((G$X - 1)^2))
G <- list(wts=G$wts, X=matrix(1, nrow(G$X), 1), XT=matrix(1, nrow(G$XT), 1), control=G$control)
init <- 10^seq(-10, 10)
f.test <- sapply(init, .f0, dat=G)
if (any(f.test != 1e20)) {
  init <- init[which.min(f.test)]
} else {
  stop("Can't find sensible starting values in [1e-10, 1e10]")
}
G$S <- matrix(0, 1, 1)
init <- evgam:::.newton_step(init, .f, .search, dat=G, control=G$control$inner)$par
replace(numeric(p), here, init)
}

## functions for point process likelihood

.f0 <- function(pars, dat) {
# negative log-likelihood 
# for point process model
out <- -sum(tcrossprod(pars, dat$X))
out <- out + sum(dat$wts * exp(tcrossprod(pars, dat$XT)[1,]))
if (!is.finite(out)) out <- 1e20
out <- min(out, 1e20)
return(out)
}

.f <- function(pars, dat, newton=FALSE) {
# negative penalised log-likelihood 
# for point process model
out <- .f0(pars, dat)
out <- out + .5 * sum(pars * crossprod(pars, dat$S)[1,])
out
}

.gH0 <- function(pars, dat) {
# gradient and Hessian of negative log-likelihood 
# for point process model
g <- -colSums(dat$X)
wLambda <- dat$wts * exp(tcrossprod(pars, dat$XT)[1,])
XTb <- dat$XT * wLambda
g <- g + colSums(XTb)
g[!is.finite(g)] <- 1e20
H <- crossprod(dat$XT, XTb)
H[!is.finite(H)] <- 1e20
list(g=g, H=H)
}

.gH <- function(pars, dat) {
# gradient and Hessian of penalised negative log-likelihood 
# for point process model
gH <- .gH0(pars, dat)
bS <- crossprod(pars, dat$S)[1,]
gH$g <- gH$g + bS
gH$H <- gH$H + dat$S
gH
}

.search <- function(pars, dat, kept, newton=TRUE) {
# Newton search direction
gH <- .gH(pars, dat)
if (newton) {
  out <- as.vector(solve(gH[[2]], gH[[1]]))
} else {
  out <- gH[[1]]
}
attr(out, "gradient") <- gH[[1]]
attr(out, "Hessian") <- gH[[2]]
out
}

## REML functions
# slightly adapted version of that from evgam

.reml0 <- function(pars, dat, beta=NULL, skipfit=FALSE) {
if (is.null(beta)) beta <- attr(pars, "beta")
sp <- exp(pars)
dat$S <- evgam:::.makeS(dat$Sd, sp)
if (!skipfit) {
fitbeta <- evgam:::.newton_step(beta, .f, .search, dat=dat, control=dat$control$inner)
if (any(abs(fitbeta$gradient) > 1)) {
it0 <- dat$control$inner$itlim
dat$control$inner$itlim <- 10
fitbeta <- evgam:::.newton_step(fitbeta$par, .f, .search, dat=dat, control=dat$control$inner, newton=FALSE, alpha0=.05)
dat$control$inner$itlim <- it0
fitbeta <- evgam:::.newton_step(fitbeta$par, .f, .search, dat=dat, control=dat$control$inner)
}
if (inherits(fitbeta, "try-error")) return(1e20)
} else {
fitbeta <- list(objective=.f(beta, dat))
fitbeta$convergence <- 0
fitbeta$gH <- .gH(beta, dat)
fitbeta$Hessian <- fitbeta$gH[[2]]
fitbeta$par <- beta
}
logdetSdata <- evgam:::.logdetS(dat$Sd, pars)
halflogdetHdata <- try(list(d0=sum(log(diag(chol(fitbeta$Hessian))))), silent=TRUE)
if (inherits(halflogdetHdata, "try-error")) return(1e20)
out <- fitbeta$objective + as.numeric(fitbeta$convergence != 0) * 1e20
out <- out + halflogdetHdata$d0 - .5 * logdetSdata$d0
out <- as.vector(out)
if (!is.finite(out)) return(1e20)
attr(out, "beta") <- fitbeta$par
attr(out, "gradient") <- fitbeta$gradient
attr(out, "Hessian") <- fitbeta$Hessian
return(out)
}

.reml1 <- function(pars, dat, H=NULL, beta=NULL) {
if (is.null(beta)) {
    beta <- attr(pars, "beta")
} else {
    attr(pars, "beta") <- beta
}
sp <- exp(pars)
spSl <- lapply(seq_along(sp), function(i) sp[i] * attr(dat$Sd, "Sl")[[i]])
S <- dat$S <- Reduce("+", spSl)
H0 <- .gH0(beta, dat)[[2]]
H <- H0 + S
iH <- MASS::ginv(H)
spSlb <- sapply(spSl, function(x) x %*% beta)
db <- crossprod(iH, spSlb)
thirdpp <- dat$wts * exp(tcrossprod(beta, dat$XT))[1,]
thirdpp <- (dat$XT * as.vector(thirdpp)) %*% db
eV <- eigen(H, symmetric=TRUE)
XTV <- t(solve(eV$vectors, t(dat$XT)))
dH <- sapply(seq_along(sp), function(i) sum(colSums(XTV * XTV * thirdpp[,i]) / eV$values))
dH <- sapply(spSl, function(x) sum(iH * x)) - dH
dH <- list(d1=dH)
dS <- evgam:::.logdetS(dat$Sd, pars, deriv=1)
d1 <- .5 * sapply(spSl, function(x) crossprod(beta, x %*% beta))
d1 <- d1 - .5 * dS$d1
d1 <- d1 + .5 * dH$d1
d1
}

## other functions

.control.ppgam <- function(i, o) {
ctrl <- list()
ctrl$inner <- list(steptol=1e-12, itlim=1e2, fntol=1e-8, gradtol=1e-4, stepmax=1e2)
ctrl$outer <- list(steptol=1e-12, itlim=1e2, fntol=1e-8, gradtol=2e-2, stepmax=3)
ctrl
}

.pivchol_rmvn <- function(n, mu, Sig) {
  R <- suppressWarnings(chol(Sig, pivot = TRUE))
  piv <- order(attr(R, "pivot"))  ## reverse pivoting index
  r <- attr(R, "rank")  ## numerical rank
  V <- R[1:r, piv]
  Y <- crossprod(V, matrix(rnorm(n * r), r))
  Y + as.vector(mu)
}

