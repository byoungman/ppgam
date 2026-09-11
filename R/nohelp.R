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

.nodes2bins <- function(x) {
  vec <- as.vector(t(x[, 1] + outer(x[, 2], .5 * c(-1, 1))))
  colMeans(matrix(rep(vec, c(2, rep(1, length(vec) - 2), 2)), 2))
}

.print.nodes <- function(x, nm) {
  cat("\n")
  cat(paste("**", nm, "**\n"))
  cat("- Nodes:\n")
  print(x[,1])
  cat("- Weights:\n")
  print(x[,2])
}

.nelder_mead_list <- function(init, f, xtol = 1e-3, ftol = 1e-4, max_iter = 500, trace = 1, ...) {
  n <- length(init)  # Dimension of the problem
  alpha <- 1         # Reflection coefficient
  gamma <- 2         # Expansion coefficient
  rho <- 0.5         # Contraction coefficient
  sigma <- 0.5       # Shrink coefficient
  
  # Initialize the simplex
  simplex <- list()
  simplex[[1]] <- init
  f_values_list <- list()
  f_values_list[[1]] <- f(init, ...)
  
  if (trace > 0) {
    cat(paste('Iteration:', 0))
    cat('\n')
    cat(paste('Value:', signif(f_values_list[[1]], 8)))
    cat('\n')
    cat(paste('Initial values: (', paste0(signif(init, 4), collapse = ', '), ')', sep = ''))
    cat('\n')
    cat(paste('Inner max |grad|:', signif(max(abs(attr(f_values_list[[1]], 'gradient'))), 4)))
    cat('\n')
    cat(paste('Inner iterations:', attr(f_values_list[[1]], 'iterations')))
    cat('\n')
    cat('\n')
  }
  
  for (i in 1:n) {
    x <- init
    x[i] <- x[i] + 0.05 * (abs(x[i]) + 1)  # Small perturbation
    simplex[[i + 1]] <- x
  }
  
  # Evaluate function at simplex points
  f_values_list <- lapply(simplex, f, ...)
  f_values <- sapply(f_values_list, as.vector)
  
  iter <- 0
  while (iter < max_iter) {
    iter <- iter + 1
    
    # Order simplex points by function values
    order_idx <- order(f_values)
    simplex <- simplex[order_idx]
    f_values <- f_values[order_idx]
    f_values_list <- f_values_list[order_idx]
    b0 <- attr(f_values_list[[1]], 'beta')
    for (i in seq_along(simplex)) attr(simplex[[i]], 'beta') <- b0
    
    # Centroid of all points except the worst
    # centroid <- colMeans(simplex[1:n, ])
    centroid <- rowMeans(do.call(cbind, simplex[1:n]))
    
    if (trace > 0) {
      cat(paste('Iteration:', iter))
      cat('\n')
      cat(paste('Value:', signif(f_values[1], 8)))
      cat('\n')
      cat(paste('Inner max |grad|:', signif(max(abs(attr(f_values_list[[1]], 'gradient'))), 4)))
      cat('\n')
      cat(paste('Inner iterations:', attr(f_values_list[[1]], 'iterations')))
      cat('\n')
      cat(paste('Centroid: (', paste0(signif(simplex[[1]], 4), collapse = ', '), ')', sep = ''))
      cat('\n')
      cat('\n')
    }
    
    # Reflection
    x_reflect <- centroid + alpha * (centroid - simplex[[n + 1]])
    f_reflect <- f(x_reflect, ...)
    
    if (f_reflect < f_values[1]) {
      # Expansion
      x_expand <- centroid + gamma * (x_reflect - centroid)
      f_expand <- f(x_expand, ...)
      if (f_expand < f_reflect) {
        simplex[[n + 1]] <- x_expand
        f_values[n + 1] <- f_expand
        f_values_list[[n + 1]] <- f_expand
      } else {
        simplex[[n + 1]] <- x_reflect
        f_values[n + 1] <- f_reflect
        f_values_list[[n + 1]] <- f_reflect
      }
    } else if (f_reflect < f_values[n]) {
      # Accept reflection
      simplex[[n + 1]] <- x_reflect
      f_values[n + 1] <- f_reflect
      f_values_list[[n + 1]] <- f_reflect
    } else {
      # Contraction
      if (f_reflect < f_values[n + 1]) {
        # Outside contraction
        x_contract <- centroid + rho * (x_reflect - centroid)
      } else {
        # Inside contraction
        x_contract <- centroid + rho * (simplex[[n + 1]] - centroid)
      }
      f_contract <- f(x_contract, ...)
      
      if (f_contract < f_values[n + 1]) {
        simplex[[n + 1]] <- x_contract
        f_values[n + 1] <- f_contract
        f_values_list[[n + 1]] <- f_contract
      } else {
        # Shrink the simplex
        for (i in 2:(n + 1)) {
          simplex[[i]] <- simplex[[1]] + sigma * (simplex[[i]] - simplex[[1]])
          f_values_list[[i]] <- f(simplex[[i]], ...)
          f_values[i] <- f_values_list[[i]]
          
        }
      }
      
    }
    
    # Check convergence
    if ((f_values[n] - f_values[1]) / abs(f_values[1]) < ftol) {
      break
    }
    if (mean(sapply(seq_along(init), function(i) diff(range(sapply(simplex, '[', i))))) < xtol) {
      break
    }
  }
  
  f1 <- f_values[1]
  attr(f1, 'beta') <- attr(simplex[[1]], 'beta')
  attr(f1, 'Hessian') <- attr(f_values_list[[1]], 'Hessian')
  list(par = simplex[[1]], objective = f1, iterations = iter, beta = attr(simplex[[1]], 'beta'))
}


## function for initial basis function coefficients

.give_beta0 <- function(G) {
p <- ncol(G$X)
here <- which.min(colSums((G$X - 1)^2))
G <- list(wts=G$wts, X=matrix(1, nrow(G$X), 1), XT=matrix(1, nrow(G$XT), 1), control=G$control)
init <- seq(-10, 10)
f.test <- sapply(init, .f0, dat=G)
if (any(f.test != 1e20)) {
  init <- init[which.min(f.test)]
} else {
  stop("Can't find sensible rate starting values in [1e-10, 1e10]")
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

.g <- function(pars, dat) {
  # gradient of penalised negative log-likelihood 
  # for point process model
  gH <- .gH0(pars, dat)
  bS <- crossprod(pars, dat$S)[1,]
  gH$g + bS
}

.H <- function(pars, dat) {
  # Hessian of penalised negative log-likelihood 
  # for point process model
  gH <- .gH0(pars, dat)
  gH$H + dat$S
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
# if (any(!is.finite(out))) {
#   browser()
#   print('fixing')
#   attr(out, "gradient") <- rep(NA, length(gH[[1]]))
# }
out
}

## REML functions
# slightly adapted version of that from evgam

.refine_rho <- function(rho0, dat) {
rho00 <- rho0
f0 <- .reml0(rho0, dat = dat)
attr(rho0, "beta") <- attr(f0, 'beta')
f1 <- 0 * as.vector(rho0)
for (i in seq_along(rho0)) {
  rho1 <- rho0
  rho1[i] <- rho0[i] + 1
  f1[i] <- .reml0(rho1, dat = dat)
}
adder <- c(0, 1)[1 + as.numeric(f1 < f0)]
other_way <- which(f1 > f0)
if (length(other_way) > 0) {
  for (i in other_way) {
    rho1 <- rho0
    rho1[i] <- rho0[i] - 1
    f1[i] <- .reml0(rho1, dat = dat)
  }
  not_adder <- which(f1[other_way] < f0)
  if (length(not_adder) > 0)
    adder[other_way[not_adder]] <- -1
}
attr(rho0, "beta") <- attr(f0, 'beta')
cond <- TRUE
it <- 0
while(cond & it < 3) {
  rho1 <- rho0 + adder
  f1 <- .reml0(rho1, dat = dat)
  if (f1 < f0) {
    rho0 <- rho1
    f0 <- f1
    attr(rho0, "beta") <- attr(f0, 'beta')
    it <- it + 1
  } else {
    cond <- FALSE
  }
}
rho0
}

.reml0 <- function(pars, dat, beta=NULL, skipfit=FALSE) {
  if (is.null(beta)) beta <- attr(pars, "beta")
  sp <- exp(pars)
  dat$S <- evgam:::.makeS(dat$Sd, sp)
  if (!skipfit) {
    fitbeta <- evgam:::.newton_step_inner(beta, .f, .search, dat = dat, control = dat$control$inner)
    # fitbeta <- try(evgam:::.newton_step_inner(beta, .f, .search, dat = dat, control = dat$control$inner), silent = TRUE)
    # if (inherits(fitbeta, 'try-error'))
    #   fitbeta <- list(par = beta, gradconv = FALSE)
    if (!fitbeta$gradconv) {
      fitbeta0 <- nlminb(fitbeta$par, .f, .g, .H, dat = dat)
      fitbeta <- evgam:::.newton_step_inner(fitbeta0$par, .f, .search, dat = dat, control = dat$control$inner)
      # fitbeta <- try(evgam:::.newton_step_inner(fitbeta0$par, .f, .search, dat = dat, control = dat$control$inner), silent = TRUE)
      # if (inherits(fitbeta, 'try-error'))
      #   return(1e20)
    }
  } else {
    fitbeta <- list(objective=.f(beta, dat))
    fitbeta$convergence <- 0
    fitbeta$gH <- .gH(beta, dat)
    fitbeta$Hessian <- fitbeta$gH[[2]]
    fitbeta$par <- beta
  }
  logdetSdata <- evgam:::.logdetS(dat$Sd, pars)
  logdetHdata <- evgam:::.d0logdetH(fitbeta, FALSE)
  # halflogdetHdata <- try(list(d0=sum(log(diag(chol(fitbeta$Hessian))))), silent=TRUE)
  # if (inherits(halflogdetHdata, "try-error")) return(1e20)
  out <- fitbeta$objective + as.numeric(fitbeta$convergence != 0) * 1e20
  # out <- out + halflogdetHdata$d0 - .5 * logdetSdata$d0
  out <- out + .5 * logdetHdata$d0 - .5 * logdetSdata$d0
  out <- as.vector(out)
  if (!is.finite(out)) return(1e20)
  attr(out, "beta") <- fitbeta$par
  attr(out, "gradient") <- fitbeta$gradient
  attr(out, "Hessian") <- fitbeta$Hessian
  attr(out, "iterations") <- fitbeta$iterations
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
eV <- eigen(H, symmetric = TRUE)
t1 <- 1 / eV$values
t1[eV$values < sqrt(.Machine$double.eps)] <- 0
iH <- crossprod(t(eV$vectors) * sqrt(t1))
spSlb <- sapply(spSl, function(x) x %*% beta)
db <- crossprod(iH, spSlb)
thirdpp <- dat$wts * exp(tcrossprod(beta, dat$XT))[1,]
thirdpp <- (dat$XT * as.vector(thirdpp)) %*% db
XTV <- dat$XT %*% eV$vectors
# XTV <- t(solve(eV$vectors, t(dat$XT)))
dH <- sapply(seq_along(sp), function(i) sum(colSums(XTV * XTV * thirdpp[,i]) / eV$values))
dH <- sapply(spSl, function(x) sum(iH * x)) - dH
dH <- list(d1=dH)
dS <- evgam:::.logdetS(dat$Sd, pars, deriv=1)
d1 <- .5 * sapply(spSl, function(x) crossprod(beta, x %*% beta))
d1 <- d1 - .5 * dS$d1
d1 <- d1 + .5 * dH$d1
d1
}

.reml1_fd <- function(pars, dat, H = NULL, beta = NULL, eps = 5e-3) {
  if (is.null(beta)) {
    beta <- attr(pars, "beta")
  } else {
    attr(pars, "beta") <- beta
  }
  f0 <- .reml0(pars, dat, H)
  b0 <- attr(f0, 'beta')
  eps <- rep(eps, length(pars))
  fu <- fl <- numeric(length(pars))
  for (i in seq_along(pars)) {
    pe <- replace(pars, i, pars[i] + eps[i])
    attr(pe, "beta") <- b0
    fu[i] <- .reml0(pe, dat, H)
  }
  for (i in seq_along(pars)) {
    pe <- replace(pars, i, pars[i] - eps[i])
    attr(pe, "beta") <- b0
    fl[i] <- .reml0(pe, dat, H)
  }
  .5 * as.vector(fu - fl)/eps
}

## other functions

.pivchol_rmvn <- function(n, mu, Sig) {
  R <- suppressWarnings(chol(Sig, pivot = TRUE))
  piv <- order(attr(R, "pivot"))  ## reverse pivoting index
  r <- attr(R, "rank")  ## numerical rank
  V <- R[1:r, piv]
  Y <- crossprod(V, matrix(rnorm(n * r), r))
  Y + as.vector(mu)
}

