#' Nakagami distribution
#'
#' The Nakagami distribution is the distribution of the square-root of a Gamma-distributed variable.
#' That is, f X ~ Gamma, then sqrt(X) ~ Nakagami.
#' 
#' It can also be used to approximate the distribution of the standard deviation of split frequencies.
#'
#' @param x, q vector of quantiles
#' @param p vector of probabilities.
#' @param n number of observations. If length(n) > 1, the length is taken to be the number required.
#' @param shape shape parameter > 0.0
#' @param spread spread parameter > 0.0
#' @param log, log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P(X <= x), otherwise P(X > x).
#' 
#' @details
#' The Nakagami distribution with shape m and scale o has density
#' f(x) = 2/gamma(m) \* (m/o)^m \* x^(2m - 1) \* exp(-m/o \* x^2)
#' 
#' @seealso \link{}
#' @export
dnakagami <- function(x, shape, spread, log = FALSE) {
  dens <- log(2) + shape * log(shape) - lgamma(shape) - shape * log(spread) +
    (2 * shape - 1) * log(x) - (shape / spread * x * x)
  if (!log) {
    dens <- exp(dens)
  }
  
  return(dens)
}

#' @rdname dnakagami
#' @export
pnakagami <- function(q, shape, spread, lower.tail = TRUE, log.p = FALSE) {
  gpar <- nakagami2gamma(shape, spread)
  stats::pgamma(
    q * q,
    shape = gpar["shape"],
    scale = gpar["scale"],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @rdname dnakagami
#' @export
qnakagami <- function(p, shape, spread, lower.tail = TRUE, log.p = FALSE) {
  gpar <- nakagami2gamma(shape, spread)
  q <- stats::qgamma(
    p,
    shape = gpar["shape"],
    scale = gpar["scale"],
    lower.tail = lower.tail,
    log.p = log.p
  )
  return(sqrt(q))
}

#' @rdname dnakagami
#' @export
rnakagami <- function(n, shape, spread) {
  gpar <- nakagami2gamma(shape, spread)
  x <- stats::rgamma(
    n,
    shape = gpar["shape"],
    scale = gpar["scale"],
  )
  return(sqrt(x))
}

#' Conversion between Nakagami and Gamma parameters
#'
#' @param shape, Nakagami or Gamma shape parameters
#' @param spread Nakagami spread parameter
#' @param scale, rate Gamma scale or rate parameters
#' 
#' @keywords internal
nakagami2gamma <- function(shape, spread, return.scale = TRUE) {
  res <- c(shape, spread / shape)
  names(res) <- c("shape", "scale")
  if (!return.scale) {
    res[2] <- 1.0 / res[2]
    names(res)[2] <- "rate"
  }
  return(res)
}

#' @rdname nakagami2gamma
#' @keywords internal
gamma2nakagami <- function(shape, rate = NA, scale = NA) {
  if (!is.numeric(rate) && !is.numeric(scale)) {
    stop("Must specify rate or scale")
  }
  if (!is.numeric(scale)) {
    scale <- 1.0 / rate
  }
  res <- c(shape, shape * scale)
  names(res) <- c("shape", "spread")
  return(res)
}

#' Moments of Nakagami distribution
#'
#' @param shape, shape parameter
#' @param spread spread parameter
#' 
#' @return mean and variance as named numeric vector
#' 
#' @seealso \link{dnakagami}
#' 
#' @export
nakagami.moms <- function(shape, spread) {
  m <- exp(
    lgamma(shape + 0.5) - lgamma(shape) + 0.5 * log(spread / shape)
  )
  v <- spread * (1 - 1 / shape * exp(
    2 * (lgamma(shape + 0.5) - lgamma(shape))
  ))
  res <- c(m, v)
  names(res) <- c("mean", "variance")
  
  return(res)
}