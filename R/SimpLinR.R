#' Simple univariate OLS
#'
#' @description Runs a simple OLS Gaussian linear regression for one numeric
#' predictor on one numeric response.
#'
#' @param x Numeric vector of predictor values.
#' @param y Numeric vector of response values.
#'
#' @return A list with coefficient estimates, standard errors, 95% confidence
#' intervals, fitted values, and residuals.
#'
#' @examples
#' x = rnorm(100); y = 1 + 2*x + rnorm(100, 0.5)
#' fit = SimpLinR(x, y)
#' str(fit)
#'
#' @export
SimpLinR = function(x, y) {
  #check if (x,y) have same length & are numeric
  if (!is.numeric(x) || !is.numeric(y)) stop("x and y must be numeric")
  if (length(x) != length(y)) stop("x and y must have the same length")

  #discard NA/NaN/Inf inputs
  legal_index = is.finite(x) & is.finite(y)
  x_legal = x[legal_index]; y_legal = y[legal_index]

  #check to ensure n>=10 & xs are not all identical
  if (length(x_legal) < 10L) stop("need at least 10 observations")
  if (stats::var(x_legal) == 0) stop("x has zero variance, cannot estimate slope")

  #perform OLS LM in Rcpp function SimpLinCpp()
  out = .SimpLinCpp(x_legal, y_legal)

  #include NAs for fitted values & residuals in the outputs if those places are NAs in the inputs
  resids_full = fitted_full = rep(NA_real_, length(x))
  resids_full[legal_index] = out$residuals; fitted_full[legal_index] = out$fitted
  out$residuals = resids_full
  out$fitted = fitted_full

  #name the predictor dynamically
  xname = deparse(substitute(x))
  names(out$coefficients) <- c("(Intercept)", xname)
  names(out[["coefficients se"]]) = c("(Intercept)", xname)
  rownames(out[["coefficients 95% CI"]]) = c("(Intercept)", xname)

  return(out)
}
