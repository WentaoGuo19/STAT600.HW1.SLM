#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export(name = ".SimpLinCpp")]]
List SimpLinCpp(NumericVector x, NumericVector y) {
  // sample size
  const int n = x.size();

  // means
  double Sx = 0.0, Sy = 0.0;
  for (int i = 0; i < n; ++i) {
    Sx += x[i];
    Sy += y[i];
  }
  const double xbar = Sx / n, ybar = Sy / n;

  // Sxx & Sxy
  double Sxx = 0.0, Sxy = 0.0;
  for (int i = 0; i < n; ++i) {
    const double dx = x[i] - xbar;
    const double dy = y[i] - ybar;
    Sxx += dx * dx;
    Sxy += dx * dy;
  }

  // coefficients estimates
  const double b1 = Sxy / Sxx;
  const double b0 = ybar - b1 * xbar;

  // fitted values & residuals
  NumericVector fitted(n), resid(n);
  for (int i = 0; i < n; ++i) {
    fitted[i] = b0 + b1 * x[i];
    resid[i]  = y[i] - fitted[i];
  }

  // error variance & coefficients' standard errors
  double RSS = 0.0;
  for (int i = 0; i < n; ++i) RSS += resid[i] * resid[i];
  const int df = n - 2;
  const double sigma2 = RSS / df;
  const double se_b1 = std::sqrt(sigma2 / Sxx);
  const double se_b0 = std::sqrt(sigma2 * (1.0 / n + (xbar * xbar) / Sxx));

  // 95% CIs
  const double tcrit = R::qt(0.975, df, /*lower_tail=*/1, /*log_p=*/0);

  NumericVector coef = NumericVector::create(
    _["(Intercept)"] = b0,
    _["x"] = b1
  );
  NumericVector se = NumericVector::create(
    _["(Intercept)"] = se_b0,
    _["x"] = se_b1
  );

  NumericMatrix ci(2, 2);
  ci(0, 0) = b0 - tcrit * se_b0; ci(0, 1) = b0 + tcrit * se_b0;
  ci(1, 0) = b1 - tcrit * se_b1; ci(1, 1) = b1 + tcrit * se_b1;
  ci.attr("dimnames") = List::create(
    CharacterVector::create("(Intercept)", "x"),
    CharacterVector::create("lower", "upper")
  );

  return List::create(
    _["coefficients"] = coef,
    _["coefficients se"] = se,
    _["coefficients 95% CI"] = ci,
    _["fitted"] = fitted,
    _["residuals"] = resid
  );
}
