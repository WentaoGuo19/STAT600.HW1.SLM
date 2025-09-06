test_that("coefficients match lm()", {
  set.seed(123)
  n = 1000
  x = rnorm(n)
  y = 1.5 + 2.0 * x + rnorm(n, sd = 0.3)

  fit_pkg = SimpLinR(x, y)
  fit_lm  = lm(y ~ x)

  # coefficients
  expect_equal(unname(coef(fit_lm)[1]), unname(fit_pkg$coefficients[1]), tolerance = 1e-8)
  expect_equal(unname(coef(fit_lm)[2]), unname(fit_pkg$coefficients[2]), tolerance = 1e-8)

  # SEs
  se_lm = sqrt(diag(vcov(fit_lm)))
  expect_equal(unname(se_lm[1]), unname(fit_pkg[["coefficients se"]][1]), tolerance = 1e-8)
  expect_equal(unname(se_lm[2]), unname(fit_pkg[["coefficients se"]][2]), tolerance = 1e-8)

  # 95% CI
  ci_lm = confint(fit_lm, level = 0.95)
  ci_pkg = fit_pkg[["coefficients 95% CI"]]
  expect_equal(unname(ci_pkg[1,1]), unname(ci_lm[1,1]), tolerance = 1e-8)
  expect_equal(unname(ci_pkg[1,2]), unname(ci_lm[1,2]), tolerance = 1e-8)
  expect_equal(unname(ci_pkg[2,1]), unname(ci_lm[2,1]), tolerance = 1e-8)
  expect_equal(unname(ci_pkg[2,2]), unname(ci_lm[2,2]), tolerance = 1e-8)
})


test_that("NA/Inf handling and lengths are preserved", {
  x = c(1:5, NA, 7, 8, Inf, 10, 11:50)
  y = 0.5 + 3 * x + rnorm(50)
  y[2] = NA
  y[9] = 0

  fit = SimpLinR(x, y)

  expect_length(fit$fitted, length(x))
  expect_length(fit$residuals, length(x))

  bad = !(is.finite(x) & is.finite(y))
  expect_true(all(is.na(fit$fitted[bad])))
  expect_true(all(is.na(fit$residuals[bad])))

  # Residuals should sum to ~0 over finite rows (with intercept)
  expect_equal(sum(fit$residuals, na.rm = TRUE), 0, tolerance = 1e-8)
})


test_that("predictor name propagates (if wrapper renames outputs)", {
  z = rnorm(50); y = 2 + 4 * z + rnorm(50, 0.2)
  fit = SimpLinR(z, y)

  expect_true("(Intercept)" %in% names(fit$coefficients))
  expect_true("z" %in% names(fit$coefficients))
  rn = rownames(fit[["coefficients 95% CI"]])
  expect_true("(Intercept)" %in% rn)
  expect_true("z" %in% rn)
})
