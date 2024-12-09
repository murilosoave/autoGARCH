#' Compute Basic Statistics and Tests
#'
#' Computes skewness, kurtosis, performs Jarque-Bera and Shapiro-Wilk tests.
#'
#' @param log_returns A numeric vector of log-returns.
#'
#' @return A list with skewness, kurtosis, jb_test, and shapiro_test results.
#' @import moments
#' @import tseries
#' @export
compute_stats <- function(log_returns) {
  log_mean <- mean(log_returns, na.rm = TRUE)
  log_sd <- sd(log_returns, na.rm = TRUE)
  log_skew <- moments::skewness(log_returns, na.rm = TRUE)
  log_kurt <- moments::kurtosis(log_returns, na.rm = TRUE)
  jb_test <- tseries::jarque.bera.test(log_returns)
  shapiro_test <- shapiro.test(log_returns)

  list(
    skewness = log_skew,
    kurtosis = log_kurt,
    jb_test = jb_test,
    shapiro_test = shapiro_test,
    mean = log_mean,
    sd = log_sd
  )
}

#' Plot Distribution Comparison
#'
#' Plots histogram of log-returns and overlays a normal density curve.
#'
#' @param log_returns A numeric vector of log-returns.
#' @param mean Mean of the log-returns.
#' @param sd Standard deviation of the log-returns.
#'
#' @export
plot_distribution <- function(log_returns, mean, sd) {
  hist(log_returns, 
       breaks = 50, 
       freq = FALSE, 
       main = "Distribuição dos Log-Retornos vs Normal",
       xlab = "Log-Retornos")
  x_vals <- seq(min(log_returns, na.rm = TRUE), max(log_returns, na.rm = TRUE), length.out = 200)
  normal_density <- dnorm(x_vals, mean = mean, sd = sd)
  lines(x_vals, normal_density, col = "red", lwd = 2)
  legend("topright", legend = c("Distribuição Empírica","Distribuição Normal"), 
         col = c("black", "red"), lty = c(NA, 1), pch = c(15, NA), lwd = c(NA,2))
}
