#' Forecast VaR Using Selected Model
#'
#' Given a fitted ARIMA-GARCH model, this function forecasts N steps ahead and computes the VaR.
#'
#' @param fit A fitted ugarchfit object.
#' @param log_returns The log-returns used in the model.
#' @param N Number of steps ahead.
#' @param alpha VaR confidence level (default 0.05).
#' @param p_t The last observed price of the asset.
#'
#' @return A list with the VaR values and price scenarios.
#' @import rugarch
#' @export
forecast_var <- function(fit, log_returns, N=10, alpha=0.05, p_t=100) {
  fcst <- rugarch::ugarchforecast(fit, n.ahead = N)
  fmean <- fitted(fcst)
  fsigma <- sigma(fcst)
  params <- coef(fit)
  nu <- params["shape"]
  q_alpha <- qt(alpha, df = nu)
  VaR_log <- fmean + fsigma * q_alpha
  cumulative_logVaR <- cumsum(VaR_log)
  price_scenario <- p_t * exp(cumulative_logVaR)
  baseline_price <- p_t * exp(cumsum(fmean))
  price_VaR <- p_t - price_scenario
  list(price_scenario=price_scenario, baseline_price=baseline_price, price_VaR=price_VaR)
}
