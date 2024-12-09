#' Prepare Log Returns
#'
#' Takes a vector of prices and computes log-returns, removing NAs.
#'
#' @param prices A numeric vector of prices.
#'
#' @return A numeric vector of log-returns.
#' @export
prepare_log_returns <- function(prices) {
  ts_prices <- ts(prices)
  log_returns <- diff(log(ts_prices))
  log_returns <- na.omit(log_returns)
  log_returns
}
