#' Objective Function for Model Selection
#'
#' This function evaluates a given ARIMA-GARCH specification, checking multiple assumptions.
#'
#' @param p_arima,q_arima,garch_p,garch_q Distribution indexes and ARIMA/GARCH parameters.
#' @param distribution_idx Index for distribution model (1=norm,2=std).
#' @return A list with evaluation metrics.
#' @import tseries
#' @import rugarch
#' @import lmtest
#' @export
objective_function <- function(p_arima, q_arima, garch_p, garch_q, distribution_idx, log_returns) {
  p_arima <- round(p_arima)
  q_arima <- round(q_arima)
  garch_p <- round(garch_p)
  garch_q <- round(garch_q)
  dist_choice <- distribution_mapping[round(distribution_idx)]

  adf_result <- tseries::adf.test(log_returns, alternative = "stationary")
  stationary_assumption <- (adf_result$p.value <= 0.05)
  d_arima <- 0
  if (p_arima < 0 || q_arima < 0) {
    return(list(
      Score = -100000, AIC = 100000,
      Stationary = stationary_assumption,
      ARIMA_No_AutoCorr = FALSE,
      ARIMA_All_Significant = FALSE,
      GARCH_Fit_Success = FALSE,
      GARCH_No_AutoCorr = TRUE,
      GARCH_No_Arch_Effects = TRUE,
      GARCH_Res_Normal = TRUE,
      GARCH_All_Significant = TRUE
    ))
  }

  arima_model <- tryCatch(
    arima(log_returns, order = c(p_arima, d_arima, q_arima), include.mean = TRUE),
    error = function(e) NULL
  )

  if (is.null(arima_model)) {
    return(list(
      Score = -100000, AIC = 100000,
      Stationary = stationary_assumption,
      ARIMA_No_AutoCorr = FALSE,
      ARIMA_All_Significant = FALSE,
      GARCH_Fit_Success = FALSE,
      GARCH_No_AutoCorr = TRUE,
      GARCH_No_Arch_Effects = TRUE,
      GARCH_Res_Normal = TRUE,
      GARCH_All_Significant = TRUE
    ))
  }

  arima_res <- residuals(arima_model)
  lb_test <- Box.test(arima_res, lag=20, type="Ljung-Box")
  arima_no_autocorr <- (lb_test$p.value >= 0.05)

  arima_coefs <- lmtest::coeftest(arima_model)
  non_significant_arima <- arima_coefs[,4] > 0.05
  arima_all_significant <- !any(non_significant_arima)

  if ((garch_p > 0 || garch_q > 0)) {
    ao <- arimaorder(arima_model)
    spec <- rugarch::ugarchspec(
      variance.model = list(model = "sGARCH", garchOrder = c(garch_p, garch_q)),
      mean.model = list(armaOrder = c(ao[1], ao[3]), include.mean = TRUE),
      distribution.model = dist_choice
    )

    fit <- tryCatch(
      rugarch::ugarchfit(spec = spec, data = log_returns),
      error = function(e) NULL
    )

    if (is.null(fit)) {
      return(list(
        Score = -100000, AIC = 100000,
        Stationary = stationary_assumption,
        ARIMA_No_AutoCorr = arima_no_autocorr,
        ARIMA_All_Significant = arima_all_significant,
        GARCH_Fit_Success = FALSE,
        GARCH_No_AutoCorr = TRUE,
        GARCH_No_Arch_Effects = TRUE,
        GARCH_Res_Normal = TRUE,
        GARCH_All_Significant = TRUE
      ))
    }

    garch_coefs <- fit@fit$matcoef
    garch_param_rows <- grepl("alpha|beta", rownames(garch_coefs), ignore.case = TRUE)
    garch_all_significant <- TRUE
    if (any(garch_param_rows)) {
      garch_param_pvalues <- garch_coefs[garch_param_rows, 4]
      garch_all_significant <- all(garch_param_pvalues <= 0.05)
    }

    garch_res <- residuals(fit, standardize=TRUE)
    lb_test_garch <- Box.test(garch_res, lag=20, type="Ljung-Box")
    garch_no_autocorr <- (lb_test_garch$p.value >= 0.05)

    arch_test_garch <- tseries::ArchTest(garch_res, lags=5)
    garch_no_arch_effects <- (arch_test_garch$p.value >= 0.05)

    jb_test <- tseries::jarque.bera.test(garch_res)
    garch_res_normal <- (jb_test$p.value >= 0.05)

    model_aic <- rugarch::infocriteria(fit)[1]
    Score <- -model_aic

    return(list(
      Score = Score, AIC = model_aic,
      Stationary = stationary_assumption,
      ARIMA_No_AutoCorr = arima_no_autocorr,
      ARIMA_All_Significant = arima_all_significant,
      GARCH_Fit_Success = TRUE,
      GARCH_No_AutoCorr = garch_no_autocorr,
      GARCH_No_Arch_Effects = garch_no_arch_effects,
      GARCH_Res_Normal = garch_res_normal,
      GARCH_All_Significant = garch_all_significant
    ))
  } else {
    jb_test_arima <- tseries::jarque.bera.test(arima_res)
    arima_res_normal <- (jb_test_arima$p.value >= 0.05)
    model_aic <- AIC(arima_model)
    Score <- -model_aic

    return(list(
      Score = Score, AIC = model_aic,
      Stationary = stationary_assumption,
      ARIMA_No_AutoCorr = arima_no_autocorr,
      ARIMA_All_Significant = arima_all_significant,
      GARCH_Fit_Success = FALSE,
      GARCH_No_AutoCorr = TRUE,
      GARCH_No_Arch_Effects = TRUE,
      GARCH_Res_Normal = TRUE,
      GARCH_All_Significant = TRUE
    ))
  }
}


#' Grid Search for ARIMA-GARCH Hyperparameters
#'
#' Performs a grid search over given parameter ranges, returning models that satisfy assumptions.
#'
#' @param log_returns Numeric vector of log-returns.
#' @param p_arima_values Range of p for ARIMA.
#' @param q_arima_values Range of q for ARIMA.
#' @param garch_p_values Range of p for GARCH.
#' @param garch_q_values Range of q for GARCH.
#' @param distribution_idx_values Distribution indexes.
#'
#' @return A data.frame with results.
#' @import rugarch
#' @export
grid_search_models <- function(log_returns, 
                               p_arima_values=0:3, 
                               q_arima_values=0:3,
                               garch_p_values=1:3,
                               garch_q_values=1:3,
                               distribution_idx_values=2:2) {
  param_grid <- expand.grid(
    p_arima = p_arima_values,
    q_arima = q_arima_values,
    garch_p = garch_p_values,
    garch_q = garch_q_values,
    distribution_idx = distribution_idx_values
  )

  best_score <- -Inf
  best_params <- NULL
  results_list <- list()

  for (i in seq_len(nrow(param_grid))) {
    p_arima <- param_grid$p_arima[i]
    q_arima <- param_grid$q_arima[i]
    garch_p <- param_grid$garch_p[i]
    garch_q <- param_grid$garch_q[i]
    distribution_idx <- param_grid$distribution_idx[i]

    out <- tryCatch(
      objective_function(p_arima, q_arima, garch_p, garch_q, distribution_idx, log_returns),
      error = function(e) {
        list(
          Score = -1e9,
          AIC = 1e9,
          Stationary = FALSE,
          ARIMA_No_AutoCorr = FALSE,
          ARIMA_All_Significant = FALSE,
          GARCH_Fit_Success = FALSE,
          GARCH_No_AutoCorr = TRUE,
          GARCH_No_Arch_Effects = TRUE,
          GARCH_Res_Normal = TRUE,
          GARCH_All_Significant = TRUE
        )
      }
    )

    result_row <- data.frame(param_grid[i, ], as.data.frame(out, stringsAsFactors=FALSE))
    results_list[[i]] <- result_row

    current_score <- out$Score
    if (!is.na(current_score) && current_score > best_score) {
      best_score <- current_score
      best_params <- result_row
    }
  }

  results_df <- do.call(rbind, results_list)
  results_df
}
