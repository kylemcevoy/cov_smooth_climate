### functions involved in fitting linear models gridbox by gridbox

#' Fit a parallel set of simple linear regressions quickly
#'
#' function takes in a response variable y and a matrix with data from many
#' gridboxes and returns a vector of t-statistics from simple linear regression
#' on each gridbox. Uses lm.fit under the hood for speed, there will be small
#' calculation differences of the t-statistic on the order of 1e-14 to 1e-16
#' between the results of lm_fast_fitter_tstat and fitting using lm().
#'
#' @param y a length n numeric vector containing the response variable
#' @param data_mat an n x p matrix whose rows are the samples at each time point
#' and the columns are data from each individual gridbox.
#' @returns a length p vector of t-statistics for the coefficients of x in the set
#' of separate simple linear regressions of y against each column of data_mat
#' in isolation. The fitted intercept terms are not returned.

fit_lm_tstat <- function(y, data_mat) {
  n <- length(y)
  p <- ncol(data_mat)
  results <- numeric(p)
  data_mat <- cbind(1, data_mat)

  for (i in 1:p) {
    fitted_lm <- lm.fit(y = y, x = data_mat[, c(1, i + 1)])
    beta_x <- fitted_lm$coefficients[2]
    MSE <- sum(fitted_lm$residuals^2) / (n - 2)
    SS_x <- (n - 1) * var(data_mat[, i + 1])
    se_beta <- sqrt(MSE / SS_x)
    results[i] <- beta_x / se_beta
  }
  return(results)
}

fast_fit_lm_tstat <- function(y, data_mat) {
  n <- length(y)
  
  varx_base <- apply(data_mat, 2, var)
  mean_x <- apply(data_mat, 2, mean)
  mean_y <- mean(y)
  
  cov_xy <- c(cov(y, data_mat))
  beta_x <- cov_xy / varx_base
  beta_0 <- mean_y - beta_x * mean_x
  
  # here sweep iterates over the columns and operates on each column with
  # the corresponding element of STATS.
  fitted <- sweep(data_mat, MARGIN = 2, STATS = beta_x, FUN = `*`)
  fitted <- sweep(fitted, MARGIN = 2, STATS = beta_0, FUN = `+`)
  
  res <- y - fitted
  MSE <- apply(res^2, MARGIN = 2, FUN = sum) / (n - 2)
  se_beta <- sqrt(MSE / ((n - 1) * varx_base))
  t_stat <- beta_x / se_beta
  
  return(t_stat)
}

#' Find the correlation between predictors and response
#'
#' This function finds the correlation between the response variable, y, and
#' each predictor variable.
#'
#' @param y length n numeric vector. The response variable of interest.
#' @param data_mat n x p numeric matrix. The p columns are the predictors with
#' n observations of each variable.
#' @returns a length p vector containing the correlations between y and the
#' columns of data_mat.

find_corr <- function(y, data_mat) {
  p <- ncol(data_mat)
  cor(data_mat, y)[, 1]
}

#' Apply Benjamini-Hochberg
#'
#' This function takes a vector of p-values and a false discovery rate control
#' level and applys the Benjamini-Hochberg (BH) algorithm to find discoveries.
#'
#' @param p_val length m vector of p-values for each of the m locations.
#' @param alpha the target level for false discovery control
#' @returns an m x 6 data.table object with columns:
#' location: a numeric sequence of location indices.
#' p_val: the p-values for each of the individual location hypothesis tests.
#' reject: whether or not the hypothesis would have been rejected at 0.05
#' significance level if there was no controlling for multiple hypotheses.
#' threshold: an intermediate threshold that needs to be calculated for the
#' BH procedure.
#' pass: another intermediate calculation step for BH.
#' BH_reject: whether the hypothesis is rejected/location is discovered after
#' application of the BH algorithm.
#'

find_BH_reject <- function(p_val, alpha) {
  m <- length(p_val)
  location <- seq_along(p_val)
  p_dt <- data.table(location, p_val)
  setkey(p_dt, p_val)

  p_dt[, threshold := (1:m) * alpha / m]
  p_dt[, pass := p_val <= threshold]

  if (any(p_dt$pass)) {
    max_cross <- max(which(p_dt$pass))
  } else {
    max_cross <- 0
  }

  p_dt[, BH_reject := (1:m) <= max_cross]
  setkey(p_dt, location)

  return(p_dt[, .(location, p_val, BH_reject)])
}

#' Estimate inference metrics on simulation output
#' 
#' This function takes the output of find_BH_reject and a vector specifying
#' which null hypotheses are true and outputs a number of metrics related
#' to the inference
#' 
#' @param p_dt An m x 6 data.table object outputted from find_BH_reject() containing
#'  the p-values and Benjamini-Hochberg rejections made at a specified level alpha
#' @param true_null a length m logical vector specifying whether the null hypothesis
#' was true at location m.
#' @returns a length 6 numeric vector containing the total number of rejections
#' made by BH for the given p-values, the number of true rejections made by BH,
#' the False Discovery Proportion, the False Non-discovery Proportion, the
#' sensitivity, and specificity.
estim_metrics <- function(p_dt, true_null) {
  m <- nrow(p_dt)

  reject <- sum(p_dt$BH_reject)
  non_reject <- m - reject

  true_reject <- sum(p_dt$BH_reject & (!true_null))
  false_reject <- sum(p_dt$BH_reject & true_null)
  true_nonreject <- sum((!p_dt$BH_reject) & true_null)
  false_nonreject <- sum((!p_dt$BH_reject) & !true_null)

  FDP <- false_reject / max(reject, 1)
  FNP <- false_nonreject / non_reject
  sens <- true_reject / sum(!true_null)
  spec <- true_nonreject / sum(true_null)

  return(c(reject = reject,
           true_reject = true_reject,
           FDP = FDP,
           FNP = FNP,
           sens = sens,
           spec = spec)
         )
}
