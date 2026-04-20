library(data.table)
library(future.apply)

yb_resample <- function(y, data_mat_std, p_vals, n_boot) {
  order_indx <- order(p_vals)
  p_ordered <- p_vals[order_indx]

  n <- length(y)
  p <- ncol(data_mat_std)
  loc_vec <- 1:p

  resample_ind <- sample(n, size = n_boot * n, replace = TRUE)
  y_star <- matrix(y[resample_ind], nrow = n)

  tstat_star <- future_apply(y_star,
    2, 
    fast_fit_lm_tstat, 
    data_mat = data_mat_std)
  
  p_star <- 2 * pt(abs(tstat_star),
                df = n - 2,
                lower.tail = FALSE)

  r_p <- round(ecdf(p_vals)(p_ordered) * p)

  R_p_star <- matrix(0, nrow = p, ncol = n_boot)

  for (i in 1:n_boot) {
   Fn <- ecdf(p_star[i, ])
    R_p_star[, i] <- round(Fn(p_ordered) * p)
  }

  upper_quant <- 1 - 0.05

  r_beta_p_star <- apply(R_p_star, 1, \(x) quantile(x, upper_quant))
  threshold <- r_p - r_beta_p_star

  prob_nonzero <- apply(R_p_star >= 1, 1, mean)
  Q_local_est <- prob_nonzero

  indices <- which(threshold >= p_ordered * p)
  if (length(indices) > 0) {
    numerator <- R_p_star[indices, ]
    correction <- r_p[indices] - (p_ordered[indices] * p)
    denom <- sweep(R_p_star[indices, ], 1, correction, "+")

    Q_star <- numerator / denom
    Q_mean <- apply(Q_star, 1, mean)

    Q_local_est[indices] <- Q_mean
  }

  reject_bool <- rep(FALSE, p)

  if (sum(Q_local_est <= 0.1) > 0) {
    max_rej_ind <- max(which(Q_local_est <= 0.1))
    ordered_loc <- loc_vec[order_indx]
    rejected_loc <- ordered_loc[1:max_rej_ind]

    reject_bool[rejected_loc] <- TRUE
  } 

  yb_reject <- data.table(location = loc_vec, 
    p_val = p_vals, 
    BH_reject = reject_bool)

  return(yb_reject)
}

run_sim_yb <- function(iter,
                        lat_lon_grid,
                        szon_dist,
                        smer_dist,
                        loc_list,
                        cov_mat,
                        beta,
                        n,
                        cov_L,
                        noise_sd,
                        n_boot,
                        verbose = FALSE) {
  
  p <- nrow(lat_lon_grid)
  
  results_mat <- matrix(0, nrow = iter * 20, ncol = 19)
  
  init_cond_mat <- matrix(c(1, 1, 1,
                            0.2 * pi, 0.5 * pi, 0.8 * pi,
                            1.5, 1.5, 1.5),
                          nrow = 3)
  
  rnd_indx <- c(replicate(iter, sample(p, size = 20, replace = FALSE)))
  
  loc_pairs_list <- lapply(loc_list,
                           \(x) find_dist_pairs(zonal_dist_mat = szon_dist,
                                            merid_dist_mat = smer_dist,
                                            loc_vec = x))
  
  
  for (i in 1:iter) {
    if (verbose) print(i)
    
    sim_data <- gen_data(n = n, L = cov_L)
    sim_data_std <- scale(sim_data)
    
    MLE_out <- future_mapply(find_MLE_mapply_wrapper,
                             loc_block = loc_list,
                             pairs_mat = loc_pairs_list,
                             MoreArgs = list(data_mat = sim_data_std,
                                             lik_func = m1.5_lik),
                             SIMPLIFY = FALSE)
    
    MLE_df <- rbindlist(MLE_out)
    MLE_mat <- data.matrix(MLE_df[, 1:3])
    
    cov_weights <- find_smooth_weights(loc_list = loc_list,
                                       zonal_sdist = szon_dist,
                                       merid_sdist = smer_dist,
                                       cov_func = matern_1.5_cov,
                                       mle_params = MLE_mat,
                                       progress = FALSE)
    
    smooth_data <- sim_data_std %*% cov_weights
    
    base_indx <- (i - 1) * 20
    
    for (k in 1:20) {
      if (verbose) print(k)
      
      true_cov <- cov_mat[rnd_indx[base_indx + k], ]
      true_null <- (true_cov == 0)
      
      y <- beta * sim_data[, rnd_indx[base_indx + k]] + rnorm(n, sd = noise_sd)
      
      base_lm_fits <- fast_fit_lm_tstat(y = y, data_mat = sim_data_std)
      smooth_lm_fits <- fast_fit_lm_tstat(y = y, data_mat = smooth_data)
      
      # two-tailed test
      base_pvals <- 2 * pt(abs(base_lm_fits),
                           df = n - 2,
                           lower.tail = FALSE)
      smooth_pvals <- 2 * pt(abs(smooth_lm_fits),
                             df = n - 2,
                             lower.tail = FALSE)
      
      base_reject <- find_BH_reject(base_pvals, alpha = 0.1)
      smooth_reject <- find_BH_reject(smooth_pvals, alpha = 0.1)
      yb_reject <- yb_resample(y, sim_data_std, base_pvals, n_boot = 1000)
      
      base_metrics <- estim_metrics(p_dt = base_reject,
                                    true_null = true_null)
      smooth_metrics <- estim_metrics(p_dt = smooth_reject,
                                      true_null = true_null)
      yb_metrics <- estim_metrics(p_dt = yb_reject,
                                  true_null = true_null)
      
      results_mat[base_indx + k, ] <- c(i,
                                        base_metrics,
                                        smooth_metrics,
                                        yb_metrics)
    }
  }
  
  results_frame <- data.table(results_mat)
  
  metric_names <- c("reject",
                    "true_reject",
                    "FDP",
                    "FNP",
                    "sens",
                    "spec")
  
  names(results_frame) <- c("iteration",
                            paste0("base_", metric_names),
                            paste0("smooth_", metric_names),
                            paste0("yb_", metric_names))
  
  return(results_frame)
}