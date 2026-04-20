library(future.apply)
library(RhpcBLASctl)
library(data.table)

find_local_MLE_par <- function(loc_block,
                               pairs_mat,
                               data_mat,
                               lik_func = m1.5_lik,
                               return_best = TRUE,
                               return_likelihood = TRUE,
                               init_cond_mat = matrix(c(1, 1, 1,
                                                        0.2 * pi, 0.5 * pi, 0.8 * pi,
                                                        1.5, 1.5, 1.5),
                                                      nrow = 3)) {
  
  data_sub <- data_mat[, loc_block]
  
  m <- nrow(init_cond_mat)
  
  phi_out <- rep(0, m)
  ani_ang_out <- rep(0, m)
  ani_rat_out <- rep(0, m)
  lik_out <- rep(0, m)
  conv_out <- rep(0, m)
  
  for (i in 1:m) {
    optim_out <- tryCatch(
      expr = optim(par = init_cond_mat[i, ],
                   fn = lik_func,
                   control = list(fnscale = -1),
                   data_mat_sub = data_sub,
                   dist_pairs_mat = pairs_mat,
                   method = "L-BFGS-B",
                   lower = c(0, 0, 1),
                   upper = c(Inf, pi, Inf)),
      error = function(cond) {
        warning("Possible error at init. cond.: ",
                i,
                ".")
        return(NA)
      }
    )
    
    if (any(is.na(optim_out))) {
      phi_out[i] <- NA
      ani_ang_out[i] <- NA
      ani_rat_out[i] <- NA
      lik_out[i] <- NA
      conv_out[i] <- NA
    } else {
      
      # can uncomment if you want to find convergence issues for all inits.
      # Generally the best of the three convergence points does not have a 
      # non-zero optim status.
      
      # if (optim_out$convergence != 0) {
      #   warning(paste0("optim convergence status was ",
      #                  optim_out$convergence,
      #                  " at location ",
      #                  loc,
      #                  ".  ",
      #                  optim_out$message)
      #   )
      # }
      
      phi_out[i] <- optim_out$par[1]
      ani_ang_out[i] <- optim_out$par[2]
      ani_rat_out[i] <- optim_out$par[3]
      lik_out[i] <- optim_out$value
      conv_out[i] <- optim_out$convergence
    }
  }
  optim_results_out <- data.frame(phi = phi_out,
                                  ani_angle = ani_ang_out,
                                  ani_ratio = ani_rat_out,
                                  convergence = conv_out)
  if (return_likelihood) {
    optim_results_out$likelihood <- lik_out
  }
  
  if (return_best) {
    max_indx <- which.max(lik_out)
    optim_results_out <- optim_results_out[max_indx, ]
  }
  
  return(optim_results_out)
}

find_MLE_mapply_wrapper <- function(loc_block,
                                    pairs_mat,
                                    data_mat,
                                    lik_func) {
  RhpcBLASctl::blas_set_num_threads(1)
  
  find_local_MLE_par(loc_block = loc_block,
                     pairs_mat = pairs_mat,
                     data_mat = data_mat,
                     lik_func = lik_func,
                     return_best = TRUE,
                     return_likelihood = FALSE,
                     init_cond_mat = matrix(c(1, 1, 1,
                                              0.2 * pi, 0.5 * pi, 0.8 * pi,
                                              1.5, 1.5, 1.5),
                                            nrow = 3))
}

run_sim_par <- function(iter,
                        lat_lon_grid,
                        szon_dist,
                        smer_dist,
                        loc_list,
                        cov_mat,
                        beta,
                        n,
                        cov_L,
                        noise_sd,
                        cor_thresh = 0,
                        verbose = FALSE) {
  
  p <- nrow(lat_lon_grid)
  
  results_mat <- matrix(0, nrow = iter * 1000, ncol = 13)
  
  init_cond_mat <- matrix(c(1, 1, 1,
                            0.2 * pi, 0.5 * pi, 0.8 * pi,
                            1.5, 1.5, 1.5),
                          nrow = 3)
  
  rnd_indx <- c(replicate(iter, sample(p, size = 1000, replace = FALSE)))
  
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
    
    base_indx <- (i - 1) * 1000
    
    for (k in 1:1000) {
      if (verbose && k %% 100 == 0) print(k)
      
      true_cov <- cov_mat[rnd_indx[base_indx + k], ]
      if (cor_thresh == 0) {
        true_null <- (true_cov == 0)
      }
      else {
        true_null <- (true_cov <= cor_thresh)
      }
      
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
      
      base_metrics <- estim_metrics(p_dt = base_reject,
                                    true_null = true_null)
      smooth_metrics <- estim_metrics(p_dt = smooth_reject,
                                      true_null = true_null)
      
      results_mat[base_indx + k, ] <- c(i,
                                        base_metrics,
                                        smooth_metrics)
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
                            paste0("smooth_", metric_names))
  
  return(results_frame)
}

run_sim_par_2loc <- function(iter,
                              lat_lon_grid,
                              szon_dist,
                              smer_dist,
                              loc_list,
                              cov_mat,
                              beta_vec,
                              n,
                              cov_L,
                              noise_sd,
                              verbose = FALSE) {
  
  p <- nrow(lat_lon_grid)
  
  results_mat <- matrix(0, nrow = iter * 1000, ncol = 13)
  
  init_cond_mat <- matrix(c(1, 1, 1,
                            0.2 * pi, 0.5 * pi, 0.8 * pi,
                            1.5, 1.5, 1.5),
                          nrow = 3)
  
  rnd_indices <- replicate(iter, 
    replicate(1000,
         sample(p, size = 2, replace = FALSE)))
  
  rnd_indices <- aperm(rnd_indices)
  
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
    
    base_indx <- (i - 1) * 1000
    
    for (k in 1:1000) {
      if (verbose && k %% 100 == 0) print(k)
      
      index_pair <- rnd_indices[iter, k, ]

      true_beta_j <- beta_vec[1] * cov_mat[index_pair[1], ] + 
        beta_vec[2] * cov_mat[index_pair[2], ]

      true_null <- (true_beta_j == 0)
      
      y <- beta_vec[1] * sim_data[, index_pair[1]] + 
        beta_vec[2] * sim_data[, index_pair[2]] +
        rnorm(n, sd = noise_sd)
      
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
      
      base_metrics <- estim_metrics(p_dt = base_reject,
                                    true_null = true_null)
      smooth_metrics <- estim_metrics(p_dt = smooth_reject,
                                      true_null = true_null)
      
      results_mat[base_indx + k, ] <- c(i,
                                        base_metrics,
                                        smooth_metrics)
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
                            paste0("smooth_", metric_names))
  
  return(results_frame)
}


run_sim_par_doubledip <- function(iter,
                                  lat_lon_grid,
                                  szon_dist,
                                  smer_dist,
                                  loc_list,
                                  cov_mat,
                                  beta,
                                  n,
                                  cov_L,
                                  noise_sd,
                                  save_dir,
                                  verbose = FALSE) {
            
  p <- nrow(lat_lon_grid)
  
  results_mat <- matrix(0, nrow = iter * 1000, ncol = 19)
  
  init_cond_mat <- matrix(c(1, 1, 1,
                            0.2 * pi, 0.5 * pi, 0.8 * pi,
                            1.5, 1.5, 1.5),
                          nrow = 3)
  
  rnd_indx <- c(replicate(iter, sample(p, size = 1000, replace = FALSE)))
  
  loc_pairs_list <- lapply(loc_list,
                           \(x) find_dist_pairs(zonal_dist_mat = szon_dist,
                                            merid_dist_mat = smer_dist,
                                            loc_vec = x))
  
  
  for (i in 1:iter) {
    if (verbose) print(i)
    
    sim_data <- gen_data(n = n, L = cov_L)
    sim_data_std <- scale(sim_data)

    sim_data2 <- gen_data(n = n, L = cov_L)
    sim_data_std2 <- scale(sim_data2)
    
    MLE_out <- future_mapply(find_MLE_mapply_wrapper,
                             loc_block = loc_list,
                             pairs_mat = loc_pairs_list,
                             MoreArgs = list(data_mat = sim_data_std,
                                             lik_func = m1.5_lik),
                             SIMPLIFY = FALSE)
    
    MLE_df <- rbindlist(MLE_out)
    MLE_mat <- data.matrix(MLE_df[, 1:3])

    MLE_out2 <- future_mapply(find_MLE_mapply_wrapper,
                          loc_block = loc_list,
                          pairs_mat = loc_pairs_list,
                          MoreArgs = list(data_mat = sim_data_std2,
                                          lik_func = m1.5_lik),
                          SIMPLIFY = FALSE)
    
    MLE_df2 <- rbindlist(MLE_out2)
    MLE_mat2 <- data.matrix(MLE_df2[, 1:3])
    
    cov_weights <- find_smooth_weights(loc_list = loc_list,
                                       zonal_sdist = szon_dist,
                                       merid_sdist = smer_dist,
                                       cov_func = matern_1.5_cov,
                                       mle_params = MLE_mat,
                                       progress = FALSE)
    
    cov_weights2 <- find_smooth_weights(loc_list = loc_list,
                                    zonal_sdist = szon_dist,
                                    merid_sdist = smer_dist,
                                    cov_func = matern_1.5_cov,
                                    mle_params = MLE_mat2,
                                    progress = FALSE)
    
    smooth_data <- sim_data_std %*% cov_weights
    smooth_data2 <- sim_data_std %*% cov_weights2
    
    base_indx <- (i - 1) * 1000
    
    for (k in 1:1000) {
      true_cov <- cov_mat[rnd_indx[base_indx + k], ]
      true_null <- (true_cov == 0)
      
      y <- beta * sim_data[, rnd_indx[base_indx + k]] + rnorm(n, sd = noise_sd)
      
      base_lm_fits <- fast_fit_lm_tstat(y = y, data_mat = sim_data_std)
      smooth_lm_fits <- fast_fit_lm_tstat(y = y, data_mat = smooth_data)
      smooth_lm_fits2 <- fast_fit_lm_tstat(y = y, data_mat = smooth_data2)
      
      # two-tailed test
      base_pvals <- 2 * pt(abs(base_lm_fits),
                           df = n - 2,
                           lower.tail = FALSE)
      smooth_pvals <- 2 * pt(abs(smooth_lm_fits),
                             df = n - 2,
                             lower.tail = FALSE)
      smooth_pvals2 <- 2 * pt(abs(smooth_lm_fits2),
                             df = n - 2,
                             lower.tail = FALSE)
      
      base_reject <- find_BH_reject(base_pvals, alpha = 0.1)
      smooth_reject <- find_BH_reject(smooth_pvals, alpha = 0.1)
      smooth_reject2 <- find_BH_reject(smooth_pvals2, alpha = 0.1)
      
      base_metrics <- estim_metrics(p_dt = base_reject,
                                    true_null = true_null)
      smooth_metrics <- estim_metrics(p_dt = smooth_reject,
                                      true_null = true_null)
      
      smooth_metrics2 <- estim_metrics(p_dt = smooth_reject2,
                                      true_null = true_null)
      
      base_reject$method <- "base"
      base_reject$true_null <- true_null

      smooth_reject$method <- "smooth"
      smooth_reject$true_null <- true_null

      smooth_reject2$method <- "iso_smooth"
      smooth_reject2$true_null <- true_null

      reject_total <- rbind(base_reject, smooth_reject, smooth_reject2)

      fname <- paste0(save_dir, "field", i, "_iter", k, "_pvals.csv")
      fwrite(reject_total, file = fname)
      
      results_mat[base_indx + k, ] <- c(i,
                                        base_metrics,
                                        smooth_metrics,
                                        smooth_metrics2)
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
                            paste0("smooth_iso_", metric_names))
  
  return(results_frame)
}
