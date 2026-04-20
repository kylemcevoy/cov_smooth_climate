# Rewritten functions for running the simulations studies. Including functions
# to build the synthetic spatial grid, draw from multivariate normal 
# distributions, and perform the simulation studies. Parallel versions of
# the simulation study functions are in parallel.R.

library(data.table)

#source("R/functions/dist_cov.R")
#source("R/functions/blockwise_MLE.R")

#' Building latitude and longitude coordinate grids
#'
#' Builds a data.table of latitude and longitude coordinate pairs for use in
#' simulating data.
#'
#' @param lat_range vector whose elements give the lower and upper bounds of
#' latitudes for the grid.
#' @param lon_range vector whose elements give the lower and upper bounds of
#' longitudes for the grid.
#' @param lat_unit gives the step size for the latitudes in the grid.
#' @param lon_unit gives the step size for the longitudes in the grid.
#' @returns a data.table whose rows specify each gridboxes' (lon, lat)
#' coordinates, as well as a location index.

build_lat_lon_grid <- function(lat_range = c(-60, 60),
                               lon_range = c(0, 358),
                               lat_unit = 2,
                               lon_unit = 2) {

  lat_vec <- seq(lat_range[1], lat_range[2], by = lat_unit)
  lon_vec <- seq(lon_range[1], lon_range[2], by = lon_unit)

  coord_dt <- data.table(expand.grid(lon = lon_vec, lat = lat_vec))
  coord_dt$location <- seq(from = 1, to = nrow(coord_dt), by = 1)

  return(coord_dt)
}

#' Building Matern covariance matrices
#'
#' The function finds isotropic Matern covariance matrices using the given
#' distances.
#'
#' @param dist_mat p x p matrix of distances between the p different locations.
#' Should be indexed by location as found in build_lat_lon_grid() or
#' extract_coord_dt().
#' @param phi the range parameter of the matern covariance
#' @param sigma2 the scaling factor for the variance.
#' @param shape parameter of the matern covariance. Can take either value 0.5
#' for matern-1/2 or 1.5 for matern-3/2
#' @returns a p x p matrix (same dimensions as dist_mat) with the specified
#' matern covariance function at each distance.

build_cov_mat <- function(dist_mat,
                          phi = 1,
                          sigma2 = 1,
                          shape = 1.5) {
  
  if (shape == 0.5) {
    cov_function <- matern_0.5_cov
  } else if (shape == 1.5) {
    cov_function <- matern_1.5_cov
  } else {
    stop("shape should be either 0.5 or 1.5")
  }

  cov_mat <- cov_function(dist_mat = dist_mat,
                          phi = phi,
                          sigma2 = sigma2)
  
  return(cov_mat)
}

#' Building Anisotropic Matern covariance matrices
#'
#' The function finds anisotropic Matern covariance matrices using the given
#' distances.
#'
#' @param zonal_sdist_mat p x p matrix of signed zonal (E-W) distances between
#' the p different locations.
#' Should be indexed by location as found in build_lat_lon_grid() or
#' extract_coord_dt().
#' @param merid_sdist_mat p x p matrix of signed meridional (N-S) distances
#'  between the p different locations.
#' Should be indexed by location as found in build_lat_lon_grid() or
#' extract_coord_dt().
#' @param phi the range parameter of the matern covariance
#' @param aniso_angle the anisotropy angle. See paper for details on the exact
#' specification.
#' @param aniso_ratio the anisotropy ratio (see paper for details).
#' @param sigma2 the scaling factor for the variance.
#' @param tau2 the nugget variance
#' @param shape parameter of the matern covariance. Can take either value 0.5
#' for matern-1/2 or 1.5 for matern-3/2
#' @returns a p x p matrix (same dimensions as dist_mat) with the specified
#' matern covariance function at each distance.

build_aniso_cov_mat <- function(zonal_sdist_mat,
                                merid_sdist_mat,
                                phi = 1,
                                aniso_angle = pi / 2,
                                aniso_ratio = 1,
                                sigma2 = 1,
                                shape = 1.5) {

  if (shape == 0.5) {
    cov_function <- matern_0.5_cov
  } else if (shape == 1.5) {
    cov_function <- matern_1.5_cov
  } else {
    stop("shape should be either 0.5 or 1.5")
  }

  aniso_mat <- find_aniso_matrix(aniso_angle = aniso_angle,
                                 aniso_ratio = aniso_ratio)

  # unravels the signed distance matrices into 2 vectors, then combines them
  # into a p^2 x 2 matrix that can be used with aniso_mat to find the distances
  # in the anisotropic space.
  zonal_sdist_vec <- c(zonal_sdist_mat)
  merid_sdist_vec <- c(merid_sdist_mat)
  dist_pairs <- cbind(zonal_sdist_vec, merid_sdist_vec)

  wdist <- sqrt(rowSums((dist_pairs %*% aniso_mat) * dist_pairs))
  wdist_mat <- matrix(wdist, nrow = p)

  cov_mat <- cov_function(dist_mat = wdist_mat,
                          phi = phi,
                          sigma2 = sigma2)
  
  return(cov_mat)
}


#' Finding square roots of covariance matrices
#'
#' @description
#'  This function is adapted from MASS::mvrnorm. It takes a covariance matrix as
#'  input and outputs a square root of that matrix. The particular square root
#'  is non-symmetric. It is formed by taking an eigendecomposition of the
#'  covariance matrix. From that decomposition we take the matrix of eigen
#'  vectors, $Q$, and multiply it with a diagonal matrix, $\Lambda^{1/2}$ whose
#'  diagonal elements are the square roots of the eigenvalues corresponding
#'  to each eigenvector. Generally unnecessary, use the Cholesky decomposition
#' (chol()) instead. However, this function can handle small negative
#' eigenvalues that can arise from imprecision.
#'
#' @param cov_mat covariance matrix -- needs to be a symmetric positive
#' semidefinite matrix to be a valid covariance matrix.
#' @param tol due to numerical instability in the eigendecomposition sometimes
#' there are small negative eigenvalues, even though the matrix should be positive
#' semidefinite. These will be replaced by 0s in the final matrix product, but
#' the function will error if the negative values are too large relative to the
#' tolerance and the largest eigenvalue.
#' @returns The square root matrix formed by the matrix of eigenvectors right
#' multiplied by square root eigenvalues matrix, $Q \Lambda^{1/2}$.

sqrt_eigen_cov <- function(cov_mat, tol = 1e-06) {
  p <- nrow(cov_mat)
  
  if (p != ncol(cov_mat)) {
    stop("covariance must be square")
  }

  eigen_decomp <- eigen(cov_mat, symmetric = TRUE)

  eigen_vec <- eigen_decomp$vectors
  eigen_vals <- eigen_decomp$values

  if (!all(eigen_vals >= -tol * abs(eigen_vals[1]))) {
    stop ("Covariance is not psd.")
  } 

  return(eigen_vec %*% diag(sqrt(pmax(eigen_vals, 0)), nrow = p))
}

#' Generate multivariate normal data
#'
#' Adapted from MASS::mvrnorm.
#' Uses either the Cholesky decomposition or the square root created by 
#' sqrt_eigen_cov() to generate multivariate normal data. The two functions are
#' adapted from rmvnorm, so that the square root/cholesky decomposition of the 
#' covariance matrix can be saved and reused to generate new data without having
#' to perform a costly matrix factorization every time.
#' @param n number of independent samples to generate from the desired
#' multivariate normal distribution.
#' @param L If Sigma is the target covariance matrix, then L must satsify:
#' L %*% t(L) = Sigma. Can be generated by taking the Cholesky decomposition
#' of Sigma with chol(). If this produces errors, the function sqrt_eigen_cov()
#' decomposes using eigen decomposition and sets small negative eigenvalues to
#' 0. This is much slower but a little more robust.
#' @param seed if a random seed is desired can be set using this argument.
#' @returns an n by p matrix whose rows are independent samples from a
#' p-dimensional multivariate normal distribution with mean 0 and covariance
#' matrix given by L %*% t(L).

gen_data <- function(n, L, seed = NULL) {
  if(!is.null(seed)) {
    set.seed(seed)
  }

  p <- ncol(L)
  # Z^T in usual dimensions.
  Z <- matrix(rnorm(p * n), nrow = p, ncol = n)
  X <- L %*% Z
  # transpose back to expected n by p.
  return(t(X))
}

#' Run simulation study
#' 
#' function to run a specific simulation study and calculate the output metrics.
#' 
#' @param iter the total number of times to generate the spatial field. In the
#' paper the iter is set to 10.
#' @param lat_lon_grid data.table containing the latitude and longitude 
#' coordinates of the locations. See build_lat_lon_grid()
#' @param t_ecdf the simulated t statistics under null.
#' @param szon_dist signed zonal distances
#' @param smer_dist signed meridional distances
#' @param loc_list list containing the nearby locations for each location.
#' @param cov_mat covariance matrix used to generate the simulated data.
#' @param beta the value of beta used to generate the response y.
#' @param delta the upper bound of the interval null
#' @param n the sample size
#' @param sqrt_cov the square root of the cov_mat calculated using eigen decomp.
#' cached for quicker multivariate normal data generation.
#' @param noise_sd the standard deviation of the normal noise used to generate y
#' @param verbose should the function report progress?
#' @returns returns a (1000 * iter) by 13 data.table containing the iteration
#' number as well as the number of rejections, number of true rejections, FDP,
#' FNP, sensitivity, and specificity for each simulation.

run_sim <- function(iter,
                    lat_lon_grid,
                    szon_dist,
                    smer_dist,
                    loc_list,
                    cov_mat,
                    beta,
                    n,
                    cov_L,
                    noise_sd,
                    verbose = FALSE) {
  
  p <- nrow(lat_lon_grid)
  
  results_mat <- matrix(NA, nrow = iter * 1000, ncol = 13)
  
  init_cond_mat <- matrix(c(1, 1, 1,
                            0.2 * pi, 0.5 * pi, 0.8 * pi,
                            1.5, 1.5, 1.5),
                          nrow = 3)
  
  MLE_df <- data.frame(phi = rep(0, p),
                       ani_angle = rep(0, p),
                       ani_ratio = rep(0, p),
                       convergence = rep(9, p))
  
  rnd_indx <- c(replicate(iter, sample(p, size = 1000, replace = FALSE)))
  
  for (i in 1:iter) {
    if (verbose) print(i)
    
    sim_data <- gen_data(n = n, L = cov_L)
    sim_data_std <- scale(sim_data)
    
    for (j in 1:p) {
      MLE_df[j, ] <- find_local_MLE(loc = j,
                                    coord_dt = lat_lon_grid,
                                    zonal_dist_mat = szon_dist,
                                    merid_dist_mat = smer_dist,
                                    data_mat_std = sim_data_std,
                                    lik_func = m1.5_lik,
                                    return_best = TRUE,
                                    init_cond = init_cond_mat,
                                    return_likelihood = FALSE)
      
      if (verbose && j %% 1000 == 0) print(j)
    }
    MLE_mat <- data.matrix(MLE_df[, 1:3])
    
    cov_weights <- find_smooth_weights(loc_list = loc_list,
                                       zonal_sdist = szon_dist,
                                       merid_sdist = smer_dist,
                                       cov_func = matern_1.5_cov,
                                       mle_params = MLE_mat,
                                       progress = FALSE)
    
    smooth_data <- sim_data_std %*% cov_weights
    
    for (k in 1:1000) {
      if (verbose && k %% 100 == 0) print(k)
      
      true_cov <- cov_mat[rnd_indx[(i - 1) * 1000 + k], ]
      true_null <- (true_cov * beta) == 0 
      
      y <- beta * sim_data[, rnd_indx[(i - 1) * 1000 + k]] + 
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

run_sim_FDRL <- function(iter,
                    lat_lon_grid,
                    szon_dist,
                    smer_dist,
                    loc_list,
                    cov_mat,
                    beta,
                    n,
                    cov_L,
                    noise_sd,
                    verbose = FALSE) {
  
  p <- nrow(lat_lon_grid)
  
  results_mat <- matrix(NA, nrow = iter * 1000, ncol = 7)
  
  rnd_indx <- c(replicate(iter, sample(p, size = 1000, replace = FALSE)))
  
  for (i in 1:iter) {
    if (verbose) print(i)
    
    sim_data <- gen_data(n = n, L = cov_L)
    sim_data_std <- scale(sim_data)
    
    base_indx <- (i - 1) * 1000

    for (k in 1:1000) {
      if (verbose && k %% 100 == 0) print(k)
      
      true_cov <-  cov_mat[rnd_indx[base_indx + k], ]
      true_null <- (true_cov == 0)
      
      y <- beta * sim_data[, rnd_indx[base_indx + k]] + rnorm(n, sd = noise_sd)
      
      base_lm_fits <- fast_fit_lm_tstat(y = y, data_mat = sim_data_std)

      # two-tailed test
      base_pvals <- 2 * pt(abs(base_lm_fits),
                           df = n - 2,
                           lower.tail = FALSE)
      
      p_star <- find_pstar(base_pvals, loc_list = loc_list)
      
      t_alpha <- find_FDRL_threshold(p_star = p_star,
                                     lambda = 0.1,
                                     alpha = 0.1)
      
      FDRL_reject <- p_star < t_alpha
      
      FDRL_reject_dt <- data.table(location = lat_lon_grid$location,
                                   p_val = p_star,
                                   BH_reject = FDRL_reject)
      
      fdrl_metrics <- estim_metrics(p_dt = FDRL_reject_dt,
                                    true_null = true_null)
      
      results_mat[(i - 1) * 1000 + k, ] <- c(i,
                                             fdrl_metrics)
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
                            paste0("fdrl_", metric_names))
  
  return(results_frame)
}


# For mean of p-values

find_pmean <- function(p_vec, loc_list) {
  p_mean <- sapply(loc_list, \(x) mean(p_vec[x]))
  return(p_mean)
}

run_sim_pmean <- function(iter,
                        lat_lon_grid,
                        szon_dist,
                        smer_dist,
                        loc_list,
                        cov_mat,
                        beta,
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
  
  rnd_indx <- c(replicate(iter, sample(p, size = 1000, replace = FALSE)))
  
  loc_pairs_list <- lapply(loc_list,
                           \(x) find_dist_pairs(zonal_dist_mat = szon_dist,
                                            merid_dist_mat = smer_dist,
                                            loc_vec = x))
  
  
  for (i in 1:iter) {
    if (verbose) print(i)
    
    sim_data <- gen_data(n = n, L = cov_L)
    sim_data_std <- scale(sim_data)

    base_indx <- (i - 1) * 1000
    
    for (k in 1:1000) {
      if (verbose && k %% 100 == 0) print(k)
      
      true_cov <- cov_mat[rnd_indx[base_indx + k], ]
      true_null <- (true_cov == 0)
      
      y <- beta * sim_data[, rnd_indx[base_indx + k]] + rnorm(n, sd = noise_sd)
      
      base_lm_fits <- fast_fit_lm_tstat(y = y, data_mat = sim_data_std)
      
      # two-tailed test
      base_pvals <- 2 * pt(abs(base_lm_fits),
                           df = n - 2,
                           lower.tail = FALSE)

      mean_pvals <- find_pmean(base_pvals, loc_list)
      
      base_reject <- find_BH_reject(base_pvals, alpha = 0.1)
      mean_reject <- find_BH_reject(mean_pvals, alpha = 0.1)
      
      base_metrics <- estim_metrics(p_dt = base_reject,
                                    true_null = true_null)
      mean_metrics <- estim_metrics(p_dt = mean_reject,
                                      true_null = true_null)
      
      results_mat[base_indx + k, ] <- c(i,
                                        base_metrics,
                                        mean_metrics)
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
                            paste0("mean_", metric_names))
  
  return(results_frame)
}

