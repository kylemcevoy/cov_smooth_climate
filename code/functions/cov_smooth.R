# MLE-based covariance weighted smoothing.

#Requires functions from the following two files to run.
#source("R/dist_cov_func.R")
#source("R/blockwise_MLE.R")

#' Finding smoothing weights from MLE of covariance function parameters
#'
#' Uses the locally fitted MLEs to calculate smoothing weights for local
#' data smoothing.
#'
#' @param loc_list list where each element is a vector containing the
#' locations for the local gridboxes.
#' @param zonal_sdist matrix of signed zonal distances between pairs.
#' @param merid_sdist matrix of signed meridional distances between pairs.
#' @param cov_func either matern_0.5_cov or matern_1.5_cov (default)
#' @param mle_params p x 3 matrix containing the locally fitted MLE of the
#' parameters. First column contains phi parameter, second contains anisotropy
#' angle, third contains anisotropy ratio.
#' @param progress if TRUE reports progress every 20 locations.
#' @returns matrix of smoothing weights by location. To output smoothed data,
#' take standardized anomaly data and multiply by smoothing weights:
#' data_mat_std %*% weight_mat.

find_smooth_weights <- function(loc_list,
                                zonal_sdist,
                                merid_sdist,
                                cov_func = matern_1.5_cov,
                                mle_params,
                                progress = TRUE) {

  p <- length(loc_list)
  weight_mat <- matrix(0, nrow = p, ncol = p)
  
  for (i in seq_along(loc_list)) {
    if (progress == TRUE && (i %% 1000 == 0)) print(i)
    
    mle_param_vec <- mle_params[i, ]
    loc_vec <- loc_list[[i]]
    m <- length(loc_vec)

    dist_pairs <- find_local_dist_pairs(zonal_dist_mat = zonal_sdist,
                                        merid_dist_mat = merid_sdist,
                                        loc = i,
                                        loc_vec = loc_vec)

    aniso_mat <- find_aniso_mat(aniso_angle = mle_param_vec[2],
                                aniso_ratio = mle_param_vec[3])
    
    #calculates Euclidean distance in anisotropic space.
    aniso_dist_vec <- sqrt(rowSums((dist_pairs %*% aniso_mat) * dist_pairs))

    cov_weights <- cov_func(dist = aniso_dist_vec,
                            sigma2 = 1,
                            phi = mle_param_vec[1])

    cov_weights <- cov_weights / sum(cov_weights)
    weight_mat[loc_vec, i] <- cov_weights
  }
  return(weight_mat)
}
