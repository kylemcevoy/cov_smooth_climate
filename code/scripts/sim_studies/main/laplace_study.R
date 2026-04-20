source("code/functions/blockwise_MLE.R")
source("code/functions/cov_smooth.R")
source("code/functions/dist_cov.R")
source("code/functions/inference.R")
source("code/functions/sim.R")
source("code/functions/parallel.R")

library(future.apply)
library(parallelly)
library(extraDistr)
library(data.table)

# necessary to prevent overhead of spinning up sandboxes for parallel sessions 
Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = FALSE)

###### Functions ######

run_sim_laplace <- function(iter,
                            lat_lon_grid,
                            szon_dist,
                            smer_dist,
                            loc_list,
                            cov_mat,
                            beta,
                            n,
                            cov_L,
                            noise_scale,
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
      true_null <- (true_cov == 0)
      
      y <- sim_data[, rnd_indx[base_indx + k]] + rlaplace(n, sigma = noise_scale)
      
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

##### Sim #####

save_dir <- "/home/data/projects/clim_smooth/sim_studies/two_sided/"

# check number of cores before running
plan(multicore, workers = 24)
# use plan sequential if you don't want to run the code in parallel.
# plan(sequential)

# selected via rolling 6 10-sided dice.
set.seed(695429)

# create lat, lon grid on which to simulate data. Choose lat. and lon. ranges
# and choose units to create resolution of grid.
lat_lon_grid <- build_lat_lon_grid(lat_range = c(-60, 60),
                                   lon_range = c(0, 358),
                                   lat_unit = 2,
                                   lon_unit = 2)

# approximate distances are found by calculating an EW distance and a NS distance
# between each gridbox. See paper/function code for details. The same distance
# is calculated in anisotropic and isotropic cases for consistency.

szon_dist <- find_signed_zonal_dist(lat_lon_grid)
smer_dist <- find_signed_merid_dist(lat_lon_grid)

#isotropic distance measure
sim_dist <- sqrt(szon_dist^2 + smer_dist^2)

taper_cov <- matern_1.5_wendland_1_cov(sim_dist, phi = 4, taper = 7)
chol_cov <- t(chol(taper_cov))

sim_data <- gen_data(50, L=chol_cov)
sim_data_std <- scale(sim_data)

loc_list <- get_loc_list(lat_lon_grid)

init_cond_mat <- matrix(c(1, 1, 1,
                          0.2 * pi, 0.5 * pi, 0.8 * pi,
                          1.5, 1.5, 1.5),
                        nrow = 3)

sim_laplace <- run_sim_laplace(iter = 10,
                               lat_lon_grid = lat_lon_grid,
                               szon_dist = szon_dist,
                               smer_dist = smer_dist,
                               loc_list = loc_list,
                               cov_mat = taper_cov,
                               beta = 1,
                               n = 50,
                               cov_L = chol_cov,
                               noise_scale = 1,
                               verbose = TRUE)

sim_laplace$sigma <- 1
sim_laplace$n <- 50
sim_laplace$phi <- 4
sim_laplace$theta <- 0
sim_laplace$psi <- 1
sim_laplace$taper <- 7

plan(sequential)

sim_laplace_path <- paste0(save_dir, "sim_n50_scale1_phi4_taper7_laplace.csv")
fwrite(sim_laplace, file = sim_laplace_path)

sim_laplace$base_FDX <- sim_laplace$base_FDP > 0.1
sim_laplace$smooth_FDX <- sim_laplace$smooth_FDP > 0.1

results_mean <- sim_laplace[, lapply(.SD, mean), by = .(sigma,
                                                        n,
                                                        phi,
                                                        theta,
                                                        psi,
                                                        taper)]

results_mean[, iteration := NULL]

results_mean <- melt(results_mean, 
                     id.vars = c("sigma", "n", "phi", "theta", "psi", "taper"),
                     measure.vars = measure(method,
                                            value.name,
                                            pattern = "^([[:alpha:]]*)_(.*)$"))

results_mean <- results_mean[order(phi)]

fwrite(results_mean, file = paste0(save_dir, "laplace_sim_study_mean.csv"))
