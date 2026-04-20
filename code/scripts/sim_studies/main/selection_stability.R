####### Selection Stability #######
library(data.table)
library(viridis)
library(colorspace)
library(ggplot2)
library(future.apply)
library(parallelly)

source("code/functions/blockwise_MLE.R")
source("code/functions/cov_smooth.R")
source("code/functions/dist_cov.R")
source("code/functions/inference.R")
source("code/functions/sim.R")
source("code/functions/parallel.R")

####### Set-up ########
# This system setting reduces overhead in creating new parallel workers.
Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = FALSE)

save_path <- '/home/data/projects/clim_smooth/sim_studies/two_sided/'

avail_cores <- availableCores()
avail_cores
# set number of workers based on available cores.
plan(multicore, workers = 16)
# use plan sequential if you don't want to run the code in parallel.
# plan(sequential)

# Seed selected via rolling 6 10-sided dice.
set.seed(704602)

n <- 50

# create lat, lon grid on which to simulate data. Choose lat. and lon. ranges
# and choose units to create resolution of grid.
lat_lon_grid <- build_lat_lon_grid(lat_range = c(-60, 60),
                                   lon_range = c(0, 358),
                                   lat_unit = 2,
                                   lon_unit = 2)

p <- nrow(lat_lon_grid)

# approximate distances are found by calculating an EW distance and a NS distance
# between each gridbox. See paper/function code for details. The same distance
# is calculated in anisotropic and isotropic cases for consistency.

szon_dist <- find_signed_zonal_dist(lat_lon_grid)
smer_dist <- find_signed_merid_dist(lat_lon_grid)

#isotropic distance measure
sim_dist <- sqrt(szon_dist^2 + smer_dist^2)

# matern-3/2 covariance with wendland-1 polynomial tapering. 
taper_cov <- matern_1.5_wendland_1_cov(sim_dist, phi = 4, taper = 7)
# transpose gives left cholesky factor
chol_cov <- t(chol(taper_cov))

loc_list <- get_loc_list(lat_lon_grid)

loc_pairs_list <- lapply(loc_list,
                         \(x) find_dist_pairs(zonal_dist_mat = szon_dist,
                                              merid_dist_mat = smer_dist,
                                              loc_vec = x))
init_cond_mat <- matrix(c(1, 1, 1,
                          0.2 * pi, 0.5 * pi, 0.8 * pi,
                          1.5, 1.5, 1.5),
                        nrow = 3)

sim_data <- gen_data(n = n, L = chol_cov)

sim_data_std <- scale(sim_data)

plan(multisession, workers = 8)
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

base_reject_loc <- matrix(NA, nrow = 1000, ncol = 10980)
smooth_reject_loc <- matrix(NA, nrow = 1000, ncol = 10980)

# picked a central location for plotting purposes.
center_loc <- 5475

for (i in 1:1000) {
  if (i %% 100 == 0) message(i)
  
  y <- sim_data[, center_loc] + rnorm(n = 50)
  base_tstat <- fast_fit_lm_tstat(y = y, data_mat = sim_data_std)
  smooth_tstat <- fast_fit_lm_tstat(y = y, data_mat = smooth_data)
  
  base_pvals <- 2 * pt(abs(base_tstat),
                       df = n - 2,
                       lower.tail = FALSE)
  smooth_pvals <- 2 * pt(abs(smooth_tstat),
                         df = n - 2,
                         lower.tail = FALSE)
  
  BH_base <- find_BH_reject(base_pvals, alpha = 0.1)
  BH_smooth <- find_BH_reject(smooth_pvals, alpha = 0.1)
  
  base_reject_loc[i, ] <- BH_base$BH_reject
  smooth_reject_loc[i, ] <- BH_smooth$BH_reject
}

base_reject_freq <- apply(base_reject_loc, 2, mean)
smooth_reject_freq <- apply(smooth_reject_loc, 2, mean)
reject_freq_diff <- smooth_reject_freq - base_reject_freq

true_cov <- taper_cov[center_loc, ]

selection_stability_dt <- data.table(location = center_loc,
                                     true_cov = true_cov,
                                     reject_freq_diff = reject_freq_diff)


fwrite(selection_stability_dt, 
       file = paste0(save_path, 'selection_stability.csv'))
