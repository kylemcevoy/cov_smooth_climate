# run FDR_L sims implementing the method of Zheng et al. (2011) See FDRL.R 
# for more details

source("code/functions/blockwise_MLE.R")
source("code/functions/cov_smooth.R")
source("code/functions/dist_cov.R")
source("code/functions/inference.R")
source("code/functions/sim.R")
source("code/functions/FDRL.R")

library(data.table)

save_dir <- "/home/data/projects/clim_smooth/sim_studies/two_sided/"

#same seed used in noise sim study.
set.seed(696145)

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

sim_data <- gen_data(50, L = chol_cov)
sim_data_std <- scale(sim_data)

loc_list <- get_loc_list(lat_lon_grid)
loc_list2 <- get_loc_list2(lat_lon_grid)

# fdrl_out1 <- run_sim_FDRL(iter = 10,
#                              lat_lon_grid = lat_lon_grid,
#                              szon_dist = szon_dist,
#                              smer_dist = smer_dist,
#                              loc_list = loc_list,
#                              cov_mat = taper_cov,
#                              cov_L = chol_cov,
#                              n = 50,
#                              beta = 1,
#                              noise_sd = 1,
#                              verbose = TRUE)

# fdrl_out1$sigma <- 1
# fdrl_out1$n <- 50
# fdrl_out1$phi <- 4
# fdrl_out1$theta <- 0
# fdrl_out1$psi <- 1
# fdrl_out1$taper <- 7

# fwrite(fdrl_out1, file = paste0(save_dir, "fdrl_sim_big_sd1.csv"))

# fdrl_out2 <- run_sim_FDRL(iter = 10,
#                              lat_lon_grid = lat_lon_grid,
#                              szon_dist = szon_dist,
#                              smer_dist = smer_dist,
#                              loc_list = loc_list,
#                              cov_mat = taper_cov,
#                              cov_L = chol_cov,
#                              n = 50,
#                              beta = 1,
#                              noise_sd = 2,
#                              verbose = TRUE)

# fdrl_out2$sigma <- 2
# fdrl_out2$n <- 50
# fdrl_out2$phi <- 4
# fdrl_out2$theta <- 0
# fdrl_out2$psi <- 1
# fdrl_out2$taper <- 7

# fwrite(fdrl_out2, file = paste0(save_dir, "fdrl_sim_big_sd2.csv"))

# fdrl_out3 <- run_sim_FDRL(iter = 10,
#                              lat_lon_grid = lat_lon_grid,
#                              szon_dist = szon_dist,
#                              smer_dist = smer_dist,
#                              loc_list = loc_list,
#                              cov_mat = taper_cov,
#                              cov_L = chol_cov,
#                              n = 50,
#                              beta = 1,
#                              noise_sd = 3,
#                              verbose = TRUE)

# fwrite(fdrl_out3, file = paste0(save_dir, "fdrl_sim_big_sd3.csv"))

# fdrl_out3$sigma <- 3
# fdrl_out3$n <- 50
# fdrl_out3$phi <- 4
# fdrl_out3$theta <- 0
# fdrl_out3$psi <- 1
# fdrl_out3$taper <- 7

# ### small grid

# fdrl_small_out1 <- run_sim_FDRL(iter = 10,
#                              lat_lon_grid = lat_lon_grid,
#                              szon_dist = szon_dist,
#                              smer_dist = smer_dist,
#                              loc_list = loc_list2,
#                              cov_mat = taper_cov,
#                              cov_L = chol_cov,
#                              n = 50,
#                              beta = 1,
#                              noise_sd = 1,
#                              verbose = TRUE)

# fdrl_small_out1$sigma <- 1
# fdrl_small_out1$n <- 50
# fdrl_small_out1$phi <- 4
# fdrl_small_out1$theta <- 0
# fdrl_small_out1$psi <- 1
# fdrl_small_out1$taper <- 7

# fwrite(fdrl_small_out1, file = paste0(save_dir, "fdrl_sim_small_sd1.csv"))

# fdrl_small_out2 <- run_sim_FDRL(iter = 10,
#                              lat_lon_grid = lat_lon_grid,
#                              szon_dist = szon_dist,
#                              smer_dist = smer_dist,
#                              loc_list = loc_list2,
#                              cov_mat = taper_cov,
#                              cov_L = chol_cov,
#                              n = 50,
#                              beta = 1,
#                              noise_sd = 2,
#                              verbose = TRUE)

# fdrl_small_out2$sigma <- 1
# fdrl_small_out2$n <- 50
# fdrl_small_out2$phi <- 4
# fdrl_small_out2$theta <- 0
# fdrl_small_out2$psi <- 1
# fdrl_small_out2$taper <- 7

# fwrite(fdrl_small_out2, file = paste0(save_dir, "fdrl_sim_small_sd2.csv"))

fdrl_small_out3 <- run_sim_FDRL(iter = 10,
                             lat_lon_grid = lat_lon_grid,
                             szon_dist = szon_dist,
                             smer_dist = smer_dist,
                             loc_list = loc_list2,
                             cov_mat = taper_cov,
                             cov_L = chol_cov,
                             n = 50,
                             beta = 1,
                             noise_sd = 3,
                             verbose = TRUE)

fdrl_small_out3$sigma <- 3
fdrl_small_out3$n <- 50
fdrl_small_out3$phi <- 4
fdrl_small_out3$theta <- 0
fdrl_small_out3$psi <- 1
fdrl_small_out3$taper <- 7
fdrl_small_out3$window <- 4

fwrite(fdrl_small_out3, file = paste0(save_dir, "fdrl_sim_small_sd3.csv"))
