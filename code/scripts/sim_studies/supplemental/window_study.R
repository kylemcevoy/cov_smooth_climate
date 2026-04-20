source("code/functions/blockwise_MLE.R")
source("code/functions/cov_smooth.R")
source("code/functions/dist_cov.R")
source("code/functions/inference.R")
source("code/functions/sim.R")
source("code/functions/parallel.R")

library(data.table)
library(future.apply)
library(parallelly)

save_path <- "/home/data/projects/clim_smooth/sim_studies/two_sided/"

# set number of workers based on available cores.
plan(multicore, workers = 16)
# use plan sequential if you don't want to run the code in parallel.
# plan(sequential)

# selected via rolling 6 10-sided dice.
set.seed(760316)

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

# Using 16 deg. by 16 deg. window instead of 20 by 20.
loc_list <- get_loc_list(lat_lon_grid, lat_margin = 8, lon_margin = 8)

init_cond_mat <- matrix(c(1, 1, 1,
                          0.2 * pi, 0.5 * pi, 0.8 * pi,
                          1.5, 1.5, 1.5),
                        nrow = 3)

sim_window <- run_sim_par(iter = 10,
                          lat_lon_grid = lat_lon_grid,
                          szon_dist = szon_dist,
                          smer_dist = smer_dist,
                          loc_list = loc_list,
                          cov_mat = taper_cov,
                          beta = 1,
                          n = 50,
                          cov_L = chol_cov,
                          noise_sd = 1,
                          verbose = TRUE)

sim_window$n <- 50
sim_window$sigma <- 1
sim_window$phi <- 4
sim_window$theta <- 0
sim_window$psi <- 1
sim_window$taper <- 7
sim_window$window <- 16

plan(sequential)

fwrite(sim_window, file = paste0(save_path, "sim_window16.csv"))

# Using 12 deg. by 12 deg. window instead of 20 by 20.
loc_list2 <- get_loc_list(lat_lon_grid, lat_margin = 6, lon_margin = 6)

sim_window2 <- run_sim_par(iter = 10,
                          lat_lon_grid = lat_lon_grid,
                          szon_dist = szon_dist,
                          smer_dist = smer_dist,
                          loc_list = loc_list2,
                          cov_mat = taper_cov,
                          beta = 1,
                          n = 50,
                          cov_L = chol_cov,
                          noise_sd = 1,
                          verbose = TRUE)

sim_window2$n <- 50
sim_window2$sigma <- 1
sim_window2$phi <- 4
sim_window2$theta <- 0
sim_window2$psi <- 1
sim_window2$taper <- 7
sim_window2$window <- 12

plan(sequential)

fwrite(sim_window2, file = paste0(save_path, "sim_window12.csv"))