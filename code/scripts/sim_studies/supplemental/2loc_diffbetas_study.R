source("code/functions/blockwise_MLE.R")
source("code/functions/cov_smooth.R")
source("code/functions/dist_cov.R")
source("code/functions/inference.R")
source("code/functions/sim.R")
source("code/functions/parallel.R")

library(future.apply)
library(parallelly)
library(data.table)

save_path <- "/home/data/projects/clim_smooth/sim_studies/two_sided/"

# set number of workers based on available cores.
plan(multicore, workers = 16)
# use plan sequential if you don't want to run the code in parallel.
# plan(sequential)

# selected via rolling 6 10-sided dice.
set.seed(939123)

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

loc_list <- get_loc_list(lat_lon_grid)

init_cond_mat <- matrix(c(1, 1, 1,
                          0.2 * pi, 0.5 * pi, 0.8 * pi,
                          1.5, 1.5, 1.5),
                        nrow = 3)

beta_vec <- c(1.5, 0.5)

sim_2loc <- run_sim_par_2loc(iter = 10,
                             lat_lon_grid = lat_lon_grid,
                             szon_dist = szon_dist,
                             smer_dist = smer_dist,
                             loc_list = loc_list,
                             cov_mat = taper_cov,
                             beta_vec = beta_vec,
                             n = 50,
                             cov_L = chol_cov,
                             noise_sd = 1,
                             verbose = TRUE)

sim_2loc$n <- 50
sim_2loc$sigma <- 1
sim_2loc$phi <- 4
sim_2loc$theta <- 0
sim_2loc$psi <- 1
sim_2loc$taper <- 7
sim_2loc$beta1 <- beta_vec[1]
sim_2loc$beta2 <- beta_vec[2]

plan(sequential)

fwrite(sim_2loc, file = paste0(save_path, "sim_2loc_beta1_5_0_5.csv"))