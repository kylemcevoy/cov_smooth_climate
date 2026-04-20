source("code/functions/blockwise_MLE.R")
source("code/functions/cov_smooth.R")
source("code/functions/dist_cov.R")
source("code/functions/inference.R")
source("code/functions/sim.R")
source("code/functions/parallel.R")

library(data.table)
library(future.apply)
library(parallelly)

save_path <- "/home/data/projects/clim_smooth/sim_studies/"

avail_cores <- availableCores()
avail_cores
# set number of workers based on available cores.
plan(multicore, workers = 24)
# use plan sequential if you don't want to run the code in parallel.
# plan(sequential)

# selected via rolling 6 10-sided dice.
set.seed(881528)

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

sim_sd1 <- run_sim_par(iter = 10,
                       lat_lon_grid = lat_lon_grid,
                       szon_dist = szon_dist,
                       smer_dist = smer_dist,
                       loc_list = loc_list,
                       cov_mat = taper_cov,
                       beta = 1,
                       n = 100,
                       cov_L = chol_cov,
                       noise_sd = 1,
                       verbose = TRUE)

sim_sd1$n <- 100
sim_sd1$sigma <- 1
sim_sd1$phi <- 4
sim_sd1$theta <- 0
sim_sd1$psi <- 1
sim_sd1$taper <- 7

plan(sequential)

fwrite(sim_sd1, file = paste0(save_path, "sim_n100_sd1_phi4_taper7.csv"))

plan(multicore, workers = 24)

sim_sd2 <- run_sim_par(iter = 10,
                       lat_lon_grid = lat_lon_grid,
                       szon_dist = szon_dist,
                       smer_dist = smer_dist,
                       loc_list = loc_list,
                       cov_mat = taper_cov,
                       beta = 1,
                       n = 100,
                       cov_L = chol_cov,
                       noise_sd = 2,
                       verbose = TRUE)

sim_sd2$n <- 100
sim_sd2$sigma <- 2
sim_sd2$phi <- 4
sim_sd2$theta <- 0
sim_sd2$psi <- 1
sim_sd2$taper <- 7

plan(sequential)

fwrite(sim_sd2, file = paste0(save_path, "sim_n100_sd2_phi4_taper7.csv"))

plan(multicore, workers = 24)

sim_sd3 <- run_sim_par(iter = 10,
                       lat_lon_grid = lat_lon_grid,
                       szon_dist = szon_dist,
                       smer_dist = smer_dist,
                       loc_list = loc_list,
                       cov_mat = taper_cov,
                       beta = 1,
                       n = 100,
                       cov_L = chol_cov,
                       noise_sd = 3,
                       verbose = TRUE)

sim_sd3$n <- 100
sim_sd3$sigma <- 3
sim_sd3$phi <- 4
sim_sd3$theta <- 0
sim_sd3$psi <- 1
sim_sd3$taper <- 7

plan(sequential)

fwrite(sim_sd3, file = paste0(save_path, "sim_n100_sd3_phi4_taper7.csv"))

sim_n100 <- rbind(sim_sd1,
                  sim_sd2,
                  sim_sd3)

sim_n100$base_FDX <- sim_n100$base_FDP > 0.1
sim_n100$smooth_FDX <- sim_n100$smooth_FDP > 0.1

n100_mean_results <- sim_n100[, lapply(.SD, mean), by = .(sigma,
                                                          n,
                                                          phi,
                                                          theta,
                                                          psi,
                                                          taper)]

n100_mean_results[, iteration := NULL]

n_sim_study_mean <- melt(n100_mean_results, 
                         id.vars = c("sigma", "n", "phi", "theta", "psi", "taper"),
                         measure.vars = measure(method,
                                                value.name,
                                                pattern = "^([[:alpha:]]*)_(.*)$"))

n_sim_study_mean <- n_sim_study_mean[order(sigma)]

fwrite(sim_n100, file = paste0(save_path, "n_sim_study.csv"))
fwrite(n_sim_study_mean, file = paste0(save_path, "n_sim_study_mean.csv"))
