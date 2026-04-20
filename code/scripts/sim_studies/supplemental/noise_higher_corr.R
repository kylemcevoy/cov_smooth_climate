source("code/functions/blockwise_MLE.R")
source("code/functions/cov_smooth.R")
source("code/functions/dist_cov.R")
source("code/functions/inference.R")
source("code/functions/sim.R")
source("code/functions/parallel.R")

library(data.table)
library(future.apply)
library(parallelly)

save_dir <- "/home/data/projects/clim_smooth/sim_studies/two_sided/"

# set number of workers based on available cores.
plan(multicore, workers = 24)
# use plan sequential if you don't want to run the code in parallel.
# plan(sequential)

# selected via rolling 6 10-sided dice.
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

sim_data <- gen_data(50, L=chol_cov)
sim_data_std <- scale(sim_data)

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
                      n = 50,
                      cov_L = chol_cov,
                      noise_sd = 1,
                      cor_thresh = 0.1,
                      verbose = TRUE)

sim_sd1$sigma <- 1
sim_sd1$n <- 50
sim_sd1$phi <- 4
sim_sd1$theta <- 0
sim_sd1$psi <- 1
sim_sd1$taper <- 7

sim_sd1_path <- paste0(save_dir, "sim_n50_sd1_phi4_taper7_corthresh01.csv")
fwrite(sim_sd1, file = sim_sd1_path)

plan(sequential)
plan(multicore, workers = 24)

sim_sd2 <- run_sim_par(iter = 10,
                       lat_lon_grid = lat_lon_grid,
                       szon_dist = szon_dist,
                       smer_dist = smer_dist,
                       loc_list = loc_list,
                       cov_mat = taper_cov,
                       beta = 1,
                       n = 50,
                       cov_L = chol_cov,
                       noise_sd = 2,
                       cor_thresh = 0.1,
                       verbose = TRUE)

sim_sd2$sigma <- 2
sim_sd2$n <- 50
sim_sd2$phi <- 4
sim_sd2$theta <- 0
sim_sd2$psi <- 1
sim_sd2$taper <- 7

sim_sd2_path <- paste0(save_dir, "sim_n50_sd2_phi4_taper7_corthresh01.csv")
fwrite(sim_sd2, file = sim_sd2_path)

plan(sequential)
plan(multicore, workers = 24)

sim_sd3 <- run_sim_par(iter = 10,
                       lat_lon_grid = lat_lon_grid,
                       szon_dist = szon_dist,
                       smer_dist = smer_dist,
                       loc_list = loc_list,
                       cov_mat = taper_cov,
                       beta = 1,
                       n = 50,
                       cov_L = chol_cov,
                       noise_sd = 3,
                       cor_thresh = 0.1,
                       verbose = TRUE)

plan(sequential)

sim_sd3$sigma <- 3
sim_sd3$n <- 50
sim_sd3$phi <- 4
sim_sd3$theta <- 0
sim_sd3$psi <- 1
sim_sd3$taper <- 7

sim_sd3_path <- paste0(save_dir, "sim_n50_sd3_phi4_taper7_corthresh01.csv")
fwrite(sim_sd3, file = sim_sd3_path)

sim_sd <- rbind(sim_sd1,
                sim_sd2,
                sim_sd3)

sim_sd$base_FDX <- sim_sd$base_FDP > 0.1
sim_sd$smooth_FDX <- sim_sd$smooth_FDP > 0.1

results_mean <- sim_sd[, lapply(.SD, mean), by = .(sigma,
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

results_mean <- results_mean[order(sigma)]

fwrite(sim_sd, file = paste0(save_dir, "noise_sim_study_corthresh01.csv"))
fwrite(results_mean, file = paste0(save_dir, "noise_sim_study_mean_corthresh01.csv"))
