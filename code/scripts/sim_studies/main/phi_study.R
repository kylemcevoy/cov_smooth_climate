source("code/functions/blockwise_MLE.R")
source("code/functions/cov_smooth.R")
source("code/functions/dist_cov.R")
source("code/functions/inference.R")
source("code/functions/sim.R")
source("code/functions/parallel.R")

library(data.table)
library(future.apply)
library(parallelly)

# necessary to prevent overhead of spinning up sandboxes for parallel sessions 
Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = FALSE)

save_dir <- "/home/data/projects/clim_smooth/sim_studies/two_sided/"

avail_cores <- availableCores()
avail_cores
# set number of workers based on available cores.
plan(multicore, workers = 24)
# use plan sequential if you don't want to run the code in parallel.
# plan(sequential)

# selected via rolling 6 10-sided dice.
set.seed(220067)

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

taper_cov_0.8 <- matern_1.5_wendland_1_cov(sim_dist, phi = 0.8, taper = 7)
chol_cov_0.8 <- t(chol(taper_cov_0.8))

taper_cov_1.3 <- matern_1.5_wendland_1_cov(sim_dist, phi = 1.3, taper = 7)
chol_cov_1.3 <- t(chol(taper_cov_1.3))

taper_cov_1.3_5t <- matern_1.5_wendland_1_cov(sim_dist, phi = 1.3, taper = 5)
chol_cov_1.3_5t <- t(chol(taper_cov_1.3_5t))

taper_cov_2 <- matern_1.5_wendland_1_cov(sim_dist, phi = 2, taper = 7)
chol_cov_2 <- t(chol(taper_cov_2))

loc_list <- get_loc_list(lat_lon_grid)

init_cond_mat <- matrix(c(1, 1, 1,
                          0.2 * pi, 0.5 * pi, 0.8 * pi,
                          1.5, 1.5, 1.5),
                        nrow = 3)

sim_cov1 <- run_sim_par(iter = 10,
                      lat_lon_grid = lat_lon_grid,
                      szon_dist = szon_dist,
                      smer_dist = smer_dist,
                      loc_list = loc_list,
                      cov_mat = taper_cov_0.8,
                      beta = 1,
                      n = 50,
                      cov_L = chol_cov_0.8,
                      noise_sd = 1,
                      verbose = TRUE)

plan(sequential)

sim_cov1$sigma <- 1
sim_cov1$n <- 50
sim_cov1$phi <- 0.8
sim_cov1$theta <- 0
sim_cov1$psi <- 1
sim_cov1$taper <- 7

sim_cov1_path <- paste0(save_dir, "sim_n50_sd1_phi08_taper7.csv")
fwrite(sim_cov1, file = sim_cov1_path)

plan(multicore, workers = 24)

sim_cov2 <- run_sim_par(iter = 10,
                       lat_lon_grid = lat_lon_grid,
                       szon_dist = szon_dist,
                       smer_dist = smer_dist,
                       loc_list = loc_list,
                       cov_mat = taper_cov_1.3,
                       beta = 1,
                       n = 50,
                       cov_L = chol_cov_1.3,
                       noise_sd = 1,
                       verbose = TRUE)

plan(sequential)

sim_cov2$sigma <- 1
sim_cov2$n <- 50
sim_cov2$phi <- 1.3
sim_cov2$theta <- 0
sim_cov2$psi <- 1
sim_cov2$taper <- 7

sim_cov2_path <- paste0(save_dir, "sim_n50_sd1_phi1.3_taper7.csv")
fwrite(sim_cov2, file = sim_cov2_path)

plan(multicore, workers = 24)

sim_cov3 <- run_sim_par(iter = 10,
                       lat_lon_grid = lat_lon_grid,
                       szon_dist = szon_dist,
                       smer_dist = smer_dist,
                       loc_list = loc_list,
                       cov_mat = taper_cov_1.3_5t,
                       beta = 1,
                       n = 50,
                       cov_L = chol_cov_1.3_5t,
                       noise_sd = 1,
                       verbose = TRUE)

plan(sequential)

sim_cov3$sigma <- 1
sim_cov3$n <- 50
sim_cov3$phi <- 1.3
sim_cov3$theta <- 0
sim_cov3$psi <- 1
sim_cov3$taper <- 5

sim_cov3_path <- paste0(save_dir, "sim_n50_sd1_phi1.3_taper5.csv")
fwrite(sim_cov3, file = sim_cov3_path)

plan(multicore, workers = 24)

sim_cov4 <- run_sim_par(iter = 10,
                        lat_lon_grid = lat_lon_grid,
                        szon_dist = szon_dist,
                        smer_dist = smer_dist,
                        loc_list = loc_list,
                        cov_mat = taper_cov_2,
                        beta = 1,
                        n = 50,
                        cov_L = chol_cov_2,
                        noise_sd = 1,
                        verbose = TRUE)

plan(sequential)

sim_cov4$sigma <- 1
sim_cov4$n <- 50
sim_cov4$phi <- 2
sim_cov4$theta <- 0
sim_cov4$psi <- 1
sim_cov4$taper <- 7

sim_cov4_path <- paste0(save_dir, "sim_n50_sd1_phi2_taper7.csv")
fwrite(sim_cov4, file = sim_cov4_path)

sim_cov <- rbind(sim_cov1,
                 sim_cov2,
                 sim_cov3,
                 sim_cov4)

sim_cov$base_FDX <- sim_cov$base_FDP > 0.1
sim_cov$smooth_FDX <- sim_cov$smooth_FDP > 0.1

results_mean <- sim_cov[, lapply(.SD, mean), by = .(sigma,
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

fwrite(sim_cov, file = paste0(save_dir, "phi_sim_study.csv"))
fwrite(results_mean, file = paste0(save_dir, "phi_sim_study_mean.csv"))
