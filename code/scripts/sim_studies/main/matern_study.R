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
set.seed(482256)

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

matern0.5_taper_cov <- matern_0.5_spherical_cov(sim_dist, phi = 4, taper = 7)
matern0.5_chol_cov <- t(chol(matern0.5_taper_cov))

matern2.5_taper_cov <- matern_2.5_wendland_2_cov(sim_dist, phi = 4, taper = 7)
matern2.5_chol_cov <- t(chol(matern2.5_taper_cov))

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
                      cov_mat = matern0.5_taper_cov,
                      beta = 1,
                      n = 50,
                      cov_L = matern0.5_chol_cov,
                      noise_sd = 1,
                      verbose = TRUE)

plan(sequential)

sim_cov1$sigma <- 1
sim_cov1$n <- 50
sim_cov1$phi <- 4
sim_cov1$theta <- 0
sim_cov1$psi <- 1
sim_cov1$taper <- 7
sim_cov1$cov <- "matern_0.5_spherical"

sim_cov1_path <- paste0(save_dir, "sim_matern0_5.csv")
fwrite(sim_cov1, file = sim_cov1_path)

plan(multicore, workers = 24)

sim_cov2 <- run_sim_par(iter = 10,
                       lat_lon_grid = lat_lon_grid,
                       szon_dist = szon_dist,
                       smer_dist = smer_dist,
                       loc_list = loc_list,
                       cov_mat = matern2.5_taper_cov,
                       beta = 1,
                       n = 50,
                       cov_L = matern2.5_chol_cov,
                       noise_sd = 1,
                       verbose = TRUE)

plan(sequential)

sim_cov2$sigma <- 1
sim_cov2$n <- 50
sim_cov2$phi <- 4
sim_cov2$theta <- 0
sim_cov2$psi <- 1
sim_cov2$taper <- 7
sim_cov2$cov <- "matern_2.5_wendland2"

sim_cov2_path <- paste0(save_dir, "sim_matern2_5.csv")
fwrite(sim_cov2, file = sim_cov2_path)

sim_cov <- rbind(sim_cov1,
                 sim_cov2)

sim_cov$base_FDX <- sim_cov$base_FDP > 0.1
sim_cov$smooth_FDX <- sim_cov$smooth_FDP > 0.1

results_mean <- sim_cov[, lapply(.SD, mean), by = .(sigma,
                                                   n,
                                                   phi,
                                                   theta,
                                                   psi,
                                                   taper,
                                                   cov)]

results_mean[, iteration := NULL]

results_mean <- melt(results_mean, 
                     id.vars = c("sigma",
                                 "n",
                                 "phi",
                                 "theta",
                                 "psi",
                                 "taper",
                                 "cov"),
                     measure.vars = measure(method,
                                            value.name,
                                            pattern = "^([[:alpha:]]*)_(.*)$"))

results_mean <- results_mean[order(cov)]

fwrite(sim_cov, file = paste0(save_dir, "cov_misspecification_sim_study.csv"))
fwrite(results_mean, file = paste0(save_dir, "cov_misspecification_study_mean.csv"))
