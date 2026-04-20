source("code/functions/blockwise_MLE.R")
source("code/functions/cov_smooth.R")
source("code/functions/dist_cov.R")
source("code/functions/inference.R")
source("code/functions/sim.R")
source("code/functions/parallel.R")

library(future.apply)
library(parallelly)
library(data.table)

# necessary to prevent overhead of spinning up sandboxes for parallel sessions 
Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = FALSE)

avail_cores <- availableCores()
avail_cores
# set number of workers based on available cores.
plan(multicore, workers = 24)
# use plan sequential if you don't want to run the code in parallel.
# plan(sequential)

save_dir <- "/home/data/projects/clim_smooth/sim_studies/two_sided/"
aniso1_save_path <- paste0(save_dir, "sim_n50_sd1_phi1.3_aniso1_taper7.csv")
aniso2_save_path <- paste0(save_dir, "sim_n50_sd1_phi1.3_aniso2_taper7.csv")

# selected via rolling 6 10-sided dice.
set.seed(166663)

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

# arbitrary choice of non-zero theta to illustrate.
theta <- pi / 4
# value of psi chosen based on pacific ocean SST variograms.
# phi equals shorter length scale, phi * sqrt(psi) equals longer length scale.
psi <- (4 / 1.3)^2

aniso_mat <- find_aniso_mat(aniso_angle = theta, aniso_ratio = psi)
aniso_mat2 <- find_aniso_mat(aniso_angle = theta, aniso_ratio = 1.5 * psi)

sim_dist_pairs <- cbind(c(szon_dist), c(smer_dist))

aniso_dist <- sqrt(rowSums((sim_dist_pairs %*% aniso_mat) * sim_dist_pairs))
aniso_dist2 <- sqrt(rowSums((sim_dist_pairs %*% aniso_mat2) * sim_dist_pairs))

aniso_cov <- matern_1.5_wendland_1_cov(aniso_dist,
                                       phi = 1.3,
                                       taper = 7,
                                       sigma2 = 1)

aniso_cov <- matrix(aniso_cov, nrow = p, ncol = p)

# The two aniso_dist calcs. incorporate the different psi parameters, 
# leading to different values for the two covariances.
aniso_cov2 <- matern_1.5_wendland_1_cov(aniso_dist2,
                                        phi = 1.3,
                                        taper = 7,
                                        sigma2 = 1)

aniso_cov2 <- matrix(aniso_cov2, nrow = p, ncol = p)

chol_cov <- t(chol(aniso_cov))
chol_cov2 <- t(chol(aniso_cov2))

loc_list <- get_loc_list(lat_lon_grid)

init_cond_mat <- matrix(c(1, 1, 1,
                          0.2 * pi, 0.5 * pi, 0.8 * pi,
                          1.5, 1.5, 1.5),
                        nrow = 3)

sim_aniso1 <- run_sim_par(iter = 10,
                          lat_lon_grid = lat_lon_grid,
                          szon_dist = szon_dist,
                          smer_dist = smer_dist,
                          loc_list = loc_list,
                          cov_mat = aniso_cov,
                          beta = 1,
                          n = 50,
                          cov_L = chol_cov,
                          noise_sd = 1,
                          verbose = TRUE)

sim_aniso1$n <- 50
sim_aniso1$sigma <- 1
sim_aniso1$phi <- 1.3
sim_aniso1$theta <- theta
sim_aniso1$psi <- psi
sim_aniso1$taper <- 7

fwrite(sim_aniso1, file = aniso1_save_path)

# sometimes worker processes don't release properly, resetting to 
# plan(sequential) before a new plan() seems to clear up some issues.
plan(sequential)
plan(multicore, workers = 24)

sim_aniso2 <- run_sim_par(iter = 10,
                          lat_lon_grid = lat_lon_grid,
                          szon_dist = szon_dist,
                          smer_dist = smer_dist,
                          loc_list = loc_list,
                          cov_mat = aniso_cov2,
                          beta = 1,
                          n = 50,
                          cov_L = chol_cov2,
                          noise_sd = 1,
                          verbose = TRUE)

sim_aniso2$n <- 50
sim_aniso2$sigma <- 1
sim_aniso2$phi <- 1.3
sim_aniso2$theta <- theta
sim_aniso2$psi <- 1.5 * psi
sim_aniso2$taper <- 7

fwrite(sim_aniso2, file = aniso2_save_path)

sim_aniso <- rbind(sim_aniso1,
                   sim_aniso2)

sim_aniso$base_FDX <- sim_aniso$base_FDP > 0.1
sim_aniso$smooth_FDX <- sim_aniso$smooth_FDP > 0.1

aniso_mean_results <- sim_aniso[, lapply(.SD, mean), by = .(sigma,
                                                          n,
                                                          phi,
                                                          theta,
                                                          psi,
                                                          taper)]

aniso_mean_results[, iteration := NULL]

aniso_sim_study_mean <- melt(aniso_mean_results, 
                             id.vars = c("sigma",
                                         "n",
                                         "phi",
                                         "theta",
                                         "psi",
                                         "taper"),
                             measure.vars = measure(method,
                                                    value.name,
                                                    pattern = "^([[:alpha:]]*)_(.*)$"))

aniso_sim_study_mean <- aniso_sim_study_mean[order(sigma)]

fwrite(sim_aniso, file = paste0(save_dir, "aniso_sim_study.csv"))
fwrite(aniso_sim_study_mean, file = paste0(save_dir, "aniso_sim_study_mean.csv"))
