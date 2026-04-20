# Script for finding the rejections for each method in the July Central US 2m 
# temperature against SST analysis in the paper.

source("code/functions/process_data.R")
source("code/functions/dist_cov.R")
source("code/functions/cov_smooth.R")
source("code/functions/blockwise_MLE.R")
source("code/functions/inference.R")

library(data.table)

data_path <- "/home/data/projects/clim_smooth/"
save_path <- "/home/data/projects/clim_smooth/"

t2m_path <- paste0(data_path, "t2m_GP_mean_jul.csv")
sst_jul_path <- paste0(data_path, "sst_detrend_anom_jul.csv")
mle_path <- paste0(data_path, "sst_jul_mle_fit.csv")

t2m_GP <- fread(t2m_path)
t2m_lm <- lm(t2m ~ time, data = t2m_GP)
t2m_anom <- t2m_lm$residuals

sst_jul <- load_obs_data(path = sst_jul_path,
                         max_const_prop = 0.25)

sst_jul_coord_dt <- extract_coord_dt(sst_jul)

sst_jul_mat_std <- extract_data_matrix(sst_jul, "sst_anoms", scale = TRUE)

sst_jul_zonal_sdist <- find_signed_zonal_dist(sst_jul_coord_dt)
sst_jul_merid_sdist <- find_signed_merid_dist(sst_jul_coord_dt)

sst_jul_mle <- fread(mle_path)
mle_mat <- data.matrix(sst_jul_mle[, 1:3])

loc_list_jul <- get_loc_list(sst_jul_coord_dt)

smoothing_weights <- find_smooth_weights(loc_list_jul,
                                         zonal_sdist = sst_jul_zonal_sdist,
                                         merid_sdist = sst_jul_merid_sdist,
                                         cov_func = matern_1.5_cov,
                                         mle_params = mle_mat)

sst_jul_mat_smooth <- sst_jul_mat_std %*% smoothing_weights

smooth_lm_fits <- fast_fit_lm_tstat(t2m_anom, sst_jul_mat_smooth)
base_lm_fits <- fast_fit_lm_tstat(t2m_anom, sst_jul_mat_std)

n <- nrow(sst_jul_mat_smooth)

smooth_pvals <- 2 * pt(abs(smooth_lm_fits), df = n - 2, lower.tail = FALSE)
base_pvals <- 2 * pt(abs(base_lm_fits), df = n - 2, lower.tail = FALSE)

smooth_reject <- find_BH_reject(smooth_pvals, alpha = 0.1)
base_reject <- find_BH_reject(base_pvals, alpha = 0.1)

smooth_reject$method <- "smooth"
base_reject$method <- "base"

reject_df <- rbind(base_reject, smooth_reject)

names(reject_df)[3] <- "reject"

no_BH <- data.table(location = 1:8195,
                    p_val = base_pvals,
                    reject = base_pvals < 0.05,
                    method = "no BH")

reject_df <- rbind(reject_df, no_BH)

fwrite(reject_df, file = paste0(save_path, "t2m_sst_rejections_jul.csv"))

