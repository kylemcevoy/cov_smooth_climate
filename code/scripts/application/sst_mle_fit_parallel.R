# Script to fit the MLEs for the Matern-3/2 model on the July SST detrended
# anomalies. 

source("code/functions/blockwise_MLE.R")
source("code/functions/process_data.R")
source("code/functions/dist_cov.R")
source("code/functions/parallel.R")

# for parallelizing code comment out all plan statements, or replace with
# plan(sequential) if running on a single CPU.
library(future.apply)
plan(multisession, workers = 16)

data_path <- "/home/data/projects/clim_smooth/"
save_path <- "/home/data/projects/clim_smooth/"

month <- 'jul'

sst_filepath <- paste0(data_path, "sst_detrend_anom_", month, ".csv")
save_filepath <- paste0(save_path, "sst_", month, "_mle_fit.csv")

# Only gridboxes with fewer than 25% of values being identical are kept.
# This arises at locations that are frozen over for a portion of the
# year.
sst_anoms <- load_obs_data(path = sst_filepath,
                           max_const_prop = 0.25)

sst_coord_dt <- extract_coord_dt(sst_anoms)
  
sst_mat_std <- extract_data_matrix(sst_anoms, "sst_anoms", scale = TRUE)

sst_zonal_sdist <- find_signed_zonal_dist(sst_coord_dt)
sst_merid_sdist <- find_signed_merid_dist(sst_coord_dt)

init_cond_mat <- matrix(c(1, 1, 1,
                          0.2 * pi, 0.5 * pi, 0.8 * pi,
                          1.5, 1.5, 1.5),
                        nrow = 3)

loc_list <- get_loc_list(sst_coord_dt)

loc_pairs_list <- lapply(loc_list,
                         \(x) find_dist_pairs(zonal_dist_mat = sst_zonal_sdist,
                                              merid_dist_mat = sst_merid_sdist,
                                              loc_vec = x))

MLE_out <- future_mapply(find_MLE_mapply_wrapper,
                         loc_block = loc_list,
                         pairs_mat = loc_pairs_list,
                         MoreArgs = list(data_mat = sst_mat_std,
                                         lik_func = m1.5_lik),
                         SIMPLIFY = FALSE)

plan(sequential)

MLE_df <- rbindlist(MLE_out)

fwrite(MLE_df, file = save_filepath)
