# Script to fit the MLEs for the Matern-3/2 model on the July SST detrended
# anomalies. Run using Rscript from the terminal.

source("code/functions/blockwise_MLE.R")
source("code/functions/process_data.R")
source("code/functions/dist_cov.R")
source("code/functions/parallel.R")

library(data.table)

# plan(sequential) runs on a single CPU, plan(multicore) or plan(multisession)
# run on multiple cpus. multicore generally preferred for memory management;
# however, it does not work with Rstudio or on Windows architectures. In those
# cases, use plan(multisession) instead.
library(future.apply)
plan(sequential)
# on Unix type OS can use plan(multicore). Do not run from Rstudio.
#plan(multicore, workers = 16)
# on Windows can use multisession, but note that the memory load will be much
# more intensive since multiple background sessions with ~4.5 GB of data will
# be spawned. 
#plan(multisession, workers = 4)

data_path <- "data/"
save_path <- "data/"

month <- 'jul'

sst_filepath <- paste0(data_path, "sst_detrend_anom_jul.csv")
save_filepath <- paste0(save_path, "sst_jul_mle_fit.csv")

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