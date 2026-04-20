source("code/functions/blockwise_MLE.R")
source("code/functions/cov_smooth.R")
source("code/functions/dist_cov.R")
source("code/functions/inference.R")
source("code/functions/sim.R")
source("code/functions/parallel.R")

library(data.table)
library(future.apply)
library(ggplot2)
library(scales)

set.seed(581751)

plan(multicore, workers = 24)

save_dir <- "/home/data/projects/clim_smooth/sim_studies/two_sided/"
plots_dir <- "/home/data/projects/clim_smooth/plots/"

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

loc_pairs_list <- lapply(loc_list,
                           \(x) find_dist_pairs(zonal_dist_mat = szon_dist,
                                            merid_dist_mat = smer_dist,
                                            loc_vec = x))

init_cond_mat <- matrix(c(1, 1, 1,
                          0.2 * pi, 0.5 * pi, 0.8 * pi,
                          1.5, 1.5, 1.5),
                        nrow = 3)


MLE_out <- future_mapply(find_MLE_mapply_wrapper,
                             loc_block = loc_list,
                             pairs_mat = loc_pairs_list,
                             MoreArgs = list(data_mat = sim_data_std,
                                             lik_func = m1.5_lik),
                             SIMPLIFY = FALSE)
    
MLE_df <- rbindlist(MLE_out)
MLE_mat <- data.matrix(MLE_df[, 1:3])

l <- sample(10980, 1)
n <- 50

y <- sim_data[, l] + rnorm(n)

cov_weights <- find_smooth_weights(loc_list = loc_list,
                                    zonal_sdist = szon_dist,
                                    merid_sdist = smer_dist,
                                    cov_func = matern_1.5_cov,
                                    mle_params = MLE_mat,
                                    progress = FALSE)

smooth_data <- sim_data_std %*% cov_weights

base_lm_fits <- fast_fit_lm_tstat(y = y, data_mat = sim_data_std)
smooth_lm_fits <- fast_fit_lm_tstat(y = y, data_mat = smooth_data)

true_cov <- taper_cov[l, ]

true_null <- true_cov == 0

null_locs <- true_cov == 0

base_pvals <- 2 * pt(abs(base_lm_fits),
                           df = n - 2,
                           lower.tail = FALSE)
smooth_pvals <- 2 * pt(abs(smooth_lm_fits),
                        df = n - 2,
                        lower.tail = FALSE)

base_reject <- find_BH_reject(base_pvals, alpha = 0.1)
smooth_reject <- find_BH_reject(smooth_pvals, alpha = 0.1)

base_metrics <- estim_metrics(p_dt = base_reject,
                              true_null = true_null)
smooth_metrics <- estim_metrics(p_dt = smooth_reject,
                                true_null = true_null)

hist(base_lm_fits)
hist(smooth_lm_fits)

null_t_stats <- data.table(base_t_stat = base_lm_fits[null_locs], smooth_t_stat = smooth_lm_fits[null_locs])

null_tstats_long <- melt(null_t_stats, measure.vars = measure(method, 
  value.name, 
  pattern = "^([[:alpha:]]*)_(.*)$"))

base_quantiles <- quantile(null_t_stats$base_t_stat, c(0.05, 0.95))
base_value_text <- paste(
  round(base_quantiles[1], 2),
  round(base_quantiles[2], 2),
  sep = ', ')

base_text <- paste0('base quant.: (', base_value_text, ')') 
base_color <- hue_pal()(2)[1]

smooth_quantiles <- quantile(null_t_stats$smooth_t_stat, c(0.05, 0.95))
smooth_value_text <- paste(
  round(smooth_quantiles[1], 2),
  round(smooth_quantiles[2], 2),
  sep = ', ')

smooth_text <- paste0('smooth quant.: (', smooth_value_text, ')') 
smooth_color <- hue_pal()(2)[2]

(null_tstat_plot <- ggplot(data = null_tstats_long, aes(x = t_stat, fill = method)) + 
  geom_histogram(position = "identity", alpha = 0.5) +
  theme_bw() +
  labs(x = 't statistics at null locations') + 
  annotate('text', x = -2, y = 925, label = base_text, size = 3.2, color = base_color) + 
  annotate('text', x = -2, y = 895, size = 3.2, label = smooth_text, color = smooth_color))

null_tstat_plot_path <- paste0(plots_dir, 'null_tstat_comp.png')
ggsave(null_tstat_plot_path,
  plot = null_tstat_plot,
  units = "in",
  width = 8,
  height = 5.5)

lapply(null_t_stats, min)
lapply(null_t_stats, max)

quantile(null_t_stats$base_t_stat, c(0.025, 0.975))
quantile(null_t_stats$smooth_t_stat, c(0.025, 0.975))

