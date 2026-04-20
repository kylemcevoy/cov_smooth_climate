# Script to generate Figure 3 showing rejections with SSTs covariates and
# synthetic response variable under multiple samplings of the noise.

library(sf)
library(ggplot2)
library(maps)
library(data.table)
library(patchwork)

source("code/functions/dist_cov.R")
source("code/functions/process_data.R")
source("code/functions/inference.R")
source("code/functions/cov_smooth.R")
source("code/functions/blockwise_MLE.R")

data_dir <- '/home/data/projects/clim_smooth/'
plot_dir <- '/home/data/projects/clim_smooth/plots/'

sst_anoms_path <- paste0(data_dir,'sst_detrend_anom_jul.csv')
sst_mle_path <- paste0(data_dir, 'sst_jul_mle_fit.csv')

center_loc <- 3500

set.seed(160326)

axis_labels_theme <- theme(axis.title.x = element_text(size = 12),
                           axis.text.x = element_text(size = 8),
                           axis.title.y = element_text(size = 12),
                           axis.text.y = element_text(size = 8),
                           legend.title = element_text(size = 10),
                           legend.text = element_text(size = 11)) 

sst_jul_mle <- fread(sst_mle_path)

sst_jul <- load_obs_data(path = sst_anoms_path,
                         max_const_prop = 0.25)


sst_jul_coord_dt <- extract_coord_dt(sst_jul)

param_dt <- cbind(sst_jul_coord_dt, sst_jul_mle)

names(param_dt)[c(5, 6)] <- c("theta", "psi")

param_dt$phi_psi <- param_dt$phi * sqrt(param_dt$psi)

land_map <- sf::st_as_sf(map(wrap = c(0, 360), plot = FALSE, fill = TRUE))

sst_jul_loc_list <- get_loc_list(sst_jul_coord_dt)

sst_zonal_sdist <- find_signed_zonal_dist(sst_jul_coord_dt)
sst_merid_sdist <- find_signed_merid_dist(sst_jul_coord_dt)

sst_mat <- extract_data_matrix(sst_jul, var = "sst_anoms", scale = TRUE)

mle_param_mat <- data.matrix(param_dt[, 4:6])

sst_cov_weights <- find_smooth_weights(loc_list = sst_jul_loc_list,
                                       zonal_sdist = sst_zonal_sdist,
                                       merid_sdist = sst_merid_sdist,
                                       cov_func = matern_1.5_cov,
                                       mle_params = mle_param_mat,
                                       progress = FALSE)

sst_jul_smooth <- sst_mat %*% sst_cov_weights

sim_rejections <- function(indx, reps) {
  p <- ncol(sst_mat)
  n <- nrow(sst_mat)
  output_data_list <- vector(mode = "list", length = reps)
  
  for (i in 1:reps) {
    y <- sst_mat[, indx] + rnorm(50)
    
    base_tstat <- fit_lm_tstat(y, data_mat = sst_mat)
    smooth_tstat <- fit_lm_tstat(y, data_mat = sst_jul_smooth)
    
    base_pvals <- 2 * pt(abs(base_tstat), df = n - 2, lower.tail = FALSE)
    smooth_pvals <- 2 * pt(abs(smooth_tstat), df = n - 2, lower.tail = FALSE)
    
    BH_base <- find_BH_reject(base_pvals, alpha = 0.1)
    BH_smooth <- find_BH_reject(smooth_pvals, alpha = 0.1)
    
    reject_base_only <- BH_base$BH_reject & !BH_smooth$BH_reject
    reject_smooth_only <- !BH_base$BH_reject & BH_smooth$BH_reject
    reject_both <- BH_base$BH_reject & BH_smooth$BH_reject
    reject_any <- BH_base$BH_reject | BH_smooth$BH_reject
    
    reject_type_vec <- vector(mode = "character", length = p)
    reject_type_vec[reject_base_only] <- "base"
    reject_type_vec[reject_smooth_only] <- "smooth"
    reject_type_vec[reject_both] <- "both"
    reject_type_vec[!reject_any] <- "none"
    reject_type_vec <- factor(reject_type_vec, levels = c("none", "base", "smooth", "both"))
    
    BH_data <- cbind(sst_jul_coord_dt,
                     rep = i, 
                     reject_type = reject_type_vec,
                     reject_any = reject_any)
    
    output_data_list[[i]] <- BH_data
  }
  
  return(output_data_list)
}

output_data <- sim_rejections(indx = center_loc, reps = 8)

BH_data <- rbindlist(output_data)

library(viridis)
color_test <- viridis(12)[c(1, 5, 9)]
names(color_test) <- c("both", "smooth", "base")
desat_color <- colorspace::desaturate(color_test)


sst_BH_example_plot <- ggplot(data = land_map) + 
  geom_sf(fill = 'gray', color = "gray", inherit.aes = FALSE) +
  coord_sf(xlim = c(0, 360), ylim = c(-60, 60), expand = FALSE) +
  geom_point(data = BH_data[reject_any == TRUE],
             aes(x = lon,
                 y = lat,
                 color = reject_type),
             alpha = 0.5,
             size = 0.5) + 
  scale_x_continuous(breaks = seq(0, 340, by = 30),
                     labels =  \(x) as.character(x)
  ) + 
  scale_y_continuous(breaks = seq(-40, 60, by = 20),
                     labels = \(x) as.character(x)) +
  geom_point(data = BH_data[location == 3500],
             aes(x = lon, y = lat),
             color = "red") + 
  labs(x = "longitude",
       y = "latitude",
       color = "rejections") +
  # scale_color_viridis_d(begin = 0.2,
  #                       end = 0.8,
  #                       option = 'D',
  #                       direction = -1) +
  scale_color_manual(values = color_test) +
 # scale_shape_manual(values = c(3, 1, 2)) + 
  facet_wrap(vars(rep), nrow = 4) + 
  guides(color = guide_legend(override.aes = list(size=4))) +
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.text = element_blank()) +
  axis_labels_theme

sst_BH_example_plot 

ggsave(file = paste0(plot_dir, 'sst_BH_example_grid.png'),
       sst_BH_example_plot)

sst_BH_example_plot_gray <- ggplot(data = land_map) + 
  geom_sf(fill = 'gray', color = "gray", inherit.aes = FALSE) +
  coord_sf(xlim = c(0, 360), ylim = c(-60, 60), expand = FALSE) +
  geom_point(data = BH_data[reject_any == TRUE],
             aes(x = lon,
                 y = lat,
                 color = reject_type),
             alpha = 0.5,
             size = 0.8) + 
  scale_x_continuous(breaks = seq(0, 340, by = 30),
                     labels =  \(x) as.character(x)
  ) + 
  scale_y_continuous(breaks = seq(-60, 60, by = 20),
                     labels = \(x) as.character(x)) +
  geom_point(data = BH_data[location == center_loc],
             aes(x = lon, y = lat),
             color = "red") + 
  labs(x = "longitude",
       y = "latitude",
       color = "rejections") +
  # scale_color_viridis_d(begin = 0.2,
  #                       end = 0.8,
  #                       option = 'D',
  #                       direction = -1) +
  scale_color_manual(values = desat_color) +
  # scale_shape_manual(values = c(3, 1, 2)) + 
  facet_wrap(vars(rep), nrow = 4) + 
  guides(color = guide_legend(override.aes = list(size=4))) 

sst_BH_example_plot_gray + theme_bw() + theme(strip.background = element_blank(),
                                         strip.text = element_blank())
