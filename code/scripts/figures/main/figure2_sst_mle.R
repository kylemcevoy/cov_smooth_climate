# Script to generate Figure 2 of the paper that shows the fitted MLE values of 
# phi, theta, and phi * sqrt(psi) for detrended July SST anomalies.

library(ggplot2)
library(maps) # world land maps
library(sf) # spatial plots
library(data.table)
library(patchwork) # spatially arrange ggplots

source("code/functions/dist_cov.R")
source("code/functions/process_data.R")
source("code/functions/inference.R")
source("code/functions/cov_smooth.R")
source("code/functions/blockwise_MLE.R")

data_dir <- "/home/data/projects/clim_smooth/"
plots_dir <- "/home/data/projects/clim_smooth/plots/"

axis_labels_theme <- theme(axis.title.x = element_text(size = 12),
                           axis.text.x = element_text(size = 11),
                           axis.title.y = element_text(size = 12),
                           axis.text.y = element_text(size = 11),
                           legend.title = element_text(size = 14),
                           legend.text = element_text(size = 11),
                           legend.key.height = unit(0.6, "cm"))

mle_path <- paste0(data_dir, "sst_jul_mle_fit.csv")
sst_path <- paste0(data_dir, "sst_detrend_anom_jul.csv")

sst_jul_mle <- fread(mle_path)

sst_jul <- load_obs_data(path = sst_path,
                         max_const_prop = 0.25)


sst_jul_coord_dt <- extract_coord_dt(sst_jul)

param_dt <- cbind(sst_jul_coord_dt, sst_jul_mle)

names(param_dt)[c(5, 6)] <- c("theta", "psi")

param_dt$phi_psi <- param_dt$phi * sqrt(param_dt$psi)

land_map <- sf::st_as_sf(map(wrap = c(0, 360), plot = FALSE, fill = TRUE))

(phi_plot <- ggplot(data = land_map) +
  geom_tile(data = param_dt, aes(x = lon, y = lat, fill = phi)) +
  geom_sf(fill = "#D2D2D2", color = "#D2D2D2", inherit.aes = FALSE) +
  scale_fill_binned(type = 'viridis',
                    name = expression(paste(~~phi, " [", 1 %.% 10^3, " km]")),
                    breaks = seq(0, 3, by = 1/2),
                    limits = c(0, 3)
  ) +
  coord_sf(xlim = c(0, 360), ylim = c(-60, 60), expand = FALSE) +
  scale_x_continuous(breaks = seq(0, 360, by = 50),
                     labels =  \(x) as.character(x)
  ) + 
  scale_y_continuous(breaks = seq(-60, 60, by = 20),
                     labels = \(x) as.character(x)) +
    labs(x = NULL, y = 'latitude')
)

(
theta_plot <- ggplot(data = land_map) +
  geom_tile(data = param_dt, aes(x = lon, y = lat, fill = theta)) +
  geom_sf(fill = "#D2D2D2", color = "#D2D2D2", inherit.aes = FALSE) +
  binned_scale(aesthetics = "fill",
                 palette = function(x) c("#6DCC57",
                                          "#1F9D87",
                                          "#440154",
                                          "#3D4988", 
                                          "#FDE725"),
                 guide = "colorsteps",
                 show.limits = FALSE,
                 name = expression(paste(theta, " [radians]")),
                 limits = c(0, pi),
                 breaks = seq(0, pi,  by = pi / 5),
                 labels = as.list(c(0,
                           expression(pi / 5),
                           expression(2 * pi / 5),
                           expression(3 * pi / 5),
                           expression(4 * pi / 5),
                           expression(pi)))
  ) +
  coord_sf(xlim = c(0, 360), ylim = c(-60, 60), expand = FALSE) +
  scale_x_continuous(breaks = seq(0, 360, by = 50),
                     labels =  \(x) as.character(x)
  ) + 
  scale_y_continuous(breaks = seq(-60, 60, by = 20),
                     labels = \(x) as.character(x)) + 
  labs(x = NULL, y = 'latitude') +
    theme_bw()
)

phipsi_plot <- ggplot(data = land_map) +
  geom_tile(data = param_dt, aes(x = lon, y = lat, fill = phi_psi)) +
  geom_sf(fill = "#D2D2D2", color = "#D2D2D2", inherit.aes = FALSE) +
  scale_fill_binned(type = "viridis",
                    name = expression(paste(phi * sqrt(psi), " [", 1 %.% 10^3, " km]")),
                    breaks = seq(0, 5, by = 1),
                    limits = c(0, 5)
  )  + 
  coord_sf(xlim = c(0, 360), ylim = c(-60, 60), expand = FALSE) +
  scale_x_continuous(breaks = seq(0, 360, by = 50),
                     labels =  \(x) as.character(x)
  ) + 
  scale_y_continuous(breaks = seq(-60, 60, by = 20),
                     labels = \(x) as.character(x)) +
  labs(x = NULL, y = 'latitude')

# aligned_plots <- AlignPlots(phi_plot, theta_plot, phipsi_plot)

# Patchwork command to stack plots
gridded_mle <- phi_plot / theta_plot / phipsi_plot + labs(x = 'longitude') & theme(legend.justification='left')

mle_plot_path <- paste0(plots_dir, 'stacked_sst_mle_test.png')

ggsave(mle_plot_path,
       gridded_mle,
       units = 'in',
       width = 8,
       height = 8)

sst_jul_loc_list <- get_loc_list(sst_jul_coord_dt)

sst1_zonal_sdist <- find_signed_zonal_dist(sst_jul_coord_dt)
sst1_merid_sdist <- find_signed_merid_dist(sst_jul_coord_dt)

sst_mat <- extract_data_matrix(sst_jul, var = "sst_anoms", scale = TRUE)

mle_param_mat <- data.matrix(param_dt[, 4:6])

sst_cov_weights <- find_smooth_weights(loc_list = sst_jul_loc_list,
                                   zonal_sdist = sst1_zonal_sdist,
                                   merid_sdist = sst1_merid_sdist,
                                   cov_func = matern_1.5_cov,
                                   mle_params = mle_param_mat,
                                   progress = FALSE)

sst_jul_smooth <- sst_mat %*% sst_cov_weights

sst_smooth_plot_dt <- cbind(sst_jul_coord_dt, 
                        unsmooth = sst_mat[33, ],
                        smooth = sst_jul_smooth[33, ])

sst_smooth_melt <- melt(sst_smooth_plot_dt,
                        id.vars = c("lon", "lat", "location"),
                        value.name = "anomaly",
                        variable.name = "smooth")

sst_smooth_melt$smooth <- factor(sst_smooth_melt$smooth,
                                 levels = c("unsmooth", "smooth"))

axis_labels_theme2 <- theme(axis.title.x = element_text(size = 12),
                           axis.text.x = element_text(size = 11),
                           axis.title.y = element_text(size = 12),
                           axis.text.y = element_text(size = 11),
                           legend.title = element_text(size = 10),
                           legend.text = element_text(size = 11),
                           legend.key.height = unit(1.4, "cm")) 

smooth_vs_unsmooth_ex_plot <- ggplot(data = land_map) +
  geom_tile(data = sst_smooth_melt,
            aes(x = lon,
                y = lat, 
                fill = anomaly)) +
  geom_sf(fill = 'gray', color = "gray", inherit.aes = FALSE) +
  coord_sf(xlim = c(0, 360), ylim = c(-60, 60), expand = FALSE) +
  scale_x_continuous(breaks = seq(0, 360, by = 50),
                     labels =  \(x) as.character(x)
  ) + 
  scale_y_continuous(breaks = seq(-60, 60, by = 20),
                     labels = \(x) as.character(x)) +
  scale_fill_continuous(name = "std. anom. [ ]",
                        type = "viridis") +
  facet_wrap(vars(smooth), nrow = 2) +
  theme_bw() +
  labs(x = "longitude", y = "latitude") + 
  axis_labels_theme2

smooth_vs_unsmooth_ex_plot

smoothing_example_path <- paste0(plots_dir, "sst_smooth_1991_example.png")

ggsave(smoothing_example_path,
       plot = smooth_vs_unsmooth_ex_plot,
       units = "in",
       width = 8, height = 5.5)
