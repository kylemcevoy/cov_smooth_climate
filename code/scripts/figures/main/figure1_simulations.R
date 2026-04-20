# Script to generate Figure 1 of the paper.
# First, histograms of the differences in number of true
# rejections for sample sizes n = 50 and n = 100 are generated.
# Second, the selection stability figure is generated that shows at a fixed
# location how rejection frequencies compare between base and smooth methods.

library(data.table)
library(viridis)
library(ggplot2)
library(patchwork)

source('code/functions/sim.R')

#### Rejection Diff. Histogram

plot_dir <- "/home/data/projects/clim_smooth/plots/"
data_dir <- "/home/data/projects/clim_smooth/sim_studies/two_sided/"

noise_study_path <- paste0(data_dir, "noise_sim_study.csv")
n_study_path <- paste0(data_dir, "n_sim_study.csv")

axis_labels_theme <- theme(axis.title.x = element_text(size = 12),
                           axis.text.x = element_text(size = 11),
                           axis.title.y = element_text(size = 12),
                           axis.text.y = element_text(size = 11),
                           legend.title = element_text(size = 11, vjust = 0.8),
                           legend.text = element_text(size = 11)) 

# n = 50
noise_study_results <- fread(noise_study_path)
# n = 100
n_study_results <- fread(n_study_path)

filter_nonzero <- function(data) {
  filter_data <- data[data$base_reject > 0 | data$smooth_reject > 0]
  filter_data$true_reject_diff <- (filter_data$smooth_true_reject -
                                     filter_data$base_true_reject)
  return(filter_data)
}

noise_sim_filtered <- filter_nonzero(noise_study_results)
n_sim_filtered <- filter_nonzero(n_study_results)

noise_sim_filtered$sigma <- factor(noise_sim_filtered$sigma, levels = c(1, 2, 3))
noise_sim_filtered$n <- factor(noise_sim_filtered$n, levels = c(50, 100))

n_sim_filtered$sigma <- factor(n_sim_filtered$sigma, levels = c(1, 2, 3))
n_sim_filtered$n <- factor(n_sim_filtered$n, levels = c(50, 100))

sim_filtered <- rbind(noise_sim_filtered, n_sim_filtered)

mean_frame <- sim_filtered[, .(mean = mean(true_reject_diff)), by = .(sigma, n)]
mean_frame$zero <- 0

facet_labels <- c("50" = "n = 50",
                  "100" = "n = 100")

reject_hist_plot <- ggplot(data = sim_filtered, aes(x = true_reject_diff,
                                  fill = sigma)) +
  geom_histogram(color = "black", binwidth = 10, position = "identity", alpha = 0.7) +
  geom_vline(data = mean_frame,
             aes(xintercept = zero),
                 color = "black",
             linetype = 2) +
  geom_vline(data = mean_frame,
             aes(xintercept = mean, color = sigma),
             linetype = 2) +
  scale_fill_viridis_d(aesthetics = c('color', 'fill'),
                       begin = 0.3,
                       end = 0.9,
                       option = 'G',
                       direction = -1) +
  labs(x = "difference in true rejections",
       color = expression(sigma[epsilon]),
       fill = expression(sigma[epsilon])) + 
  facet_wrap(vars(n), nrow = 2, labeller = as_labeller(facet_labels)) + 
  coord_cartesian(xlim = c(-140, 350)) + 
  theme_bw() +
  axis_labels_theme +
  guides(fill = guide_legend(position = 'bottom', label.position = 'bottom'),
         color = guide_legend(position = 'bottom', label.position = 'bottom'))



##### Selection Stability Plot ##### 

###### Functions #######
# this function finds convex hulls of sets of points. Used to create hulls for
# regions above a certain covariance in the plot.
find_chull <- function(coord_dt, condition) {
  map_back <- which(condition)
  convex_hull <- chull(coord_dt[condition == TRUE, ])
  convex_hull <- c(convex_hull, convex_hull[1])
  orig_indices <- map_back[convex_hull]
  chull_orig_data <- coord_dt[orig_indices, ]
  
  return(chull_orig_data)
}


##### Data ######

data_path <- paste0('/home/data/projects/clim_smooth/sim_studies/two_sided/',
                    'selection_stability.csv')

plot_path <- '/home/data/projects/clim_smooth/plots/selection_stability.png'
combined_path <- '/home/data/projects/clim_smooth/plots/combined_fig1.png'

selection_stability_dt <- fread(data_path)

# location variable stores center location for the simulation.
center_loc <- selection_stability_dt$location[1]

true_cov <- selection_stability_dt$true_cov

lat_lon_grid <- build_lat_lon_grid(lat_range = c(-60, 60),
                                   lon_range = c(0, 358),
                                   lat_unit = 2,
                                   lon_unit = 2)

selection_stability_dt <- cbind(selection_stability_dt[, -c("location")],
                                lat_lon_grid)

##### Plotting ######

beta <- 1
test_chull <- find_chull(lat_lon_grid, !(true_cov * beta < 1e-3))
test_chull2 <- find_chull(lat_lon_grid, !(true_cov * beta < 1e-2))
test_chull3 <- find_chull(lat_lon_grid, !(true_cov * beta < 1e-1))

test_chull$corr <- 1e-3
test_chull2$corr <- 1e-2
test_chull3$corr <- 1e-1

test_chull_comb <- rbind(test_chull, test_chull2, test_chull3)

axis_labels_theme <- theme(axis.title.x = element_text(size = 13),
                           axis.text.x = element_text(size = 11),
                           axis.title.y = element_text(size = 13),
                           axis.text.y = element_text(size = 11),
                           legend.title = element_text(size = 12),
                           legend.text = element_text(size = 11),
                           legend.key.width = unit(1, "cm")) 


(selection_stability_plot_virid <- ggplot(data = selection_stability_dt,
       aes(x = lon, y = lat)) + 
  geom_tile(aes(fill = reject_freq_diff)) +
  geom_point(data = lat_lon_grid[center_loc, ], color = "red", alpha = 0.3) +
  # scale_fill_stepsn(colors = rev(diverging_hcl(20,
  #                                              h = c(180, 50),
  #                                              c = 80,
  #                                              l = c(20, 95),
  #                                              power = c(0.7, 1.3))),
  #                   n.breaks = 8,
  #                   limits = c(-0.25, 0.25),
  #                   name = "reject freq. diff.") +
  scale_fill_viridis_b(limits = c(-0.05, 0.3),
                       breaks = seq(-0.05, 0.3, by = 0.05),
                       name = "reject freq. diff.") +
  geom_path(data = test_chull, color = "red") +
  geom_path(data = test_chull2, color = "purple") +
  geom_path(data = test_chull3, color = "blue")+ 
  coord_sf(xlim = c(96, 200), ylim = c(-50, 50)) + 
  labs(x = "longitude", y = "latitude") +
  theme_bw() +
  axis_labels_theme + 
  guides(fill = guide_bins(
    title.position = "top",
    title.hjust = 0.5,
    position = "bottom")) +
    theme(legend.axis.line = element_blank()))

combined_figure_1 <- ((reject_hist_plot + selection_stability_plot_virid) + 
                        plot_annotation(tag_levels = 'a', tag_suffix = ')'))

ggsave(plot_path, selection_stability_plot_virid, width = 4, height = 5.1, units = 'in')

ggsave(combined_path, combined_figure_1, width = 8, height = 5.1, units = 'in')

selection_stability_plot_gray <- ggplot(data = selection_stability_dt,
                                   aes(x = lon, y = lat)) + 
  geom_tile(aes(fill = reject_freq_diff)) +
  geom_point(data = lat_lon_grid[center_loc, ], color = "red", alpha = 0.3) +
  scale_fill_stepsn(colors = colorspace::desaturate(rev(diverging_hcl(20,
                                               h = c(180, 50),
                                               c = 80,
                                               l = c(20, 95),
                                               power = c(0.7, 1.3)))),
                    n.breaks = 8,
                    limits = c(-0.25, 0.25),
                    name = "reject freq. diff.") +
  geom_path(data = test_chull, color = "red") +
  geom_path(data = test_chull2, color = "purple") +
  geom_path(data = test_chull3, color = "blue")+ 
  coord_sf(xlim = c(96, 200), ylim = c(-50, 50)) + 
  labs(x = "longitude", y = "latitude") +
  theme_bw() +
  axis_labels_theme


  
