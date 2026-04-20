# old version of the rejection histogram.

library(ggplot2)
library(data.table)

axis_labels_theme <- theme(axis.title.x = element_text(size = 12),
                           axis.text.x = element_text(size = 11),
                           axis.title.y = element_text(size = 12),
                           axis.text.y = element_text(size = 11),
                           legend.title = element_text(size = 11, vjust = 0.8),
                           legend.text = element_text(size = 11)) 

plot_dir <- "/home/data/projects/clim_smooth/plots/"
data_dir <- "/home/data/projects/clim_smooth/sim_studies/two_sided/"

noise_study_path <- paste0(data_dir, "noise_sim_study.csv")
n_study_path <- paste0(data_dir, "n_sim_study.csv")

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

reject_hist_plot

ggsave(paste0(plot_dir, 'rej_hist_plot.png'), reject_hist_plot,
       width = 4,
       height = 5.1,
       units = 'in')

sum(sim_filtered[n == 50]$true_reject_diff > 450)
# 5
sum(sim_filtered[n == 100]$true_reject_diff > 450)
# 4

sum(sim_filtered[n == 50]$true_reject_diff < -175)
# 0
sum(sim_filtered[n == 100]$true_reject_diff < -175)
# 0 

ggsave(filename = paste0(plot_dir, "reject_hist_plot.png"), 
       plot = reject_hist_plot,
       units = "in",
       width = 6,
       height = 6)

gray_palette <- hcl.colors(n = 3, palette = "Light Grays")
names(gray_palette) <- c("1", "2", "3")

reject_hist_plot_gray <- ggplot(data = sim_filtered, aes(x = true_reject_diff,
                                                    fill = sigma)) +
  geom_histogram(color = "black", binwidth = 10, position = "identity", alpha = 0.7) +
  geom_vline(data = mean_frame,
             aes(xintercept = zero),
             color = "black",
             linetype = 2) +
  scale_color_manual(values = gray_palette, aesthetics = c("color", "fill")) + 
  geom_vline(data = mean_frame,
             aes(xintercept = mean, color = sigma),
             linetype = 2) +
  labs(x = "difference in # of true rejections per sim. (smoothed - base)",
       color = expression(sigma[epsilon]),
       fill = expression(sigma[epsilon])) + 
  facet_wrap(vars(n), nrow = 2, labeller = as_labeller(facet_labels)) + 
  coord_cartesian(xlim = c(-175, 450)) + 
  theme_bw() +
  axis_labels_theme

reject_hist_plot_gray

ggsave(filename = paste0(plot_dir, "reject_hist_plot_gray.png"), 
       plot = reject_hist_plot_gray,
       units = "in",
       width = 8,
       height = 6)
