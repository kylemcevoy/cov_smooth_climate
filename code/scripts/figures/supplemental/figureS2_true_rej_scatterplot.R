# Script for generating histograms of the differences in number of true
# rejections for sample sizes n = 50 and n = 100.

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

noise_study_results$sigma <- factor(noise_study_results$sigma, levels = c(1, 2, 3))
noise_study_results$n <- factor(noise_study_results$n, levels = c(50, 100))

n_study_results$sigma <- factor(n_study_results$sigma, levels = c(1, 2, 3))
n_study_results$n <- factor(n_study_results$n, levels = c(50, 100))

sim_filtered <- rbind(noise_study_results, n_study_results)

facet_labels <- c("50" = "n = 50",
                  "100" = "n = 100")

reject_scatter <- ggplot(data = sim_filtered, 
  aes(x = base_true_reject, y = smooth_true_reject, color = sigma)) + 
  scale_fill_viridis_d(aesthetics = c('color'),
                       begin = 0.3,
                       end = 0.9,
                       option = 'G',
                       direction = -1) +
  facet_wrap("n") +
  geom_point(alpha = 0.7) + 
  geom_abline(slope = 1, color = "red", linetype = 2) + 
  theme_bw() +
  labs(x = "true rejections for base method",
       y = "true rejections for smooth method",
       color = expression(sigma[epsilon])) + 
  facet_wrap(vars(n), nrow = 1, labeller = as_labeller(facet_labels)) + 
  theme_bw() +
  coord_fixed(xlim = c(0, 1275), ylim = c(0, 1275)) +
  axis_labels_theme 

reject_scatter

ggsave(filename = paste0(plot_dir, "true_rej_scatterplot.png"), 
       plot = reject_scatter,
       units = "in",
       width = 10,
       height = 6)
