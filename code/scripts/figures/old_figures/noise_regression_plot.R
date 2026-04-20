# Noise Regression Figure
# Old figure, showing how regressing noise against SST looks with no FDR
# control. Saving code in case I need it.

source('code/functions/inference.R')

library(ggplot2)
library(sf)
library(maps)
library(pals)

axis_labels_theme <- theme(axis.title.x = element_text(size = 14),
                           axis.text.x = element_text(size = 12),
                           axis.title.y = element_text(size = 14),
                           axis.text.y = element_text(size = 12),
                           legend.title = element_text(size = 14),
                           legend.text = element_text(size = 12),
                           legend.key.height = unit(0.67, "cm"),
                           legend.key.width = unit(0.75, "cm")
) 

set.seed(3914)

y <- rnorm(n = 50, mean = 0, sd = 1)
t_stat_vec <- fit_lm_tstat(y, sst_mat_sub)

p_vals <- pt(t_stat_vec, df = 48, lower.tail = FALSE)

rejects <- p_vals < 0.05

land_map <- sf::st_as_sf(map(wrap = c(0, 360), plot = FALSE, fill = TRUE))

noise_reg <- cbind(sst_jan_coord_dt, t_stat = t_stat_vec, reject = rejects)


noise_reg_plot <- ggplot(data = noise_reg, aes(x = lon, y = lat)) + 
  geom_tile(aes(fill = t_stat)) + 
  geom_sf(data = land_map,
          fill = 'gray',
          color = "gray",
          inherit.aes = FALSE) +
  coord_sf(xlim = c(0, 360), ylim = c(-60, 60), expand = FALSE) +
  scale_fill_viridis_c(name = "t stat.") + 
  geom_point(data = noise_reg[reject == TRUE],
             color = "red",
             size = 0.2,
             shape = 3) + 
  scale_x_continuous(breaks = seq(0, 340, by = 30),
                     labels =  \(x) as.character(x)
  ) + 
  scale_y_continuous(breaks = seq(-60, 60, by = 20),
                     labels = \(x) as.character(x)) +
  theme_bw() +
  labs(x = "longitude",
       y = "latitude") +
  axis_labels_theme

noise_reg_plot

ggsave("images/noise_regression.png", width = 8, height = 3, units = "in")
