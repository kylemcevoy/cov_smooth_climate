# Script to generate Figure 4 of the paper showing the results in the 
# regression of Central US temperatures in July on July SST anomalies.

## Data Manipulation
library(data.table)
## Plotting Packages
library(maps)
library(patchwork) # stacking ggplots
library(ggplot2)

source("code/functions/process_data.R")

data_dir <- '/home/data/projects/clim_smooth/'
plot_dir <- '/home/data/projects/clim_smooth/plots/'

land_map <- sf::st_as_sf(map(wrap = c(0, 360), plot = FALSE, fill = TRUE))

t2m_rejections <- fread(paste0(data_dir, 't2m_sst_rejections_jul.csv'))

t2m_rejections$method <- factor(t2m_rejections$method,
                                levels = c('no BH', 'base', 'smooth'),
                                labels = c('no FDR control', 'base', 'smooth'))

sst_jul <- load_obs_data(path = paste0(data_dir, "sst_detrend_anom_jul.csv"))
sst_jul_coord_dt <- extract_coord_dt(sst_jul)

t2m_rejections <- t2m_rejections[sst_jul_coord_dt, on = "location"]

no_fdr_rejections <- t2m_rejections[method == 'no FDR control']
base_rejections <- t2m_rejections[method == 'base']
smooth_rejections <- t2m_rejections[method == 'smooth']

no_fdr_plot <- ggplot(data = land_map) +
  geom_sf(fill = "gray",
          color = "gray",
          inherit.aes = FALSE) +
  geom_point(data = no_fdr_rejections[reject == TRUE],
             size = 0.4,
             color = '#5A5A5A',
             aes(x = lon, y = lat)) +
  coord_sf(xlim = c(0, 360),
           ylim = c(-60, 60),
           expand = FALSE) +
  scale_x_continuous(breaks = seq(0, 360, by = 50),
                     labels =  \(x) as.character(x)) +
  scale_y_continuous(breaks = seq(-60, 60, by = 20),
                     labels = \(x) as.character(x)) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  theme_bw() +
  labs(title = 'a) no FDR control') + 
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        title = element_text(size = 8))

base_plot <- ggplot(data = land_map) +
  geom_sf(fill = "gray",
          color = "gray",
          inherit.aes = FALSE) +
  geom_point(data = base_rejections[reject == TRUE],
             size = 0.4,
             color = '#5A5A5A',
             aes(x = lon, y = lat)) +
  coord_sf(xlim = c(0, 360),
           ylim = c(-60, 60),
           expand = FALSE) +
  scale_x_continuous(breaks = seq(0, 360, by = 50),
                     labels =  \(x) as.character(x)) +
  scale_y_continuous(breaks = seq(-60, 60, by = 20),
                     labels = \(x) as.character(x)) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  theme_bw() +
  labs(title = 'b) base') + 
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        title = element_text(size = 8))

smooth_plot <- ggplot(data = land_map) +
  geom_sf(fill = "gray",
          color = "gray",
          inherit.aes = FALSE) +
  geom_point(data = smooth_rejections[reject == TRUE],
             size = 0.4,
             color = '#5A5A5A',
             aes(x = lon, y = lat)) +
  coord_sf(xlim = c(0, 360),
           ylim = c(-60, 60),
           expand = FALSE) +
  scale_x_continuous(breaks = seq(0, 360, by = 50),
                     labels =  \(x) as.character(x)) +
  scale_y_continuous(breaks = seq(-60, 60, by = 20),
                     labels = \(x) as.character(x)) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  labs(title='c) smooth') + 
  theme_bw() +
  theme(title = element_text(size = 8))

patchwork_stacked_plot <- (no_fdr_plot / base_plot / smooth_plot)

ggsave(filename = paste0(plot_dir, 'july_GP_t2m_against_sst_stacked_patch.png'),
       plot = patchwork_stacked_plot)

