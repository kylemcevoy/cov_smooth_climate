# Script to generate variograms that were used to estimate SST length-scale
# parameters that were used in the simulated data.
library(data.table)
library(ggplot2)
library(viridis)

source("code/functions/dist_cov.R")
source("code/functions/process_data.R")

data_path <- "/home/data/projects/clim_smooth/sst_detrend_anom_jul.csv"
plot_dir <- "/home/data/projects/clim_smooth/plots/"

# Functions ---------------------------------------------------------------

matern_0.5_taper_vario <- function(dist, sigma2, phi, taper) {
  1 - matern_0.5_spherical_cov(dist = dist,
                               sigma2 = sigma2,
                               phi = phi,
                               taper = taper)
}

matern_1.5_taper_vario <- function(dist, sigma2, phi, taper){
  1 - matern_1.5_wendland_1_cov(dist = dist,
                                sigma2 = sigma2,
                                phi = phi,
                                taper = taper)
}

matern_2.5_taper_vario <- function(dist, sigma2, phi, taper){
  1 - matern_2.5_wendland_2_cov(dist = dist,
                                sigma2 = sigma2,
                                phi = phi,
                                taper = taper)
}

calc_zonal_bin_vario <- function(coord_dt,
                                 band_lat,
                                 zonal_dist,
                                 sst_mat,
                                 bins = 0:10) {
  zon_band_locs <- coord_dt[lat == band_lat]$location
  zon_band_dist <- zonal_dist[zon_band_locs, zon_band_locs]
  
  n_bins <- length(bins) - 1
  bin_counts <- rep(0, n_bins)
  mean_bin_vario <- rep(0, n_bins)
  for (i in 1:n_bins) {
    tmp_dist <- 0
    for (j in 1:length(zon_band_locs)) {
      in_bin_locs <-
        (abs(zon_band_dist[j,]) > bins[i]) &
        (abs(zon_band_dist[j,]) < bins[i + 1])
      bin_counts[i] <- bin_counts[i] + sum(in_bin_locs)
      if (sum(in_bin_locs) == 0) {
        next
      } else {
        bin_locs <- zon_band_locs[in_bin_locs]
        for (k in 1:length(bin_locs)) {
          tmp_dist <- tmp_dist + mean((sst_mat[, zon_band_locs[j]] -
                                         sst_mat[, bin_locs[k]]) ^
                                        2)
        }
      }
    }
    mean_bin_vario[i] <- tmp_dist
  }
  return(mean_bin_vario / (2 * bin_counts))
}

calc_merid_bin_vario <- function(coord_dt,
                                 band_lon,
                                 merid_dist,
                                 sst_mat,
                                 bins = 0:10) {
  merid_band_locs <- coord_dt[lon == band_lon]$location
  merid_band_dist <- merid_dist[merid_band_locs, merid_band_locs]
  
  n_bins <- length(bins) - 1
  bin_counts <- rep(0, n_bins)
  mean_bin_vario <- rep(0, n_bins)
  for (i in 1:n_bins) {
    tmp_dist <- 0
    for (j in 1:length(merid_band_locs)) {
      in_bin_locs <-
        (abs(merid_band_dist[j,]) > bins[i]) &
        (abs(merid_band_dist[j,]) < bins[i + 1])
      bin_counts[i] <- bin_counts[i] + sum(in_bin_locs)
      if (sum(in_bin_locs) == 0) {
        next
      } else {
        bin_locs <- merid_band_locs[in_bin_locs]
        for (k in 1:length(bin_locs)) {
          tmp_dist <- tmp_dist + mean((sst_mat[, merid_band_locs[j]] -
                                         sst_mat[, bin_locs[k]]) ^
                                        2)
        }
      }
    }
    mean_bin_vario[i] <- tmp_dist
  }
  return(mean_bin_vario / (2 * bin_counts))
}

# Load Data ---------------------------------------------------------------

sst_jul <- load_obs_data(path = data_path,
                         max_const_prop = 0.25)

sst_jul_coord_dt <- extract_coord_dt(sst_jul)

sst_zonal_sdist <- find_signed_zonal_dist(sst_jul_coord_dt)
sst_merid_sdist <- find_signed_merid_dist(sst_jul_coord_dt)

sst_mat <- extract_data_matrix(sst_jul,
                               var = "sst_anoms",
                               scale = TRUE)


# Fit Variograms ----------------------------------------------------------

zonal_bands <- 10 * seq(-6, 6, by=1)
zonal_vario_list <- list(length(zonal_bands))
zonal_bins <- 0:10
zonal_bin_midpoints <- zonal_bins[-length(zonal_bins)] + diff(zonal_bins) / 2

merid_bands <- 10 * seq(16, 24, by=1)
merid_vario_list <- list(length(merid_bands))
merid_bins <- 0:10
merid_bin_midpoints <- merid_bins[-length(merid_bins)] + diff(merid_bins) / 2

for (i in seq_along(zonal_bands)) {
  zonal_vario_list[[i]] <- calc_zonal_bin_vario(coord_dt = sst_jul_coord_dt,
                                                band_lat = zonal_bands[i],
                                                zonal_dist = sst_zonal_sdist,
                                                sst_mat = sst_mat,
                                                bins = zonal_bins)
  }

zonal_vario_mat <- t(matrix(unlist(zonal_vario_list), ncol = 13))
zonal_vario_dt <- data.table(zonal_vario_mat)
names(zonal_vario_dt) <- as.character(zonal_bin_midpoints)
zonal_vario_dt$lat_band <- zonal_bands

zonal_vario_melted <- melt(zonal_vario_dt,
                           id.vars = "lat_band",
                           variable.name = "distance",
                           variable.factor = FALSE,
                           value.name = "variogram")

zonal_vario_melted$distance <- as.numeric(zonal_vario_melted$distance)

for (i in seq_along(merid_bands)) {
  merid_vario_list[[i]] <- calc_merid_bin_vario(coord_dt = sst_jul_coord_dt,
                                                band_lon = merid_bands[i],
                                                merid_dist = sst_merid_sdist,
                                                sst_mat = sst_mat,
                                                bins = merid_bins)
}

merid_vario_mat <- t(matrix(unlist(merid_vario_list), ncol = 9))
merid_vario_dt <- data.table(merid_vario_mat)
names(merid_vario_dt) <- as.character(merid_bin_midpoints)
merid_vario_dt$lon_band <- merid_bands

merid_vario_melted <- melt(merid_vario_dt,
                           id.vars = "lon_band",
                           variable.name = "distance",
                           variable.factor = FALSE,
                           value.name = "variogram")

merid_vario_melted$distance <- as.numeric(merid_vario_melted$distance)


# Plotting ----------------------------------------------------------------
latitude_colors <- viridis(5)
names(latitude_colors) <- c('-30', '-20', '20', '30', 'Cov. Func.')
  
low_lat <- ggplot(data = zonal_vario_melted[abs(lat_band) > 15 & abs(lat_band) < 35],
       aes(x = distance,
           y = variogram)) + 
  geom_point(aes(color = factor(lat_band))) + 
  geom_line(aes(color = factor(lat_band)), alpha = 0.5) + 
  geom_function(fun = matern_1.5_taper_vario,
                args = list(phi = 4, taper = 7, sigma2 = 1),
                linetype = 'dashed',
                aes(color = 'Cov. Func.')) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  scale_color_manual(values = latitude_colors) + 
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  scale_y_continuous(breaks = seq(0, 1.2, by = 0.2)) + 
  coord_cartesian(xlim = c(0, 10)) + 
  labs(color = "latitude band") + 
  theme_bw()

ggsave(low_lat, filename = paste0(plot_dir, 'lowlat_variogram.png'))

midlat_colors <- viridis(5)
names(midlat_colors) <- c('-50', '-40', '40', '50', 'Cov. Func.')

mid_lat <- ggplot(data = zonal_vario_melted[abs(lat_band) > 35 & abs(lat_band) < 55],
       aes(x = distance,
           y = variogram)) + 
  geom_point(aes(color = factor(lat_band))) + 
  geom_line(aes(color = factor(lat_band)), alpha = 0.5) + 
  geom_function(fun = matern_1.5_taper_vario, 
                args = list(phi = 4, taper = 7, sigma2 = 1),
                aes(color = 'Cov. Func.'),
                linetype = 'dashed') +
  scale_color_manual(values = midlat_colors) + 
  geom_hline(yintercept = 1, linetype = 'dashed') +
  scale_x_continuous(breaks = seq(0, 12, by = 2)) +
  scale_y_continuous(breaks = seq(0, 1.2, by = 0.2)) + 
  coord_cartesian(xlim = c(0, 10)) + 
  labs(color = "latitude band") + 
  theme_bw()

ggsave(mid_lat, filename = paste0(plot_dir, 'midlat_variogram.png'))

longitude_colors <- viridis::viridis(6)
names(longitude_colors) <- c('190', '200', '210', '220', 'Taper at 5', 'Taper at 7')

longitude <- ggplot(data = merid_vario_melted[lon_band > 185 & lon_band < 230],
       aes(x = distance,
           y = variogram)) + 
  geom_point(aes(color = factor(lon_band))) + 
  geom_line(aes(color = factor(lon_band)), alpha = 0.5) + 
  geom_function(fun = matern_1.5_taper_vario, 
                args = list(phi = 1.3, taper = 7, sigma2 = 1),
                aes(color = 'Taper at 7'),
                linetype = 'dashed') +
  geom_function(fun = matern_1.5_taper_vario, args = list(phi = 1.3, taper = 5, sigma2 = 1),
                aes(color = 'Taper at 5'),
                linetype = 'dashed') +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  scale_color_manual(values = longitude_colors) + 
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  scale_y_continuous(breaks = seq(0, 1.2, by = 0.2)) + 
  coord_cartesian(xlim = c(0, 10)) + 
  labs(color = "longitude band",
       x = 'distance [10^3 km]') + 
  theme_bw()

ggsave(longitude, filename = paste0(plot_dir, 'longitude_variogram.png'))


# phi * sqrt(psi) should equal the longer length scale and phi should equal the
# shorter length scale.
psi = (4 / 1.3)^2
psi

##### Matern 2.5 ######

ggplot(data = zonal_vario_melted[abs(lat_band) > 35 & abs(lat_band) < 55],
       aes(x = distance,
           y = variogram)) + 
  geom_point(aes(color = factor(lat_band))) + 
  geom_line(aes(color = factor(lat_band)), alpha = 0.5) + 
  geom_function(fun = matern_2.5_taper_vario, args = list(phi = 4, taper = 7, sigma2 = 1),
                color = 'red', linetype = 'dashed') +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  scale_x_continuous(breaks = seq(0, 12, by = 2)) +
  scale_y_continuous(breaks = seq(0, 1.2, by = 0.2)) + 
  coord_cartesian(xlim = c(0, 10)) + 
  labs(color = "latitude band") + 
  theme_bw()

ggplot(data = zonal_vario_melted[abs(lat_band) > 15 & abs(lat_band) < 35],
       aes(x = distance,
           y = variogram)) + 
  geom_point(aes(color = factor(lat_band))) + 
  geom_line(aes(color = factor(lat_band)), alpha = 0.5) + 
  geom_function(fun = matern_1.5_taper_vario, args = list(phi = 4, taper = 7, sigma2 = 1),
                color = 'red', linetype = 'dashed') +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  scale_y_continuous(breaks = seq(0, 1.2, by = 0.2)) + 
  coord_cartesian(xlim = c(0, 10)) + 
  labs(color = "latitude band") + 
  theme_bw()

##### Matern 0.5 ######

ggplot(data = zonal_vario_melted[abs(lat_band) > 35 & abs(lat_band) < 55],
       aes(x = distance,
           y = variogram)) + 
  geom_point(aes(color = factor(lat_band))) + 
  geom_line(aes(color = factor(lat_band)), alpha = 0.5) + 
  geom_function(fun = matern_0.5_taper_vario, args = list(phi = 4, taper = 7, sigma2 = 1),
                color = 'red', linetype = 'dashed') +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  scale_x_continuous(breaks = seq(0, 12, by = 2)) +
  scale_y_continuous(breaks = seq(0, 1.2, by = 0.2)) + 
  coord_cartesian(xlim = c(0, 10)) + 
  labs(color = "latitude band") + 
  theme_bw()

ggplot(data = zonal_vario_melted[abs(lat_band) > 15 & abs(lat_band) < 35],
       aes(x = distance,
           y = variogram)) + 
  geom_point(aes(color = factor(lat_band))) + 
  geom_line(aes(color = factor(lat_band)), alpha = 0.5) + 
  geom_function(fun = matern_0.5_taper_vario, args = list(phi = 4, taper = 7, sigma2 = 1),
                color = 'red', linetype = 'dashed') +
  geom_function(fun = matern_1.5_taper_vario, args = list(phi = 4, taper = 7, sigma2 = 1),
                color = 'gray', linetype = 'dashed') +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  scale_y_continuous(breaks = seq(0, 1.2, by = 0.2)) + 
  coord_cartesian(xlim = c(0, 10)) + 
  labs(color = "latitude band") + 
  theme_bw()

