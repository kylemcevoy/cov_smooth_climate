## functions for calculating distance and covariance matrices as well as
## processing covariance matrices for use in simulation functions

#### Finding Distances and Anisotropy Matrix ####

#' Find signed meridional distances
#'
#' function to find signed meridional distances between (lon, lat) coordinate pairs,
#' i.e. find the signed distance between the two pairs in just the N-S, latitudinal
#' direction. Used for anisotropic covariances where the relative signs of meridional and
#' zonal distances matter. Measured on a spherical approximation of earth with radius
#' 6378 km.
#'
#' @param coord_dt data.frame or data.table containing longitude and latitude
#' pairs for the data grid.
#' @returns matrix of signed meridional distances in units of 1 * 10^3 km.

find_signed_merid_dist <- function(coord_dt) {
  n <- nrow(coord_dt)
  lat_dist_unit <- 6378 * (pi / 180) 
  
  merid_dist_mat <- t(sapply(1:n, \(x) (coord_dt$lat[x] - coord_dt$lat)))
  merid_dist_mat <- (merid_dist_mat) * lat_dist_unit
  
  return(merid_dist_mat / 1000)
}

#' Find signed zonal distances
#'
#' function to find signed zonal distances between (lon, lat) coordinate pairs,
#' i.e. find the distance between the two pairs in just the E-W, longitudinal
#' direction. The distance is calculated along a sphere with radius 6378 km
#' using a cosine formula. The midpoint of the two input latitudes is used to
#' set the latitude at which the spherical distance is calculated. The sign of
#' the distance from the second point to the first is preserved for use in
#' anistropric covariance functions.
#'
#' @param coord_dt data.frame or data.table containing longitude and latitude
#' pairs for the data grid.
#' @returns matrix of signed zonal distances in 1e3 km.
#'

find_signed_zonal_dist <- function(coord_dt) {
  n <- nrow(coord_dt)

  lon_vec <- sort(unique(coord_dt$lon))
  lat_vec <- sort(unique(coord_dt$lat))

  lon_step <- lon_vec[2] - lon_vec[1]
  lat_length <- length(lat_vec)

  mid_lats <- apply(combn(lat_vec, 2), MARGIN = 2, FUN = mean)
  mid_lats <- sort(unique(c(mid_lats, lat_vec)))

  zonal_dist_mat <- matrix(0, n, n)

  lon_deg_step <- 6378 * pi / 180 * cos(pi / 180 * mid_lats)
  lon_seq <- seq(0, 180, by = lon_step)

  stored_dists <- outer(lon_deg_step, lon_seq)
  zonal_dist_mat <- matrix(0, n, n)

  for (i in 1:n) {
    tmp_lat <- coord_dt$lat[i]
    tmp_lon <- coord_dt$lon[i]

    tmp_mid_lats <- (tmp_lat + coord_dt$lat) / 2

    tmp_dist <- abs(tmp_lon - coord_dt$lon)
    lon_dists <- pmin(tmp_dist, 360 - tmp_dist)

    tmp_signs <- sign(tmp_lon - coord_dt$lon)
    sign_flip <- -2 * ((360 - tmp_dist) < tmp_dist) + 1
    tmp_signs <- tmp_signs * sign_flip

    lon_indx <- match(lon_dists, lon_seq)
    lat_indx <- match(tmp_mid_lats, mid_lats)

    zonal_dist_mat[i, ] <- tmp_signs * stored_dists[cbind(lat_indx, lon_indx)]
  }
  return(zonal_dist_mat / 1000)
}

#' Find anisotropy matrix
#'
#' given the anisotropy angle and anisotropy ratio calculates an anistropy matrix
#' for use in weighting the inner products of distances. The matrix is the product
#' of $R^T \Lambda R$ where R is a rotation matrix with angle aniso_angle (note
#' that R is the transpose of the usual formula for a rotation matrix with angle
#' $\theta$), and $\Lambda$ is a diagonal matrix with diagonal elements
#' (1, 1 / aniso_ratio).
#'
#' @param aniso_angle the angle of the rotation matrix. The angle given is the
#' negative of the angle in the usual rotation matrix parameterization (this is
#' equivalent to transposing the rotation matrix). This is done so that the angle
#' is the angle to rotate the original distance pairs to where they align with the
#' eigenvector directions.
#' @param aniso_ratio should be larger than 1. Gives the scaling factor with which
#' to divide the distances in the direction of the aniso_angle direction + $pi/2$.
#' Note that this means that the longer correlation length will be in the direction
#' of aniso_angle + $pi / 2$.
#' @returns a 2 x 2 matrix $R^T \Lambda R$ that rotates and scales distances for
#' use in covariance functions.

find_aniso_mat <- function(aniso_angle, aniso_ratio) {
  rotation_mat <- matrix(c(cos(aniso_angle),
                           -sin(aniso_angle),
                           sin(aniso_angle),
                           cos(aniso_angle)),
                         nrow = 2)
  scale_mat <- diag(c(1, 1 / aniso_ratio))

  return(t(rotation_mat) %*% scale_mat %*% rotation_mat)
}

##### Matern Covariance Functions #####

#' Calculate matern-1/2 covariance matrix
#'
#' Calculates Matern-1/2 covariance matrix on matrix of distances.
#' $\sigma^2 \exp(-dist / \phi)$
#'
#' @param dist matrix of distances between each grid box.
#' @param sigma2 the covariance scale of the function.
#' @param phi the correlation length scale parameter.
#' @returns the matern-1/2 covariance matrix.

matern_0.5_cov <- function(dist, sigma2 = 1, phi = 1) {
  sigma2 * exp(-dist / phi)
}

#' Calculate matern-3/2 covariance function
#'
#' Calculates Matern-3/2 covariance on vector or matrix of distances.
#' $\sigma^2 (1 + (\sqrt{3} dist / \phi)) \exp(-(\sqrt{3} dist / \phi)) + \tau^2$
#'
#' @param dist vector or matrix of distances between each grid box.
#' @param sigma2 the variance scale of the function.
#' @param phi the correlation length scale parameter.
#' @returns the matern-3/2 covariance function evaluated on the input distances.

matern_1.5_cov <- function(dist, sigma2 = 1, phi = 1) {
  sigma2 * (1 + (sqrt(3) * dist / phi)) * exp(-sqrt(3) * dist / phi)
}

#' Calculate matern-5/2 covariance function
#'
#' Calculates Matern-5/2 covariance on vector or matrix of distances.
#' $\sigma^2 (1 + (\sqrt{3} dist / \phi)) \exp(-(\sqrt{3} dist / \phi)) + \tau^2$
#'
#' @param dist vector or matrix of distances between each grid box.
#' @param sigma2 the variance scale of the function.
#' @param phi the correlation length scale parameter.
#' @returns the matern-5/2 covariance function evaluated on the input distances.
#' 
matern_2.5_cov <- function(dist, sigma2 = 1, phi = 1) {
  sigma2 * (1 + (sqrt(5) * dist / phi) + ((5 * dist^2) / (3 * phi^2))) * exp(-sqrt(5) * dist / phi)
}

spherical_taper <- function(dist, taper) {
  positive_part <- pmax(0, (1 - (dist / taper)))
  (positive_part)^2 * (1 + dist / (2 * taper))
}

wendland_1 <- function(dist, taper) {
  positive_part <- pmax(0, (1 - (dist / taper)))
  (positive_part^4) * (1 + 4 * (dist / taper))
}

wendland_2 <- function(dist, taper) {
  positive_part <- pmax(0, (1 - (dist / taper)))
  (positive_part^6) * (1 + 6 * (dist / taper) + (35 / 3) * (dist^2 / taper^2))
}

matern_0.5_spherical_cov <- function(dist, sigma2 = 1, phi = 1, taper) {
  cov_mat <- matern_0.5_cov(dist, sigma2 = sigma2, phi = phi)
  taper_mat <- spherical_taper(dist, taper = taper)
  cov_mat * taper_mat
}

matern_1.5_wendland_1_cov <- function(dist, sigma2 = 1, phi = 1, taper) {
  cov_mat <- matern_1.5_cov(dist, sigma2 = sigma2, phi = phi)
  taper_mat <- wendland_1(dist, taper=taper)
  cov_mat * taper_mat
}

matern_2.5_wendland_2_cov <- function(dist, sigma2 = 1, phi = 1, taper) {
  cov_mat <- matern_2.5_cov(dist, sigma2 = sigma2, phi = phi)
  taper_mat <- wendland_2(dist, taper=taper)
  cov_mat * taper_mat
}
