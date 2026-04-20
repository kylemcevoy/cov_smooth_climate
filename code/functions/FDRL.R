# implementing method of Zhang et al. (2011). Multiple Testing via FDR_L for 
# large scale imaging data.

get_loc_list2 <- function(coord_dt) {
  n <- nrow(coord_dt)
  loc_list <- vector(mode = "list", length = n)
  
  for (i in 1:n) {
    center_lon <- coord_dt[i]$lon
    center_lat <- coord_dt[i]$lat
    
    if (center_lon == 358) {
      subset_dt <- coord_dt[(lon == 356 & lat == center_lat) | 
                            (lon == 0 & lat == center_lat) | 
                            (lon == 358 & (abs(lat - center_lat) <= 2))]
    } else if (center_lon == 0) {
      subset_dt <- coord_dt[(lon == 358 & lat == center_lat) | 
                            (lon == 2 & lat == center_lat) |
                            (lon == 0 & (abs(lat - center_lat) <= 2))] 
    } else {
      subset_dt <- coord_dt[((lon == center_lon) & (abs(lat - center_lat) <= 2)) |
                            ((abs(lon - center_lon) <= 2) & (lat == center_lat))]
    }
    loc_list[[i]] <- subset_dt$location
  }
  
  return(loc_list)
}

find_pstar <- function(p_vec, loc_list) {
  p_star <- sapply(loc_list, \(x) median(p_vec[x]))
  return(p_star)
}

find_pmean <- function(p_vec, loc_list) {
  p_mean <- sapply(loc_list, \(x) mean(p_vec[x]))
  return(p_mean)
}

find_Gstar <- function(p_star, t) {
  G_denom <- 2 * sum(p_star > 0.5) + sum(p_star == 0.5)
  if (t <= 0.5) {
    G_star <- sum(p_star >= (1 - t)) / G_denom
  } else {
    G_star <- 1 - ((sum(p_star > t)) / G_denom)
  }
  return(G_star)
}

find_FDRL <- function(p_star, lambda, t) {
  R_star <- sum(p_star <= t)
  W_star <- sum(p_star > lambda)
  G_star_t <- find_Gstar(p_star, t)
  G_star_lambda <- find_Gstar(p_star, lambda)
  FDRL <- (W_star * G_star_t) / (max(R_star, 1) * (1 - G_star_lambda))
  if (t == 0) FDRL <- 0
  return(FDRL)
}

find_FDRL_threshold <- function(p_star, lambda, alpha) {
  fixed_FDRL <- \(x) find_FDRL(p_star = p_star, lambda = lambda, t = x)
  
  test_half <- seq(2 * alpha, 1, by = 0.1)
  FDRL_half <- sapply(test_half, fixed_FDRL)
  
  if (all(FDRL_half > (1.5 * alpha))) {
    seq1 <- seq(0, 2 * alpha, by = 1e-4)
  } else {
    seq1 <- seq(0, 1, by = 1e-4)
  }
  
  FDRL1 <- sapply(seq1, fixed_FDRL)
  bounded1 <- FDRL1 < alpha
  if (any(bounded1)) {
    largest_indx <- which(bounded1)[length(which(bounded1))]
    t_alpha <- seq1[largest_indx]
  } else {
    t_alpha <- 0
  }
  
  seq2 <- seq(t_alpha, t_alpha + 1e-4, by = 1e-6)
  FDRL2 <- sapply(seq2, fixed_FDRL)
  bounded2 <- FDRL2 < alpha
  if (any(bounded2)) {
    largest_indx <- which(bounded2)[length(which(bounded2))]
    t_alpha <- seq2[largest_indx]
  } else {
    t_alpha <- 0
  }
  
  return(t_alpha)
}




