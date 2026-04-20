### Functions for reading and processing observational data for the project.
### Main variables used are 500 mb geopotential heights and SST.
### Pre-processed anomaly csv files for january and july should be the inputs.

## developed using data.table version 1.14.9

library(data.table)
#' Load SST data
#' 
#' function to load in and pre-process the STT csv files for use in the analysis.
#' The expected .csv files are for a time series of lat/lon gridded data for a 
#' particular month.
#' Each gridbox has been centered at 0. So the data is in anomalies. No trend
#' has been removed.
#' @param path_to_data file path to the directory containing the data file.
#' @param max_const_prop For SST some areas freeze over during winter which is
#' represented in the original dataset by a limiting value of -1.7 degrees C.
#' After centering, the values will change but this will result in many repeated
#' values in certain gridboxes. If the number of constant values exceeds the
#' max_const_prop threshold the gridbox is discarded from the analysis.
#' @returns a data.table with variables year, month, lon, lat, location, and var.

load_obs_data <- function(path, max_const_prop = 0.25) {
  ## specifying colClasses makes this step slower because
  ## it skips lazy loading of the columns.
  sst_class_vector <- c("integer", "Date", "double", "double", "double")

  obs_df <- fread(path,
                  sep = ",",
                  header = TRUE,
                  colClasses = sst_class_vector,
                  drop = 1)
    
  obs_df[, year := data.table::year(time)]
  obs_df[, month := data.table::month(time)]
  n_yr <- length(unique(obs_df$year))

  #the following lines removes all locations with constant data (sea ice?)
  #only a problem with january SST
  obs_df[, length_unique := length(unique(sst_anoms)), by = c("lon", "lat")]
  obs_df <- obs_df[(length_unique / n_yr)  >= (1 - max_const_prop)]
  obs_df[, length_unique := NULL]
  # the next line removes seas that are isolated from open ocean:
  # gulf of finland, caspian sea, black sea, hudson bay, great lakes. These
  # regions did not fit well using the MLE method and are not the focus of
  # the analysis.

  obs_df <- obs_df[!((lat >= 50 & lon >= 16 & lon <= 30) |
                     (lon >= 30 & lon <= 60 & lat >= 38) |
                     (lon >= 260 & lon <= 284 & lat >= 40))]


  setkey(obs_df, lon, lat, year)
  n_loc <- nrow(obs_df) / n_yr
  obs_df[, location :=  rep(1:n_loc, each = n_yr)]
  obs_df <- obs_df[, c("year", "month", "lon", "lat", "location", "sst_anoms")]

  return(obs_df)
}

extract_coord_dt <- function(obs_df) {
  first_year <- unique(obs_df$year)[1]
  year_sub <- obs_df[year == first_year]
  column_sub <- year_sub[, c("lon", "lat", "location")]
  return(column_sub)
}

extract_data_matrix <- function(obs_df, var, scale = TRUE) {
  ## expect data.table input to be sorted as outputted by load_obs_data()
  n_yr <- length(unique(obs_df$year))
  data_mat <- matrix(obs_df[[var]],
                     nrow = n_yr)
  if (scale) {
    data_mat <- scale(data_mat)
  }
  
  return(data_mat)
}