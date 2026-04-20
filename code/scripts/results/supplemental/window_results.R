library(data.table)
library(knitr)

options(knitr.kable.NA = "-")

sim_dir <- "/home/data/projects/clim_smooth/sim_studies/two_sided/"
window12_path <- paste0(sim_dir, "sim_window12.csv")
window16_path <- paste0(sim_dir, "sim_window16.csv")

process_study <- function(sim_study) {
  study_mean <- sim_study[, lapply(.SD, mean), by = .(sigma,
                                                   n,
                                                   phi,
                                                   theta,
                                                   psi,
                                                   taper)]
  
  study_mean[, iteration := NULL]

  study_mean <- melt(study_mean, 
                     id.vars = c("sigma", "n", "phi", "theta", "psi", "taper"),
                     measure.vars = measure(method,
                                            value.name,
                                            pattern = "^([[:alpha:]]*)_(.*)$"))

  study_mean[, spec := NULL]
  study_mean[, FNP := NULL]
  study_mean[, reject := NULL]

  study_mean <- study_mean[order(sigma)]

  study_sd <- sim_study[, lapply(.SD, sd), 
  by = .(sigma,
         n,
         phi,
         theta,
         psi,
         taper)]

  study_sd[, iteration := NULL]

  study_sd <- melt(study_sd, 
                     id.vars = c("sigma", "n", "phi", "theta", "psi", "taper"),
                     measure.vars = measure(method,
                                            value.name,
                                            pattern = "^([[:alpha:]]*)_(.*)$"))


  study_sd[, spec := NULL]
  study_sd[, FNP := NULL]
  study_sd[, reject := NULL]

  study_sd$FDX <- study_sd$FDX / sqrt(10000)

  names(study_sd)[8:10] <- paste(names(study_sd)[8:10], "sd", sep="_")
  names(study_sd)[11] <- "FDX_se"

  study_sd <- study_sd[order(sigma)]

  joined_table <- study_mean[study_sd,
                             on = .(sigma, n, phi, theta, psi, taper, method)]
  
  return(joined_table)
}

window_results12 <- fread(window12_path)
window_results12$base_FDX <- window_results12$base_FDP > 0.1
window_results12$smooth_FDX <- window_results12$smooth_FDP > 0.1

window_results16 <- fread(window16_path)
window_results16$base_FDX <- window_results16$base_FDP > 0.1
window_results16$smooth_FDX <- window_results16$smooth_FDP > 0.1

processed_window12 <- process_study(window_results12)
processed_window16 <- process_study(window_results16)

processed_window12$sqrt_psi <- sqrt(processed_window12$psi)
processed_window16$sqrt_psi <- sqrt(processed_window16$psi)

processed_window12$window <- 12
processed_window16$window <- 16

processed_window <- rbind(processed_window12, processed_window16)

processed_window <- processed_window[, .(method, 
                                   window,
                                   sigma,
                                   n,
                                   phi,
                                   theta,
                                   sqrt_psi,
                                   taper,
                                   FDP,
                                   FDP_sd,
                                   FDX,
                                   FDX_se,
                                   true_reject,
                                   true_reject_sd,
                                   sens, 
                                   sens_sd)]

window_results_table <- processed_window

perc_delta_sens_vals <- ((window_results_table[method == "smooth"]$sens - 
                          window_results_table[method == "base"]$sens) / 
                          window_results_table[method == "base"]$sens) * 100

row_num <- nrow(window_results_table)
perc_delta_sens <- rep(NA, row_num)
perc_delta_sens[seq_len(row_num) %% 2 == 0] <- perc_delta_sens_vals

window_results_table$perc_delta_sens <- perc_delta_sens

window_results_table

window_results_table_1half <- window_results_table[, .(method, 
       window,
       sigma,
       n,
       phi,
       theta,
       sqrt_psi,
       taper,
       FDP,
       FDP_sd,
       FDX,
       FDX_se)]

window_results_table_2half <- window_results_table[, .(method,
       window,
       true_reject,
       true_reject_sd,
       sens,
       sens_sd,
       perc_delta_sens)]

kable(window_results_table_1half,
      format = "latex",
      digits = c(NA, 1, 1, 2, 2, 0, 2, 1, 3, 3, 3, 3),
      col.names = c("method",
                    "window",
                    "$\\sigma_\\varepsilon$",
                    "$n$",
                    "$\\phi$",
                    "$\\theta$",
                    "$\\sqrt{\\psi}$",
                    "$\\gamma$",
                    "FDP",
                    "s.d.",
                    "FDX",
                    "s.e."
      ),
      linesep = c("", "", "", "", "", "", "\\addlinespace",
                  "", "", "", "", "", "\\addlinespace",
                  "", "", "", "", "", "", "", "\\addlinespace",
                  "", "", "", ""),
      booktabs = TRUE,
      escape = FALSE
)

kable(window_results_table_2half,
      format = "latex",
      digits = c(NA, 1, 1, 2, 3, 3, 1),
      col.names = c("method",
                    "window",
                    "true rej.",
                    "s.d.",
                    "sens.",
                    "s.d.",
                    "$\\%\\Delta$sens."
      ),
      linesep = c("", "", "", "", "", "\\addlinespace",
                  "", "", "", "", "", "\\addlinespace",
                  "", "", "", "", "", "", "", "\\addlinespace",
                  "", "", "", ""),
      booktabs = TRUE,
      escape = FALSE
)
