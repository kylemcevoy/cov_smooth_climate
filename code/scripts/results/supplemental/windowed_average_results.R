library(data.table)
library(knitr)

save_dir <- "/home/data/projects/clim_smooth/sim_studies/two_sided/"

sim_sd1_path <- paste0(save_dir, "sim_pmean_n50_sd1_phi4_taper7.csv")
sim_sd2_path <- paste0(save_dir, "sim_pmean_n50_sd2_phi4_taper7.csv")
sim_sd3_path <- paste0(save_dir, "sim_pmean_n50_sd3_phi4_taper7.csv")

sim_sd1_small_path <- paste0(save_dir, "sim_pmean_small_n50_sd1_phi4_taper7.csv")
sim_sd2_small_path <- paste0(save_dir, "sim_pmean_small_n50_sd2_phi4_taper7.csv")
sim_sd3_small_path <- paste0(save_dir, "sim_pmean_small_n50_sd3_phi4_taper7.csv")

sim_sd1 <- fread(sim_sd1_path)
sim_sd2 <- fread(sim_sd2_path)
sim_sd3 <- fread(sim_sd3_path)

sim_sd1_small <- fread(sim_sd1_small_path)
sim_sd2_small <- fread(sim_sd2_small_path)
sim_sd3_small <- fread(sim_sd3_small_path)

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

  return(study_mean)
}

sim_study_small <- rbind(sim_sd1_small, sim_sd2_small, sim_sd3_small)
sim_mean_small <- process_study(sim_study_small)

sim_mean_small$window <- 4

sim_study_large <- rbind(sim_sd1, sim_sd2, sim_sd3)
sim_mean_large <- process_study(sim_study_large)

sim_mean_large$window <- 10

sim_mean_large
sim_mean_small

pmean_table <- rbind(sim_mean_small, sim_mean_large)

short_table <- pmean_table[, .(method, window, sigma, true_reject, FDP, sens)]

kable(short_table,
      format = "latex",
      digits = c(NA, 1, 1, 1, 3, 3),
      col.names = c("method", 
                    "window (deg.)",
                    "$\\sigma_\\varepsilon$",
                    "true rej.",
                    "FDP",
                    "sens."
      ),
      linesep = c("", "", "", "", "", "\\addlinespace",
                  "", "", "", "", "", "\\addlinespace"),
      booktabs = TRUE,
      escape = FALSE
)