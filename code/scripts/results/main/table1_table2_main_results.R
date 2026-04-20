## Tables 1 and 2 of the results section
library(data.table)
library(knitr)

sim_study_dir <- "/home/data/projects/clim_smooth/sim_studies/two_sided/"
noise_path <- paste0(sim_study_dir, "noise_sim_study.csv")
n_path <- paste0(sim_study_dir, "n_sim_study.csv")
phi_path <- paste0(sim_study_dir, "phi_sim_study.csv")
aniso_path <- paste0(sim_study_dir, "aniso_sim_study.csv")

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

noise_sim_study <- fread(noise_path)
n_sim_study <- fread(n_path)
phi_sim_study <- fread(phi_path)
aniso_sim_study <- fread(aniso_path)

noise_table <- process_study(noise_sim_study)
n_table <- process_study(n_sim_study)
phi_table <- process_study(phi_sim_study)
aniso_table <- process_study(aniso_sim_study)

results_table <- rbind(noise_table,
                       n_table,
                       phi_table,
                       aniso_table)

results_table$sqrt_psi <- sqrt(results_table$psi)

results_table <- results_table[, .(method,
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


output_results_table <- results_table[order(sqrt_psi, theta, -phi, taper)]

perc_delta_sens_vals <- ((output_results_table[method == "smooth"]$sens - 
                          output_results_table[method == "base"]$sens) / 
                          output_results_table[method == "base"]$sens) * 100

row_num <- nrow(output_results_table)
perc_delta_sens <- rep(NA, row_num)
perc_delta_sens[seq_len(row_num) %% 2 == 0] <- perc_delta_sens_vals

output_results_table$perc_delta_sens <- perc_delta_sens

setcolorder(output_results_table,
            neworder = "perc_delta_sens",
            after = "sens_sd")

fwrite(output_results_table, 
       file = paste0(sim_study_dir, "final_results_table.csv"))

options(knitr.kable.NA = "-")

output_results_table_1half <- output_results_table[, .(method,
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

output_results_table_2half <- output_results_table[, .(method,
       sigma,
       n,
       phi,
       theta,
       sqrt_psi,
       taper,
       true_reject,
       true_reject_sd,
       sens,
       sens_sd,
       perc_delta_sens)]

kable(output_results_table_1half,
      format = "latex",
      digits = c(NA, 1, 2, 2, 0, 2, 1, 3, 3, 3, 3),
      col.names = c("method",
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
      linesep = c("", "", "", "", "", "\\addlinespace",
                  "", "", "", "", "", "\\addlinespace",
                  "", "", "", "", "", "", "", "\\addlinespace",
                  "", "", "", ""),
      booktabs = TRUE,
      escape = FALSE
)

kable(output_results_table_2half,
      format = "latex",
      digits = c(NA, 1, 2, 2, 0, 2, 1, 1, 1, 3, 3, 3),
      col.names = c("method",
                    "$\\sigma_\\varepsilon$",
                    "$n$",
                    "$\\phi$",
                    "$\\theta$",
                    "$\\sqrt{\\psi}$",
                    "$\\gamma$",
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

