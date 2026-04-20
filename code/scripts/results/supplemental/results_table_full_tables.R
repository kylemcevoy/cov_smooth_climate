## Old version of tables in results section before they were split in half.
library(data.table)
library(knitr)

sim_study_dir <- "/home/data/projects/clim_smooth/sim_studies/two_sided/"
noise_path <- paste0(sim_study_dir, "noise_sim_study_mean.csv")
n_path <- paste0(sim_study_dir, "n_sim_study_mean.csv")
phi_path <- paste0(sim_study_dir, "phi_sim_study_mean.csv")
aniso_path <- paste0(sim_study_dir, "aniso_sim_study_mean.csv")

noise_sim_study_mean <- fread(noise_path)
n_sim_study_mean <- fread(n_path)
phi_sim_study_mean <- fread(phi_path)
aniso_sim_study_mean <- fread(aniso_path)

results_table <- rbind(noise_sim_study_mean,
                       n_sim_study_mean,
                       phi_sim_study_mean,
                       aniso_sim_study_mean)

results_table[, FNP := NULL]
results_table[, reject := NULL]

results_table$sqrt_psi <- sqrt(results_table$psi)

results_table <- results_table[, .(method,
                                   sigma,
                                   n,
                                   phi,
                                   theta,
                                   sqrt_psi,
                                   taper,
                                   FDP,
                                   FDX,
                                   true_reject,
                                   sens)]


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
            after = "sens")

fwrite(output_results_table, 
       file = paste0(sim_study_dir, "final_results_table.csv"))

options(knitr.kable.NA = "-")

kable(output_results_table,
      format = "latex",
      digits = c(NA, 1, 2, 2, 0, 2, 1, 3, 3, 1, 3, 1),
      col.names = c("method",
                    "$\\sigma_\\varepsilon$",
                    "$n$",
                    "$\\phi$",
                    "$\\theta$",
                    "$\\sqrt{\\psi}$",
                    "$\\gamma$",
                    "FDP",
                    "FDX",
                    "true rej.",
                    "sens.",
                    "$\\%\\Delta$sens."
      ),
      linesep = c("", "", "", "", "", "\\addlinespace",
                  "", "", "", "", "", "\\addlinespace",
                  "", "", "", "", "", "", "", "\\addlinespace",
                  "", "", "", ""),
      booktabs = TRUE,
      escape = FALSE
)



sens_table <- results_table[, .(method,
                                sigma,
                                n,
                                phi,
                                theta,
                                sqrt_psi,
                                taper,
                                sens)]

wide_sens_table <- dcast(sens_table,
                         sigma + n + phi + theta + psi + taper ~ method,
                         value.var = "sens")

wide_sens_table[, perc_delta_sens := 100 * ((smooth - base) / base)]

sorted_wide_sens_table <- wide_sens_table[order(theta, psi, -phi, n, sigma)]

fwrite(sorted_wide_sens_table, 
       file = paste0(sim_study_dir, "sens_results_table.csv"))

kable(sorted_wide_sens_table,
      format = "latex",
      digits = c(1, 2, 2, 0, 2, 1, 4, 4, 2),
      col.names = c("$\\sigma_\\varepsilon$",
                    "$n$",
                    "$\\phi$",
                    "$\\theta$",
                    "$\\psi$",
                    "taper",
                    "base",
                    "smooth",
                    "\\% $\\Delta$ sens."
                    ),
      linesep = c("", "", "\\addlinespace",
                  "", "", "\\addlinespace",
                  "", "", "", "\\addlinespace",
                  "", ""),
      booktabs = TRUE,
      escape = FALSE
      )


