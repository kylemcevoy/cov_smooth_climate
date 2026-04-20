# Results for Table 3 of the paper on misspecification.

library(data.table)
library(knitr)

sim_study_dir <- "/home/data/projects/clim_smooth/sim_studies/two_sided/"

# covariance misspec. results

cov_misspec_results_out <- fread(paste0(sim_study_dir, 
                                        "cov_misspecification_sim_study.csv"))


cov_misspec_results_out$misspecification <- rep(c("Mat\\'ern-1/2", "Mat\\'ern-5/2"),
                                                each = 10000)

cov_misspec_results_out$cov <- NULL

cov_study_mean <- cov_misspec_results_out[, lapply(.SD, mean), by = .(sigma,
                                                                  n,
                                                                  phi,
                                                                  theta,
                                                                  psi,
                                                                  taper,
                                                                  misspecification)]


cov_study_mean[, iteration := NULL]

cov_study_mean <- melt(cov_study_mean, 
                  id.vars = c("sigma", "n", "phi", "theta", "psi", "taper", "misspecification"),
                  measure.vars = measure(method,
                                          value.name,
                                          pattern = "^([[:alpha:]]*)_(.*)$"))

cov_study_mean[, `:=`(n = NULL,
                      theta = NULL,
                     phi = NULL,
                     psi = NULL,
                     sigma = NULL,
                     taper = NULL,
                     FNP = NULL,
                     reject = NULL,
                     spec = NULL)]

cov_study_mean <- cov_study_mean[order(misspecification)]

cov_study_sd <- cov_misspec_results_out[, lapply(.SD, sd), 
by = .(sigma,
      n,
      phi,
      theta,
      psi,
      taper,
      misspecification)]

cov_study_sd[, iteration := NULL]

cov_study_sd <- melt(cov_study_sd, 
                  id.vars = c("sigma", "n", "phi", "theta", "psi", "taper", "misspecification"),
                  measure.vars = measure(method,
                                          value.name,
                                          pattern = "^([[:alpha:]]*)_(.*)$"))


cov_study_sd[, `:=`(n = NULL,
                      theta = NULL,
                     phi = NULL,
                     psi = NULL,
                     sigma = NULL,
                     taper = NULL,
                     FNP = NULL,
                     reject = NULL,
                     spec = NULL)]

cov_study_sd$FDX <- cov_study_sd$FDX / sqrt(10000)

names(cov_study_sd)[3:5] <- paste(names(cov_study_sd)[3:5], "sd", sep="_")
names(cov_study_sd)[6] <- "FDX_se"

cov_study_sd <- cov_study_sd[order(misspecification)]

joined_table <- cov_study_mean[cov_study_sd,
                        on = .(misspecification, method)]


# Laplace results

laplace_results_out <- fread(paste0(sim_study_dir,
                                    "sim_n50_scale1_phi4_taper7_laplace.csv"))


laplace_results_out$misspecification <- "Laplace"

laplace_results_out$base_FDX <- laplace_results_out$base_FDP > 0.1
laplace_results_out$smooth_FDX <- laplace_results_out$smooth_FDP > 0.1


laplace_study_mean <- laplace_results_out[, lapply(.SD, mean), by = .(sigma,
                                                                  n,
                                                                  phi,
                                                                  theta,
                                                                  psi,
                                                                  taper,
                                                                  misspecification)]


laplace_study_mean[, iteration := NULL]

laplace_study_mean <- melt(laplace_study_mean, 
                  id.vars = c("sigma", "n", "phi", "theta", "psi", "taper", "misspecification"),
                  measure.vars = measure(method,
                                          value.name,
                                          pattern = "^([[:alpha:]]*)_(.*)$"))

laplace_study_mean[, `:=`(n = NULL,
                      theta = NULL,
                     phi = NULL,
                     psi = NULL,
                     sigma = NULL,
                     taper = NULL,
                     FNP = NULL,
                     reject = NULL,
                     spec = NULL)]

laplace_study_mean <- laplace_study_mean[order(misspecification)]

laplace_study_sd <- laplace_results_out[, lapply(.SD, sd), 
by = .(sigma,
      n,
      phi,
      theta,
      psi,
      taper,
      misspecification)]

laplace_study_sd[, iteration := NULL]

laplace_study_sd <- melt(laplace_study_sd, 
                  id.vars = c("sigma", "n", "phi", "theta", "psi", "taper", "misspecification"),
                  measure.vars = measure(method,
                                          value.name,
                                          pattern = "^([[:alpha:]]*)_(.*)$"))


laplace_study_sd[, `:=`(n = NULL,
                      theta = NULL,
                     phi = NULL,
                     psi = NULL,
                     sigma = NULL,
                     taper = NULL,
                     FNP = NULL,
                     reject = NULL,
                     spec = NULL)]

laplace_study_sd$FDX <- laplace_study_sd$FDX / sqrt(10000)

names(laplace_study_sd)[3:5] <- paste(names(laplace_study_sd)[3:5], "sd", sep="_")
names(laplace_study_sd)[6] <- "FDX_se"

laplace_study_sd <- laplace_study_sd[order(misspecification)]

laplace_joined_table <- laplace_study_mean[laplace_study_sd,
                        on = .(misspecification, method)]


misspec_df <- rbind(joined_table, laplace_joined_table)



misspec_df <- misspec_df[, .(misspecification,
                             method,
                             FDP,
                             FDP_sd,
                             FDX,
                             FDX_se,
                             true_reject,
                             true_reject_sd,
                             sens, 
                             sens_sd)]

perc_delta_sens_vals <- ((misspec_df[method == "smooth"]$sens - 
                            misspec_df[method == "base"]$sens) / 
                           misspec_df[method == "base"]$sens) * 100

row_num <- nrow(misspec_df)
perc_delta_sens <- rep(NA, row_num)
perc_delta_sens[seq_len(row_num) %% 2 == 0] <- perc_delta_sens_vals

misspec_df$perc_delta_sens <- perc_delta_sens

misspec_df

fwrite(misspec_df, 
       file = paste0(sim_study_dir, "misspec_results_table2.csv"))

options(knitr.kable.NA = "-")

kable(misspec_df,
      format = "latex",
      digits = c(NA, NA, 3, 3, 3, 3, 1, 1, 3, 3, 1),
      col.names = c("misspec.",
                    "method",
                    "FDP",
                    "s.d.",
                    "FDX",
                    "s.e.",
                    "true rej.",
                    "s.d.",
                    "sens.",
                    "s.d.",
                    "$\\%\\Delta$sens."
      ),
      linesep = c("", "\\addlinespace"),
      booktabs = TRUE,
      escape = FALSE
)
