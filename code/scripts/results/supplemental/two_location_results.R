## Results for Supplemental Table 

library(data.table)
library(knitr)

options(knitr.kable.NA = "-")

sim_dir <- "/home/data/projects/clim_smooth/sim_studies/two_sided/"

sim_2loc <- fread(paste0(sim_dir, 'sim_2loc.csv'))
sim_2loc_diff_betas <- fread(paste0(sim_dir, 'sim_2loc_beta1_5_0_5.csv'))

sim_2loc$base_FDX <- sim_2loc$base_FDP > 0.1
sim_2loc$smooth_FDX <- sim_2loc$smooth_FDP > 0.1

sim_2loc_proc <- process_study(sim_2loc)
sim_2loc_proc$beta1 <- 1
sim_2loc_proc$beta2 <- 1


sim_2loc_diff_betas$base_FDX <- sim_2loc_diff_betas$base_FDP > 0.1
sim_2loc_diff_betas$smooth_FDX <- sim_2loc_diff_betas$smooth_FDP > 0.1

sim_2loc_diff_betas_proc <- process_study(sim_2loc_diff_betas)
sim_2loc_diff_betas_proc$beta1 <- 1.5
sim_2loc_diff_betas_proc$beta2 <- 0.5

sim_2loc_proc
sim_2loc_diff_betas_proc

sim_2loc_mean <- rbind(sim_2loc_proc, sim_2loc_diff_betas_proc)

perc_delta_sens_vals <- ((sim_2loc_mean[method == "smooth"]$sens - 
                          sim_2loc_mean[method == "base"]$sens) / 
                          sim_2loc_mean[method == "base"]$sens) * 100

row_num <- nrow(sim_2loc_mean)
perc_delta_sens <- rep(NA, row_num)
perc_delta_sens[seq_len(row_num) %% 2 == 0] <- perc_delta_sens_vals

sim_2loc_mean$perc_delta_sens <- perc_delta_sens

sim_2loc_mean

sim_2loc_table <- sim_2loc_mean[, .(beta1,
       beta2,
        method,
         FDP, 
         FDP_sd,
          FDX,
           FDX_se,
            true_reject,
             true_reject_sd,
              sens,
               sens_sd,
                perc_delta_sens)]

kable(sim_2loc_table,
      format = "latex",
      digits = c(1, 1, NA, 3, 3, 3, 3, 1, 1, 3, 3, 1),
      col.names = c("$\\beta_1$", 
                    "$\\beta_2$",
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
      linesep = c("", "", "\\addlinespace",
                  "", "", ""),
      booktabs = TRUE,
      escape = FALSE
)

sim_study_dir <- "/home/data/projects/clim_smooth/sim_studies/two_sided/"
noise_path <- paste0(sim_study_dir, "noise_sim_study.csv")

noise_sim_study <- fread(noise_path)
noise_sd1 <- noise_sim_study[sigma == 1]

noise_sd1_table <- process_study(noise_sd1)

(sim_2loc_proc$sens - noise_sd1_table$sens) / noise_sd1_table$sens

sim_2loc_proc$sens
sim_2loc_proc
noise_sd1_table
