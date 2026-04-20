#Code for supplemental Figure 3.

library(data.table)
library(ggplot2)
library(knitr)
library(patchwork)
library(stringr)

data_dir <- '/home/data/projects/clim_smooth/'
doubledip_dir <- paste0(data_dir, 'sim_studies/double_dipping/p_values/')

doubledip_big_dir <- paste0(data_dir, 'sim_studies/double_dipping/')

plot_dir <- "/home/data/projects/clim_smooth/plots/"

set.seed(609010)

fields <- 1:10
# randomly select one iteration per field
rand_iters <- sample(10, size = 10, replace = TRUE)

iter_list <- vector(10, mode = "list")
for (i in 1:10) {
  filename <- paste0("field",
  fields[i],
  "_iter",
  rand_iters[i],
  "_pvals.csv")
  iter_list[[i]] <- fread(paste0(doubledip_dir, filename))
}

difference_list <- vector(10, mode = "list")
for (i in 1:10) {
   p_val_differences <- iter_list[[i]][method == 'smooth']$p_val - 
    iter_list[[i]][method == 'iso_smooth']$p_val
  difference_list[[i]] <- data.frame(p_val_difference = p_val_differences)
}

lapply(difference_list, \(x) summary(x$p_val_difference))

### Histograms

plot_list <- vector(10, mode = "list")
for (i in 1:10) {
  plot_list[[i]] <- ggplot(data = difference_list[[i]], aes(x = p_val_difference)) + 
  geom_histogram(aes(y = after_stat(density)),
 bins = 30,
  fill = "darkorchid",
   color = "black") +
  labs(x = "smoothing method p-value difference (one-sample - two-sample)",
  title = 'b)') + 
   theme_bw()
}

p_val_diff_path <- paste0(plot_dir, 'double_dipping_hist.png')

### 9 is typical of the general pattern

ggsave(p_val_diff_path,
  plot = plot_list[[9]],
  units = "in",
  width = 8,
  height = 5.5)

### Scatterplots

two_pvals_list <- vector(10, mode = "list")
for (i in 1:10) {
   p_val_1sample <- iter_list[[i]][method == 'smooth']$p_val
   p_val_2sample <- iter_list[[i]][method == 'iso_smooth']$p_val
  two_pvals_list[[i]] <- data.frame(p_val_1sample = p_val_1sample,
                                     p_val_2sample = p_val_2sample)
}

scatterplot_list <- vector(10, mode = "list")
for (i in 1:10) {
  scatterplot_list[[i]] <- ggplot(data = two_pvals_list[[i]], 
      aes(x = p_val_1sample, y = p_val_2sample)) + 
  geom_point(size=0.5, alpha=0.6) +
  labs(title = 'a)',
       x = "smoothing method p-value (one-sample)",
       y = 'smoothing method p-value (two-sample)') + 
   theme_bw()
}

stacked_scatter_hist <- scatterplot_list[[9]] / plot_list[[9]]

p_val_comp_path <- paste0(plot_dir, 'double_dipping_scatter_hist_combined.png')

ggsave(p_val_comp_path,
  plot = stacked_scatter_hist,
  units = "in",
  width = 8,
  height = 8)

doubledip_output <- fread(paste0(doubledip_big_dir, 'sim_n50_sd1_phi4_taper7_doubledip_big.csv'))

doubledip_mean  <- doubledip_output[, lapply(.SD, mean), by = .(sigma,
                                                   n,
                                                   phi,
                                                   theta,
                                                   psi,
                                                   taper)]

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

  study_sd$FDX <- study_sd$FDX / sqrt(100)

  names(study_sd)[8:10] <- paste(names(study_sd)[8:10], "sd", sep="_")
  names(study_sd)[11] <- "FDX_se"

  study_sd <- study_sd[order(sigma)]

  joined_table <- study_mean[study_sd,
                             on = .(sigma, n, phi, theta, psi, taper, method)]
  
  return(joined_table)
}

names(doubledip_output) <- str_replace(names(doubledip_output), "smooth_iso", "doublesmooth")

doubledip_output$base_FDX <- doubledip_output$base_FDP > 0.1
doubledip_output$smooth_FDX <- doubledip_output$smooth_FDP > 0.1
doubledip_output$doublesmooth_FDX <- doubledip_output$doublesmooth_FDP > 0.1

doubledip_table <- process_study(doubledip_output)

doubledip_table$sqrt_psi <- sqrt(doubledip_table$psi)

doubledip_table <- doubledip_table[, .(method,
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


doubledip_table <- doubledip_table[order(sqrt_psi, theta, -phi, taper)]

doubledip_table_1half <- doubledip_table[, .(method,
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

doubledip_table_2half <- doubledip_table[, .(method,
       sigma,
       n,
       phi,
       theta,
       sqrt_psi,
       taper,
       true_reject,
       true_reject_sd,
       sens,
       sens_sd)]

doubledip_table_1half

kable(doubledip_table_1half,
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

kable(doubledip_table_2half,
      format = "latex",
      digits = c(NA, 1, 2, 2, 0, 2, 1, 1, 1, 3, 3),
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
                    "s.d."
      ),
      linesep = c("", "", "", "", "", "\\addlinespace",
                  "", "", "", "", "", "\\addlinespace",
                  "", "", "", "", "", "", "", "\\addlinespace",
                  "", "", ""),
      booktabs = TRUE,
      escape = FALSE
)

sim_dir <- '/home/data/projects/clim_smooth/sim_studies/two_sided/'
cor01_path <- paste0(sim_dir, 'sim_n50_sd1_phi4_taper7_corthresh01.csv')
cor02_path <- paste0(sim_dir, 'sim_n50_sd2_phi4_taper7_corthresh01.csv')
cor03_path <- paste0(sim_dir, 'sim_n50_sd3_phi4_taper7_corthresh01.csv')

cor01_sd1 <- fread(cor01_path)
cor01_sd2 <- fread(cor02_path)
cor01_sd3 <- fread(cor03_path)

cor01_study <- rbind(cor01_sd1, cor01_sd2, cor01_sd3)

cor01_study$base_FDX <- cor01_study$base_FDP > 0.1
cor01_study$smooth_FDX <- cor01_study$smooth_FDP > 0.1

cor01_table <- process_study(cor01_study)

cor01_table$sqrt_psi <- sqrt(cor01_table$psi)

cor01_sd1$base_FDX <- cor01_sd1$base_FDP > 0.1
cor01_sd1$smooth_FDX <- cor01_sd1$smooth_FDP > 0.1

cor01_table2 <- process_study(cor01_sd1)

cor01_table2$sqrt_psi <- sqrt(cor01_table$psi)

cor01_sd1$base_FDX <- cor01_sd1$base_FDP > 0.1
cor01_sd1$smooth_FDX <- cor01_sd1$smooth_FDP > 0.1

cor01_table <- process_study(cor01_study)

cor01_table$sqrt_psi <- sqrt(cor01_table$psi)
cor01_sd1$base_FDX <- cor01_sd1$base_FDP > 0.1
cor01_sd1$smooth_FDX <- cor01_sd1$smooth_FDP > 0.1

cor01_table <- process_study(cor01_sd1)

cor01_table$sqrt_psi <- sqrt(cor01_table$psi)

cor01_table <- cor01_table[, .(method,
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


cor01_table <- cor01_table[order(sqrt_psi, theta, -phi, taper)]

perc_delta_sens_vals <- ((cor01_table[method == "smooth"]$sens - 
                          cor01_table[method == "base"]$sens) / 
                          cor01_table[method == "base"]$sens) * 100

row_num <- nrow(cor01_table)
perc_delta_sens <- rep(NA, row_num)
perc_delta_sens[seq_len(row_num) %% 2 == 0] <- perc_delta_sens_vals

cor01_table$perc_delta_sens <- perc_delta_sens

setcolorder(cor01_table,
            neworder = "perc_delta_sens",
            after = "sens_sd")

options(knitr.kable.NA = "-")

cor01_table_1half <- cor01_table[, .(method,
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

cor01_table_2half <- cor01_table[, .(method,
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

kable(cor01_table_1half,
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

kable(cor01_table_2half,
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

noise_table <- process_study(noise_sim_study)

cor01_table$sens / noise_table$sens


base_sens_2loc1 <- 0.069
smooth_sens_2loc1 <- 0.089

base_sens <- 0.074
smooth_sens <- 0.093

base_sens_2loc2 <- 0.063
smooth_sens_2loc2 <- 0.080

1 - base_sens_2loc1 / base_sens

# 6.8% decline

1 - smooth_sens_2loc1 / smooth_sens

# 4.3% decline
# 
#  

1 - base_sens_2loc2 / base_sens

# 14.9% decline

1 - smooth_sens_2loc2 / smooth_sens

# 13.9% decline
# 
#  

base_fdp_2loc1 <- 0.050
smooth_fdp_2loc1 <- 0.049

base_fdp <- 0.065
smooth_fdp <- 0.061

base_fdp_2loc2 <- 0.053
smooth_fdp_2loc2 <- 0.052

1 - base_fdp_2loc1 / base_fdp

# 23% decline

1 - smooth_fdp_2loc1 / smooth_fdp

# 19.7% decline
# 
#  

1 - base_fdp_2loc2 / base_fdp

# 18.5% decline

1 - smooth_fdp_2loc2 / smooth_fdp

# 14.8% decline
# 




save_dir <- "/home/data/projects/clim_smooth/sim_studies/two_sided/"
sim_sd1_path <- paste0(save_dir, "sim_yb_resample_n50_sd1_phi4_taper7.csv")

sim_sd1 <- fread(sim_sd1_path)

sim_yb_table <- process_study(sim_sd1)

sim_yb_table <- sim_yb_table[c(1, 3, 2), ]

sim_yb_table

sim_yb_table_limited <- sim_yb_table[, .(method, true_reject, FDP, sens)]

sim_yb_table_limited

kable(sim_yb_table_limited,
      format = "latex",
      digits = c(NA, 1, 3, 3),
      col.names = c("method",
                     "true rej.",
                    "FDP",
                    "sens."
      ),
      linesep = c("", "", "", "", ""),
      booktabs = TRUE,
      escape = FALSE
)

