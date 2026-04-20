# Results for the FDR_L method of Zhang et al. 2011. 
# See paper for full citation.

library(data.table)
library(knitr)

options(knitr.kable.NA = "-")

sim_dir <- "/home/data/projects/clim_smooth/sim_studies/two_sided/"

fdrl_big1 <- fread(paste0(sim_dir, "fdrl_sim_big_sd1.csv"))
fdrl_big2 <- fread(paste0(sim_dir, "fdrl_sim_big_sd2.csv"))
fdrl_big3 <- fread(paste0(sim_dir, "fdrl_sim_big_sd3.csv"))

fdrl_small1 <- fread(paste0(sim_dir, "fdrl_sim_small_sd1.csv"))
fdrl_small2 <- fread(paste0(sim_dir, "fdrl_sim_small_sd2.csv"))
fdrl_small3 <- fread(paste0(sim_dir, "fdrl_sim_small_sd3.csv"))

fdrl_small3$window <- fdrl_small2$window
fdrl_small1$sigma <- 1
fdrl_small2$sigma <- 2
fdrl_small3$sigma <- 3

fdrl_big2$sigma <- 2

fdrl_big3$sigma <- 3
fdrl_big3$n <- fdrl_big2$n
fdrl_big3$phi <- fdrl_big2$phi
fdrl_big3$theta <- fdrl_big2$theta
fdrl_big3$psi <- fdrl_big2$psi
fdrl_big3$taper <- fdrl_big2$taper

fdrl_small <- rbind(fdrl_small1, fdrl_small2, fdrl_small3)
fdrl_big <- rbind(fdrl_big1, fdrl_big2, fdrl_big3)

fdrl_small$window <- 4
fdrl_big$window <- 20

fdrl <- rbind(fdrl_small, fdrl_big)

fdrl_mean <- fdrl[, lapply(.SD, mean), by = .(sigma,
                                                   n,
                                                   phi,
                                                   theta,
                                                   psi,
                                                   taper,
                                                  window)]

fdrl_mean_table <- fdrl_mean[, .(window, sigma, fdrl_FDP, fdrl_true_reject, fdrl_sens)]

fdrl_mean_table


kable(fdrl_mean_table,
      format = "latex",
      digits = c(1, 1, 3, 1, 3),
      col.names = c("window",
                    "$\\sigma_\\varepsilon$",
                    "FDP",
                    "true rej.",
                    "sens."
      ),
      linesep = c("", "", "\\addlinespace",
                  "", "", ""),
      booktabs = TRUE,
      escape = FALSE
)