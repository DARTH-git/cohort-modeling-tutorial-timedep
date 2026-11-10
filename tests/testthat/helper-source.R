# Load the packages and source the repository's functions so that the tests
# exercise exactly the code used by the analysis scripts.
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

source(test_path("..", "..", "R", "Functions.R"))

# The two models define functions with the same names (decision_model(),
# calculate_ce_out(), generate_psa_params(), ...), so each is sourced into its
# own environment. Their parent is the global environment, where the shared
# functions of "R/Functions.R" (e.g., gen_wcc()) live.
env_sim <- new.env(parent = globalenv())
sys.source(test_path("..", "..", "R", "Functions_cSTM_time_dep_simulation.R"),
           envir = env_sim)
env_sr <- new.env(parent = globalenv())
sys.source(test_path("..", "..", "R", "Functions_cSTM_time_dep_state_residence.R"),
           envir = env_sr)

# Age-specific background mortality (same as in the analysis scripts)
n_age_init <- 25
n_age_max  <- 100
v_r_HDage  <- read.csv(test_path("..", "..", "data", "LifeTable_USA_Mx_2015.csv")) %>%
  dplyr::filter(Age >= n_age_init & Age < n_age_max) %>%
  dplyr::pull(Total)

# Parameters shared by both models (same as in the analysis scripts)
l_params_common <- list(
  v_r_HDage = v_r_HDage,
  r_HS1 = 0.15, r_S1H = 0.5, hr_S1 = 3, hr_S2 = 10, hr_S1S2_trtB = 0.6,
  c_H = 2000, c_S1 = 4000, c_S2 = 15000, c_D = 0,
  c_trtA = 12000, c_trtB = 13000,
  u_H = 1, u_S1 = 0.75, u_S2 = 0.5, u_D = 0, u_trtA = 0.95,
  du_HS1 = 0.01, ic_HS1 = 1000, ic_D = 2000,
  n_age_init = n_age_init, n_age_max = n_age_max,
  d_c = 0.03, d_e = 0.03,
  cycle_length = 1
)
# Simulation-time-dependent model: constant S1 -> S2 rate
l_params_sim <- c(l_params_common, list(r_S1S2 = 0.105))
# State-residence-dependent model: Weibull S1 -> S2 hazard
l_params_sr  <- c(l_params_common, list(r_S1S2_scale = 0.08, r_S1S2_shape = 1.1))
