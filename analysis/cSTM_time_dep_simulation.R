# Appendix code to Time-dependent cSTMs in R: Simulation-time dependency ----

#* This code forms the basis for the state-transition model of the tutorial: 
#* 'A Tutorial on Time-Dependent Cohort State-Transition Models in R using a 
#* Cost-Effectiveness Analysis Example' 
#* Authors: 
#* - Fernando Alarid-Escudero <fernando.alarid@cide.edu>
#* - Eline Krijkamp
#* - Eva A. Enns
#* - Alan Yang
#* - M.G. Myriam Hunink
#* - Petros Pechlivanoglou
#* - Hawre Jalal
#* Please cite the article when using this code
#* 
#* To program this tutorial we used:
#* R version 4.0.5 (2021-03-31)
#* Platform: 64-bit operating system, x64-based processor
#* Running under: Mac OS 12.2.1
#* RStudio: Version 1.4.1717 2009-2021 RStudio, Inc

#******************************************************************************#
# Description ----

#* This code implements a simulation-time-dependent Sick-Sicker cSTM model to 
#* conduct a CEA of four strategies:
#* - Standard of Care (SoC): best available care for the patients with the 
#*   disease. This scenario reflects the natural history of the disease 
#*   progression.
#* - Strategy A: treatment A is given to patients in the Sick and Sicker states, 
#*   but only improves the quality of life of those in the Sick state.
#* - Strategy B: treatment B is given to all sick patients and reduces disease 
#*   progression from the Sick to Sicker state.
#* - Strategy AB: This strategy combines treatment A and treatment B. The disease 
#*   progression is reduced, and individuals in the Sick state have an improved 
#*   quality of life.

#******************************************************************************#
# Initial setup ---- 
rm(list = ls())    # remove any variables in R's memory 

## Install required packages ----
# install.packages("dplyr")      # to manipulate data
# install.packages("tidyr")      # to manipulate data
# install.packages("reshape2")   # to manipulate data
# install.packages("ggplot2")    # to visualize data
# install.packages("ggrepel")    # to visualize data
# install.packages("gridExtra")  # to visualize data
# install.packages("ellipse")    # to visualize data
# install.packages("scales")     # for dollar signs and commas
# install.packages(patchwork)    # for combining ggplot2 figures
# install.packages("dampack")    # for CEA and calculate ICERs
# install.packages("devtools")   # to install packages from GitHub
# devtools::install_github("DARTH-git/darthtools") # to install darthtools from GitHub using devtools
# install.packages("doParallel") # to handle parallel processing

## Load packages ----
library(dplyr)
library(tidyr)
library(reshape2)   # For melting data
library(ggplot2)    # For plotting
library(ggrepel)    # For plotting
library(gridExtra)  # For plotting
library(ellipse)    # For plotting
library(scales)     # For dollar signs and commas
library(patchwork)  # For combining ggplot2 figures
# library(dampack)  # Uncomment to use CEA and PSA visualization functionality from dampack instead of the functions included in this repository
# library(darthtools) # Uncomment to use WCC, parameter transformation, and matrix checks from darthtools instead of the functions included in this repository
# library(doParallel) # For running PSA in parallel

## Load supplementary functions ----
source("R/Functions.R")

# Model input ----
## General setup ----
cycle_length <- 1 # cycle length equal to one year (use 1/12 for monthly)
n_age_init <- 25  # age at baseline
n_age_max  <- 100 # maximum age of follow up
n_cycles <- (n_age_max - n_age_init)/cycle_length # time horizon, number of cycles
#* Age labels 
v_age_names <- paste(rep(n_age_init:(n_age_max-1), each = 1/cycle_length), 
                     1:(1/cycle_length), 
                     sep = ".")
#* the 4 health states of the model:
v_names_states <- c("H",  # Healthy (H)
                    "S1", # Sick (S1)
                    "S2", # Sicker (S2)
                    "D")  # Dead (D)
                                           
n_states <- length(v_names_states)   # number of health states 

### Discounting factors ----
d_c <- 0.03 # annual discount rate for costs 
d_e <- 0.03 # annual discount rate for QALYs

### Strategies ----
v_names_str <- c("Standard of care",      # store the strategy names
                 "Strategy A", 
                 "Strategy B",
                 "Strategy AB") 
n_str       <- length(v_names_str)        # number of strategies

## Within-cycle correction (WCC) using Simpson's 1/3 rule ----
v_wcc <- gen_wcc(n_cycles = n_cycles,  # Function included in "R/Functions.R". The latest version can be found in `darthtools` package
                 method = "Simpson1/3") # vector of wcc

### Transition rates (annual), and hazard ratios (HRs) ----
r_HS1  <- 0.15  # constant annual rate of becoming Sick when Healthy
r_S1H  <- 0.5   # constant annual rate of becoming Healthy when Sick
r_S1S2 <- 0.105 # constant annual rate of becoming Sicker when Sick
hr_S1  <- 3     # hazard ratio of death in Sick vs Healthy 
hr_S2  <- 10    # hazard ratio of death in Sicker vs Healthy 

### Effectiveness of treatment B ----
hr_S1S2_trtB <- 0.6  # hazard ratio of becoming Sicker when Sick under treatment B

## Age-dependent mortality rates ----
lt_usa_2015 <- read.csv("data/LifeTable_USA_Mx_2015.csv")
#* Extract age-specific all-cause mortality for ages in model time horizon
v_r_mort_by_age <- lt_usa_2015 %>% 
  dplyr::filter(Age >= n_age_init & Age < n_age_max) %>%
  dplyr::select(Total) %>%
  as.matrix()

### State rewards ----
#### Costs ----
c_H    <- 2000  # annual cost of being Healthy
c_S1   <- 4000  # annual cost of being Sick
c_S2   <- 15000 # annual cost of being Sicker
c_D    <- 0     # annual cost of being dead
c_trtA <- 12000 # annual cost of receiving treatment A
c_trtB <- 13000 # annual cost of receiving treatment B
#### Utilities ----
u_H    <- 1     # annual utility of being Healthy
u_S1   <- 0.75  # annual utility of being Sick
u_S2   <- 0.5   # annual utility of being Sicker
u_D    <- 0     # annual utility of being dead
u_trtA <- 0.95  # annual utility when receiving treatment A

### Transition rewards ----
du_HS1 <- 0.01  # disutility when transitioning from Healthy to Sick
ic_HS1 <- 1000  # increase in cost when transitioning from Healthy to Sick
ic_D   <- 2000  # increase in cost when dying

### Discount weight for costs and effects ----
v_dwc  <- 1 / ((1 + (d_e * cycle_length)) ^ (0:n_cycles))
v_dwe  <- 1 / ((1 + (d_c * cycle_length)) ^ (0:n_cycles))

# Process model inputs ----
## Age-specific transition rates to the Dead state for all cycles ----
v_r_HDage  <- rep(v_r_mort_by_age, each = 1/cycle_length)
#* Name age-specific mortality vector 
names(v_r_HDage) <- v_age_names

#* compute mortality rates
v_r_S1Dage <- v_r_HDage * hr_S1 # Age-specific mortality rate in the Sick state 
v_r_S2Dage <- v_r_HDage * hr_S2 # Age-specific mortality rate in the Sicker state 
#* transform rates to probabilities adjusting by cycle length
#* Function included in "R/Functions.R". The latest version can be found in `darthtools` package
p_HS1  <- rate_to_prob(r = r_HS1, t = cycle_length) # constant annual probability of becoming Sick when Healthy conditional on surviving 
p_S1H  <- rate_to_prob(r = r_S1H, t = cycle_length) # constant annual probability of becoming Healthy when Sick conditional on surviving
p_S1S2 <- rate_to_prob(r = r_S1S2, t = cycle_length)# constant annual probability of becoming Sicker when Sick conditional on surviving
v_p_HDage  <- rate_to_prob(v_r_HDage, t = cycle_length)  # Age-specific mortality risk in the Healthy state 
v_p_S1Dage <- rate_to_prob(v_r_S1Dage, t = cycle_length) # Age-specific mortality risk in the Sick state
v_p_S2Dage <- rate_to_prob(v_r_S2Dage, t = cycle_length) # Age-specific mortality risk in the Sicker state

## Annual transition probability of becoming Sicker when Sick for treatment B ----
#* Apply hazard ratio to rate to obtain transition rate of becoming Sicker when 
#* Sick for treatment B
r_S1S2_trtB <- r_S1S2 * hr_S1S2_trtB
#* Transform rate to probability to become Sicker when Sick under treatment B 
#* adjusting by cycle length conditional on surviving
#* (Function included in "R/Functions.R". The latest version can be found in 
#* `darthtools` package)
p_S1S2_trtB <- rate_to_prob(r = r_S1S2_trtB, t = cycle_length)

# Construct state-transition models ----
## Initial state vector ----
#* All starting healthy
v_m_init <- c(H = 1, S1 = 0, S2 = 0, D = 0) # initial state vector

## Initialize cohort traces ----
### Initialize cohort trace under SoC ----
m_M_SoC <- matrix(NA, 
              nrow = (n_cycles + 1), ncol = n_states, 
              dimnames = list(0:n_cycles, v_names_states))
#* Store the initial state vector in the first row of the cohort trace
m_M_SoC[1, ] <- v_m_init

### Initialize cohort trace for strategies A, B, and AB ----
#* Structure and initial states are the same as for SoC
m_M_strA  <- m_M_SoC # Strategy A
m_M_strB  <- m_M_SoC # Strategy B
m_M_strAB <- m_M_SoC # Strategy AB

## Create transition probability arrays for strategy SoC ----
### Initialize transition probability array for strategy SoC ----
#* All transitions to a non-death state are assumed to be conditional on survival
a_P_SoC <- array(0,
                 dim  = c(n_states, n_states, n_cycles),
                 dimnames = list(v_names_states, 
                                 v_names_states, 
                                 0:(n_cycles - 1)))
### Fill in array
## From H
a_P_SoC["H", "H", ]   <- (1 - v_p_HDage) * (1 - p_HS1)
a_P_SoC["H", "S1", ]  <- (1 - v_p_HDage) * p_HS1
a_P_SoC["H", "D", ]   <- v_p_HDage
## From S1
a_P_SoC["S1", "H", ]  <- (1 - v_p_S1Dage) * p_S1H
a_P_SoC["S1", "S1", ] <- (1 - v_p_S1Dage) * (1 - (p_S1H + p_S1S2))
a_P_SoC["S1", "S2", ] <- (1 - v_p_S1Dage) * p_S1S2
a_P_SoC["S1", "D", ]  <- v_p_S1Dage
## From S2
a_P_SoC["S2", "S2", ] <- 1 - v_p_S2Dage
a_P_SoC["S2", "D", ]  <- v_p_S2Dage
## From D
a_P_SoC["D", "D", ]   <- 1

### Initialize transition probability array for strategy A as a copy of SoC's ----
a_P_strA <- a_P_SoC

### Initialize transition probability array for strategy B ----
a_P_strB <- a_P_SoC
#* Update only transition probabilities from S1 involving p_S1S2
a_P_strB["S1", "S1", ] <- (1 - v_p_S1Dage) * (1 - (p_S1H + p_S1S2_trtB))
a_P_strB["S1", "S2", ] <- (1 - v_p_S1Dage) * p_S1S2_trtB

### Initialize transition probability array for strategy AB as a copy of B's ----
a_P_strAB <- a_P_strB

## Check if transition probability arrays are valid ----
#* Functions included in "R/Functions.R". The latest version can be found in `darthtools` package
### Check that transition probabilities are [0, 1] ----
check_transition_probability(a_P_SoC, verbose = TRUE)
check_transition_probability(a_P_strA, verbose = TRUE)
check_transition_probability(a_P_strB, verbose = TRUE)
check_transition_probability(a_P_strAB, verbose = TRUE)
### Check that all rows for each slice of the array sum to 1 ----
check_sum_of_transition_array(a_P_SoC, n_states = n_states, n_cycles = n_cycles, verbose = TRUE)
check_sum_of_transition_array(a_P_strA, n_states = n_states, n_cycles = n_cycles, verbose = TRUE)
check_sum_of_transition_array(a_P_strB, n_states = n_states, n_cycles = n_cycles, verbose = TRUE)
check_sum_of_transition_array(a_P_strAB, n_states = n_states, n_cycles = n_cycles, verbose = TRUE)

## Create transition dynamics arrays ----
#* These arrays will capture transitions from each state to another over time 
### Initialize transition dynamics array for strategy SoC ----
a_A_SoC <- array(0,
             dim      = c(n_states, n_states, n_cycles + 1),
             dimnames = list(v_names_states, v_names_states, 0:n_cycles))
#* Set first slice of a_A_SoC with the initial state vector in its diagonal
diag(a_A_SoC[, , 1]) <- v_m_init
### Initialize transition-dynamics array for strategies A, B, and AB ----
#* Structure and initial states are the same as for SoC
a_A_strA  <- a_A_SoC
a_A_strB  <- a_A_SoC
a_A_strAB <- a_A_SoC

#  Run Markov model ----
#* Iterative solution of age-dependent cSTM
for(t in 1:n_cycles){
  ## Fill in cohort trace
  # For SoC
  m_M_SoC[t + 1, ]  <- m_M_SoC[t, ]  %*% a_P_SoC[, , t]
  # For strategy A
  m_M_strA[t + 1, ] <- m_M_strA[t, ] %*% a_P_strA[, , t]
  # For strategy B 
  m_M_strB[t + 1, ] <- m_M_strB[t, ] %*% a_P_strB[, , t]
  # For strategy ZB 
  m_M_strAB[t + 1, ] <- m_M_strAB[t, ] %*% a_P_strAB[, , t]
  
  ## Fill in transition-dynamics array
  # For SoC
  a_A_SoC[, , t + 1]  <- diag(m_M_SoC[t, ]) %*% a_P_SoC[, , t]
  # For strategy A
  a_A_strA[, , t + 1] <- diag(m_M_strA[t, ]) %*% a_P_strA[, , t]
  # For strategy B
  a_A_strB[, , t + 1] <- diag(m_M_strB[t, ]) %*% a_P_strB[, , t]
  # For strategy AB
  a_A_strAB[, , t + 1] <- diag(m_M_strAB[t, ]) %*% a_P_strAB[, , t]
}

## Store the cohort traces in a list ----
l_m_M <- list(SoC =  m_M_SoC,
              A   =  m_M_strA,
              B   =  m_M_strB,
              AB  =  m_M_strAB)
names(l_m_M) <- v_names_str

## Store the transition dynamics array for each strategy in a list ----
l_a_A <- list(SoC =  a_A_SoC,
              A   =  a_A_strA,
              B   =  a_A_strB,
              AB  =  a_A_strAB)
names(l_a_A) <- v_names_str

# Plot Outputs ----
#* (Functions included in "R/Functions.R"; depends on the `ggplot2` package)

## Plot the cohort trace for strategy SoC ----
plot_trace(m_M_SoC)
## Plot the cohort trace for all strategies ----
plot_trace_strategy(l_m_M)

## Plot the epidemiology outcomes ----
### Survival ----
survival_plot <- plot_surv(l_m_M, v_names_death_states = "D") +
  theme(legend.position = "bottom")
survival_plot
### Prevalence ----
prevalence_S1_plot   <- plot_prevalence(l_m_M, 
                                        v_names_sick_states = c("S1"), 
                                        v_names_dead_states = "D")  +
                        theme(legend.position = "")
prevalence_S1_plot
prevalence_S2_plot   <- plot_prevalence(l_m_M, 
                                        v_names_sick_states = c("S2"), 
                                        v_names_dead_states = "D")  +
                        theme(legend.position = "")
prevalence_S2_plot
prevalence_S1S2_plot <- plot_prevalence(l_m_M, 
                                        v_names_sick_states = c("S1", "S2"), 
                                        v_names_dead_states = "D") +
                        theme(legend.position = "")
prevalence_S1S2_plot
prop_sicker_plot     <- plot_proportion_sicker(l_m_M, 
                                               v_names_sick_states = c("S1", "S2"), 
                                               v_names_sicker_states = c("S2")) +
                        theme(legend.position = "bottom")
prop_sicker_plot

## Combine plots ----
gridExtra::grid.arrange(survival_plot, 
                        prevalence_S1_plot, 
                        prevalence_S2_plot, 
                        prevalence_S1S2_plot, 
                        prop_sicker_plot, 
                        ncol = 1, heights = c(0.75, 0.75, 0.75, 0.75, 1))

# State Rewards ----
## Scale by the cycle length ----
#* Vector of state utilities under strategy SoC
v_u_SoC    <- c(H  = u_H, 
                S1 = u_S1, 
                S2 = u_S2, 
                D  = u_D) * cycle_length
#* Vector of state costs under strategy SoC
v_c_SoC    <- c(H  = c_H, 
                S1 = c_S1,
                S2 = c_S2, 
                D  = c_D) * cycle_length
#* Vector of state utilities under strategy A
v_u_strA   <- c(H  = u_H, 
                S1 = u_trtA, 
                S2 = u_S2, 
                D  = u_D) * cycle_length
#* Vector of state costs under strategy A
v_c_strA   <- c(H  = c_H, 
                S1 = c_S1 + c_trtA,
                S2 = c_S2 + c_trtA, 
                D  = c_D) * cycle_length
#* Vector of state utilities under strategy B
v_u_strB   <- c(H  = u_H, 
                S1 = u_S1, 
                S2 = u_S2, 
                D  = u_D) * cycle_length
#* Vector of state costs under strategy B
v_c_strB   <- c(H  = c_H, 
                S1 = c_S1 + c_trtB, 
                S2 = c_S2 + c_trtB, 
                D  = c_D) * cycle_length
#* Vector of state utilities under strategy AB
v_u_strAB  <- c(H  = u_H, 
                S1 = u_trtA, 
                S2 = u_S2, 
                D  = u_D) * cycle_length
#* Vector of state costs under strategy AB
v_c_strAB  <- c(H  = c_H, 
                S1 = c_S1 + (c_trtA + c_trtB), 
                S2 = c_S2 + (c_trtA + c_trtB), 
                D  = c_D) * cycle_length

## Store state rewards ----
#* Store the vectors of state utilities for each strategy in a list 
l_u <- list(SoC = v_u_SoC,
            A   = v_u_strA,
            B   = v_u_strB,
            AB  = v_u_strAB)
#* Store the vectors of state cost for each strategy in a list 
l_c <- list(SoC =  v_c_SoC,
            A   =  v_c_strA,
            B   =  v_c_strB,
            AB  =  v_c_strAB)

#* assign strategy names to matching items in the lists
names(l_u) <- names(l_c) <- v_names_str

# Compute expected outcomes ----
#* Create empty vectors to store total utilities and costs 
v_tot_qaly <- v_tot_cost <- vector(mode = "numeric", length = n_str)
names(v_tot_qaly) <- names(v_tot_cost) <- v_names_str

## Loop through each strategy and calculate total utilities and costs ----
for (i in 1:n_str) { # i <- 1
  v_u_str <- l_u[[i]]   # select the vector of state utilities for the i-th strategy
  v_c_str <- l_c[[i]]   # select the vector of state costs for the i-th strategy
  a_A_str <- l_a_A[[i]] # select the transition array for the i-th strategy
  
  ##* Array of state rewards 
  #* Create transition matrices of state utilities and state costs for the i-th strategy 
  m_u_str   <- matrix(v_u_str, nrow = n_states, ncol = n_states, byrow = T)
  m_c_str   <- matrix(v_c_str, nrow = n_states, ncol = n_states, byrow = T)
  #* Expand the transition matrix of state utilities across cycles to form a transition array of state utilities
  a_R_u_str <- array(m_u_str, 
                     dim      = c(n_states, n_states, n_cycles + 1),
                     dimnames = list(v_names_states, v_names_states, 0:n_cycles))
  # Expand the transition matrix of state costs across cycles to form a transition array of state costs
  a_R_c_str <- array(m_c_str, 
                     dim      = c(n_states, n_states, n_cycles + 1),
                     dimnames = list(v_names_states, v_names_states, 0:n_cycles))
  
  ##* Apply transition rewards
  #* Apply disutility due to transition from H to S1
  a_R_u_str["H", "S1", ]      <- a_R_u_str["H", "S1", ]       - du_HS1
  #* Add transition cost per cycle due to transition from H to S1
  a_R_c_str["H", "S1", ]      <- a_R_c_str["H", "S1", ]       + ic_HS1
  #* Add transition cost  per cycle of dying from all non-dead states
  a_R_c_str[-n_states, "D", ] <- a_R_c_str[-n_states, "D", ] + ic_D
  
  ###* Expected QALYs and costs for all transitions per cycle
  #* QALYs = life years x QoL
  #* Note: all parameters are annual in our example. In case your own case example is different make sure you correctly apply.
  a_Y_c_str <- a_A_str * a_R_c_str
  a_Y_u_str <- a_A_str * a_R_u_str 
  
  ###* Expected QALYs and costs per cycle
  ##* Vector of QALYs and costs
  v_qaly_str <- apply(a_Y_u_str, 3, sum) # sum the proportion of the cohort across transitions 
  v_cost_str <- apply(a_Y_c_str, 3, sum) # sum the proportion of the cohort across transitions
  
  ##* Discounted total expected QALYs and Costs per strategy and apply within-cycle correction if applicable
  #* QALYs
  v_tot_qaly[i] <- t(v_qaly_str) %*% (v_dwe * v_wcc)
  #* Costs
  v_tot_cost[i] <- t(v_cost_str) %*% (v_dwc * v_wcc)
}

# Cost-effectiveness analysis (CEA) ----
## Incremental cost-effectiveness ratios (ICERs) ----
#* Function included in "R/Functions.R"; depends on the `dplyr` package
#* The latest version can be found in `dampack` package
df_cea <- calculate_icers(cost       = v_tot_cost, 
                          effect     = v_tot_qaly,
                          strategies = v_names_str)
df_cea

## CEA table in proper format ----
table_cea <- format_table_cea(df_cea) # Function included in "R/Functions.R"; depends on the `scales` package
table_cea

## CEA frontier -----
#* Function included in "R/Functions.R"; depends on the `ggplot2`  and `ggrepel` packages.
#* The latest version can be found in `dampack` package
plot(df_cea, label = "all", txtsize = 16) +
  expand_limits(x = max(table_cea$QALYs) + 0.1) +
  theme(legend.position = c(0.8, 0.2))

#******************************************************************************#
# Probabilistic Sensitivity Analysis (PSA) -----
## Load model, CEA and PSA functions ----
source('R/Functions_cSTM_time_dep_simulation.R')
source('R/Functions.R')

## List of input parameters -----
l_params_all <- list(
  # Transition probabilities (per cycle), hazard ratios
  v_r_HDage   = v_r_HDage, # constant rate of dying when Healthy (all-cause mortality)
  r_HS1       = 0.15,  # constant annual rate of becoming Sick when Healthy conditional on surviving
  r_S1H       = 0.5,   # constant annual rate of becoming Healthy when Sick conditional on surviving
  r_S1S2      = 0.105, # constant annual rate of becoming Sicker when Sick conditional on surviving
  hr_S1       = 3,     # hazard ratio of death in Sick vs Healthy 
  hr_S2       = 10,    # hazard ratio of death in Sicker vs Healthy 
  # Effectiveness of treatment B 
  hr_S1S2_trtB = 0.6,  # hazard ratio of becoming Sicker when Sick under treatment B
  ## State rewards
  # Costs
  c_H    = 2000,  # cost of remaining one cycle in Healthy 
  c_S1   = 4000,  # cost of remaining one cycle in Sick 
  c_S2   = 15000, # cost of remaining one cycle in Sicker 
  c_D    = 0,     # cost of being dead (per cycle)
  c_trtA = 12000, # cost of treatment A
  c_trtB = 13000, # cost of treatment B
  # Utilities
  u_H    = 1,     # utility when Healthy 
  u_S1   = 0.75,  # utility when Sick 
  u_S2   = 0.5,   # utility when Sicker
  u_D    = 0,     # utility when Dead 
  u_trtA = 0.95,  # utility when being treated with A
  ## Transition rewards
  du_HS1 = 0.01,  # disutility when transitioning from Healthy to Sick
  ic_HS1 = 1000,  # increase in cost when transitioning from Healthy to Sick
  ic_D   = 2000,   # increase in cost when dying
  # Initial and maximum ages
  n_age_init = 25,
  n_age_max = 100,
  # Discount rates
  d_c = 0.03, # annual discount rate for costs 
  d_e = 0.03, # annual discount rate for QALYs,
  # Cycle length
  cycle_length = 1
)

#* Store the parameter names into a vector
v_names_params <- names(l_params_all)

## Test functions to generate CE outcomes and PSA dataset ----
#* Test function to compute CE outcomes
calculate_ce_out(l_params_all) # Function included in "R/Functions_cSTM_time_dep_simulation.R"

#* Test function to generate PSA input dataset
generate_psa_params(10) # Function included in "R/Functions_cSTM_time_dep_simulation.R"

## Generate PSA dataset ----
#* Number of simulations
n_sim <- 1000

#* Generate PSA input dataset
df_psa_input <- generate_psa_params(n_sim = n_sim)
#* First six observations
head(df_psa_input)

### Histogram of parameters ----
ggplot(melt(df_psa_input, variable.name = "Parameter"), aes(x = value)) +
  facet_wrap(~Parameter, scales = "free") +
  geom_histogram(aes(y = ..density..)) +
  ylab("") +
  theme_bw(base_size = 16) + 
  theme(axis.text = element_text(size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank()) 

## Run PSA ----
#* Initialize data.frames with PSA output 
#* data.frame of costs
df_c <- as.data.frame(matrix(0, 
                             nrow = n_sim,
                             ncol = n_str))
colnames(df_c) <- v_names_str
#* data.frame of effectiveness
df_e <- as.data.frame(matrix(0, 
                             nrow = n_sim,
                             ncol = n_str))
colnames(df_e) <- v_names_str

#* data.frame of survival
m_surv <- matrix(NA, nrow = n_sim, ncol = (n_cycles + 1), 
                 dimnames = list(1:n_sim, 0:n_cycles))
df_surv <- data.frame(Outcome = "Survival",
                      m_surv, check.names = "FALSE")
#* data.frame of life expectancy
df_le <- data.frame(Outcome = "Life Expectancy",
                    LE = m_surv[, 1])
#* data.frame of prevalence of S1
m_prev <- matrix(NA, nrow = n_sim, ncol = (n_cycles + 1), 
                 dimnames = list(1:n_sim, 0:n_cycles))
df_prevS1 <- data.frame(States = "S1",
                        m_prev, check.names = "FALSE")
#* data.frame of prevalence of S2
df_prevS2 <- data.frame(States = "S2",
                        m_prev, check.names = "FALSE")

#* data.frame of prevalence of S1 & S2
df_prevS1S2 <- data.frame(States = "S1 + S2",
                          m_prev, check.names = "FALSE")

#* Conduct probabilistic sensitivity analysis
#* Run Markov model on each parameter set of PSA input dataset
n_time_init_psa_series <- Sys.time()
for (i in 1:n_sim) { # i <- 1
  l_psa_input <- update_param_list(l_params_all, df_psa_input[i,])
  # Economics Measures
  l_out_ce_temp  <- calculate_ce_out(l_psa_input)
  df_c[i, ]  <- l_out_ce_temp$Cost  
  df_e[i, ]  <- l_out_ce_temp$Effect
  # Epidemiological Measures
  l_out_epi_temp <- generate_epi_measures_SoC(l_psa_input)
  df_surv[i, -1] <- l_out_epi_temp$S
  df_le[i, -1] <- l_out_epi_temp$LE
  df_prevS1[i, -1] <- l_out_epi_temp$PrevS1
  df_prevS2[i, -1] <- l_out_epi_temp$PrevS2
  df_prevS1S2[i, -1] <- l_out_epi_temp$PrevS1S2
  # Display simulation progress
  if (i/(n_sim/100) == round(i/(n_sim/100), 0)) { # display progress every 5%
    cat('\r', paste(i/n_sim * 100, "% done", sep = " "))
  }
}
n_time_end_psa_series <- Sys.time()
n_time_total_psa_series <- n_time_end_psa_series - n_time_init_psa_series
print(paste0("PSA with ", scales::comma(n_sim), " simulations run in series in ", 
             round(n_time_total_psa_series, 2), " ", 
             units(n_time_total_psa_series)))

## Run Markov model on each parameter set of PSA input dataset in parallel
# Uncomment next section to run in parallel
## Get OS
os <- get_os()

no_cores <- parallel::detectCores() - 1

print(paste0("Parallelized PSA on ", os, " using ", no_cores, " cores."))

n_time_init_psa_parallel <- Sys.time()

# ## Run parallelized PSA based on OS
# if (os == "osx") {
#   # Initialize cluster object
#   cl <- parallel::makeForkCluster(no_cores)
#   # Register clusters
#   doParallel::registerDoParallel(cl)
#   # Run parallelized PSA
#   df_ce <- foreach::foreach(i = 1:n_sim, .combine = rbind) %dopar% {
#     l_out_temp <- calculate_ce_out(df_psa_input[i, ])
#     df_ce <- c(l_out_temp$Cost, l_out_temp$Effect)
#   }
#   # Extract costs and effects from the PSA dataset
#   df_c <- df_ce[, 1:n_str]
#   df_e <- df_ce[, (n_str + 1):(2*n_str)]
#   # Register end time of parallelized PSA
#   n_time_end_psa_parallel <- Sys.time()
# }
# if (os == "windows") {
#   # Initialize cluster object
#   cl <- parallel::makeCluster(no_cores)
#   # Register clusters
#   doParallel::registerDoParallel(cl)
#   opts <- list(attachExportEnv = TRUE)
#   # Run parallelized PSA
#   df_ce <- foreach::foreach(i = 1:n_samp, .combine = rbind,
#                             .export = ls(globalenv()),
#                             .packages = c("dampack"),
#                             .options.snow = opts) %dopar% {
#                               l_out_temp <- calculate_ce_out(df_psa_input[i, ])
#                               df_ce <- c(l_out_temp$Cost, l_out_temp$Effect)
#                             }
#   # Extract costs and effects from the PSA dataset
#   df_c <- df_ce[, 1:n_str]
#   df_e <- df_ce[, (n_str + 1):(2*n_str)]
#   # Register end time of parallelized PSA
#   n_time_end_psa_parallel <- Sys.time()
# }
# if (os == "linux") {
#   # Initialize cluster object
#   cl <- parallel::makeCluster(no_cores)
#   # Register clusters
#   doParallel::registerDoMC(cl)
#   # Run parallelized PSA
#   df_ce <- foreach::foreach(i = 1:n_sim, .combine = rbind) %dopar% {
#     l_out_temp <- calculate_ce_out(df_psa_input[i, ])
#     df_ce <- c(l_out_temp$Cost, l_out_temp$Effect)
#   }
#   # Extract costs and effects from the PSA dataset
#   df_c <- df_ce[, 1:n_str]
#   df_e <- df_ce[, (n_str + 1):(2*n_str)]
#   # Register end time of parallelized PSA
#   n_time_end_psa_parallel <- Sys.time()
# }
# # Stop clusters
# stopCluster(cl)
# n_time_total_psa_parallel <- n_time_end_psa_parallel - n_time_init_psa_parallel
# print(paste0("PSA with ", scales::comma(n_sim), " simulations run in series in ",
#              round(n_time_total_psa_parallel, 2), " ",
#              units(n_time_total_psa_parallel)))

## Visualize PSA results for CEA ----
### Create PSA object ----
#* Function included in "R/Functions.R" The latest version can be found in `dampack` package
l_psa <- make_psa_obj(cost          = df_c, 
                      effectiveness = df_e, 
                      parameters    = df_psa_input, 
                      strategies    = v_names_str)
l_psa$strategies <- v_names_str
colnames(l_psa$effectiveness) <- v_names_str
colnames(l_psa$cost) <- v_names_str

#* Vector with willingness-to-pay (WTP) thresholds.
v_wtp <- seq(0, 200000, by = 5000)

### Cost-Effectiveness Scatter plot ----
txtsize <- 13
#* Function included in "R/Functions.R"; depends on `tidyr` and `ellipse` packages.
#* The latest version can be found in `dampack` package
gg_scattter <- plot.psa(l_psa, txtsize = txtsize) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  scale_y_continuous("Cost (Thousand $)", 
                     breaks = number_ticks(10),
                     labels = function(x) x/1000) +
  xlab("Effectiveness (QALYs)") +
  guides(col = guide_legend(nrow = 2)) +
  theme(legend.position = "bottom")
gg_scattter

### Incremental cost-effectiveness ratios (ICERs) with probabilistic output ----
#* Compute expected costs and effects for each strategy from the PSA
#* Function included in "R/Functions.R". The latest version can be found in `dampack` package
df_out_ce_psa <- summary(l_psa)

#* Function included in "R/Functions.R"; depends on the `dplyr` package
#* The latest version can be found in `dampack` package
df_cea_psa <- calculate_icers(cost       = df_out_ce_psa$meanCost, 
                              effect     = df_out_ce_psa$meanEffect,
                              strategies = df_out_ce_psa$Strategy)
df_cea_psa

### Plot cost-effectiveness frontier with probabilistic output ----
#* Function included in "R/Functions.R"; depends on the `ggplot2`  and `ggrepel` packages.
#* The latest version can be found in `dampack` package
plot.icers(df_cea_psa, label = "all", txtsize = txtsize) +
  expand_limits(x = max(table_cea$QALYs) + 0.1) +
  theme(legend.position = c(0.8, 0.2))

### Cost-effectiveness acceptability curves (CEACs) and frontier (CEAF) ---
#* Functions included in "R/Functions.R". The latest versions can be found in `dampack` package
ceac_obj <- ceac(wtp = v_wtp, psa = l_psa)
#* Regions of highest probability of cost-effectiveness for each strategy
summary(ceac_obj)
#* CEAC & CEAF plot
gg_ceac <- plot.ceac(ceac_obj, txtsize = txtsize, xlim = c(0, NA), n_x_ticks = 14) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  theme(legend.position = c(0.8, 0.48))
gg_ceac

### Expected Loss Curves (ELCs) ----
#* Function included in "R/Functions.R".The latest version can be found in `dampack` package
elc_obj <- calc_exp_loss(wtp = v_wtp, psa = l_psa)
elc_obj

#* ELC plot
gg_elc <- plot.exp_loss(elc_obj, log_y = FALSE, 
               txtsize = txtsize, xlim = c(0, NA), n_x_ticks = 14,
               col = "full") +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  # geom_point(aes(shape = as.name("Strategy"))) +
  scale_y_continuous("Expected Loss (Thousand $)", 
                     breaks = number_ticks(10),
                     labels = function(x) x/1000) +
  theme(legend.position = c(0.4, 0.7),)
gg_elc

### Expected value of perfect information (EVPI) ----
#* Function included in "R/Functions.R". The latest version can be found in `dampack` package
evpi <- calc_evpi(wtp = v_wtp, psa = l_psa)
#* EVPI plot
gg_evpi <- plot.evpi(evpi, effect_units = "QALY", 
                     txtsize = txtsize, xlim = c(0, NA), n_x_ticks = 14) +
  scale_y_continuous("EVPI (Thousand $)", 
                     breaks = number_ticks(10),
                     labels = function(x) x/1000)
gg_evpi

### Combine all figures into one ----
patched_cea <- (gg_scattter +  gg_ceac + plot_layout(guides = "keep"))/(gg_elc + gg_evpi)
gg_psa_plots <- patched_cea + 
  plot_annotation(tag_levels = 'A')
gg_psa_plots

## Visualize PSA results for Epidemiological Measures ----
### Wrangle PSA output ----
#* Combine prevalence measures 
df_prev <- dplyr::bind_rows(df_prevS1,
                            df_prevS2,
                            df_prevS1S2)
#* Transform to long format
df_surv_lng <- reshape2::melt(df_surv, 
                              id.vars = c("Outcome"), 
                              # value.name = "Survival",
                              variable.name = "Time")
df_prev_lng <- reshape2::melt(df_prev, 
                              id.vars = c("States"), 
                              # value.name = "Proportion",
                              variable.name = "Time")
#* Compute posterior-predicted 95% CI
df_surv_summ <- data_summary(df_surv_lng, varname = "value",
                             groupnames = c("Outcome", "Time"))
df_le_summ <- data_summary(df_le, varname = "LE",
                             groupnames = c("Outcome"))
df_prev_summ <- data_summary(df_prev_lng, varname = "value",
                             groupnames = c("States", "Time"))
df_prev_summ$States <- ordered(df_prev_summ$States,
                               levels = c("S1", "S2", "S1 + S2"))

### Plot epidemiological measures ---
txtsize_epi <- 16
#### Survival ---
gg_surv_psa <- ggplot(df_surv_summ, aes(x = as.numeric(Time), y = value, 
                                        ymin = lb, ymax = ub)) +
  geom_line() +
  geom_ribbon(alpha = 0.4) +
  scale_x_continuous(breaks = number_ticks(8)) + 
  xlab("Cycle") +
  ylab("Proportion alive") +
  theme_bw(base_size = txtsize_epi) +
  theme()
gg_surv_psa

#### Life Expectancy ---
gg_le_psa <- ggplot(df_le, aes(x = LE)) +
  geom_density(color = "darkblue", fill = "lightblue") +
  scale_x_continuous(breaks = number_ticks(8)) + 
  xlab("Life expectancy") +
  ylab("") +
  theme_bw(base_size = txtsize_epi) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks = element_blank())
gg_le_psa

#### Prevalence ---
gg_prev_psa <- ggplot(df_prev_summ, aes(x = as.numeric(Time), y = value, 
                                        ymin = lb, ymax = ub,
                                        color = States, linetype = States, 
                                        fill = States)) +
  geom_line() +
  geom_ribbon(alpha = 0.4, color = NA) +
  scale_x_continuous(breaks = number_ticks(8)) + 
  scale_y_continuous(breaks = number_ticks(8),
                     labels = scales::percent_format(accuracy = 1)) + 
  scale_color_discrete(name = "Health state", l = 50) +
  scale_linetype(name = "Health state") +
  scale_fill_discrete(name = "Health state", l = 50) +
  xlab("Cycle") +
  ylab("Prevalence (%)") +
  theme_bw(base_size = 16) +
  theme(legend.position = c(0.83, 0.83))
gg_prev_psa

### Combine all figures into one ----
patched_epi <- (gg_surv_psa / gg_le_psa) | gg_prev_psa
gg_psa_epi_plots <- patched_epi + 
  plot_annotation(tag_levels = 'A')
gg_psa_epi_plots