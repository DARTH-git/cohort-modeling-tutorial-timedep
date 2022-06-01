#------------------------------------------------------------------------------#
####                         Decision Model                                 ####
#------------------------------------------------------------------------------#
#' Decision Model
#'
#' \code{decision_model} implements the decision model used.
#'
#' @param l_params_all List with all parameters of decision model
#' @param verbose Logical variable to indicate print out of messages
#' @return The transition probability array and the cohort trace matrix.
#' @export
decision_model <- function(l_params_all, verbose = FALSE) {
  with(as.list(l_params_all), {
    ########################### Process model inputs ###########################
    ##* Number of cycles
    n_cycles <- (n_age_max - n_age_init)/cycle_length # time horizon, number of cycles
    ## Age-specific transition probabilities to the Dead state
    # compute mortality rates
    v_r_S1Dage <- v_r_HDage * hr_S1        # Age-specific mortality rate in the Sick state 
    v_r_S2Dage <- v_r_HDage * hr_S2        # Age-specific mortality rate in the Sicker state 
    #* transform rates to probabilities adjusting by cycle length
    #* Function included in "R/Functions.R". The latest version can be found in `darthtools` package
    p_HS1  <- rate_to_prob(r = r_HS1, t = cycle_length) # constant annual probability of becoming Sick when Healthy conditional on surviving 
    p_S1H  <- rate_to_prob(r = r_S1H, t = cycle_length) # constant annual probability of becoming Healthy when Sick conditional on surviving
    p_S1S2 <- rate_to_prob(r = r_S1S2, t = cycle_length)# constant annual probability of becoming Sicker when Sick conditional on surviving
    v_p_HDage  <- rate_to_prob(v_r_HDage, t = cycle_length)  # Age-specific mortality risk in the Healthy state 
    v_p_S1Dage <- rate_to_prob(v_r_S1Dage, t = cycle_length) # Age-specific mortality risk in the Sick state
    v_p_S2Dage <- rate_to_prob(v_r_S2Dage, t = cycle_length) # Age-specific mortality risk in the Sicker state
    
    ##* Annual transition probability of becoming Sicker when Sick for treatment B
    #* Apply hazard ratio to rate to obtain transition rate of becoming Sicker when 
    #* Sick for treatment B
    r_S1S2_trtB <- r_S1S2 * hr_S1S2_trtB
    #* Transform rate to probability to become Sicker when Sick under treatment B 
    #* adjusting by cycle length conditional on surviving
    #* (Function included in "R/Functions.R". The latest version can be found in 
    #* `darthtools` package)
    p_S1S2_trtB <- rate_to_prob(r = r_S1S2_trtB, t = cycle_length)
    
    ###################### Construct state-transition models ###################
    ##* Initial state vector
    #* All starting healthy
    v_m_init <- c(H = 1, S1 = 0, S2 = 0, D = 0) # initial state vector
    #* Number of health states 
    n_states    <- length(v_m_init)
    #* Health state names
    v_names_states <- names(v_m_init)
    
    ##* Initialize cohort trace for age-dependent (ad) cSTM for strategies SoC and A
    m_M_SoC <- matrix(0, 
                      nrow     = (n_cycles + 1), ncol = n_states, 
                      dimnames = list(0:n_cycles, v_names_states))
    #* Store the initial state vector in the first row of the cohort trace
    m_M_SoC[1, ] <- v_m_init
    ##* Initialize cohort trace for strategies A, B, and AB
    #* Structure and initial states are the same as for SoC
    m_M_strA  <- m_M_SoC # Strategy A
    m_M_strB  <- m_M_SoC # Strategy B
    m_M_strAB <- m_M_SoC # Strategy AB
    
    ##* Initialize transition dynamics arrays which will capture transitions 
    ##* from each state to another over time for strategy SoC
    a_A_SoC <- array(0,
                     dim      = c(n_states, n_states, n_cycles + 1),
                     dimnames = list(v_names_states, v_names_states, 0:n_cycles))
    #* Set first slice of a_A_SoC with the initial state vector in its diagonal
    diag(a_A_SoC[, , 1]) <- v_m_init
    #* Initialize transition-dynamics array for strategies A, B, and AB
    #* Structure and initial states are the same as for SoC
    a_A_strA  <- a_A_SoC
    a_A_strB  <- a_A_SoC
    a_A_strAB <- a_A_SoC
    
    ##* Create transition arrays
    #* Initialize 3-D array
    a_P_SoC <- array(0, dim   = c(n_states, n_states, n_cycles),
                 dimnames = list(v_names_states, v_names_states, 0:(n_cycles - 1)))
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
    
    ## Initialize transition probability matrix for strategy A as a copy of SoC's
    a_P_strA <- a_P_SoC
    
    ### For strategies B and AB
    ## Initialize transition probability array for strategies B and AB
    a_P_strB <- a_P_SoC
    ## Only need to update the probabilities involving the transition from Sick to Sicker, p_S1S2
    # From S1
    a_P_strB["S1", "S1", ] <- (1 - v_p_S1Dage) * (1 - (p_S1H + p_S1S2_trtB))
    a_P_strB["S1", "S2", ] <- (1 - v_p_S1Dage) * p_S1S2_trtB
    ## Initialize transition probability matrix for strategy AB as a copy of B's
    a_P_strAB <- a_P_strB
    
    ### Check if transition probability matrices are valid
    ## Check that transition probabilities are [0, 1]
    check_transition_probability(a_P_SoC, verbose = TRUE)
    check_transition_probability(a_P_strA, verbose = TRUE)
    check_transition_probability(a_P_strB, verbose = TRUE)
    check_transition_probability(a_P_strAB, verbose = TRUE)
    ### Check that all rows for each slice of the array sum to 1
    check_sum_of_transition_array(a_P_SoC, n_states = n_states, n_cycles = n_cycles, verbose = TRUE)
    check_sum_of_transition_array(a_P_strA, n_states = n_states, n_cycles = n_cycles, verbose = TRUE)
    check_sum_of_transition_array(a_P_strB, n_states = n_states, n_cycles = n_cycles, verbose = TRUE)
    check_sum_of_transition_array(a_P_strAB, n_states = n_states, n_cycles = n_cycles, verbose = TRUE)
    
    #### Run Markov model ####
    ## Iterative solution of age-dependent cSTM
    for (t in 1:n_cycles) {
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
    
    ## Store the cohort traces in a list
    l_m_M <- list(m_M_SoC,
                  m_M_strA,
                  m_M_strB,
                  m_M_strAB)
    names(l_m_M) <- v_names_str
    
    ## Store the transition array for each strategy in a list
    l_a_A <- list(a_A_SoC,
                  a_A_strA,
                  a_A_strB,
                  a_A_strAB)
    names(l_m_M) <- v_names_str
    
    ########################################## RETURN OUTPUT  ##########################################
    out <- list(l_m_M = l_m_M,
                l_a_A = l_a_A)
    
    return(out)
  }
  )
}

#------------------------------------------------------------------------------#
####              Calculate cost-effectiveness outcomes                     ####
#------------------------------------------------------------------------------#
#' Calculate cost-effectiveness outcomes
#'
#' \code{calculate_ce_out} calculates costs and effects for a given vector of 
#' parameters using a simulation model.
#' @param l_params_all List with all parameters of decision model
#' @param n_wtp Willingness-to-pay threshold to compute net benefits
#' @return A data frame with discounted costs, effectiveness and NMB.
#' @export
calculate_ce_out <- function(l_params_all, n_wtp = 100000){ # User defined
  with(as.list(l_params_all), {
    
    ### Run decision model to get transition dynamics array
    model <- decision_model(l_params_all = l_params_all)
    l_a_A <- model[["l_a_A"]]
    
    #### State Rewards ####
    ## Vector of state utilities under strategy SoC
    v_u_SoC    <- c(H  = u_H, 
                    S1 = u_S1, 
                    S2 = u_S2, 
                    D  = u_D) * cycle_length
    ## Vector of state costs under strategy SoC
    v_c_SoC    <- c(H  = c_H, 
                    S1 = c_S1,
                    S2 = c_S2, 
                    D  = c_D) * cycle_length
    ## Vector of state utilities under strategy A
    v_u_strA   <- c(H  = u_H, 
                    S1 = u_trtA, 
                    S2 = u_S2, 
                    D  = u_D) * cycle_length
    ## Vector of state costs under strategy A
    v_c_strA   <- c(H  = c_H, 
                    S1 = c_S1 + c_trtA,
                    S2 = c_S2 + c_trtA, 
                    D  = c_D) * cycle_length
    ## Vector of state utilities under strategy B
    v_u_strB   <- c(H  = u_H, 
                    S1 = u_S1, 
                    S2 = u_S2, 
                    D  = u_D) * cycle_length
    ## Vector of state costs under strategy B
    v_c_strB   <- c(H  = c_H, 
                    S1 = c_S1 + c_trtB, 
                    S2 = c_S2 + c_trtB, 
                    D  = c_D) * cycle_length
    ## Vector of state utilities under strategy AB
    v_u_strAB  <- c(H  = u_H, 
                    S1 = u_trtA, 
                    S2 = u_S2, 
                    D  = u_D) * cycle_length
    ## Vector of state costs under strategy AB
    v_c_strAB  <- c(H  = c_H, 
                    S1 = c_S1 + (c_trtA + c_trtB), 
                    S2 = c_S2 + (c_trtA + c_trtB), 
                    D  = c_D) * cycle_length
    
    ## Store the vectors of state utilities for each strategy in a list 
    l_u   <- list(SoC = v_u_SoC,
                  A  = v_u_strA,
                  B  = v_u_strB,
                  AB = v_u_strAB)
    ## Store the vectors of state cost for each strategy in a list 
    l_c   <- list(SoC = v_c_SoC,
                  A  = v_c_strA,
                  B  = v_c_strB,
                  AB = v_c_strAB)
    
    # assign strategy names to matching items in the lists
    names(l_u) <- names(l_c) <- names(l_a_A) <- v_names_str
    
    ## create empty vectors to store total utilities and costs 
    v_tot_qaly <- v_tot_cost <- vector(mode = "numeric", length = n_str)
    names(v_tot_qaly) <- names(v_tot_cost) <- v_names_str
    
    ## Number of cycles
    n_cycles <- (n_age_max - n_age_init)/cycle_length # time horizon, number of cycles
    
    ## Discount weight for costs and effects
    v_dwc  <- 1 / ((1 + d_e * cycle_length) ^ (0:n_cycles))
    v_dwe  <- 1 / ((1 + d_c * cycle_length) ^ (0:n_cycles))
    
    ## Within-cycle correction (WCC) using Simpson's 1/3 rule
    v_wcc <- darthtools::gen_wcc(n_cycles = n_cycles, 
                                 method = "Simpson1/3") # vector of wcc
    
    #### Loop through each strategy and calculate total utilities and costs ####
    for (i in 1:n_str) { # i <- 1
      v_u_str <- l_u[[i]]   # select the vector of state utilities for the i-th strategy
      v_c_str <- l_c[[i]]   # select the vector of state costs for the i-th strategy
      a_A_str <- l_a_A[[i]] # select the transition array for the i-th strategy
      
      #### Array of state rewards ####
      # Create transition matrices of state utilities and state costs for the i-th strategy 
      m_u_str   <- matrix(v_u_str, nrow = n_states, ncol = n_states, byrow = T)
      m_c_str   <- matrix(v_c_str, nrow = n_states, ncol = n_states, byrow = T)
      # Expand the transition matrix of state utilities across cycles to form a transition array of state utilities
      a_R_u_str <- array(m_u_str, 
                         dim      = c(n_states, n_states, n_cycles + 1),
                         dimnames = list(v_names_states, v_names_states, 0:n_cycles))
      # Expand the transition matrix of state costs across cycles to form a transition array of state costs
      a_R_c_str <- array(m_c_str, 
                         dim      = c(n_states, n_states, n_cycles + 1),
                         dimnames = list(v_names_states, v_names_states, 0:n_cycles))
      
      #### Apply transition rewards ####  
      # Apply disutility due to transition from H to S1
      a_R_u_str["H", "S1", ]      <- a_R_u_str["H", "S1", ]       - du_HS1
      # Add transition cost per cycle due to transition from H to S1
      a_R_c_str["H", "S1", ]      <- a_R_c_str["H", "S1", ]       + ic_HS1
      # Add transition cost  per cycle of dying from all non-dead states
      a_R_c_str[-n_states, "D", ] <- a_R_c_str[-n_states, "D", ] + ic_D
      
      #### Expected QALYs and Costs for all transitions per cycle ####
      # QALYs = life years x QoL
      # Note: all parameters are annual in our example. In case your own case example is different make sure you correctly apply .
      a_Y_c_str <- a_A_str * a_R_c_str
      a_Y_u_str <- a_A_str * a_R_u_str 
      
      #### Expected QALYs and Costs per cycle ####
      ## Vector of QALYs and Costs
      v_qaly_str <- apply(a_Y_u_str, 3, sum) # sum the proportion of the cohort across transitions 
      v_cost_str <- apply(a_Y_c_str, 3, sum) # sum the proportion of the cohort across transitions
      
      #### Discounted total expected QALYs and Costs per strategy and apply half-cycle correction if applicable ####
      ## QALYs
      v_tot_qaly[i] <- t(v_qaly_str) %*% (v_dwe * v_wcc)
      ## Costs
      v_tot_cost[i] <- t(v_cost_str) %*% (v_dwc * v_wcc)
    }
    
    ## Vector with discounted net monetary benefits (NMB)
    v_nmb <- v_tot_qaly * n_wtp - v_tot_cost
    
    ## data.frame with discounted costs, effectiveness and NMB
    df_ce <- data.frame(Strategy = v_names_str,
                        Cost     = v_tot_cost,
                        Effect   = v_tot_qaly,
                        NMB      = v_nmb)
    
    return(df_ce)
  }
  )
}

#------------------------------------------------------------------------------#
####                Generate Epidemiological Measures                       ####
#------------------------------------------------------------------------------#
generate_epi_measures_SoC <- function(l_params_all){ # User defined
  with(as.list(l_params_all), {
    ### Run decision model to get cohort trace and transition dynamics array
    model <- decision_model(l_params_all = l_params_all)
    m_M_SoC <- model$l_m_M$`Standard of care`
    ## Survival curve
    v_S_SoC <- rowSums(m_M_SoC[, -which(v_names_states == "D")])
    ## Life expectancy
    le <- sum(v_S_SoC)
    ## Prevalence
    # Prevalence of Sick
    v_prev_S1_SoC   <- m_M_SoC[, "S1"] / v_S_SoC
    # Prevalence of Sicker
    v_prev_S2_SoC   <- m_M_SoC[, "S2"] / v_S_SoC
    # Prevalence of Sick and Sicker
    v_prev_S1S2_SoC <- rowSums(m_M_SoC[, c("S1", "S2")])/v_S_SoC
    l_out_epi <- list(S  = v_S_SoC,
                      LE = le,
                      PrevS1   = v_prev_S1_SoC,
                      PrevS2   = v_prev_S2_SoC,
                      PrevS1S2 = v_prev_S1S2_SoC)
    return(l_out_epi)
  }
 )
}

#------------------------------------------------------------------------------#
####             Generate a PSA input parameter dataset                     ####
#------------------------------------------------------------------------------#
#' Generate parameter sets for the probabilistic sensitivity analysis (PSA)
#'
#' \code{generate_psa_params} generates a PSA dataset of the parameters of the 
#' cost-effectiveness analysis.
#' @param n_sim Number of parameter sets for the PSA dataset
#' @param seed Seed for the random number generation
#' @return A data.frame with a PSA dataset of he parameters of the 
#' cost-effectiveness analysis
#' @export
generate_psa_params <- function(n_sim = 1000, seed = 071818){
  set.seed(seed) # set a seed to be able to reproduce the same results
  df_psa <- data.frame(
    # Transition probabilities (per cycle)
    r_HS1    = rgamma(n_sim, shape = 30, rate = 170 + 30), # constant rate of becoming Sick when Healthy conditional on surviving
    r_S1H    = rgamma(n_sim, shape = 60, rate = 60 + 60),  # constant rate of becoming Healthy when Sick conditional on surviving
    r_S1S2   = rgamma(n_sim, shape = 84, rate = 716 + 84), # constant rate of becoming Sicker when Sick conditional on surviving
    hr_S1    = rlnorm(n_sim, meanlog = log(3), sdlog = 0.01),  # hazard ratio of death in Sick vs healthy
    hr_S2    = rlnorm(n_sim, meanlog = log(10), sdlog = 0.02), # hazard ratio of death in Sicker vs healthy 
    
    # Effectiveness of treatment B 
    hr_S1S2_trtB = rlnorm(n_sim, meanlog = log(0.6), sdlog = 0.02), # hazard ratio of becoming Sicker when Sick under treatment B
    
    ## State rewards
    # Costs
    c_H    = rgamma(n_sim, shape = 100,   scale = 20),   # cost of remaining one cycle in state H
    c_S1   = rgamma(n_sim, shape = 177.8, scale = 22.5), # cost of remaining one cycle in state S1
    c_S2   = rgamma(n_sim, shape = 225,   scale = 66.7), # cost of remaining one cycle in state S2
    c_trtA = rgamma(n_sim, shape = 73.5, scale = 163.3), # cost of treatment A (per cycle)
    c_trtB = rgamma(n_sim, shape = 86.2, scale = 150.8), # cost of treatment B (per cycle)
    c_D    = 0,                                          # cost of being in the death state
    
    # Utilities
    u_H    = rbeta(n_sim, shape1 = 200, shape2 = 3),     # utility when healthy
    u_S1   = rbeta(n_sim, shape1 = 130, shape2 = 45),    # utility when sick
    u_S2   = rbeta(n_sim, shape1 = 230, shape2 = 230),  # utility when sicker
    u_D    = 0,                                          # utility when dead
    u_trtA = rbeta(n_sim, shape1 = 300, shape2 = 15),    # utility when being treated
    
    ## Transition rewards
    du_HS1 = rbeta(n_sim, shape1 = 11,  shape2 = 1088),  # disutility when transitioning from Healthy to Sick
    ic_HS1 = rgamma(n_sim, shape = 25,  scale = 40),     # increase in cost when transitioning from Healthy to Sick
    ic_D   = rgamma(n_sim, shape = 100, scale = 20)      # increase in cost when dying
  )
  return(df_psa)
}