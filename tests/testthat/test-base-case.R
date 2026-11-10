# Regression tests for the base-case cost-effectiveness analyses of the
# simulation-time-dependent and state-residence-time-dependent Sick-Sicker
# models with Simpson's 1/3 (+ 3/8) within-cycle correction (v0.3.0).
# Reference values were produced by analysis/cSTM_time_dep_simulation.R and
# analysis/cSTM_time_dep_state_residence.R after the gen_wcc() fix; the PSA
# wrappers (calculate_ce_out()) reproduce them exactly.

v_names_str <- c("Standard of care", "Strategy A", "Strategy B", "Strategy AB")

test_that("simulation-time model: base-case costs and QALYs match reference values", {
  df_ce <- env_sim$calculate_ce_out(l_params_sim)
  expect_equal(df_ce$Strategy, v_names_str)
  expect_equal(df_ce$Cost,
               c(117126.245546, 219922.676317, 203706.603639, 297850.022748),
               tolerance = 1e-9)
  expect_equal(df_ce$Effect,
               c(19.1992859658, 19.9623547451, 20.5188822358, 21.4230220996),
               tolerance = 1e-9)
})

test_that("state-residence model: base-case costs and QALYs match reference values", {
  df_ce <- env_sr$calculate_ce_out(l_params_sr)
  expect_equal(df_ce$Strategy, v_names_str)
  expect_equal(df_ce$Cost,
               c(103871.786480, 199197.153147, 187435.837556, 275959.928325),
               tolerance = 1e-9)
  expect_equal(df_ce$Effect,
               c(20.3247916838, 21.2083091619, 21.4234469130, 22.4239838857),
               tolerance = 1e-9)
})

test_that("base-case ICERs and dominance match reference values", {
  for (m in list(list(env = env_sim, par = l_params_sim,
                      icer = c(65611.2479688, 104124.8405021)),
                 list(env = env_sr,  par = l_params_sr,
                      icer = c(76060.3043256, 88476.5812566)))) {
    df_ce  <- m$env$calculate_ce_out(m$par)
    df_cea <- suppressWarnings(
      calculate_icers(df_ce$Cost, df_ce$Effect, df_ce$Strategy)
    )
    status <- setNames(df_cea$Status, df_cea$Strategy)
    expect_equal(unname(status[c("Standard of care", "Strategy B",
                                 "Strategy AB", "Strategy A")]),
                 c("ND", "ND", "ND", "D"))
    icer <- setNames(df_cea$ICER, df_cea$Strategy)
    expect_equal(unname(icer[c("Strategy B", "Strategy AB")]), m$icer,
                 tolerance = 1e-9)
  }
})

test_that("PSA wrappers apply the corrected within-cycle correction", {
  # With all utilities set to 1 and no transition disutility, total discounted
  # QALYs equal the discounted, WCC-corrected life years of the cohort trace.
  n_cycles <- 75
  v_dw     <- 1 / (1.03^(0:n_cycles))
  l_par <- modifyList(l_params_sim,
                      list(u_S1 = 1, u_S2 = 1, u_trtA = 1, du_HS1 = 0))
  m_M   <- env_sim$decision_model(l_par)$l_m_M[[1]]
  expect_equal(env_sim$calculate_ce_out(l_par)$Effect[1],
               sum((1 - m_M[, "D"]) * v_dw * gen_wcc(n_cycles)))
  l_par <- modifyList(l_params_sr,
                      list(u_S1 = 1, u_S2 = 1, u_trtA = 1, du_HS1 = 0))
  m_M   <- env_sr$decision_model(l_par)$l_m_M[[1]]
  expect_equal(env_sr$calculate_ce_out(l_par)$Effect[1],
               sum((1 - m_M[, "D"]) * v_dw * gen_wcc(n_cycles)))
})

test_that("PSA wrappers do not depend on objects in the global environment", {
  # calculate_ce_out() used to read v_wcc, v_dwc and v_dwe (state-residence
  # model) and the strategy and state names (both models) from the global
  # environment. Run them in a clean environment.
  for (nm in c("v_wcc", "v_dwc", "v_dwe", "v_names_str", "n_str",
               "v_names_states", "n_states")) {
    expect_false(exists(nm, envir = globalenv(), inherits = FALSE), label = nm)
  }
  expect_no_error(env_sim$calculate_ce_out(l_params_sim))
  expect_no_error(env_sr$calculate_ce_out(l_params_sr))
})

test_that("cost and QALY discount rates are not swapped", {
  # Changing only the cost discount rate must not change QALYs
  for (m in list(list(env = env_sim, par = l_params_sim),
                 list(env = env_sr,  par = l_params_sr))) {
    df_base <- m$env$calculate_ce_out(m$par)
    df_dc   <- m$env$calculate_ce_out(modifyList(m$par, list(d_c = 0.05)))
    expect_equal(df_dc$Effect, df_base$Effect)
    expect_true(all(df_dc$Cost < df_base$Cost))
  }
})

test_that("transition array check tolerates floating-point rounding", {
  # Valid rows (each sums to 1 up to rounding) whose total over 3 states is
  # 3 - 4.4e-16, so an exact `==` comparison with n_states fails
  v_p <- c(0.54756458822449339, 0.21743485789087841, 0.23500055388462818)
  a_P <- array(rep(v_p, each = 3), dim = c(3, 3, 2))
  expect_false(sum(rowSums(a_P[, , 1])) == 3)
  expect_no_error(check_sum_of_transition_array(a_P, n_states = 3, n_cycles = 2,
                                                verbose = FALSE))
  # An invalid array is still rejected
  a_P[1, 1, 2] <- 0.9
  expect_error(check_sum_of_transition_array(a_P, n_states = 3, n_cycles = 2,
                                             verbose = FALSE),
               "not a valid transition array")
})

test_that("state-residence decision model accepts every PSA parameter set", {
  # Regression test: PSA sample 220 used to fail the (exact) transition array
  # check because of floating-point rounding, which stopped the PSA.
  df_psa <- env_sr$generate_psa_params(n_sim = 1000, seed = 071818)
  l_par  <- update_param_list(l_params_sr, df_psa[220, ])
  expect_no_error(invisible(capture.output(env_sr$decision_model(l_par))))
})
