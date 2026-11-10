# The PSA distributions should be centred on the base-case values, so that
# the parameter table, the base-case code and the PSA agree. (The published
# parameter table listed c_trtB as $12,000 instead of the $13,000 used in the
# text, the code and the results; see manuscript/README.md.)

test_that("PSA parameter means match the base-case values", {
  for (m in list(list(env = env_sim, par = l_params_sim),
                 list(env = env_sr,  par = l_params_sr))) {
    df_psa <- m$env$generate_psa_params(n_sim = 10000, seed = 071818)
    expect_true(all(names(df_psa) %in% names(m$par)))
    for (par in names(df_psa)) {
      # 2% relative tolerance: utilities follow beta distributions whose means
      # differ slightly from the base case (e.g., u_H has mean 200/203 vs 1)
      expect_equal(mean(df_psa[[par]]), m$par[[par]],
                   tolerance = 0.02, label = paste("mean of", par))
    }
  }
})

test_that("cost of treatment B is $13,000 in the base case and the PSA", {
  expect_equal(l_params_common$c_trtB, 13000)
  # Gamma mean = shape * scale
  expect_equal(86.2 * 150.8, 13000, tolerance = 1e-3)
})
