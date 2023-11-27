bayesian_routine <- function(f_r1,
                             f_r2,
                             lower,
                             upper,
                             r1_iter = 2 * 10^5,
                             r2_iter = 10^5) {
  # Bayesian core function
  #   f_r1 <func>: function to apply for round 1
  #   f_r1 <func>: function to apply for round 2
  #   lower <num vec>: lower limit for parameter ranges
  #   upper <num vec>: upper limit for parameter ranges
  #   r1_iter <int>: number of iterations for round 1
  #   r2_iter <int>: number of iterations for round 2
  #----------------------------------------------------------------------------#
  # Round 1
  #----------------------------------------------------------------------------#
  # Define the round 1 likelihood function
  likelihood_r1 <- function(params) {
    residuals <- f_r1(params)

    ll_pt <- dnorm(residuals$pt_residuals,
      sd = params[length(params)],
      log = TRUE
    )

    ll_rf <- dnorm(residuals$rf_residuals,
      sd = params[length(params)],
      log = TRUE
    )

    ll <- c(ll_pt, ll_rf)

    return(sum(ll))
  }

  # Round 1 prior, simple uniform distribution
  r1_prior <- createUniformPrior(
    lower = unlist(lower),
    upper = unlist(upper)
  )

  # Create Bayesian set up
  bayesian_setup_r1 <- createBayesianSetup(
    likelihood_r1,
    r1_prior,
    names = names(lower)
  )

  # Set values for MCMC sampler
  settings <- list(iterations = r1_iter)

  out_r1 <- runMCMC(
    bayesianSetup = bayesian_setup_r1,
    sampler = "DEzs",
    settings = settings
  )

  #----------------------------------------------------------------------------#
  # Analyse round 1 posterior
  #----------------------------------------------------------------------------#
  # Posterior distribution of round 1
  r1_post <- getSample(out_r1) %>%
    as_tibble()

  # Get the highest density interval (HDI) of round 1, credible interval is 0.95
  r1_post_hdr <- r1_post %>%
    pivot_longer(everything()) %>%
    mutate(name = factor(name, colnames(r1_post))) %>%
    group_by(name) %>%
    summarise(hdr = list(HDInterval::hdi(value))) %>%
    unnest_wider(hdr)

  # Filter R1 posterior distribution by HDI
  r1_post_filt <- r1_post %>%
    mutate(n = seq_len(n())) %>%
    pivot_longer(-n) %>%
    left_join(r1_post_hdr) %>%
    filter(value > lower, value < upper) %>%
    pivot_wider(id_cols = n) %>%
    select(-n) %>%
    drop_na()

  # Get the mean and sd of the filtered r1 posterior
  r1_post_filt_msd <- r1_post_filt %>%
    pivot_longer(everything()) %>%
    mutate(name = factor(name, levels = names(lower))) %>%
    group_by(name) %>%
    summarise(
      mean = mean(value),
      sd = sd(value)
    )

  # r1_post_filt_density <- r1_post_filt_msd %>%
  #   group_by(name) %>%
  #   nest() %>%
  #   mutate(data = purrr::map(data, function(p) {
  #     tibble(n = rnorm(nrow(prior_sample), p$mean, p$sd))
  #   })) %>%
  #   unnest(data)

  #----------------------------------------------------------------------------#
  # Round 2
  #----------------------------------------------------------------------------#
  # Define the likelihood function
  likelihood_r2 <- function(params) {
    residuals <- f_r2(params)

    ll <- dnorm(residuals,
      sd = params[length(params)],
      log = TRUE
    )

    return(sum(ll))
  }

  # Informative prior for round 2 is truncated normal prior of the filtered
  # posterior from round 1
  r2_prior <- createTruncatedNormalPrior(
    r1_post_filt_msd$mean,
    r1_post_filt_msd$sd * 2, # double the sd to allow for more flexible fitting
    lower = unlist(lower),
    upper = unlist(upper)
  )

  # Create Bayesian set up round 2
  bayesian_setup_r2 <- createBayesianSetup(
    likelihood_r2,
    r2_prior,
    names = names(lower),
  )

  # Set values for MCMC sampler
  settings <- list(iterations = r2_iter)

  out_r2 <- runMCMC(
    bayesianSetup = bayesian_setup_r2,
    sampler = "DEzs",
    settings = settings
  )

  return(out_r2)
}
