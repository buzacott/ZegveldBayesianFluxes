ch4_bayesian_routine <- function(ch4_model_data,
                                 pt_split = 0.7,
                                 rf_split = 0.3) {
  # Define functions here, unfortunately functions can't be passed directly
  # to the BayesianTools sampler

  #----------------------------------------------------------------------------#
  # Define data split
  #----------------------------------------------------------------------------#
  ch4_pt_data <- ch4_model_data[ch4_model_data$Contribution >= pt_split, ]
  ch4_rf_data <- ch4_model_data[ch4_model_data$Contribution < rf_split, ]

  #----------------------------------------------------------------------------#
  # Define parameter limits
  #----------------------------------------------------------------------------#
  ch4_lower <- c(
    pt_a = 0,
    pt_b = 0,
    pt_k = 0,
    rf_a = 0,
    rf_b = 0,
    rf_k = 0,
    sd_err = sd(ch4_model_data$FCH4) * 0.02
  )

  ch4_upper <- c(
    pt_a = 0.5,
    pt_b = 0.5,
    pt_k = 1,
    rf_a = 0.5,
    rf_b = 0.5,
    rf_k = 1,
    sd_err = sd(ch4_model_data$FCH4) * 2
  )

  #----------------------------------------------------------------------------#
  # Functions
  #----------------------------------------------------------------------------#
  # CH4 round 1 model
  ch4_model_r1 <- function(params) {
    pt_fch4 <- fch4_temp_wl(
      ch4_pt_data$TS_PT,
      ch4_pt_data$WL_PT,
      params[1], params[2], params[3]
    ) # pt_a, pt_b, pt_k
    rf_fch4 <- fch4_temp_wl(
      ch4_rf_data$TS_RF,
      ch4_rf_data$WL_RF,
      params[4], params[5], params[6]
    ) # rf_a, rf_b, rf_k

    pt_residuals <- ch4_pt_data$FCH4 - pt_fch4
    rf_residuals <- ch4_rf_data$FCH4 - rf_fch4

    return(list(
      pt_residuals = pt_residuals,
      rf_residuals = rf_residuals
    ))
  }

  ch4_model_r2 <- function(params) {
    pt_fch4 <- fch4_temp_wl(
      ch4_model_data$TS_PT,
      ch4_model_data$WL_PT,
      params[1], params[2], params[3]
    ) # pt_a, pt_b, pt_k
    rf_fch4 <- fch4_temp_wl(
      ch4_model_data$TS_RF,
      ch4_model_data$WL_RF,
      params[4], params[5], params[6]
    ) # rf_a, rf_b, rf_k

    fch4_pred <- pt_fch4 * ch4_model_data$Contribution +
      rf_fch4 * (1 - ch4_model_data$Contribution)

    residuals <- ch4_model_data$FCH4 - fch4_pred

    return(residuals)
  }

  # Run routine ----
  out <- bayesian_routine(
    f_r1 = ch4_model_r1,
    f_r2 = ch4_model_r2,
    lower = ch4_lower,
    upper = ch4_upper
  )

  return(out)
}
