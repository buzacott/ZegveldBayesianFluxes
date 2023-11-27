co2_bayesian_routine <- function(co2_model_data,
                                 pt_split = 0.7,
                                 rf_split = 0.3) {
  # Define functions here, unfortunately functions can't be passed to the
  # BayesianTools sampler

  #----------------------------------------------------------------------------#
  # Define data split ----
  #----------------------------------------------------------------------------#
  co2_pt_data <- co2_model_data[co2_model_data$Contribution >= pt_split, ]
  co2_rf_data <- co2_model_data[co2_model_data$Contribution < rf_split, ]

  #----------------------------------------------------------------------------#
  # Parameter limits ----
  #----------------------------------------------------------------------------#
  co2_lower <- c(
    pt_r_ref = 0,
    pt_e_0 = 50,
    pt_alpha = 0,
    pt_beta = 0.001,
    pt_alpha_amplitude = 0,
    pt_beta_amplitude = 0,
    pt_phase_shift = 0,
    rf_r_ref = 0,
    rf_e_0 = 50,
    rf_alpha = 0,
    rf_beta = 0,
    rf_alpha_amplitude = 0,
    rf_beta_amplitude = 0,
    rf_phase_shift = 0,
    nee_sd = sd(co2_model_data$NEE) * 0.02
  )

  co2_upper <- c(
    pt_r_ref = max(co2_pt_data[co2_pt_data$night, ]$NEE),
    pt_e_0 = 400,
    pt_alpha = 0.22,
    pt_beta = 250,
    pt_alpha_amplitude = 0.11,
    pt_beta_amplitude = 50,
    pt_phase_shift = 180,
    rf_r_ref = max(co2_rf_data[co2_rf_data$night, ]$NEE),
    rf_e_0 = 400,
    rf_alpha = 0.22,
    rf_beta = 250,
    rf_alpha_amplitude = 0.11,
    rf_beta_amplitude = 50,
    rf_phase_shift = 180,
    nee_sd = sd(co2_model_data$NEE) * 2
  )

  #----------------------------------------------------------------------------#
  # Functions ----
  #----------------------------------------------------------------------------#
  ## CO2 round 1 model ----
  co2_model_r1 <- function(params) {
    pt_r_ref <- params[1]
    pt_e_0 <- params[2]
    pt_alpha <- params[3]
    pt_beta <- params[4]
    pt_alpha_amplitude <- params[5]
    pt_beta_amplitude <- params[6]
    pt_phase_shift <- params[7]
    rf_r_ref <- params[8]
    rf_e_0 <- params[9]
    rf_alpha <- params[10]
    rf_beta <- params[11]
    rf_alpha_amplitude <- params[12]
    rf_beta_amplitude <- params[13]
    rf_phase_shift <- params[14]

    pt_reco <- respiration_lloyd_taylor(
      r_ref = pt_r_ref,
      e_0 = pt_e_0,
      t = co2_pt_data$TS_PT
    )
    pt_gpp <- gpp_sin(
      alpha = pt_alpha,
      beta = pt_beta,
      Rg = co2_pt_data$Rg,
      doy = co2_pt_data$doy,
      phase_shift = pt_phase_shift,
      alpha_amplitude = pt_alpha_amplitude,
      beta_amplitude = pt_beta_amplitude
    )
    pt_nee <- pt_gpp + pt_reco

    rf_reco <- respiration_lloyd_taylor(
      r_ref = rf_r_ref,
      e_0 = rf_e_0,
      t = co2_rf_data$TS_RF
    )
    rf_gpp <- gpp_sin(
      alpha = rf_alpha,
      beta = rf_beta,
      Rg = co2_rf_data$Rg,
      doy = co2_rf_data$doy,
      phase_shift = rf_phase_shift,
      alpha_amplitude = rf_alpha_amplitude,
      beta_amplitude = rf_beta_amplitude
    )
    rf_nee <- rf_gpp + rf_reco

    # Residuals
    pt_residuals <- (co2_pt_data$NEE - pt_nee)
    rf_residuals <- (co2_rf_data$NEE - rf_nee)

    return(list(
      pt_residuals = pt_residuals,
      rf_residuals = rf_residuals
    ))
  }

  ## CO2 round 2 model ----
  co2_model_r2 <- function(params) {
    pt_r_ref <- params[1]
    pt_e_0 <- params[2]
    pt_alpha <- params[3]
    pt_beta <- params[4]
    pt_alpha_amplitude <- params[5]
    pt_beta_amplitude <- params[6]
    pt_phase_shift <- params[7]
    rf_r_ref <- params[8]
    rf_e_0 <- params[9]
    rf_alpha <- params[10]
    rf_beta <- params[11]
    rf_alpha_amplitude <- params[12]
    rf_beta_amplitude <- params[13]
    rf_phase_shift <- params[14]

    pt_reco <- respiration_lloyd_taylor(
      r_ref = pt_r_ref,
      e_0 = pt_e_0,
      t = co2_model_data$TS_PT
    )
    pt_gpp <- gpp_sin(
      alpha = pt_alpha,
      beta = pt_beta,
      Rg = co2_model_data$Rg,
      doy = co2_model_data$doy,
      phase_shift = pt_phase_shift,
      alpha_amplitude = pt_alpha_amplitude,
      beta_amplitude = pt_beta_amplitude
    )
    pt_nee <- pt_gpp + pt_reco

    rf_reco <- respiration_lloyd_taylor(
      r_ref = rf_r_ref,
      e_0 = rf_e_0,
      t = co2_model_data$TS_RF
    )
    rf_gpp <- gpp_sin(
      alpha = rf_alpha,
      beta = rf_beta,
      Rg = co2_model_data$Rg,
      doy = co2_model_data$doy,
      phase_shift = rf_phase_shift,
      alpha_amplitude = rf_alpha_amplitude,
      beta_amplitude = rf_beta_amplitude
    )
    rf_nee <- rf_gpp + rf_reco

    nee_sim <- pt_nee * co2_model_data$Contribution +
      rf_nee * (1 - co2_model_data$Contribution)
    residuals <- nee_sim - co2_model_data$NEE

    return(residuals)
  }

  # Run routine ----
  out <- bayesian_routine(
    f_r1 = co2_model_r1,
    f_r2 = co2_model_r2,
    lower = co2_lower,
    upper = co2_upper
  )

  return(out)
}
