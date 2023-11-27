# CO2 modelling functions ----
gpp <- function(alpha, beta, Rg) {
  # GPP light response curve
  #   alpha <num>: light use efficiency parameter
  #   beta <num>: GPP max parameter
  #   Rg <num vec>: global radiation vector
  return(-(alpha * beta * Rg) / (alpha * Rg + beta))
}

gpp_sin <- function(
    alpha,
    beta,
    Rg,
    doy,
    phase_shift,
    alpha_amplitude,
    beta_amplitude,
    alpha_min = 0,
    beta_min = 0.001) {
  # GPP predictions with annual sine modifications to alpha/beta
  #   alpha <num>: light use efficiency parameter
  #   beta <num>: GPP max parameter
  #   Rg <num vec>: global radiation vector, or PPFD_IN
  #   doy <num vec>: day of the year vector
  #   phase_shift <num>: offset for sine curve (phi parameter)
  #   alpha_amplitude, beta_amplitude <num>: sine offsets for alpha and beta
  #   alpha_min <num>: force minimum value for alpha (prevent GPP emissions)
  #   beta_min <num>: force minimum value for beta for stability
  sin_term <- (sin(2 * pi * (doy - phase_shift) / 365))
  alpha_sin <- alpha + sin_term * alpha_amplitude
  beta_sin <- beta + sin_term * beta_amplitude
  alpha_sin[alpha_sin < alpha_min] <- alpha_min
  beta_sin[beta_sin < beta_min] <- beta_min

  gpp <- gpp(alpha_sin, beta_sin, Rg)

  return(gpp)
}

respiration_lloyd_taylor <- function(
    r_ref,
    e_0,
    t,
    t_ref = 15,
    t_0 = -46.02) {
  # Lloyd-Taylor soil respiration function
  #   r_ref: Respiration rate at reference temperature
  #   e_0: Activation energy parameter
  #   t: soil temperature (K)
  #   t_ref: reference temperature (K) (10C)
  #   t_0: temperature as fitted by LloydTaylor (1994)

  r_eco <- r_ref * exp(e_0 * (1 / (t_ref - t_0) - 1 / (t - t_0)))

  return(r_eco)
}

# CH4 modelling functions ----
fch4_temp <- function(temp, a, b) {
  # Exponential temperature CH4 response function
  #   temp <num vec>: vector of temperature
  #   a, b <num>: power law coefficients
  return(a * exp(b * temp))
}

fch4_temp_wl <- function(temp, wl, a, b, k) {
  # Exponential temperature CH4 response function with sigmoid water level func
  #   temp <num vec>: vector of temperature
  #   wl <num vec>: vector of water levels
  #   a, b <num>: power law coefficients
  #   k <num> decay constant for wl sigmoid curve
  return(a * exp(b * temp) * 1 / (1 + exp(-wl))^k)
}
