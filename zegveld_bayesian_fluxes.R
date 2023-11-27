#------------------------------------------------------------------------------#
# Main script that runs the modelling procedure as described in the manuscript:
# Buzacott et al. A Bayesian inference approach to determine experimental
# Typha latifolia paludiculture greenhouse gas exchange measured
# with eddy covariance
#------------------------------------------------------------------------------#

library(dplyr)
library(tibble)
library(readr)
library(ggplot2)
library(lubridate)
library(tidyr)
library(stringr)
library(randomForest)
library(data.table)
library(BayesianTools)

set.seed(123)

#------------------------------------------------------------------------------#
source("src/bayesian_core.R")
source("src/modelling_functions.R")
source("src/co2_bayesian_routine.R")
source("src/ch4_bayesian_routine.R")

#------------------------------------------------------------------------------#
# Create necessary directories
make_dir <- function(d) {
  if (!dir.exists(d)) dir.create(d)
}
make_dir("output")
make_dir("figures")
make_dir("results")

#------------------------------------------------------------------------------#
# Load data ----
df <- read_csv("data/2021_2022_zegveld.csv")

#------------------------------------------------------------------------------#
# Gapfill contribution ----
rf_gf_fname <- file.path("output", "cb_gapfilled_rf.rds")
if (!file.exists(rf_gf_fname)) {
  # Select data to train RF
  rf_cb_data <- df %>%
    select(WD_F, WS_F, doy, Contribution) %>%
    drop_na()

  # Split into train/test
  idx <- sort(sample(seq_len(nrow(rf_cb_data)), size = nrow(rf_cb_data) * 0.7))
  train <- rf_cb_data[idx, ]
  test <- rf_cb_data[-idx, ]

  # Train model
  rf_cb <- randomForest(Contribution ~ ., data = train)

  # Predict
  rf_cb_test <- predict(rf_cb, newdata = test)

  # Plot
  p1 <- ggplot() +
    geom_point(aes(train$Contribution, rf_cb$predicted), size = 0.5) +
    geom_abline(col = "red") +
    ylim(0, 1) +
    xlim(0, 1) +
    labs(
      x = "Obs",
      y = "Sim",
      title = bquote(
        paste(
          Train ~ R^2: .(round(cor(train$Contribution, rf_cb$predicted)^2, 2))
        )
      )
    )
  p2 <- ggplot() +
    geom_point(aes(test$Contribution, rf_cb_test), size = 0.5) +
    geom_abline(col = "red") +
    ylim(0, 1) +
    xlim(0, 1) +
    labs(
      x = "Obs",
      y = "Sim",
      title = bquote(
        paste(
          Test ~ R^2: .(round(cor(test$Contribution, rf_cb_test)^2, 2))
        )
      )
    )

  ggsave(file.path("figures", "cb_gf_rf.png"),
         egg::ggarrange(p1, p2),
         width = 4,
         height = 6,
         dpi = 300)

  # Save predictions
  cb_pred <- predict(rf_cb, newdata = df)
  write_rds(cb_pred, rf_gf_fname)
} else {
  cb_pred <- read_rds(rf_gf_fname)
}

df <- df %>%
  mutate(Contribution_GF = if_else(is.na(Contribution), cb_pred, Contribution))

write_rds(df, file.path("output", "df.rds"), compress = "gz")

#------------------------------------------------------------------------------#
# Years to model ----
years <- 2021:2022

#------------------------------------------------------------------------------#
# CO2 modelling ----
co2_model_data <- df %>%
  select(datetime,
         year,
         season,
         doy,
         night,
         Contribution,
         Rg = PPFD_IN_F, # Use photosynthetically active radiation for lrc
         TS_PT,
         TS_RF,
         NEE) %>%
  drop_na() %>%
  as.data.table()

write_rds(co2_model_data,
          file.path("output", "co2_model_data.rds"),
          compress = "gz")

for (i in seq_along(years)) {
  co2_r_fname <- file.path("output", paste0("co2_r", years[i], ".rds"))
  if (!file.exists(co2_r_fname)) {
    # Run routine for selected year
    co2_ry <- co2_bayesian_routine(co2_model_data %>% filter(year == years[i]))
    write_rds(co2_ry, co2_r_fname, compress = "gz")
  }
}

#------------------------------------------------------------------------------#
# CH4 modelling ----
ch4_model_data <- df %>%
  select(datetime,
         year,
         season,
         doy,
         night,
         Contribution,
         FCH4,
         TS_PT,
         TS_RF,
         WL_PT,
         WL_RF) %>%
  drop_na() %>%
  as.data.table()

write_rds(ch4_model_data,
          file.path("output", "ch4_model_data.rds"),
          compress = "gz")

for (i in seq_along(years)) {
  ch4_r_fname <- file.path("output", paste0("ch4_r", years[i], ".rds"))
  if (!file.exists(ch4_r_fname)) {
    # Run routine for selected year
    ch4_ry <- ch4_bayesian_routine(ch4_model_data %>% filter(year == years[i]))
    write_rds(ch4_ry, ch4_r_fname, compress = "gz")
  }
}
#------------------------------------------------------------------------------#
