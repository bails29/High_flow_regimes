###############################################################################
# 1. Generate synthetic daily data for one hydrological year
###############################################################################
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

source("./functions.R") # load in my functions

set.seed(123)

# Hydrological year: 1 Oct 2020 – 30 Sep 2021 (365 days)
hydro_start <- as.Date("2020-10-01")
N_days      <- 365

daily_synth <- data.table(
  date = hydro_start + 0:(N_days - 1),
  dohy = 1:N_days
)

# Fundamental annual frequency
omega <- 2 * pi / N_days

# Reference "spring" peak day for phase anchoring
peak_day  <- 200
base_flow <- 5

# zero-based time index for harmonics
tt <- daily_synth$dohy - 1

# ---------------------------------------------------------------------------
# DAILY regime: unimodal, mostly A1 but A2 non-negligible
# ---------------------------------------------------------------------------
A1_d  <- 7      # dominant annual cycle
A2_d  <- 3      # non-trivial semiannual structure
A3_d  <- 0.3    # tiny third harmonic (fine structure only)

phi1_d <- 0
phi2_d <- -pi / 6
phi3_d <-  pi / 6

seasonal_daily <- base_flow +
  A1_d * sin(1 * omega * (tt - (peak_day - 1)) + phi1_d) +
  A2_d * sin(2 * omega * (tt - (peak_day - 1)) + phi2_d) +
  A3_d * sin(3 * omega * (tt - (peak_day - 1)) + phi3_d)

# small noise so that A1 + A2 still dominate after smoothing
sd_sig_d   <- sd(seasonal_daily)
noise_sd_d <- 0.15 * sd_sig_d   # adjust factor to get R2 in high 90s

daily_synth[, mean_q := pmax(
  0,
  seasonal_daily + rnorm(.N, mean = 0, sd = noise_sd_d)
)]


# ---------------------------------------------------------------------------
# HOURLY regime: clearly bimodal, A2 very strong
# ---------------------------------------------------------------------------
A1_h  <- 3      # smaller annual component
A2_h  <- 8      # large semiannual component -> strong bimodality
A3_h  <- 0.3    # tiny third harmonic

phi1_h <- 0
phi2_h <-  pi / 3
phi3_h <- -pi / 4

seasonal_hourly <- base_flow +
  A1_h * sin(1 * omega * (tt - (peak_day - 1)) + phi1_h) +
  A2_h * sin(2 * omega * (tt - (peak_day - 1)) + phi2_h) +
  A3_h * sin(3 * omega * (tt - (peak_day - 1)) + phi3_h)

sd_sig_h   <- sd(seasonal_hourly)
noise_sd_h <- 0.15 * sd_sig_h   # again, noise small vs signal

daily_synth[, max_q := pmax(
  0,
  seasonal_hourly + rnorm(.N, mean = 0, sd = noise_sd_h)
)]

# NOTE: We drop the extra spike events (peak_days) here,
# because they introduce non-harmonic local structure that
# harms the R^2 for A1 + A2 after smoothing.

# Quick check plot of the raw synthetic data (before POT + smoothing)
ggplot(daily_synth) +
  geom_line(aes(x = dohy, y = mean_q, colour = "Daily mean")) +
  geom_line(aes(x = dohy, y = max_q,  colour = "Hourly max")) +
  scale_colour_manual(values = c("Daily mean" = "red", "Hourly max" = "black")) +
  labs(x = "Day of hydrologic year", y = "Flow", colour = "") +
  theme_bw()

###############################################################################
# 2. Run get_smoothed_curves() on the synthetic data
###############################################################################

smoothed_curves <- get_smoothed_curves(
  daily       = daily_synth,
  threshold   = 0.8,
  spar        = 0.7,
  smoothing   = 365,
  window_days = 10
)

head(smoothed_curves)

# Plot smoothed curves
ggplot(smoothed_curves)+
  geom_line(aes(x = dohy, y = smoothed_daily), color = "red")+
  geom_line(aes(x = dohy, y = smoothed_hourly), color = "black")+
  labs(y = "Smoothed regime", x = "Day of the hydrologic year")+
  theme_bw()
  

###############################################################################
# 3. Use estimate_harmonic_seasonality() on smoothed curves
###############################################################################

# Seasonality of declustered & smoothed daily POT frequency
seasonality_daily <- estimate_harmonic_seasonality(
  dt       = smoothed_curves,
  dohy_col = "dohy",
  val_col  = "smoothed_daily"
)

# Seasonality of declustered & smoothed hourly-max POT frequency
seasonality_hourly <- estimate_harmonic_seasonality(
  dt       = smoothed_curves,
  dohy_col = "dohy",
  val_col  = "smoothed_hourly"
)

seasonality_daily
seasonality_hourly


 # Extract R² values for plot annotation
r2_df <- data.frame(
  type      = c("Uni-modal", "Uni-modal", "Uni-modal",
                "Bi-modal",  "Bi-modal",  "Bi-modal"),
  component = c("A1", "A1 + A2", "A1 + A2 + A3",
                "A1", "A1 + A2", "A1 + A2 + A3"),
  r2        = c(
    seasonality_daily$r2_first_only,
    seasonality_daily$r2_first_second,
    seasonality_daily$r2_first_second_third,
    seasonality_hourly$r2_first_only,
    seasonality_hourly$r2_first_second,
    seasonality_hourly$r2_first_second_third
  )
)


# Offsets for stacking labels (A1 top, A1+A2 below, A1+A2+A3 lowest)
r2_df$y_offset <- rep(c(0.006, 0.008, 0.01), 2)

#plot reconstructed curves 
# ---------------------------------------------------------------------------
# 3. Reconstruct curves with A1, A1+A2, A1+A2+A3
# ---------------------------------------------------------------------------

N  <- nrow(smoothed_curves)
w1 <- 2 * pi / N
tt <- smoothed_curves$dohy - 1  # zero-based, same convention as your function

# Add daily fits
smoothed_curves <- smoothed_curves %>%
  mutate(
    # Daily: A1 only
    fit_daily_h1 =
      seasonality_daily$mean_level +
      seasonality_daily$a1 * cos(1 * w1 * tt) +
      seasonality_daily$b1 * sin(1 * w1 * tt),
    
    # Daily: A1 + A2
    fit_daily_h12 =
      seasonality_daily$mean_level +
      seasonality_daily$a1 * cos(1 * w1 * tt) +
      seasonality_daily$b1 * sin(1 * w1 * tt) +
      seasonality_daily$a2 * cos(2 * w1 * tt) +
      seasonality_daily$b2 * sin(2 * w1 * tt),
    
    # Daily: A1 + A2 + A3
    fit_daily_h123 =
      seasonality_daily$mean_level +
      seasonality_daily$a1 * cos(1 * w1 * tt) +
      seasonality_daily$b1 * sin(1 * w1 * tt) +
      seasonality_daily$a2 * cos(2 * w1 * tt) +
      seasonality_daily$b2 * sin(2 * w1 * tt) +
      seasonality_daily$a3 * cos(3 * w1 * tt) +
      seasonality_daily$b3 * sin(3 * w1 * tt)
  )

# Add hourly fits
smoothed_curves <- smoothed_curves %>%
  mutate(
    # Hourly: A1 only
    fit_hourly_h1 =
      seasonality_hourly$mean_level +
      seasonality_hourly$a1 * cos(1 * w1 * tt) +
      seasonality_hourly$b1 * sin(1 * w1 * tt),
    
    # Hourly: A1 + A2
    fit_hourly_h12 =
      seasonality_hourly$mean_level +
      seasonality_hourly$a1 * cos(1 * w1 * tt) +
      seasonality_hourly$b1 * sin(1 * w1 * tt) +
      seasonality_hourly$a2 * cos(2 * w1 * tt) +
      seasonality_hourly$b2 * sin(2 * w1 * tt),
    
    # Hourly: A1 + A2 + A3
    fit_hourly_h123 =
      seasonality_hourly$mean_level +
      seasonality_hourly$a1 * cos(1 * w1 * tt) +
      seasonality_hourly$b1 * sin(1 * w1 * tt) +
      seasonality_hourly$a2 * cos(2 * w1 * tt) +
      seasonality_hourly$b2 * sin(2 * w1 * tt) +
      seasonality_hourly$a3 * cos(3 * w1 * tt) +
      seasonality_hourly$b3 * sin(3 * w1 * tt)
  )

# ---------------------------------------------------------------------------
# 4. Long format for facetted plot
# ---------------------------------------------------------------------------

plot_long <- smoothed_curves %>%
  pivot_longer(
    cols = c(
      smoothed_daily, smoothed_hourly,
      fit_daily_h1, fit_daily_h12, fit_daily_h123,
      fit_hourly_h1, fit_hourly_h12, fit_hourly_h123
    ),
    names_to  = "series",
    values_to = "value"
  ) %>%
  mutate(
    type = dplyr::case_when(
      grepl("daily",  series) ~ "B. Uni-modal",
      grepl("hourly", series) ~ "A. Bi-modal",
      TRUE                    ~ NA_character_
    ),
    component = dplyr::case_when(
      series %in% c("smoothed_daily", "smoothed_hourly") ~ "Smoothed regime",
      grepl("_h1$",   series) ~ "A1",
      grepl("_h12$",  series) ~ "A1 + A2",
      grepl("_h123$", series) ~ "A1 + A2 + A3",
      TRUE                    ~ NA_character_
    )
  ) %>%
  filter(!is.na(type), !is.na(component))

# ---------------------------------------------------------------------------
# 5. Facetted plot: daily vs hourly with all harmonic reconstructions
# ---------------------------------------------------------------------------

example <- ggplot(plot_long, aes(x = dohy, y = value,
                      colour = component,
                      linetype = component)) +
  geom_line(size = 1) +
  facet_wrap(~ type, ncol = 1, scales = "free_y") +
  ylim(0, 0.015)+
  scale_color_manual(values = c('#cb181d','#fb6a4a','#fcae91','black'))+
  scale_linetype_manual(values = c("dashed", "dashed", "dashed", "solid"))+
  labs(
    x = "Day of the hydrologic year",
    y = "Regime value",
    colour   = "Fourier Reconstruction",
    linetype = "Fourier Reconstruction"
  ) +
  theme_bw()+
  theme(legend.position = "bottom", 
        legend.title.position = "top", 
        legend.title = element_text(hjust = 0.5)) +
  theme_bw()

ggsave("./SI_harmonic_amp.png", example, units = "mm", height = 150, width = 150)

##############################################################################
#classification, just as a test. 

classified <- classify_seasonality(
  peak_day_fit123  = seasonality_hourly$peak_day_fit123,
  peak2_day_fit123 = seasonality_hourly$peak2_day_fit123,
  A1 = seasonality_daily$amp1,
  A2 = seasonality_daily$amp2, 
  min_peak_separation = 30, 
  A2A1_weak = 0.5
)
classified
