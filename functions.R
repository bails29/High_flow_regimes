#' Build Smoothed Hydrological Event Frequency Curves with Declustering
#'
#' Given a daily time series for a single catchment, this function identifies
#' high-flow exceedance events based on user-defined thresholds, declusters
#' them to avoid multiple counts of the same event, and then constructs
#' **smoothed 365-day frequency curves** for both mean daily flow and daily
#' maximum flow.
#'
#' The workflow is:
#' \enumerate{
#'   \item Compute POT-style thresholds from the daily series:
#'         \code{mean_q} and \code{max_q} (e.g., 0.99 quantile).
#'   \item Extract exceedance series for daily mean and daily maximum flows.
#'   \item Decluster exceedances using \code{\link{decluster_events}}, keeping
#'         at most one event within a \code{window_days}-day window.
#'   \item Count declustered events per hydrological day (\code{dohy}).
#'   \item Convert raw counts to relative frequencies (scaled to sum to 1).
#'   \item Wrap the series at the start and end of the year and apply spline
#'         smoothing to obtain smooth 365-point curves.
#' }
#'
#' The output provides both raw and smoothed daily/hourly (daily max) event
#' frequencies on a hydrological day-of-year axis and a reconstructed
#' (arbitrary) calendar date axis.
#'
#' @param daily `data.frame` or `data.table` with at least the columns:
#'   \itemize{
#'     \item \code{date}: Date or POSIXct time stamp.
#'     \item \code{dohy}: integer day-of-hydrological-year (1–365).
#'     \item \code{mean_q}: daily mean discharge (or other daily metric).
#'     \item \code{max_q}: daily maximum discharge (or other daily metric).
#'   }
#'   It is assumed that \code{daily} contains one full hydrological year.
#' @param threshold Numeric in (0, 1). Quantile used to define exceedance
#'   thresholds for \code{mean_q} and \code{max_q} (e.g., 0.99 for top 1\%).
#' @param spar Numeric smoothing parameter passed to \code{\link{smooth.spline}}.
#'   Larger values yield smoother curves. I recommend 0.7 as a starting point.
#' @param smoothing Integer. Number of days at each end of the year used for
#'   wrapping before smoothing (e.g., 365 to wrap a full year).
#' @param window_days Integer. Declustering window length (in days) used in
#'   \code{\link{decluster_events}}. Within a moving window of this width,
#'   only the largest exceedance is retained. suggest 10 days as a starting point if working with floods.
#'
#' @details
#' The declustering step reduces dependence among exceedance events by ensuring
#' that closely spaced peaks are not counted multiple times. Daily and hourly
#' (daily max) series are declustered independently but on the same time axis.
#'
#' Event frequencies are scaled separately for daily and hourly exceedances so
#' that the sums of \code{freq_daily_scaled} and \code{freq_hourly_scaled} over
#' the 365 days each equal 1 (if any events exist).
#'
#' Smoothing is performed on an extended domain where the first and last
#' \code{smoothing} days are wrapped around the year to reduce edge effects.
#' Negative smoothed values, if any, are truncated to 0.
#'
#' The function currently assumes a 365-day hydrological year. If your setup
#' uses 366 days or a different calendar, adjust \code{dohy} and wrapping logic
#' accordingly.
#'
#' @return A `tibble`/`data.frame` with 365 rows and the following columns:
#'   \itemize{
#'     \item \code{dohy}: day-of-hydrological-year (1–365).
#'     \item \code{day}: corresponding (nominal) calendar day (1–365).
#'     \item \code{freq_daily}, \code{freq_hourly}: declustered event counts
#'           per \code{dohy}.
#'     \item \code{freq_daily_scaled}, \code{freq_hourly_scaled}: counts scaled
#'           to sum to 1 across the year.
#'     \item \code{smoothed_daily}, \code{smoothed_hourly}: spline-smoothed
#'           frequency curves on \code{dohy}.
#'     \item \code{hydro_date}: synthetic hydrological date constructed from
#'           \code{dohy} assuming a start at 1 October 2020. This is an arbitrary date added here for plotting purposes
#'   }
#'
#' @seealso \code{\link{decluster_events}} for the declustering logic.
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' library(tidyr)
#'
#' # daily must contain: date, dohy, mean_q, max_q
#' smoothed <- get_smoothed_curves(
#'   daily       = daily_data,
#'   threshold   = 0.99,
#'   spar        = 0.7,
#'   smoothing   = 365,
#'   window_days = 10
#' )
#'
#' head(smoothed)
#' }
#'
#' @export
#' 

get_smoothed_curves <- function(daily,
                                threshold = 0.99,
                                spar = 0.7,
                                smoothing = 365, 
                                window_days = 10) {
  
  # 1) set thresholds once for this daily table
  thr_daily  <- quantile(daily$mean_q, threshold, na.rm = TRUE)
  thr_hourly <- quantile(daily$max_q,  threshold, na.rm = TRUE)
  
  # 2) build exceedance tables (these will be declustered)
  #    important: your decluster_events() expects columns named `date` and a value col called `q`
  exc_daily <- daily %>%
    dplyr::filter(mean_q >= thr_daily) %>%
    dplyr::select(date, dohy, q = mean_q)
  
  exc_hourly <- daily %>%
    dplyr::filter(max_q >= thr_hourly) %>%
    dplyr::select(date, dohy, q = max_q)
  
  # 3) DECLUSTER here using your function defined above
  #    (10-day window like your colleague’s example — change if needed)
  exc_daily_decl  <- decluster_events(exc_daily,
                                      time_col = "date",
                                      value_col = "q",
                                      window_days = window_days)
  exc_hourly_decl <- decluster_events(exc_hourly,
                                      time_col = "date",
                                      value_col = "q",
                                      window_days = window_days)
  
  # 4) count declustered events per hydrological day
  #    we now have at most 1 event per 10 days (for each series separately)
  flow_sum <- dplyr::full_join(
    exc_daily_decl  %>% dplyr::count(dohy, name = "freq_daily"),
    exc_hourly_decl %>% dplyr::count(dohy, name = "freq_hourly"),
    by = "dohy"
  ) %>%
    dplyr::mutate(
      freq_daily  = tidyr::replace_na(freq_daily, 0),
      freq_hourly = tidyr::replace_na(freq_hourly, 0)
    )
  
  # 5) make sure we have 1..365 in DOHY, even if no events
  flow_sum <- flow_sum %>%
    tidyr::complete(dohy = 1:365,
                    fill = list(freq_daily = 0, freq_hourly = 0))
  
  # reconstruct calendar day from dohy 
  flow_sum <- flow_sum %>%
    dplyr::mutate(
      day = dohy + 274,
      day = ifelse(day > 365, day - 365, day)
    )
  
  # 6) scale
  total_daily  <- sum(flow_sum$freq_daily)
  total_hourly <- sum(flow_sum$freq_hourly)
  
  flow_sum <- flow_sum %>%
    dplyr::mutate(
      freq_daily_scaled  = if (total_daily  > 0) freq_daily  / total_daily  else 0,
      freq_hourly_scaled = if (total_hourly > 0) freq_hourly / total_hourly else 0
    )
  
  # 7) wrap for smoothing (same as your version)
  end_part <- flow_sum %>%
    dplyr::filter(dohy %in% (365 - smoothing + 1):365) %>%
    dplyr::mutate(dohy = dohy - 365)
  
  start_part <- flow_sum %>%
    dplyr::filter(dohy %in% 1:smoothing) %>%
    dplyr::mutate(dohy = dohy + 365)
  
  ext <- dplyr::bind_rows(end_part, flow_sum, start_part) %>%
    dplyr::arrange(dohy)
  
  # 8) smooth daily + hourly
  spl_daily  <- smooth.spline(ext$dohy, ext$freq_daily_scaled,  spar = spar)
  spl_hourly <- smooth.spline(ext$dohy, ext$freq_hourly_scaled, spar = spar)
  
  ext$smoothed_daily  <- predict(spl_daily,  ext$dohy)$y
  ext$smoothed_hourly <- predict(spl_hourly, ext$dohy)$y
  
  ext$smoothed_daily[ext$smoothed_daily < 0]   <- 0
  ext$smoothed_hourly[ext$smoothed_hourly < 0] <- 0
  
  # 9) keep only real 1..365 again
  out <- ext %>%
    dplyr::filter(dohy >= 1, dohy <= 365) %>%
    dplyr::arrange(dohy)
  
  # 10) add hydro_date here so every gauge df is self-contained
  hydro_start <- as.Date("2020-10-01")
  out <- out %>%
    dplyr::mutate(hydro_date = hydro_start + dohy - 1)
  
  out
}


#' Decluster Peak-Over-Threshold (POT) Events on a Time Axis
#'
#' Implements a simple **time-window declustering** scheme for exceedance
#' (POT) events. Given a data frame of exceedances (already filtered above a
#' threshold), the function retains at most one event within any
#' \code{window_days}-day window, keeping the event with the larger value in
#' \code{value_col}.
#'
#' This is useful for hydrological and environmental extremes where multiple
#' nearby peaks are likely part of the same physical event and should not be
#' counted separately in frequency analyses.
#'
#' @param df `data.frame` or `data.table` containing exceedance events.
#'   Must already be filtered so that each row represents an exceedance.
#' @param time_col Character string giving the name of the time column (e.g.
#'   `"date"`). The column should be `Date` or `POSIXct`-like and sortable.
#' @param value_col Character string giving the name of the variable whose
#'   magnitude defines which event to keep when two events fall within the
#'   same declustering window (e.g. discharge or precipitation).
#' @param window_days Integer. Length of the declustering window in days.
#'   If two consecutive events are separated by \code{<= window_days}, only
#'   the larger one (in \code{value_col}) is retained.
#'
#' @details
#' The function first orders \code{df} by \code{time_col} and then walks
#' through the events sequentially. When two consecutive events are closer
#' than or equal to \code{window_days} apart, the event with the larger
#' \code{value_col} is kept and the smaller one is discarded. The process
#' is repeated until all events satisfy the minimum spacing constraint.
#'
#' This is a simple, greedy declustering approach suitable for many POT
#' applications. For more complex storm tracking or multi-variable criteria,
#' users may wish to adapt the logic.
#'
#' @return A filtered `data.frame` consisting of declustered exceedance events,
#'   with the same columns as the input \code{df} (no additional columns).
#'
#' @examples
#' # df_exc has columns: date (Date), q (flow), dohy, etc.
#' df_decl <- decluster_events(
#'   df         = df_exc,
#'   time_col   = "date",
#'   value_col  = "q",
#'   window_days = 10
#' )
#'
#' nrow(df_exc)     # original number of exceedances
#' nrow(df_decl)    # reduced after declustering
#'
#' @export


# THIS FUNCITON IS USED INSIDE THE PREVIOUS ONE AS A HELPER
decluster_events <- function(df, time_col = "date", value_col = "q", window_days = 10) {
  # df must already be filtered to exceedances
  df <- df[order(df[[time_col]]), , drop = FALSE]
  df$keep <- TRUE
  
  i <- 2
  while (i <= nrow(df)) {
    # are two events closer than window_days?
    if (as.numeric(df[[time_col]][i] - df[[time_col]][i - 1]) <= window_days) {
      # drop the smaller one
      if (df[[value_col]][i] > df[[value_col]][i - 1]) {
        df$keep[i - 1] <- FALSE
      } else {
        df$keep[i] <- FALSE
      }
      df <- df[df$keep, , drop = FALSE]
    } else {
      i <- i + 1
    }
  }
  df$keep <- NULL
  df
}




#' Estimate Seasonal Harmonics and Peak Timing from a Daily Time Series
#'
#' Fits a truncated discrete Fourier series (first three harmonics) to a seasonal
#' hydrological or environmental time series and extracts robust seasonality
#' metrics — including amplitudes, peak timing, variance explained, and indicators
#' of potential bimodality.
#'
#' The model assumes one full hydrological year (length \eqn{N}, not required to be
#' 365) and expresses the fitted curve as:
#' \deqn{
#' x(t) = a_0 +
#'        a_1 \cos(\omega t) + b_1 \sin(\omega t) +
#'        a_2 \cos(2 \omega t) + b_2 \sin(2 \omega t) +
#'        a_3 \cos(3 \omega t) + b_3 \sin(3 \omega t),
#' }
#' where \eqn{\omega = 2\pi/N}.
#'
#' ## Returned metrics include:
#' **Seasonality strength**
#' * `amp1`, `amp2`, `amp3` — absolute amplitudes of the 1st–3rd harmonics  
#' * `rel_amp1`, `rel_amp2`, `rel_amp3` — amplitudes scaled by mean → strength of seasonal variation  
#'
#' **Peak timing**
#' * `peak_day_h1`, `peak_day_h2`, `peak_day_h3` — harmonic peak days *if each harmonic acted alone*  
#' * `peak_day_fit123` — dominant peak timing of the full 3-harmonic fit  
#' * `peak2_day_fit123` — timing of a secondary peak if present; `NA` if unimodal  
#'
#' **Prominence filtering**
#' The secondary peak is only reported if it rises at least  
#' `peak2_prom_frac * mean_level` above the mean. Otherwise it is ignored as
#' hydrologically insignificant (default: 20% above the mean).
#'
#' **Explained variance**
#' * `r2_first_only`           — fraction of variance explained by the 1st harmonic (annual cycle)  
#' * `r2_first_second`         — fraction explained by 1st + 2nd harmonics  
#' * `r2_first_second_third`   — fraction explained by 1st–3rd harmonics  
#'
#' **Bimodality indicator**
#' * `bimodality_ratio = amp2 / amp1` — strength of semiannual structure  
#'   (large values suggest possible bimodality when combined with a meaningful 2nd peak)
#'
#' **Fourier coefficients**
#' * (`a1`, `b1`, `a2`, `b2`, `a3`, `b3`) — raw harmonic components if further diagnostics are needed
#'
#' @param dt `data.table` or `data.frame` containing one full-year seasonal cycle.
#' @param dohy_col Name of column with day-of-hydrological-year index (1..N).
#' @param val_col  Name of column containing daily values to analyze.
#' @param peak2_prom_frac Minimum required prominence (relative to mean) for a
#'   secondary peak to be considered hydrologically meaningful. Default: `0.2`.
#'
#' @details
#' Harmonic phases determine individual component peak timing; peak timing of the
#' **full** seasonal curve is determined directly from the fitted reconstruction.
#' Prominence filtering ensures that secondary peaks reflect physically distinct
#' seasons rather than small shoulders of a single peak.
#'
#' @return A one-row `data.table` of seasonal harmonic metrics.
#'
#' @examples
#' sm <- data.table(dohy = 1:365, smoothed_daily = runif(365))
#' estimate_harmonic_seasonality(sm)
#'
#' @export

estimate_harmonic_seasonality <- function(dt,
                                          dohy_col = "dohy",
                                          val_col  = "smoothed_daily",
                                          peak2_prom_frac = 0.05) {
  
  x <- dt[[val_col]]
  t <- dt[[dohy_col]]
  
  N <- length(x)
  if (length(t) != N) stop("dohy and smoothed_daily must have same length")
  
  w1 <- 2 * pi / N
  a0 <- mean(x, na.rm = TRUE)
  
  # ---- Fourier coefficients ----
  get_coeffs <- function(k) {
    tt <- (t - 1)
    list(
      a = (2 / N) * sum(x * cos(k * w1 * tt)),
      b = (2 / N) * sum(x * sin(k * w1 * tt))
    )
  }
  c1 <- get_coeffs(1)
  c2 <- get_coeffs(2)
  c3 <- get_coeffs(3)
  
  A1 <- sqrt(c1$a^2 + c1$b^2)
  A2 <- sqrt(c2$a^2 + c2$b^2)
  A3 <- sqrt(c3$a^2 + c3$b^2)
  
  phi1 <- atan2(c1$b, c1$a)
  phi2 <- atan2(c2$b, c2$a)
  phi3 <- atan2(c3$b, c3$a)
  
  peak1 <- ((phi1 / (w1 * 1) + 1 - 1) %% N) + 1
  peak2 <- ((phi2 / (w1 * 2) + 1 - 1) %% N) + 1
  peak3 <- ((phi3 / (w1 * 3) + 1 - 1) %% N) + 1
  
  # ---- Reconstruct fits ----
  tt <- (t - 1)
  
  fit1 <- a0 +
    c1$a * cos(1 * w1 * tt) + c1$b * sin(1 * w1 * tt)
  
  fit12 <- a0 +
    c1$a * cos(1 * w1 * tt) + c1$b * sin(1 * w1 * tt) +
    c2$a * cos(2 * w1 * tt) + c2$b * sin(2 * w1 * tt)
  
  fit123 <- a0 +
    c1$a * cos(1 * w1 * tt) + c1$b * sin(1 * w1 * tt) +
    c2$a * cos(2 * w1 * tt) + c2$b * sin(2 * w1 * tt) +
    c3$a * cos(3 * w1 * tt) + c3$b * sin(3 * w1 * tt)
  
  # ---- Variance explained (R^2) ----
  xvar <- var(x)
  if (is.na(xvar) || xvar <= 0) {
    r2_1   <- NA_real_
    r2_12  <- NA_real_
    r2_123 <- NA_real_
  } else {
    r2_1   <- 1 - var(x - fit1)   / xvar
    r2_12  <- 1 - var(x - fit12)  / xvar
    r2_123 <- 1 - var(x - fit123) / xvar
  }
  
  relA1 <- if (a0 > 0) A1 / a0 else NA
  relA2 <- if (a0 > 0) A2 / a0 else NA
  relA3 <- if (a0 > 0) A3 / a0 else NA
  bimodality_ratio <- if (A1 > 0) A2 / A1 else NA
  
  # ---- Peak detection on reconstructed 3-harmonic fit ----
  n <- length(fit123)
  idx_peak1   <- which.max(fit123)
  peak_day_fit123 <- t[idx_peak1]
  peak1_val   <- fit123[idx_peak1]
  
  # Identify all candidate peaks (circular)
  v_ext <- c(fit123[n], fit123, fit123[1])
  cand_idx <- which(
    v_ext[2:(n+1)] > v_ext[1:n] &
      v_ext[2:(n+1)] > v_ext[3:(n+2)]
  )
  
  # Remove primary peak
  cand_idx <- setdiff(cand_idx, idx_peak1)
  
  if (length(cand_idx) == 0) {
    peak2_day_fit123 <- NA_real_
  } else {
    cand_vals <- fit123[cand_idx]
    ord <- order(cand_vals, decreasing = TRUE)
    idx_peak2 <- cand_idx[ord[1]]
    peak2_val <- fit123[idx_peak2]
    
    # Prominence filter relative to mean
    prominence_threshold <- a0 * peak2_prom_frac
    
    if ((peak2_val - a0) >= prominence_threshold) {
      peak2_day_fit123 <- t[idx_peak2]
    } else {
      peak2_day_fit123 <- NA_real_
    }
  }
  
  # ---- Return all metrics ----
  out <- data.table(
    N_days               = N,
    mean_level           = a0,
    amp1                 = A1,
    amp2                 = A2,
    amp3                 = A3,
    rel_amp1             = relA1,
    rel_amp2             = relA2,
    rel_amp3             = relA3,
    peak_day_h1          = peak1,
    peak_day_h2          = peak2,
    peak_day_h3          = peak3,
    peak_day_fit123      = peak_day_fit123,
    peak2_day_fit123     = peak2_day_fit123,
    r2_first_only        = r2_1,
    r2_first_second      = r2_12,
    r2_first_second_third = r2_123,
    bimodality_ratio     = bimodality_ratio,
    a1                   = c1$a,
    b1                   = c1$b,
    a2                   = c2$a,
    b2                   = c2$b,
    a3                   = c3$a,
    b3                   = c3$b
  )
  return(out)
}


#' Classify seasonal regime based on harmonic peak timing and amplitudes
#'
#' Uses the primary and secondary peak timing from the reconstructed 3-harmonic
#' fit (peak_day_fit123 and peak2_day_fit123) along with the harmonic ratio
#' A2/A1 to classify seasonality as:
#' * "unimodal"
#' * "weakly_bimodal"
#' * "bimodal"
#'
#' @param peak_day_fit123 Day-of-year of the dominant peak from fit123
#' @param peak2_day_fit123 Day-of-year of the secondary peak (or NA)
#' @param A1 Amplitude of 1st harmonic
#' @param A2 Amplitude of 2nd harmonic
#' @param min_peak_separation Minimum separation (days) for peaks to be treated as distinct (default 60)
#' @param A2A1_weak Threshold on A2/A1 to distinguish weak from clear bimodality (default 0.7)
#'
#' @return One of "unimodal", "weakly_bimodal", "bimodal"
#'
classify_seasonality <- function(peak_day_fit123,
                                 peak2_day_fit123,
                                 A1, A2,
                                 min_peak_separation = 60,
                                 A2A1_weak = 0.7) {
  
  # Calculate harmonic ratio
  A2A1 <- if (A1 > 0) A2 / A1 else 0
  
  # If no second peak -> unimodal
  if (is.na(peak2_day_fit123)) {
    return("unimodal")
  }
  
  # Compute peak separation with wrap-around correction
  sep <- abs(peak_day_fit123 - peak2_day_fit123)
  sep <- min(sep, 365 - sep)
  
  # If peaks too close -> unimodal (just shoulders)
  if (sep < min_peak_separation) {
    return("unimodal")
  }
  
  # Two distinct hydrologically meaningful peaks
  # Strength of bimodality determined by harmonic ratio size
  if (A2A1 >= A2A1_weak) {
    return("bimodal")
  } else {
    return("weakly_bimodal")
  }
}

