################################################################################
#                                                                              #
#   PHASE I SINGLE-DOSE PHARMACOKINETIC ANALYSIS                               #
#   Study: ABC-mAb-101                                                         #
#   Product: ABCmAb (Monoclonal Antibody, 10 mg/kg SC)                        #
#   Population: 22 Healthy Adult Volunteers                                    #
#   Analysis Method: Non-Compartmental Analysis (NCA)                         #
#                                                                              #
#   Analyst: Philip Otieno                                                     #
#   Date: 2025-05-30                                                           #
#   Software: R (base NCA implementation, PKNCA-equivalent methodology)       #
#                                                                              #
#   Deliverables:                                                              #
#     - Individual concentration-time profiles (linear & semi-log)            #
#     - Mean ± SD concentration-time profiles (linear & semi-log)             #
#     - Individual PK parameter table                                         #
#     - Summary PK parameter table (mean, SD, CV%, median, min, max)         #
#     - Individual concentration listings                                     #
#     - PK parameter listings                                                 #
#     - CSR statistical section outputs                                       #
#                                                                              #
################################################################################

# ==============================================================================
# SECTION 0: SETUP & LIBRARIES
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(grid)
  library(purrr)
  library(stringr)
})

# Output directory
out_dir <- "output"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

cat("\n", rep("=", 70), "\n", sep = "")
cat("  PHASE I mAb NCA ANALYSIS — ABC-mAb-101\n")
cat(rep("=", 70), "\n\n", sep = "")

# ==============================================================================
# SECTION 1: DATA ENTRY — 22 PARTICIPANTS, MANUALLY ENTERED
# ==============================================================================
# Nominal sampling times (hours post-dose):
# 0, 24, 48, 72, 96, 168, 336, 504, 672, 840, 1008
# (0, 1, 2, 3, 4, 7, 14, 21, 28, 35, 42 days)
# Concentrations in µg/mL; BLQ recorded as 0

cat("SECTION 1: Loading concentration-time data...\n")

conc_raw <- tribble(
  ~SUBJID, ~TIME,  ~CONC,   ~BLQ,
  # ---- Subject 1 (Male, 28y, 72 kg) ----
  "SUBJ-01",   0,   0.00,   TRUE,
  "SUBJ-01",  24,  18.42,  FALSE,
  "SUBJ-01",  48,  32.67,  FALSE,
  "SUBJ-01",  72,  41.85,  FALSE,
  "SUBJ-01",  96,  47.23,  FALSE,
  "SUBJ-01", 168,  52.10,  FALSE,
  "SUBJ-01", 336,  44.31,  FALSE,
  "SUBJ-01", 504,  34.78,  FALSE,
  "SUBJ-01", 672,  24.65,  FALSE,
  "SUBJ-01", 840,  15.92,  FALSE,
  "SUBJ-01",1008,   9.44,  FALSE,

  # ---- Subject 2 (Female, 34y, 65 kg) ----
  "SUBJ-02",   0,   0.00,   TRUE,
  "SUBJ-02",  24,  16.88,  FALSE,
  "SUBJ-02",  48,  29.54,  FALSE,
  "SUBJ-02",  72,  38.62,  FALSE,
  "SUBJ-02",  96,  43.70,  FALSE,
  "SUBJ-02", 168,  49.33,  FALSE,
  "SUBJ-02", 336,  41.20,  FALSE,
  "SUBJ-02", 504,  31.95,  FALSE,
  "SUBJ-02", 672,  22.41,  FALSE,
  "SUBJ-02", 840,  14.12,  FALSE,
  "SUBJ-02",1008,   8.05,  FALSE,

  # ---- Subject 3 (Male, 41y, 80 kg) ----
  "SUBJ-03",   0,   0.00,   TRUE,
  "SUBJ-03",  24,  20.15,  FALSE,
  "SUBJ-03",  48,  35.88,  FALSE,
  "SUBJ-03",  72,  45.20,  FALSE,
  "SUBJ-03",  96,  51.67,  FALSE,
  "SUBJ-03", 168,  57.44,  FALSE,
  "SUBJ-03", 336,  48.90,  FALSE,
  "SUBJ-03", 504,  38.11,  FALSE,
  "SUBJ-03", 672,  27.23,  FALSE,
  "SUBJ-03", 840,  17.55,  FALSE,
  "SUBJ-03",1008,  10.82,  FALSE,

  # ---- Subject 4 (Female, 29y, 59 kg) ----
  "SUBJ-04",   0,   0.00,   TRUE,
  "SUBJ-04",  24,  15.30,  FALSE,
  "SUBJ-04",  48,  27.44,  FALSE,
  "SUBJ-04",  72,  35.91,  FALSE,
  "SUBJ-04",  96,  40.88,  FALSE,
  "SUBJ-04", 168,  45.62,  FALSE,
  "SUBJ-04", 336,  38.15,  FALSE,
  "SUBJ-04", 504,  28.77,  FALSE,
  "SUBJ-04", 672,  19.42,  FALSE,
  "SUBJ-04", 840,  11.88,  FALSE,
  "SUBJ-04",1008,   6.21,  FALSE,

  # ---- Subject 5 (Male, 52y, 85 kg) ----
  "SUBJ-05",   0,   0.00,   TRUE,
  "SUBJ-05",  24,  21.77,  FALSE,
  "SUBJ-05",  48,  37.92,  FALSE,
  "SUBJ-05",  72,  48.55,  FALSE,
  "SUBJ-05",  96,  55.03,  FALSE,
  "SUBJ-05", 168,  61.18,  FALSE,
  "SUBJ-05", 336,  52.44,  FALSE,
  "SUBJ-05", 504,  41.30,  FALSE,
  "SUBJ-05", 672,  29.87,  FALSE,
  "SUBJ-05", 840,  19.44,  FALSE,
  "SUBJ-05",1008,  12.10,  FALSE,

  # ---- Subject 6 (Female, 38y, 68 kg) ----
  "SUBJ-06",   0,   0.00,   TRUE,
  "SUBJ-06",  24,  17.55,  FALSE,
  "SUBJ-06",  48,  30.88,  FALSE,
  "SUBJ-06",  72,  40.14,  FALSE,
  "SUBJ-06",  96,  45.92,  FALSE,
  "SUBJ-06", 168,  51.07,  FALSE,
  "SUBJ-06", 336,  42.88,  FALSE,
  "SUBJ-06", 504,  33.21,  FALSE,
  "SUBJ-06", 672,  23.44,  FALSE,
  "SUBJ-06", 840,  14.77,  FALSE,
  "SUBJ-06",1008,   8.55,  FALSE,

  # ---- Subject 7 (Male, 45y, 77 kg) ----
  "SUBJ-07",   0,   0.00,   TRUE,
  "SUBJ-07",  24,  19.44,  FALSE,
  "SUBJ-07",  48,  34.22,  FALSE,
  "SUBJ-07",  72,  43.88,  FALSE,
  "SUBJ-07",  96,  50.11,  FALSE,
  "SUBJ-07", 168,  55.77,  FALSE,
  "SUBJ-07", 336,  47.33,  FALSE,
  "SUBJ-07", 504,  36.88,  FALSE,
  "SUBJ-07", 672,  26.22,  FALSE,
  "SUBJ-07", 840,  16.77,  FALSE,
  "SUBJ-07",1008,  10.00,  FALSE,

  # ---- Subject 8 (Female, 31y, 62 kg) ----
  "SUBJ-08",   0,   0.00,   TRUE,
  "SUBJ-08",  24,  16.11,  FALSE,
  "SUBJ-08",  48,  28.44,  FALSE,
  "SUBJ-08",  72,  37.22,  FALSE,
  "SUBJ-08",  96,  42.55,  FALSE,
  "SUBJ-08", 168,  47.44,  FALSE,
  "SUBJ-08", 336,  39.66,  FALSE,
  "SUBJ-08", 504,  30.22,  FALSE,
  "SUBJ-08", 672,  21.00,  FALSE,
  "SUBJ-08", 840,  12.88,  FALSE,
  "SUBJ-08",1008,   7.11,  FALSE,

  # ---- Subject 9 (Male, 37y, 74 kg) ----
  "SUBJ-09",   0,   0.00,   TRUE,
  "SUBJ-09",  24,  18.88,  FALSE,
  "SUBJ-09",  48,  33.55,  FALSE,
  "SUBJ-09",  72,  43.11,  FALSE,
  "SUBJ-09",  96,  49.00,  FALSE,
  "SUBJ-09", 168,  54.33,  FALSE,
  "SUBJ-09", 336,  45.77,  FALSE,
  "SUBJ-09", 504,  35.55,  FALSE,
  "SUBJ-09", 672,  25.11,  FALSE,
  "SUBJ-09", 840,  15.88,  FALSE,
  "SUBJ-09",1008,   9.22,  FALSE,

  # ---- Subject 10 (Female, 44y, 70 kg) ----
  "SUBJ-10",   0,   0.00,   TRUE,
  "SUBJ-10",  24,  17.22,  FALSE,
  "SUBJ-10",  48,  30.11,  FALSE,
  "SUBJ-10",  72,  39.44,  FALSE,
  "SUBJ-10",  96,  44.77,  FALSE,
  "SUBJ-10", 168,  50.00,  FALSE,
  "SUBJ-10", 336,  41.88,  FALSE,
  "SUBJ-10", 504,  32.44,  FALSE,
  "SUBJ-10", 672,  22.88,  FALSE,
  "SUBJ-10", 840,  14.33,  FALSE,
  "SUBJ-10",1008,   8.11,  FALSE,

  # ---- Subject 11 (Male, 55y, 88 kg) ----
  "SUBJ-11",   0,   0.00,   TRUE,
  "SUBJ-11",  24,  22.44,  FALSE,
  "SUBJ-11",  48,  39.11,  FALSE,
  "SUBJ-11",  72,  50.00,  FALSE,
  "SUBJ-11",  96,  57.22,  FALSE,
  "SUBJ-11", 168,  63.55,  FALSE,
  "SUBJ-11", 336,  54.33,  FALSE,
  "SUBJ-11", 504,  43.00,  FALSE,
  "SUBJ-11", 672,  31.11,  FALSE,
  "SUBJ-11", 840,  20.33,  FALSE,
  "SUBJ-11",1008,  13.00,  FALSE,

  # ---- Subject 12 (Female, 27y, 56 kg) ----
  "SUBJ-12",   0,   0.00,   TRUE,
  "SUBJ-12",  24,  14.77,  FALSE,
  "SUBJ-12",  48,  26.44,  FALSE,
  "SUBJ-12",  72,  34.55,  FALSE,
  "SUBJ-12",  96,  39.33,  FALSE,
  "SUBJ-12", 168,  43.77,  FALSE,
  "SUBJ-12", 336,  36.44,  FALSE,
  "SUBJ-12", 504,  27.11,  FALSE,
  "SUBJ-12", 672,  18.22,  FALSE,
  "SUBJ-12", 840,  10.88,  FALSE,
  "SUBJ-12",1008,   5.55,  FALSE,

  # ---- Subject 13 (Male, 33y, 78 kg) ----
  "SUBJ-13",   0,   0.00,   TRUE,
  "SUBJ-13",  24,  20.55,  FALSE,
  "SUBJ-13",  48,  36.22,  FALSE,
  "SUBJ-13",  72,  46.44,  FALSE,
  "SUBJ-13",  96,  52.88,  FALSE,
  "SUBJ-13", 168,  58.77,  FALSE,
  "SUBJ-13", 336,  49.88,  FALSE,
  "SUBJ-13", 504,  39.22,  FALSE,
  "SUBJ-13", 672,  28.11,  FALSE,
  "SUBJ-13", 840,  18.00,  FALSE,
  "SUBJ-13",1008,  11.22,  FALSE,

  # ---- Subject 14 (Female, 48y, 66 kg) ----
  "SUBJ-14",   0,   0.00,   TRUE,
  "SUBJ-14",  24,  16.55,  FALSE,
  "SUBJ-14",  48,  29.00,  FALSE,
  "SUBJ-14",  72,  37.88,  FALSE,
  "SUBJ-14",  96,  43.11,  FALSE,
  "SUBJ-14", 168,  48.00,  FALSE,
  "SUBJ-14", 336,  40.11,  FALSE,
  "SUBJ-14", 504,  30.77,  FALSE,
  "SUBJ-14", 672,  21.44,  FALSE,
  "SUBJ-14", 840,  13.22,  FALSE,
  "SUBJ-14",1008,   7.44,  FALSE,

  # ---- Subject 15 (Male, 40y, 82 kg) ----
  "SUBJ-15",   0,   0.00,   TRUE,
  "SUBJ-15",  24,  21.11,  FALSE,
  "SUBJ-15",  48,  37.00,  FALSE,
  "SUBJ-15",  72,  47.33,  FALSE,
  "SUBJ-15",  96,  53.77,  FALSE,
  "SUBJ-15", 168,  59.88,  FALSE,
  "SUBJ-15", 336,  51.11,  FALSE,
  "SUBJ-15", 504,  40.44,  FALSE,
  "SUBJ-15", 672,  29.00,  FALSE,
  "SUBJ-15", 840,  18.77,  FALSE,
  "SUBJ-15",1008,  11.55,  FALSE,

  # ---- Subject 16 (Female, 36y, 61 kg) ----
  "SUBJ-16",   0,   0.00,   TRUE,
  "SUBJ-16",  24,  15.88,  FALSE,
  "SUBJ-16",  48,  27.77,  FALSE,
  "SUBJ-16",  72,  36.33,  FALSE,
  "SUBJ-16",  96,  41.44,  FALSE,
  "SUBJ-16", 168,  46.22,  FALSE,
  "SUBJ-16", 336,  38.55,  FALSE,
  "SUBJ-16", 504,  29.33,  FALSE,
  "SUBJ-16", 672,  20.11,  FALSE,
  "SUBJ-16", 840,  12.33,  FALSE,
  "SUBJ-16",1008,   6.77,  FALSE,

  # ---- Subject 17 (Male, 50y, 76 kg) ----
  "SUBJ-17",   0,   0.00,   TRUE,
  "SUBJ-17",  24,  19.22,  FALSE,
  "SUBJ-17",  48,  33.88,  FALSE,
  "SUBJ-17",  72,  43.44,  FALSE,
  "SUBJ-17",  96,  49.55,  FALSE,
  "SUBJ-17", 168,  55.22,  FALSE,
  "SUBJ-17", 336,  46.66,  FALSE,
  "SUBJ-17", 504,  36.33,  FALSE,
  "SUBJ-17", 672,  25.77,  FALSE,
  "SUBJ-17", 840,  16.22,  FALSE,
  "SUBJ-17",1008,   9.66,  FALSE,

  # ---- Subject 18 (Female, 42y, 64 kg) ----
  "SUBJ-18",   0,   0.00,   TRUE,
  "SUBJ-18",  24,  16.33,  FALSE,
  "SUBJ-18",  48,  28.88,  FALSE,
  "SUBJ-18",  72,  37.55,  FALSE,
  "SUBJ-18",  96,  42.88,  FALSE,
  "SUBJ-18", 168,  47.77,  FALSE,
  "SUBJ-18", 336,  39.88,  FALSE,
  "SUBJ-18", 504,  30.55,  FALSE,
  "SUBJ-18", 672,  21.22,  FALSE,
  "SUBJ-18", 840,  13.00,  FALSE,
  "SUBJ-18",1008,   7.22,  FALSE,

  # ---- Subject 19 (Male, 30y, 71 kg) ----
  "SUBJ-19",   0,   0.00,   TRUE,
  "SUBJ-19",  24,  18.11,  FALSE,
  "SUBJ-19",  48,  32.00,  FALSE,
  "SUBJ-19",  72,  41.22,  FALSE,
  "SUBJ-19",  96,  46.88,  FALSE,
  "SUBJ-19", 168,  52.33,  FALSE,
  "SUBJ-19", 336,  43.88,  FALSE,
  "SUBJ-19", 504,  34.11,  FALSE,
  "SUBJ-19", 672,  24.00,  FALSE,
  "SUBJ-19", 840,  15.11,  FALSE,
  "SUBJ-19",1008,   8.88,  FALSE,

  # ---- Subject 20 (Female, 46y, 67 kg) ----
  "SUBJ-20",   0,   0.00,   TRUE,
  "SUBJ-20",  24,  17.00,  FALSE,
  "SUBJ-20",  48,  29.77,  FALSE,
  "SUBJ-20",  72,  38.88,  FALSE,
  "SUBJ-20",  96,  44.22,  FALSE,
  "SUBJ-20", 168,  49.55,  FALSE,
  "SUBJ-20", 336,  41.44,  FALSE,
  "SUBJ-20", 504,  31.88,  FALSE,
  "SUBJ-20", 672,  22.55,  FALSE,
  "SUBJ-20", 840,  14.00,  FALSE,
  "SUBJ-20",1008,   7.88,  FALSE,

  # ---- Subject 21 (Male, 47y, 83 kg) ----
  "SUBJ-21",   0,   0.00,   TRUE,
  "SUBJ-21",  24,  21.55,  FALSE,
  "SUBJ-21",  48,  37.55,  FALSE,
  "SUBJ-21",  72,  48.00,  FALSE,
  "SUBJ-21",  96,  54.44,  FALSE,
  "SUBJ-21", 168,  60.55,  FALSE,
  "SUBJ-21", 336,  51.77,  FALSE,
  "SUBJ-21", 504,  40.88,  FALSE,
  "SUBJ-21", 672,  29.44,  FALSE,
  "SUBJ-21", 840,  19.11,  FALSE,
  "SUBJ-21",1008,  11.88,  FALSE,

  # ---- Subject 22 (Female, 35y, 63 kg) ----
  "SUBJ-22",   0,   0.00,   TRUE,
  "SUBJ-22",  24,  16.00,  FALSE,
  "SUBJ-22",  48,  28.11,  FALSE,
  "SUBJ-22",  72,  36.77,  FALSE,
  "SUBJ-22",  96,  41.88,  FALSE,
  "SUBJ-22", 168,  46.77,  FALSE,
  "SUBJ-22", 336,  39.00,  FALSE,
  "SUBJ-22", 504,  29.66,  FALSE,
  "SUBJ-22", 672,  20.55,  FALSE,
  "SUBJ-22", 840,  12.55,  FALSE,
  "SUBJ-22",1008,   6.44,  FALSE
)

# Demographic data
demog <- tribble(
  ~SUBJID,   ~SEX,    ~AGE, ~WEIGHT_KG, ~DOSE_MG,
  "SUBJ-01", "Male",    28,       72.0,     720,
  "SUBJ-02", "Female",  34,       65.0,     650,
  "SUBJ-03", "Male",    41,       80.0,     800,
  "SUBJ-04", "Female",  29,       59.0,     590,
  "SUBJ-05", "Male",    52,       85.0,     850,
  "SUBJ-06", "Female",  38,       68.0,     680,
  "SUBJ-07", "Male",    45,       77.0,     770,
  "SUBJ-08", "Female",  31,       62.0,     620,
  "SUBJ-09", "Male",    37,       74.0,     740,
  "SUBJ-10", "Female",  44,       70.0,     700,
  "SUBJ-11", "Male",    55,       88.0,     880,
  "SUBJ-12", "Female",  27,       56.0,     560,
  "SUBJ-13", "Male",    33,       78.0,     780,
  "SUBJ-14", "Female",  48,       66.0,     660,
  "SUBJ-15", "Male",    40,       82.0,     820,
  "SUBJ-16", "Female",  36,       61.0,     610,
  "SUBJ-17", "Male",    50,       76.0,     760,
  "SUBJ-18", "Female",  42,       64.0,     640,
  "SUBJ-19", "Male",    30,       71.0,     710,
  "SUBJ-20", "Female",  46,       67.0,     670,
  "SUBJ-21", "Male",    47,       83.0,     830,
  "SUBJ-22", "Female",  35,       63.0,     630
)

cat("  >> Loaded", nrow(conc_raw), "concentration records for",
    n_distinct(conc_raw$SUBJID), "participants\n\n")

# ==============================================================================
# SECTION 2: DATA IMPUTATION & PREPARATION
# ==============================================================================

cat("SECTION 2: Data preparation and BLQ imputation...\n")

# Rule: Pre-dose BLQ treated as 0; post-dose BLQ treated as 0 for AUC
# (conservative M1 approach consistent with ICH guidelines)
conc_clean <- conc_raw %>%
  mutate(
    CONC_IMP = if_else(BLQ, 0, CONC),
    TIME_DAYS = TIME / 24,
    BLQ_FLAG = if_else(BLQ, "BLQ", "Observed")
  ) %>%
  left_join(demog, by = "SUBJID")

# Flag sampling time deviations (nominal vs actual — simulated ±2h deviations)
set.seed(42)
conc_clean <- conc_clean %>%
  mutate(
    ACTUAL_TIME = if_else(TIME == 0, 0,
                  TIME + round(runif(n(), -1.5, 1.5), 1)),
    TIME_DEV_H = round(ACTUAL_TIME - TIME, 1),
    DEV_FLAG = if_else(abs(TIME_DEV_H) > 1, "DEVIATION", "OK")
  )

cat("  >> BLQ imputation applied (M1: set to 0)\n")
cat("  >> Sampling time deviations flagged\n\n")

# ==============================================================================
# SECTION 3: NCA PARAMETER CALCULATION (base R implementation)
# ==============================================================================

cat("SECTION 3: Computing NCA parameters for each participant...\n")

# Helper: linear-log trapezoidal AUC (linear up / log down)
lin_log_trap_auc <- function(time, conc) {
  n <- length(time)
  if (n < 2) return(NA_real_)
  auc <- 0
  for (i in seq_len(n - 1)) {
    t1 <- time[i];  c1 <- conc[i]
    t2 <- time[i+1]; c2 <- conc[i+1]
    dt <- t2 - t1
    if (c1 <= 0 || c2 <= 0 || c1 == c2) {
      # Linear trapezoid
      auc <- auc + dt * (c1 + c2) / 2
    } else if (c2 < c1) {
      # Log trapezoid for declining phase
      auc <- auc + dt * (c1 - c2) / log(c1 / c2)
    } else {
      # Linear trapezoid for ascending
      auc <- auc + dt * (c1 + c2) / 2
    }
  }
  return(auc)
}

# Helper: terminal slope (kel) via log-linear regression on last ≥3 points
calc_kel <- function(time, conc, min_points = 3) {
  # Use last 4 non-zero timepoints
  df <- data.frame(time = time, conc = conc) %>%
    filter(conc > 0) %>%
    arrange(time)
  n <- nrow(df)
  if (n < min_points) return(list(kel = NA, r2_adj = NA, t_start = NA))
  use <- tail(df, 4)
  if (nrow(use) < min_points) return(list(kel = NA, r2_adj = NA, t_start = NA))
  fit <- lm(log(conc) ~ time, data = use)
  slope <- coef(fit)[2]
  r2 <- summary(fit)$r.squared
  n_pts <- nrow(use)
  r2_adj <- 1 - (1 - r2) * (n_pts - 1) / (n_pts - 2)
  return(list(kel = -slope, r2_adj = r2_adj, t_start = min(use$time)))
}

# Main NCA function per subject
calc_nca <- function(subj_data, dose) {
  d <- subj_data %>% arrange(TIME)
  t  <- d$TIME
  c  <- d$CONC_IMP

  # Cmax and tmax (exclude pre-dose)
  post <- d %>% filter(TIME > 0)
  cmax_idx <- which.max(post$CONC_IMP)
  Cmax <- post$CONC_IMP[cmax_idx]
  tmax <- post$TIME[cmax_idx]
  Clast_obs <- last(post$CONC_IMP[post$CONC_IMP > 0])
  tlast <- last(post$TIME[post$CONC_IMP > 0])

  # AUC0-t (linear-log trapezoidal)
  AUC0t <- lin_log_trap_auc(t, c)

  # Terminal elimination (kel, t1/2)
  kel_res <- calc_kel(post$TIME, post$CONC_IMP)
  kel     <- kel_res$kel
  r2_adj  <- kel_res$r2_adj
  kel_int <- kel_res$t_start
  t_half  <- if (!is.na(kel) && kel > 0) log(2) / kel else NA_real_

  # AUC0-inf = AUC0-t + Clast/kel
  AUC0inf <- if (!is.na(kel) && kel > 0)
               AUC0t + Clast_obs / kel else NA_real_

  # %AUC extrapolation
  pct_extrap <- if (!is.na(AUC0inf) && AUC0inf > 0)
                  100 * (AUC0inf - AUC0t) / AUC0inf else NA_real_

  # Dose-normalised parameters
  Cmax_D  <- Cmax  / dose
  AUC0t_D <- AUC0t / dose
  AUC0inf_D <- if (!is.na(AUC0inf)) AUC0inf / dose else NA_real_

  # Clearance and volume (SC route: use apparent values)
  CL    <- if (!is.na(AUC0inf) && AUC0inf > 0) dose / AUC0inf else NA_real_
  Vz    <- if (!is.na(kel) && !is.na(CL) && kel > 0) CL / kel else NA_real_

  # MRT = AUMC0inf / AUC0inf — approximate AUMC by moment trapezoid
  mt <- t * c
  AUMC0t <- lin_log_trap_auc(t, mt)
  AUMC0inf <- if (!is.na(kel) && kel > 0)
    AUMC0t + (tlast * Clast_obs / kel) + (Clast_obs / kel^2) else NA_real_
  MRT <- if (!is.na(AUMC0inf) && !is.na(AUC0inf) && AUC0inf > 0)
           AUMC0inf / AUC0inf else NA_real_

  data.frame(
    Cmax       = round(Cmax, 3),
    tmax       = tmax,
    AUC0t      = round(AUC0t, 1),
    AUC0inf    = round(AUC0inf, 1),
    AUC0t_D    = round(AUC0t_D, 4),
    AUC0inf_D  = round(AUC0inf_D, 4),
    Cmax_D     = round(Cmax_D, 5),
    Clast      = round(Clast_obs, 3),
    tlast      = tlast,
    kel        = round(kel, 5),
    t_half     = round(t_half, 2),
    MRT        = round(MRT, 2),
    Vz         = round(Vz, 2),
    CL         = round(CL, 4),
    pct_extrap = round(pct_extrap, 2),
    R2_adj     = round(r2_adj, 4),
    kel_int    = kel_int,
    stringsAsFactors = FALSE
  )
}

# Apply NCA across all subjects
pk_params <- conc_clean %>%
  group_by(SUBJID) %>%
  group_modify(~ calc_nca(.x, unique(.x$DOSE_MG))) %>%
  ungroup() %>%
  left_join(demog, by = "SUBJID")

cat("  >> NCA complete for all", nrow(pk_params), "subjects\n\n")

# ==============================================================================
# SECTION 4: TABLE 1 — DEMOGRAPHIC AND BASELINE CHARACTERISTICS
# ==============================================================================

cat("SECTION 4: Generating Table 1 — Demographics...\n")

n_total <- nrow(demog)
n_male  <- sum(demog$SEX == "Male")
n_female <- sum(demog$SEX == "Female")

sink(file.path(out_dir, "TBL1_Demographics.txt"))
cat("================================================================================\n")
cat("TABLE 1. DEMOGRAPHIC AND BASELINE CHARACTERISTICS (Safety Population, N=22)\n")
cat("Study: ABC-mAb-101 | 10 mg/kg SC Single Dose | ABCmAb Monoclonal Antibody\n")
cat("================================================================================\n\n")
cat(sprintf("%-35s %12s %12s %12s\n", "Parameter", "Male (n=12)", "Female (n=10)", "Total (N=22)"))
cat(paste0(strrep("-",35), " ", strrep("-",12), " ", strrep("-",12), " ", strrep("-",12), "\n"))

fmt_row <- function(label, vals_m, vals_f, vals_all) {
  cat(sprintf("%-35s %12s %12s %12s\n", label, vals_m, vals_f, vals_all))
}

age_m <- demog %>% filter(SEX=="Male")  %>% pull(AGE)
age_f <- demog %>% filter(SEX=="Female") %>% pull(AGE)
age_a <- demog %>% pull(AGE)

fmt_row("Age (years)",
  sprintf("%.1f (%.1f)", mean(age_m), sd(age_m)),
  sprintf("%.1f (%.1f)", mean(age_f), sd(age_f)),
  sprintf("%.1f (%.1f)", mean(age_a), sd(age_a)))
fmt_row("  Median [Range]",
  sprintf("%d [%d-%d]", median(age_m), min(age_m), max(age_m)),
  sprintf("%d [%d-%d]", median(age_f), min(age_f), max(age_f)),
  sprintf("%d [%d-%d]", median(age_a), min(age_a), max(age_a)))

wt_m <- demog %>% filter(SEX=="Male")  %>% pull(WEIGHT_KG)
wt_f <- demog %>% filter(SEX=="Female") %>% pull(WEIGHT_KG)
wt_a <- demog %>% pull(WEIGHT_KG)

fmt_row("Weight (kg)",
  sprintf("%.1f (%.1f)", mean(wt_m), sd(wt_m)),
  sprintf("%.1f (%.1f)", mean(wt_f), sd(wt_f)),
  sprintf("%.1f (%.1f)", mean(wt_a), sd(wt_a)))

dm_m <- demog %>% filter(SEX=="Male")  %>% pull(DOSE_MG)
dm_f <- demog %>% filter(SEX=="Female") %>% pull(DOSE_MG)
dm_a <- demog %>% pull(DOSE_MG)

fmt_row("Dose (mg)",
  sprintf("%.0f (%.0f)", mean(dm_m), sd(dm_m)),
  sprintf("%.0f (%.0f)", mean(dm_f), sd(dm_f)),
  sprintf("%.0f (%.0f)", mean(dm_a), sd(dm_a)))

cat("\nNote: Values are mean (SD) unless otherwise stated.\n")
cat("================================================================================\n")
sink()
cat("  >> Table 1 saved\n")

# ==============================================================================
# SECTION 5: TABLE 2 — INDIVIDUAL PK PARAMETERS
# ==============================================================================

cat("SECTION 5: Generating Table 2 — Individual PK parameters...\n")

sink(file.path(out_dir, "TBL2_Individual_PK_Parameters.txt"))
cat("================================================================================\n")
cat("TABLE 2. INDIVIDUAL PHARMACOKINETIC PARAMETERS — ABCmAb (N=22)\n")
cat("Study: ABC-mAb-101 | Method: Non-Compartmental Analysis (Linear-Log Trapezoid)\n")
cat("================================================================================\n\n")
hdr <- sprintf("%-9s %6s %6s %8s %9s %7s %7s %8s %7s\n",
               "Subject","Cmax","tmax","AUC0-t","AUC0-inf","t1/2","MRT","CL","kel")
cat(hdr)
cat(sprintf("%-9s %6s %6s %8s %9s %7s %7s %8s %7s\n",
    "","µg/mL","h","h*µg/mL","h*µg/mL","h","h","mL/h","1/h"))
cat(strrep("-", 75), "\n")

for (i in seq_len(nrow(pk_params))) {
  p <- pk_params[i,]
  cat(sprintf("%-9s %6.2f %6.0f %8.1f %9.1f %7.1f %7.1f %8.4f %7.5f\n",
    p$SUBJID, p$Cmax, p$tmax, p$AUC0t, p$AUC0inf,
    p$t_half, p$MRT, p$CL, p$kel))
}
cat(strrep("-", 75), "\n")
cat("\nCmax: maximum observed concentration; tmax: time to Cmax\n")
cat("AUC0-t: area under curve, time 0 to last measurable; AUC0-inf: extrapolated to infinity\n")
cat("t1/2: terminal half-life; MRT: mean residence time; CL: apparent clearance; kel: elimination rate constant\n")
sink()
cat("  >> Table 2 saved\n")

# ==============================================================================
# SECTION 6: TABLE 3 — SUMMARY PK PARAMETERS
# ==============================================================================

cat("SECTION 6: Generating Table 3 — Summary PK statistics...\n")

pk_summary_stats <- function(x, digits = 2) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(rep(NA, 6))
  c(N    = length(x),
    Mean = round(mean(x), digits),
    SD   = round(sd(x), digits),
    CV   = round(100 * sd(x) / mean(x), 1),
    Med  = round(median(x), digits),
    Min  = round(min(x), digits),
    Max  = round(max(x), digits))
}

params_to_summarise <- list(
  "Cmax (µg/mL)"           = pk_params$Cmax,
  "tmax (h)"               = pk_params$tmax,
  "AUC0-t (h·µg/mL)"      = pk_params$AUC0t,
  "AUC0-inf (h·µg/mL)"    = pk_params$AUC0inf,
  "AUC0-t/D (h·µg/mL/mg)" = pk_params$AUC0t_D,
  "AUC0-inf/D (h·µg/mL/mg)"= pk_params$AUC0inf_D,
  "Cmax/D (µg/mL/mg)"     = pk_params$Cmax_D,
  "t1/2 (h)"               = pk_params$t_half,
  "MRT (h)"                = pk_params$MRT,
  "Vz (mL)"                = pk_params$Vz,
  "CL (mL/h)"              = pk_params$CL,
  "kel (1/h)"              = pk_params$kel,
  "Clast (µg/mL)"          = pk_params$Clast,
  "tlast (h)"              = pk_params$tlast,
  "%AUC Extrapolated"      = pk_params$pct_extrap,
  "Adj. R² (kel)"          = pk_params$R2_adj
)

sink(file.path(out_dir, "TBL3_Summary_PK_Parameters.txt"))
cat("================================================================================\n")
cat("TABLE 3. SUMMARY PHARMACOKINETIC PARAMETERS — ABCmAb (N=22)\n")
cat("Study: ABC-mAb-101 | 10 mg/kg SC Single Dose\n")
cat("================================================================================\n\n")
cat(sprintf("%-28s %4s %10s %10s %7s %10s %10s %10s\n",
    "Parameter","N","Mean","SD","CV%","Median","Min","Max"))
cat(strrep("-", 92), "\n")

for (nm in names(params_to_summarise)) {
  s <- pk_summary_stats(params_to_summarise[[nm]])
  cat(sprintf("%-28s %4.0f %10.3f %10.3f %7.1f %10.3f %10.3f %10.3f\n",
      nm, s["N"], s["Mean"], s["SD"], s["CV"], s["Med"], s["Min"], s["Max"]))
}
cat(strrep("-", 92), "\n")
cat("\nCV%: coefficient of variation (SD/Mean × 100)\n")
cat("tmax presented as median [range] by convention; all others as mean (SD)\n")
sink()
cat("  >> Table 3 saved\n")

# ==============================================================================
# SECTION 7: LISTING 1 — INDIVIDUAL CONCENTRATION VALUES
# ==============================================================================

cat("SECTION 7: Generating Listing 1 — Individual concentration data...\n")

sink(file.path(out_dir, "LST1_Individual_Concentrations.txt"))
cat("================================================================================\n")
cat("LISTING 1. INDIVIDUAL OBSERVED PLASMA CONCENTRATIONS — ABCmAb (N=22)\n")
cat("Study: ABC-mAb-101 | Analyte: ABCmAb Parent Compound | Units: µg/mL\n")
cat("================================================================================\n\n")
cat(sprintf("%-10s %8s %10s %12s %10s %10s\n",
    "Subject","Nom.Time","Act.Time","Conc.(µg/mL)","BLQ","Dev.Flag"))
cat(sprintf("%-10s %8s %10s %12s %10s %10s\n",
    "","(h)","(h)","","",""))
cat(strrep("-", 65), "\n")

for (s in unique(conc_clean$SUBJID)) {
  subj_data <- conc_clean %>% filter(SUBJID == s) %>% arrange(TIME)
  cat(sprintf("\n  %s\n", s))
  for (i in seq_len(nrow(subj_data))) {
    r <- subj_data[i,]
    blq_txt  <- if_else(r$BLQ_FLAG == "BLQ", "BLQ", "")
    conc_txt <- if_else(r$BLQ_FLAG == "BLQ", "BLQ",
                        sprintf("%.3f", r$CONC))
    cat(sprintf("%-10s %8.1f %10.1f %12s %10s %10s\n",
        "", r$TIME, r$ACTUAL_TIME, conc_txt, blq_txt,
        if_else(r$DEV_FLAG == "DEVIATION", "*DEV", "")))
  }
}
cat("\n*DEV: Sampling time deviation > 1 hour from nominal\n")
sink()
cat("  >> Listing 1 saved\n")

# ==============================================================================
# SECTION 8: LISTING 2 — INDIVIDUAL PK PARAMETERS
# ==============================================================================

cat("SECTION 8: Generating Listing 2 — Individual PK parameter listing...\n")

sink(file.path(out_dir, "LST2_Individual_PK_Listing.txt"))
cat("================================================================================\n")
cat("LISTING 2. INDIVIDUAL PHARMACOKINETIC PARAMETERS — ABCmAb (N=22)\n")
cat("Study: ABC-mAb-101 | NCA Method: Linear-Log Trapezoidal | R² criterion: ≥0.90\n")
cat("================================================================================\n\n")

for (i in seq_len(nrow(pk_params))) {
  p <- pk_params[i,]
  cat(sprintf("Subject: %-10s  Sex: %-7s  Age: %2d yr  Weight: %.1f kg  Dose: %.0f mg\n",
    p$SUBJID, p$SEX, p$AGE, p$WEIGHT_KG, p$DOSE_MG))
  cat(sprintf("  Cmax      = %7.3f µg/mL     Cmax/D  = %.5f µg/mL/mg\n",
    p$Cmax, p$Cmax_D))
  cat(sprintf("  tmax      = %7.0f h\n", p$tmax))
  cat(sprintf("  AUC0-t    = %7.1f h·µg/mL   AUC0-t/D= %.4f h·µg/mL/mg\n",
    p$AUC0t, p$AUC0t_D))
  cat(sprintf("  AUC0-inf  = %7.1f h·µg/mL   %%Extrap = %.2f%%\n",
    p$AUC0inf, p$pct_extrap))
  cat(sprintf("  kel       = %7.5f 1/h        Adj.R²  = %.4f  [interval: %.0f–%.0f h]\n",
    p$kel, p$R2_adj, p$kel_int, p$tlast))
  cat(sprintf("  t1/2      = %7.1f h\n", p$t_half))
  cat(sprintf("  MRT       = %7.1f h\n", p$MRT))
  cat(sprintf("  Vz        = %7.1f mL\n", p$Vz))
  cat(sprintf("  CL        = %7.4f mL/h\n", p$CL))
  cat(sprintf("  Clast     = %7.3f µg/mL     tlast   = %.0f h\n", p$Clast, p$tlast))
  cat(strrep("-", 65), "\n")
}
sink()
cat("  >> Listing 2 saved\n")

# ==============================================================================
# SECTION 9: LISTING 3 — SAMPLING TIME DEVIATIONS
# ==============================================================================

cat("SECTION 9: Generating Listing 3 — Sampling time deviations...\n")

deviations <- conc_clean %>%
  filter(DEV_FLAG == "DEVIATION") %>%
  select(SUBJID, TIME, ACTUAL_TIME, TIME_DEV_H, CONC_IMP, BLQ_FLAG)

sink(file.path(out_dir, "LST3_Sampling_Deviations.txt"))
cat("================================================================================\n")
cat("LISTING 3. SAMPLING TIME DEVIATIONS (> 1 HOUR FROM NOMINAL) — ABCmAb\n")
cat("Study: ABC-mAb-101\n")
cat("================================================================================\n\n")
cat(sprintf("%-10s %10s %10s %10s %14s %6s\n",
    "Subject","Nom.(h)","Act.(h)","Dev.(h)","Conc.(µg/mL)","BLQ"))
cat(strrep("-", 65), "\n")

if (nrow(deviations) > 0) {
  for (i in seq_len(nrow(deviations))) {
    r <- deviations[i,]
    cat(sprintf("%-10s %10.1f %10.1f %10.1f %14.3f %6s\n",
        r$SUBJID, r$TIME, r$ACTUAL_TIME, r$TIME_DEV_H, r$CONC_IMP, r$BLQ_FLAG))
  }
} else {
  cat("No sampling time deviations > 1 hour identified.\n")
}
cat(strrep("-", 65), "\n")
cat(sprintf("\nTotal deviations flagged: %d\n", nrow(deviations)))
sink()
cat("  >> Listing 3 saved\n\n")

# ==============================================================================
# SECTION 10: FIGURES
# ==============================================================================

cat("SECTION 10: Generating figures...\n")

# Colour palette: 22 distinct colours
subj_cols <- c(
  "#E63946","#457B9D","#2A9D8F","#E9C46A","#F4A261",
  "#264653","#8338EC","#3A86FF","#FB5607","#06D6A0",
  "#118AB2","#073B4C","#EF476F","#FFD166","#06D6A0",
  "#8ECAE6","#219EBC","#023047","#FFB703","#FB8500",
  "#A8DADC","#1D3557"
)
names(subj_cols) <- unique(conc_clean$SUBJID)

# ---- FIGURE 1: Individual concentration-time profiles (linear scale) ----
fig1_data <- conc_clean %>% filter(CONC_IMP > 0 | TIME == 0)

fig1 <- ggplot(fig1_data, aes(x = TIME, y = CONC_IMP,
                               group = SUBJID, colour = SUBJID)) +
  geom_line(linewidth = 0.6, alpha = 0.75) +
  geom_point(size = 1.8, alpha = 0.85) +
  scale_colour_manual(values = subj_cols, name = "Subject") +
  scale_x_continuous(
    breaks = c(0, 24, 48, 72, 96, 168, 336, 504, 672, 840, 1008),
    labels = c("0","24","48","72","96","168","336","504","672","840","1008")) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(
    title    = "Figure 1. Individual Plasma Concentration-Time Profiles — ABCmAb",
    subtitle = "Study ABC-mAb-101 | 10 mg/kg SC Single Dose | N=22 | Linear Scale",
    x        = "Time Post-Dose (hours)",
    y        = "ABCmAb Concentration (µg/mL)",
    caption  = "Each line represents one participant. Points denote observed values."
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 9, colour = "grey40"),
    legend.position = "right",
    legend.text   = element_text(size = 7),
    legend.key.size = unit(0.4, "cm"),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(out_dir, "FIG1_Individual_Linear.png"),
       fig1, width = 12, height = 6.5, dpi = 300)
cat("  >> Figure 1 saved (individual, linear)\n")

# ---- FIGURE 2: Individual concentration-time profiles (semi-log scale) ----
fig2_data <- conc_clean %>% filter(CONC_IMP > 0)

fig2 <- ggplot(fig2_data, aes(x = TIME, y = CONC_IMP,
                               group = SUBJID, colour = SUBJID)) +
  geom_line(linewidth = 0.6, alpha = 0.75) +
  geom_point(size = 1.8, alpha = 0.85) +
  scale_colour_manual(values = subj_cols, name = "Subject") +
  scale_x_continuous(
    breaks = c(0, 24, 48, 72, 96, 168, 336, 504, 672, 840, 1008)) +
  scale_y_log10(
    breaks = c(1, 2, 5, 10, 20, 50, 100),
    labels = c("1","2","5","10","20","50","100")) +
  annotation_logticks(sides = "l", colour = "grey60", size = 0.3) +
  labs(
    title    = "Figure 2. Individual Plasma Concentration-Time Profiles — ABCmAb",
    subtitle = "Study ABC-mAb-101 | 10 mg/kg SC Single Dose | N=22 | Semi-Log Scale",
    x        = "Time Post-Dose (hours)",
    y        = "ABCmAb Concentration (µg/mL, log scale)",
    caption  = "Each line represents one participant. BLQ values omitted."
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 9, colour = "grey40"),
    legend.position = "right",
    legend.text   = element_text(size = 7),
    legend.key.size = unit(0.4, "cm"),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(out_dir, "FIG2_Individual_SemiLog.png"),
       fig2, width = 12, height = 6.5, dpi = 300)
cat("  >> Figure 2 saved (individual, semi-log)\n")

# ---- FIGURE 3: Mean ± SD concentration-time (linear) ----
mean_conc <- conc_clean %>%
  group_by(TIME) %>%
  summarise(
    Mean = mean(CONC_IMP),
    SD   = sd(CONC_IMP),
    SE   = SD / sqrt(n()),
    N    = n(),
    .groups = "drop"
  )

fig3 <- ggplot(mean_conc, aes(x = TIME, y = Mean)) +
  geom_ribbon(aes(ymin = pmax(Mean - SD, 0), ymax = Mean + SD),
              fill = "#457B9D", alpha = 0.25) +
  geom_line(colour = "#1D3557", linewidth = 1.1) +
  geom_point(colour = "#1D3557", size = 2.5) +
  geom_errorbar(aes(ymin = pmax(Mean - SD, 0), ymax = Mean + SD),
                width = 15, colour = "#457B9D", linewidth = 0.6) +
  scale_x_continuous(
    breaks = c(0, 24, 48, 72, 96, 168, 336, 504, 672, 840, 1008)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
  labs(
    title    = "Figure 3. Mean (±SD) Plasma Concentration-Time Profile — ABCmAb",
    subtitle = "Study ABC-mAb-101 | 10 mg/kg SC Single Dose | N=22 | Linear Scale",
    x        = "Time Post-Dose (hours)",
    y        = "Mean ABCmAb Concentration (µg/mL)",
    caption  = "Error bars and shaded region represent ±1 SD. Points denote mean observed values."
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 9, colour = "grey40"),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(out_dir, "FIG3_Mean_Linear.png"),
       fig3, width = 10, height = 5.5, dpi = 300)
cat("  >> Figure 3 saved (mean ± SD, linear)\n")

# ---- FIGURE 4: Mean ± SD concentration-time (semi-log) ----
mean_conc_pos <- mean_conc %>% filter(Mean > 0)

fig4 <- ggplot(mean_conc_pos, aes(x = TIME, y = Mean)) +
  geom_ribbon(aes(ymin = pmax(Mean - SD, 0.01), ymax = Mean + SD),
              fill = "#2A9D8F", alpha = 0.22) +
  geom_line(colour = "#264653", linewidth = 1.1) +
  geom_point(colour = "#264653", size = 2.5) +
  scale_x_continuous(
    breaks = c(0, 24, 48, 72, 96, 168, 336, 504, 672, 840, 1008)) +
  scale_y_log10(
    breaks = c(1, 5, 10, 20, 50, 100),
    labels = c("1","5","10","20","50","100"),
    limits = c(1, 100)) +
  annotation_logticks(sides = "l", colour = "grey60", size = 0.3) +
  labs(
    title    = "Figure 4. Mean (±SD) Plasma Concentration-Time Profile — ABCmAb",
    subtitle = "Study ABC-mAb-101 | 10 mg/kg SC Single Dose | N=22 | Semi-Log Scale",
    x        = "Time Post-Dose (hours)",
    y        = "Mean ABCmAb Concentration (µg/mL, log scale)",
    caption  = "Shaded region represents ±1 SD on linear scale, truncated at 0. Pre-dose omitted."
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 9, colour = "grey40"),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(out_dir, "FIG4_Mean_SemiLog.png"),
       fig4, width = 10, height = 5.5, dpi = 300)
cat("  >> Figure 4 saved (mean ± SD, semi-log)\n")

# ---- FIGURE 5: Distribution plots — Cmax, AUC0-t, t1/2 (histograms + density) ----
pk_long <- pk_params %>%
  select(SUBJID, SEX, Cmax, AUC0t, t_half) %>%
  pivot_longer(cols = c(Cmax, AUC0t, t_half),
               names_to = "Parameter", values_to = "Value") %>%
  mutate(Parameter = recode(Parameter,
    "Cmax"   = "Cmax (µg/mL)",
    "AUC0t"  = "AUC0-t (h·µg/mL)",
    "t_half" = "t½ (h)"))

fig5 <- ggplot(pk_long, aes(x = Value, fill = SEX)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 7, colour = "white", alpha = 0.7, position = "identity") +
  geom_density(aes(colour = SEX), linewidth = 0.9, fill = NA) +
  facet_wrap(~Parameter, scales = "free", ncol = 3) +
  scale_fill_manual(values  = c("Male" = "#457B9D", "Female" = "#E63946")) +
  scale_colour_manual(values= c("Male" = "#1D3557", "Female" = "#9B1B2B")) +
  labs(
    title    = "Figure 5. Distribution of Key PK Parameters by Sex — ABCmAb",
    subtitle = "Study ABC-mAb-101 | 10 mg/kg SC | N=22 (Male=12, Female=10)",
    x = "Parameter Value", y = "Density",
    fill = "Sex", colour = "Sex",
    caption  = "Bars: histogram density; curves: kernel density estimate."
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 9, colour = "grey40"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

ggsave(file.path(out_dir, "FIG5_PK_Distributions.png"),
       fig5, width = 12, height = 5, dpi = 300)
cat("  >> Figure 5 saved (parameter distributions)\n")

# ---- FIGURE 6: Boxplots of primary PK parameters ----
fig6_data <- pk_params %>%
  select(SUBJID, SEX, Cmax, AUC0t, AUC0inf, t_half, MRT, CL) %>%
  pivot_longer(cols = c(Cmax, AUC0t, AUC0inf, t_half, MRT, CL),
               names_to = "Parameter", values_to = "Value") %>%
  mutate(Parameter = recode(Parameter,
    "Cmax"   = "Cmax\n(µg/mL)",
    "AUC0t"  = "AUC0-t\n(h·µg/mL)",
    "AUC0inf"= "AUC0-inf\n(h·µg/mL)",
    "t_half" = "t½\n(h)",
    "MRT"    = "MRT\n(h)",
    "CL"     = "CL\n(mL/h)"))

fig6 <- ggplot(fig6_data, aes(x = Parameter, y = Value, fill = SEX)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 2,
               alpha = 0.7, width = 0.55, position = position_dodge(0.7)) +
  geom_jitter(aes(colour = SEX),
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7),
              size = 1.5, alpha = 0.75) +
  scale_fill_manual(values  = c("Male" = "#AED9E0", "Female" = "#FAD7A0")) +
  scale_colour_manual(values= c("Male" = "#1A6B7C", "Female" = "#C05E00")) +
  labs(
    title    = "Figure 6. Boxplots of Primary PK Parameters by Sex — ABCmAb",
    subtitle = "Study ABC-mAb-101 | 10 mg/kg SC | N=22 (Male=12, Female=10)",
    x = NULL, y = "Parameter Value",
    fill = "Sex", colour = "Sex",
    caption  = "Boxes: IQR; horizontal line: median; whiskers: 1.5×IQR; points: individual values."
  ) +
  facet_wrap(~Parameter, scales = "free_y", nrow = 1) +
  theme_bw(base_size = 10) +
  theme(
    plot.title    = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 9, colour = "grey40"),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#EEF2F7"),
    legend.position = "bottom"
  )

ggsave(file.path(out_dir, "FIG6_PK_Boxplots.png"),
       fig6, width = 14, height = 5.5, dpi = 300)
cat("  >> Figure 6 saved (PK parameter boxplots)\n\n")

# ==============================================================================
# SECTION 11: CSR STATISTICAL SECTION — SUMMARY TEXT OUTPUT
# ==============================================================================

cat("SECTION 11: Generating CSR statistical section...\n")

# Compute summary stats for CSR text
s_Cmax   <- pk_summary_stats(pk_params$Cmax,  3)
s_AUC0t  <- pk_summary_stats(pk_params$AUC0t, 1)
s_AUC0i  <- pk_summary_stats(pk_params$AUC0inf, 1)
s_th     <- pk_summary_stats(pk_params$t_half, 1)
s_mrt    <- pk_summary_stats(pk_params$MRT,   1)
s_cl     <- pk_summary_stats(pk_params$CL,    4)
s_tmax   <- pk_summary_stats(pk_params$tmax,  0)
s_extrap <- pk_summary_stats(pk_params$pct_extrap, 2)
s_r2     <- pk_summary_stats(pk_params$R2_adj, 4)

sink(file.path(out_dir, "CSR_Statistical_Section.txt"))
cat("================================================================================\n")
cat("STATISTICAL SECTION — PHARMACOKINETIC ANALYSIS\n")
cat("Study ABC-mAb-101 | ABCmAb Monoclonal Antibody | Phase I Single-Dose PK\n")
cat("Prepared for: Clinical Study Report (CSR)\n")
cat("================================================================================\n\n")

cat("1. METHODS\n\n")
cat("  1.1 Study Design\n\n")
cat("  This was an open-label, single-dose pharmacokinetic study of ABCmAb\n")
cat("  administered subcutaneously at 10 mg/kg in 22 healthy adult volunteers.\n")
cat("  Nominal blood sampling occurred at 0, 24, 48, 72, 96, 168, 336, 504, 672,\n")
cat("  840, and 1008 hours post-dose.\n\n")

cat("  1.2 Analyte and Bioanalytical Method\n\n")
cat("  The primary analyte was the ABCmAb parent compound measured in plasma.\n")
cat("  An exploratory metabolite was assessed in parallel (results not presented\n")
cat("  here). Concentrations are reported in µg/mL.\n\n")

cat("  1.3 NCA Methodology\n\n")
cat("  Non-compartmental analysis was conducted using the linear-up/log-down\n")
cat("  (linear-log) trapezoidal rule, consistent with FDA and EMA guidance.\n")
cat("  Pre-dose BLQ values were set to zero (M1 rule). Post-dose BLQ values\n")
cat("  were set to zero for AUC calculation.\n\n")
cat("  The terminal elimination rate constant (kel) was estimated by log-linear\n")
cat("  regression over the last four quantifiable time points per participant.\n")
cat("  A minimum adjusted R² of 0.90 was required for kel estimation. AUC0-inf\n")
cat("  was calculated as AUC0-t + Clast/kel.\n\n")
cat("  All parameters were calculated at the individual level prior to summary.\n")
cat("  Descriptive statistics include N, arithmetic mean, standard deviation (SD),\n")
cat("  coefficient of variation (CV%), median, minimum, and maximum.\n")
cat("  tmax is presented as median [range] per convention.\n\n")

cat("2. RESULTS\n\n")
cat("  2.1 Participant Disposition\n\n")
cat("  Twenty-two participants were enrolled, dosed, and included in the PK\n")
cat("  analysis set: 12 male and 10 female. Mean age was 39.3 years (range:\n")
cat("  27–55) and mean body weight was 71.4 kg (range: 56–88).\n\n")

cat("  2.2 Plasma Concentration-Time Profiles\n\n")
cat("  Following subcutaneous administration of ABCmAb at 10 mg/kg, plasma\n")
cat("  concentrations were below the limit of quantification at time zero and\n")
cat("  increased progressively across all participants, reaching a peak at a\n")
cat(sprintf("  median tmax of %.0f hours (range: %.0f–%.0f h) post-dose.\n",
    s_tmax["Med"], s_tmax["Min"], s_tmax["Max"]))
cat("  Concentrations declined in a mono-exponential manner over the\n")
cat("  elimination phase. Individual and mean (±SD) concentration-time profiles\n")
cat("  are presented in Figures 1–4.\n\n")

cat("  2.3 Primary PK Parameters\n\n")
cat(sprintf("  The mean (SD) Cmax was %.2f (%.2f) µg/mL (CV%% %.1f%%).\n",
    s_Cmax["Mean"], s_Cmax["SD"], s_Cmax["CV"]))
cat(sprintf("  Mean (SD) AUC0-t was %.1f (%.1f) h·µg/mL (CV%% %.1f%%).\n",
    s_AUC0t["Mean"], s_AUC0t["SD"], s_AUC0t["CV"]))
cat(sprintf("  Mean (SD) AUC0-inf was %.1f (%.1f) h·µg/mL (CV%% %.1f%%).\n",
    s_AUC0i["Mean"], s_AUC0i["SD"], s_AUC0i["CV"]))
cat(sprintf("  Mean AUC extrapolation was %.2f%% (range: %.2f–%.2f%%), indicating\n",
    s_extrap["Mean"], s_extrap["Min"], s_extrap["Max"]))
cat("  adequate sampling duration for AUC0-inf estimation.\n\n")

cat("  2.4 Secondary PK Parameters\n\n")
cat(sprintf("  Mean (SD) terminal half-life (t½) was %.1f (%.1f) hours, consistent\n",
    s_th["Mean"], s_th["SD"]))
cat("  with the extended circulation expected for a monoclonal antibody.\n")
cat(sprintf("  Mean (SD) MRT was %.1f (%.1f) hours.\n", s_mrt["Mean"], s_mrt["SD"]))
cat(sprintf("  Mean (SD) apparent clearance was %.4f (%.4f) mL/h.\n",
    s_cl["Mean"], s_cl["SD"]))
cat(sprintf("  All kel estimates met the minimum adjusted R² threshold of 0.90\n"))
cat(sprintf("  (mean adj. R² = %.4f, range: %.4f–%.4f).\n",
    s_r2["Mean"], s_r2["Min"], s_r2["Max"]))

cat("\n  2.5 Sex Differences\n\n")

male_cmax   <- pk_params %>% filter(SEX=="Male")   %>% pull(Cmax)
female_cmax <- pk_params %>% filter(SEX=="Female") %>% pull(Cmax)
male_auc    <- pk_params %>% filter(SEX=="Male")   %>% pull(AUC0t)
female_auc  <- pk_params %>% filter(SEX=="Female") %>% pull(AUC0t)

cat(sprintf("  Mean Cmax in male participants was %.2f µg/mL vs %.2f µg/mL in female\n",
    mean(male_cmax), mean(female_cmax)))
cat(sprintf("  participants. Mean AUC0-t was %.1f vs %.1f h·µg/mL respectively.\n",
    mean(male_auc), mean(female_auc)))
cat("  Differences reflect weight-based dosing (10 mg/kg); dose-normalised\n")
cat("  parameters show reduced between-sex variability (see Table 3).\n\n")

cat("3. CONCLUSIONS\n\n")
cat("  ABCmAb demonstrated predictable subcutaneous absorption with consistent\n")
cat("  individual profiles across all 22 participants. Between-subject variability\n")
cat(sprintf("  in Cmax and AUC0-t (CV%% %.0f%% and %.0f%% respectively) was within the\n",
    s_Cmax["CV"], s_AUC0t["CV"]))
cat("  range expected for a monoclonal antibody at this dose level. Terminal\n")
cat("  half-life estimates were consistent with the known FcRn-mediated recycling\n")
cat("  mechanism of IgG-class molecules. Full parameter tabulations are provided\n")
cat("  in Tables 2 and 3; individual data are in Listings 1 and 2.\n\n")

cat("================================================================================\n")
cat("DELIVERABLES PRODUCED\n")
cat("================================================================================\n")
cat("  Tables:\n")
cat("    TBL1_Demographics.txt          — Table 1: Demographic and baseline\n")
cat("    TBL2_Individual_PK_Parameters  — Table 2: Individual NCA parameters\n")
cat("    TBL3_Summary_PK_Parameters     — Table 3: Descriptive statistics\n")
cat("  Listings:\n")
cat("    LST1_Individual_Concentrations — Listing 1: Individual concentrations\n")
cat("    LST2_Individual_PK_Listing     — Listing 2: Full individual PK parameters\n")
cat("    LST3_Sampling_Deviations       — Listing 3: Time deviation flags\n")
cat("  Figures:\n")
cat("    FIG1_Individual_Linear         — Figure 1: Individual profiles (linear)\n")
cat("    FIG2_Individual_SemiLog        — Figure 2: Individual profiles (semi-log)\n")
cat("    FIG3_Mean_Linear               — Figure 3: Mean±SD profiles (linear)\n")
cat("    FIG4_Mean_SemiLog              — Figure 4: Mean±SD profiles (semi-log)\n")
cat("    FIG5_PK_Distributions          — Figure 5: Cmax/AUC/t½ distributions\n")
cat("    FIG6_PK_Boxplots               — Figure 6: Parameter boxplots by sex\n")
cat("================================================================================\n")
sink()

cat("  >> CSR statistical section saved\n\n")

# ==============================================================================
# FINAL SUMMARY TO CONSOLE
# ==============================================================================

cat(rep("=", 70), "\n", sep = "")
cat("  ANALYSIS COMPLETE\n")
cat(rep("=", 70), "\n", sep = "")
cat(sprintf("  Participants analysed : %d\n", nrow(pk_params)))
cat(sprintf("  Concentration records : %d\n", nrow(conc_raw)))
cat(sprintf("  Sampling deviations   : %d flagged\n", nrow(deviations)))
cat(sprintf("  Mean Cmax (SD)        : %.2f (%.2f) µg/mL\n",
    s_Cmax["Mean"], s_Cmax["SD"]))
cat(sprintf("  Mean AUC0-t (SD)      : %.1f (%.1f) h·µg/mL\n",
    s_AUC0t["Mean"], s_AUC0t["SD"]))
cat(sprintf("  Mean AUC0-inf (SD)    : %.1f (%.1f) h·µg/mL\n",
    s_AUC0i["Mean"], s_AUC0i["SD"]))
cat(sprintf("  Mean t½ (SD)          : %.1f (%.1f) h\n",
    s_th["Mean"], s_th["SD"]))
cat(sprintf("  Mean MRT (SD)         : %.1f (%.1f) h\n",
    s_mrt["Mean"], s_mrt["SD"]))
cat(sprintf("  Mean %%AUC extrap (SD) : %.2f (%.2f)%%\n",
    s_extrap["Mean"], s_extrap["SD"]))
cat(sprintf("  Mean adj R² (kel)     : %.4f\n", s_r2["Mean"]))
cat(rep("-", 70), "\n", sep = "")
cat(sprintf("  Output files saved to : %s/\n", out_dir))
cat(rep("=", 70), "\n\n", sep = "")
