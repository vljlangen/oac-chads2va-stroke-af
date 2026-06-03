# =============================================================================
# OAC and CHA2DS2-VA score — analysis code (Cox models)
# =============================================================================
#
# Manuscript: Does the association between oral anticoagulation and stroke risk
# differ by CHA2DS2-VA score in patients with atrial fibrillation? A nationwide
# cohort study.
#
# Workflow:
#   1. Load baseline cohort and prescription-period (OAC) data
#   2. Define follow-up and ischemic stroke outcome
#   3. Build counting-process data with time-dependent OAC exposure
#   4. Fit Cox models (any OAC; separate VKA and DOAC) with CHA2DS2-VA interaction
#   5. Export stratum-specific hazard ratios, rate tables, NNT, and figures
#
# Required files in the working directory (rename your local exports to match):
#   cohort_baseline.sav          — one row per patient; baseline covariates and dates
#   oac_prescription_periods.rds — OAC purchase periods (ATC, start/end per period)
#
# These generic names avoid publishing internal file naming in public repositories.
# Map your origin export column names in SECTION 2 if they differ from the defaults
# listed there.
#
# R packages: data.table, survival, haven, lubridate, ggplot2, gt;
#             showtext, sysfonts (optional; for figure fonts — needs internet once)
#
# Run: source("OAC_CHA2DS2VA.R")
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Packages and paths
# -----------------------------------------------------------------------------

library(data.table)
library(survival)
library(haven)
library(lubridate)
library(ggplot2)
library(gt)

# Origin data filenames (place renamed copies in the working directory)
COHORT_BASELINE_FILE     <- "cohort_baseline.sav"
OAC_PRESCRIPTION_FILE    <- "oac_prescription_periods.rds"

# Study calendar limits
STUDY_END_DATE   <- as.Date("2018-12-31")
STUDY_START_DATE <- as.Date("2007-01-01")

# -----------------------------------------------------------------------------
# 0B. Figure style (publication)
# -----------------------------------------------------------------------------

plot_font   <- "sans"
plot_width  <- 7
plot_height <- 5

if (requireNamespace("showtext", quietly = TRUE) &&
    requireNamespace("sysfonts", quietly = TRUE)) {
  tryCatch({
    library(showtext)
    library(sysfonts)
    font_add_google("Rosario", "rosario")
    showtext_auto()
    plot_font <- "rosario"
  }, error = function(e) {
    warning(
      "Could not load Google font 'Rosario'; using default sans. ",
      conditionMessage(e),
      call. = FALSE
    )
  })
}

colour_any_oac <- "#003f5c"
colour_vka     <- "#003f5c"
colour_doac    <- "#bc5090"
colour_ref     <- "grey35"
colour_grid    <- "grey88"
colour_text    <- "grey20"
colour_annot   <- "grey25"

x_lab_chadsva <- expression("CHA"["2"] * "DS"["2"] * "-VA score")
y_lab_hr      <- "Adjusted hazard ratio"
ref_label     <- "No OAC as reference"

common_y_limits <- c(0.2, 1.2)
common_y_breaks <- seq(0.2, 1.2, by = 0.2)

position_drug <- position_dodge(width = 0.25)

standard_pdf_save <- function(filename, plot) {
  ggsave(
    filename = filename,
    plot = plot,
    width = plot_width,
    height = plot_height,
    device = cairo_pdf
  )
}

theme_publication <- function() {
  theme_minimal(base_family = plot_font, base_size = 14) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      axis.title = element_text(face = "bold", size = 13),
      axis.text = element_text(size = 12, colour = colour_text),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(colour = colour_grid, linewidth = 0.4),
      legend.title = element_blank(),
      legend.text = element_text(size = 12, face = "bold"),
      legend.background = element_rect(fill = "white", colour = "grey85", linewidth = 0.3),
      plot.margin = margin(12, 18, 12, 12)
    )
}

add_reference_line <- function(max_x) {
  list(
    geom_hline(
      yintercept = 1,
      linetype = "dashed",
      linewidth = 0.6,
      colour = colour_ref
    ),
    annotate(
      "text",
      x = max_x,
      y = 1.03,
      label = ref_label,
      hjust = 1,
      vjust = 0,
      family = plot_font,
      size = 4.2,
      colour = colour_annot
    )
  )
}

format_p <- function(p) {
  fifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3))
}

# -----------------------------------------------------------------------------
# 1. Load origin data
# -----------------------------------------------------------------------------

if (!file.exists(COHORT_BASELINE_FILE)) {
  stop(
    "Baseline cohort file not found: ", COHORT_BASELINE_FILE,
    ". Save your origin export under this name in the working directory.",
    call. = FALSE
  )
}
if (!file.exists(OAC_PRESCRIPTION_FILE)) {
  stop(
    "OAC period file not found: ", OAC_PRESCRIPTION_FILE,
    ". Save your origin prescription-period export under this name.",
    call. = FALSE
  )
}

cohort_dt <- as.data.table(read_sav(COHORT_BASELINE_FILE))
oac_periods_dt <- as.data.table(readRDS(OAC_PRESCRIPTION_FILE))

# -----------------------------------------------------------------------------
# 2. Standardise column names (edit old= if your export uses different labels)
# -----------------------------------------------------------------------------

# Baseline cohort — expected origin columns -> analysis names
cohort_name_map <- c(
  SID = "patient_id",
  CohortEntryDate = "cohort_entry_date",
  kuolpvmSPSSdate = "death_date",
  ISfirstdateAfterCohortfirsteverandrecurrent = "first_stroke_date",
  SukupuoliBin = "sex_binary",
  CHADSVAbaseline = "chads2va_baseline",
  BleedingsBOAC = "prior_bleeding_flag",
  AbnormalRenalFunctionBOAC = "prior_renal_flag",
  AbnormalLiverFunctionBOAC = "prior_liver_flag",
  AlcoholBOAC = "prior_alcohol_flag",
  CancerBeforeOrAtCohort = "prior_cancer_flag",
  DementiaBOAC = "prior_dementia_flag",
  PsychiatricDiseaseBOAC = "prior_psychiatric_flag",
  Incometertiles = "income_tertile_code"
)

present <- intersect(names(cohort_name_map), names(cohort_dt))
if (length(present) > 0L) {
  setnames(
    cohort_dt,
    old = present,
    new = unname(cohort_name_map[present])
  )
}

# OAC prescription periods — expected origin columns -> analysis names
oac_name_map <- c(
  SID = "patient_id",
  dup_start = "period_start_date",
  dup_end = "period_end_date"
)

present_oac <- intersect(names(oac_name_map), names(oac_periods_dt))
if (length(present_oac) > 0L) {
  setnames(
    oac_periods_dt,
    old = present_oac,
    new = unname(oac_name_map[present_oac])
  )
}

# -----------------------------------------------------------------------------
# 3. Dates and patient identifiers
# -----------------------------------------------------------------------------

cohort_dt[, cohort_entry_date := as.Date(cohort_entry_date)]
cohort_dt[, first_stroke_date := as.Date(first_stroke_date)]
cohort_dt[, death_date := as.Date(death_date)]

oac_periods_dt[, period_start_date := as.Date(period_start_date)]
oac_periods_dt[, period_end_date := as.Date(period_end_date)]

cohort_dt[, patient_id := as.character(patient_id)]
oac_periods_dt[, patient_id := as.character(patient_id)]

# -----------------------------------------------------------------------------
# 4. Follow-up and stroke outcome
# -----------------------------------------------------------------------------
# Follow-up starts at the later of cohort entry and the study start date.
# Follow-up ends at the earliest of death, first stroke, or fixed study end.
# Stroke event = 1 when follow-up ends on the first stroke date.

cohort_dt[, cohort_start_date := pmax(cohort_entry_date, STUDY_START_DATE)]
cohort_dt <- cohort_dt[cohort_start_date < STUDY_END_DATE]

cohort_dt[, followup_end_date := pmin(
  death_date,
  first_stroke_date,
  STUDY_END_DATE,
  na.rm = TRUE
)]

cohort_dt[, stroke_event := 0L]
cohort_dt[
  !is.na(first_stroke_date) & first_stroke_date == followup_end_date,
  stroke_event := 1L
]

cohort_dt <- cohort_dt[
  !is.na(cohort_start_date) &
    !is.na(followup_end_date) &
    followup_end_date > cohort_start_date
]

cohort_dt[, followup_days := as.integer(followup_end_date - cohort_start_date)]
cohort_dt[, calendar_year := year(cohort_start_date)]

# -----------------------------------------------------------------------------
# 5. Time-dependent any-OAC exposure (counting process)
# -----------------------------------------------------------------------------
# Prescription periods are clipped to each patient's follow-up window, then
# converted to day offsets from cohort start for survival::tmerge.

oac_periods_dt <- merge(
  oac_periods_dt,
  cohort_dt[, .(patient_id, cohort_start_date, followup_end_date)],
  by = "patient_id"
)

oac_periods_dt <- oac_periods_dt[
  period_end_date > cohort_start_date & period_start_date < followup_end_date
]

oac_periods_dt[period_start_date < cohort_start_date, period_start_date := cohort_start_date]
oac_periods_dt[period_end_date > followup_end_date, period_end_date := followup_end_date]

oac_periods_dt[, period_start_day := as.integer(period_start_date - cohort_start_date)]
oac_periods_dt[, period_end_day := as.integer(period_end_date - cohort_start_date)]
oac_periods_dt <- oac_periods_dt[period_end_day > period_start_day]

counting_base <- cohort_dt[, .(patient_id, followup_days, stroke_event)]

cox_any_oac <- tmerge(
  data1 = counting_base,
  data2 = counting_base,
  id = patient_id,
  tstart = 0,
  tstop = followup_days,
  event = event(followup_days, stroke_event)
)

oac_status_changes <- rbindlist(list(
  oac_periods_dt[, .(patient_id, change_day = period_start_day, on_oac = 1L)],
  oac_periods_dt[, .(patient_id, change_day = period_end_day, on_oac = 0L)]
))

oac_status_changes <- merge(
  oac_status_changes,
  counting_base[, .(patient_id, followup_days)],
  by = "patient_id"
)

# Interior changes only (tmerge handles start at 0 separately)
oac_status_changes <- unique(
  oac_status_changes[change_day > 0 & change_day < followup_days]
)

cox_any_oac <- tmerge(
  data1 = cox_any_oac,
  data2 = oac_status_changes,
  id = patient_id,
  on_oac = tdc(change_day, on_oac)
)

cox_any_oac <- as.data.table(cox_any_oac)
cox_any_oac[is.na(on_oac), on_oac := 0L]

# -----------------------------------------------------------------------------
# 6. Baseline covariates merged onto counting-process rows
# -----------------------------------------------------------------------------

baseline_covariates <- cohort_dt[, .(
  patient_id,
  calendar_year = factor(calendar_year),
  chads2va_baseline,
  sex = factor(sex_binary, levels = 0:1, labels = c("female", "male")),
  prior_bleeding = factor(prior_bleeding_flag, levels = 0:1, labels = c("no", "yes")),
  renal = factor(prior_renal_flag, levels = 0:1, labels = c("no", "yes")),
  liver = factor(prior_liver_flag, levels = 0:1, labels = c("no", "yes")),
  alcohol = factor(prior_alcohol_flag, levels = 0:1, labels = c("no", "yes")),
  cancer = factor(prior_cancer_flag, levels = 0:1, labels = c("no", "yes")),
  dementia = factor(prior_dementia_flag, levels = 0:1, labels = c("no", "yes")),
  psychiatric = factor(prior_psychiatric_flag, levels = 0:1, labels = c("no", "yes")),
  income = factor(income_tertile_code, levels = 1:3, labels = c("low", "mid", "high"))
)]

baseline_covariates[, sex := relevel(sex, ref = "male")]

cox_any_oac <- merge(cox_any_oac, baseline_covariates, by = "patient_id", all.x = TRUE)
cox_any_oac <- as.data.table(cox_any_oac)

# Analysis restricted to CHA2DS2-VA score >= 1 (score 0 excluded by design)
cox_any_oac <- cox_any_oac[!is.na(chads2va_baseline) & chads2va_baseline >= 1]

cox_any_oac[, chads2va_baseline := factor(
  chads2va_baseline,
  levels = sort(unique(chads2va_baseline))
)]

chads2va_ref_level <- levels(cox_any_oac$chads2va_baseline)[1]
cox_any_oac[, chads2va_baseline := relevel(chads2va_baseline, ref = chads2va_ref_level)]

# -----------------------------------------------------------------------------
# 7. Cox models — any OAC vs no OAC
# -----------------------------------------------------------------------------

cox_any_oac_adj <- coxph(
  Surv(tstart, tstop, event) ~
    on_oac * chads2va_baseline +
    calendar_year +
    sex +
    prior_bleeding +
    renal +
    liver +
    alcohol +
    cancer +
    dementia +
    psychiatric +
    income,
  data = cox_any_oac
)

summary(cox_any_oac_adj)

cox_any_oac_adj_additive <- coxph(
  Surv(tstart, tstop, event) ~
    on_oac +
    chads2va_baseline +
    calendar_year +
    sex +
    prior_bleeding +
    renal +
    liver +
    alcohol +
    cancer +
    dementia +
    psychiatric +
    income,
  data = cox_any_oac
)

any_oac_interaction_test <- anova(
  cox_any_oac_adj_additive,
  cox_any_oac_adj,
  test = "LRT"
)

print(any_oac_interaction_test)

# -----------------------------------------------------------------------------
# 8. Stratum-specific hazard ratios (within each CHA2DS2-VA score)
# -----------------------------------------------------------------------------
# Contrasts: OAC vs no OAC at each score; reference score stratum uses main OAC effect only.

extract_any_oac_hazard_ratios <- function(model, model_data, model_label = "Adjusted") {
  coefs <- coef(model)
  vc <- vcov(model)
  chads_levels <- levels(model_data$chads2va_baseline)
  ref_score <- chads_levels[1]

  rbindlist(lapply(chads_levels, function(score_level) {
    if (score_level == ref_score) {
      beta <- coefs["on_oac"]
      var <- vc["on_oac", "on_oac"]
    } else {
      interaction_term <- paste0("on_oac:chads2va_baseline", score_level)
      if (!interaction_term %in% names(coefs)) {
        interaction_term <- paste0("chads2va_baseline", score_level, ":on_oac")
      }
      beta <- coefs["on_oac"] + coefs[interaction_term]
      var <- vc["on_oac", "on_oac"] +
        vc[interaction_term, interaction_term] +
        2 * vc["on_oac", interaction_term]
    }

    se <- sqrt(var)

    data.table(
      chads2va_score = as.numeric(as.character(score_level)),
      exposure = "Any OAC",
      model = model_label,
      hazard_ratio = exp(beta),
      ci_low = exp(beta - 1.96 * se),
      ci_high = exp(beta + 1.96 * se)
    )
  }))
}

extract_drug_class_hazard_ratios <- function(model, model_data, model_label = "Adjusted") {
  coefs <- coef(model)
  vc <- vcov(model)
  chads_levels <- levels(model_data$chads2va_baseline)
  ref_score <- chads_levels[1]

  rbindlist(lapply(chads_levels, function(score_level) {
    rbindlist(lapply(c("VKA", "DOAC"), function(drug_class) {
      main_term <- paste0("oac_drug_class", drug_class)

      if (score_level == ref_score) {
        beta <- coefs[main_term]
        var <- vc[main_term, main_term]
      } else {
        interaction_term <- paste0(
          "oac_drug_class", drug_class, ":chads2va_baseline", score_level
        )
        if (!interaction_term %in% names(coefs)) {
          interaction_term <- paste0(
            "chads2va_baseline", score_level, ":oac_drug_class", drug_class
          )
        }

        beta <- coefs[main_term] + coefs[interaction_term]
        var <- vc[main_term, main_term] +
          vc[interaction_term, interaction_term] +
          2 * vc[main_term, interaction_term]
      }

      se <- sqrt(var)

      data.table(
        chads2va_score = as.numeric(as.character(score_level)),
        exposure = drug_class,
        model = model_label,
        hazard_ratio = exp(beta),
        ci_low = exp(beta - 1.96 * se),
        ci_high = exp(beta + 1.96 * se)
      )
    }))
  }))
}

hr_any_oac_adj <- extract_any_oac_hazard_ratios(cox_any_oac_adj, cox_any_oac, "Adjusted")
print(hr_any_oac_adj)

# -----------------------------------------------------------------------------
# 9. Figure — any OAC hazard ratios by CHA2DS2-VA score
# -----------------------------------------------------------------------------

plot_any_oac <- ggplot(hr_any_oac_adj, aes(x = chads2va_score, y = hazard_ratio)) +
  add_reference_line(max(hr_any_oac_adj$chads2va_score, na.rm = TRUE)) +
  geom_errorbar(
    aes(ymin = ci_low, ymax = ci_high),
    width = 0.12,
    linewidth = 0.7,
    colour = colour_any_oac
  ) +
  geom_line(linewidth = 0.9, colour = colour_any_oac) +
  geom_point(
    size = 3.4,
    shape = 21,
    fill = "white",
    colour = colour_any_oac,
    stroke = 1.1
  ) +
  scale_x_continuous(breaks = sort(unique(hr_any_oac_adj$chads2va_score))) +
  scale_y_continuous(
    limits = common_y_limits,
    breaks = common_y_breaks,
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = x_lab_chadsva, y = y_lab_hr) +
  theme_publication() +
  theme(legend.position = "none")

print(plot_any_oac)
standard_pdf_save("HR_by_CHA2DS2_VA_any_OAC.pdf", plot_any_oac)

# -----------------------------------------------------------------------------
# 10. Crude stroke rates per 100 person-years
# -----------------------------------------------------------------------------

rate_data <- copy(cox_any_oac)
rate_data[, chads2va_score := as.numeric(as.character(chads2va_baseline))]
rate_data[, person_years := (tstop - tstart) / 365.25]

rate_data <- rate_data[
  !is.na(chads2va_score) &
    chads2va_score >= 1 &
    !is.na(on_oac) &
    person_years > 0
]

stroke_rates_by_score <- rate_data[
  ,
  .(
    n_patients = uniqueN(patient_id),
    person_years = sum(person_years, na.rm = TRUE),
    stroke_events = sum(event, na.rm = TRUE),
    rate_per_100_py = 100 * sum(event, na.rm = TRUE) / sum(person_years, na.rm = TRUE)
  ),
  by = .(chads2va_score, on_oac)
]

stroke_rates_by_score[, oac_status := fifelse(
  on_oac == 1L,
  "OAC period",
  "No OAC period"
)]

setorder(stroke_rates_by_score, chads2va_score, on_oac)
print(stroke_rates_by_score)

stroke_rates_table <- stroke_rates_by_score[
  ,
  .(
    chads2va_score,
    oac_status,
    n_patients,
    person_years = round(person_years, 1),
    stroke_events,
    rate_per_100_py = round(rate_per_100_py, 2)
  )
]

stroke_rates_html <- gt(stroke_rates_table) |>
  tab_header(
    title = "Stroke incidence rates by CHA2DS2-VA score and OAC status",
    subtitle = "Rates per 100 person-years"
  ) |>
  cols_label(
    chads2va_score = "CHA2DS2-VA score",
    oac_status = "OAC status",
    n_patients = "Patients",
    person_years = "Person-years",
    stroke_events = "Stroke events",
    rate_per_100_py = "Stroke rate per 100 PY"
  )

gtsave(stroke_rates_html, "stroke_rates_by_chadsva_oac_status.html")

# -----------------------------------------------------------------------------
# 11. Rate-based 10-year number needed to treat (NNT)
# -----------------------------------------------------------------------------
# Uses crude off-OAC rate and adjusted HR to estimate treated rate and ARD;
# NNTB = benefit, NNH = harm depending on sign of rate difference.

hr_any_oac_for_nnt <- hr_any_oac_adj[
  ,
  .(
    chads2va_score,
    hr_adjusted = hazard_ratio,
    hr_adjusted_low = ci_low,
    hr_adjusted_high = ci_high
  )
]

nnt_10y_table <- merge(
  stroke_rates_by_score[
    oac_status == "No OAC period",
    .(
      chads2va_score,
      untreated_patients = n_patients,
      untreated_person_years = person_years,
      untreated_stroke_events = stroke_events,
      untreated_rate_per_100_py = rate_per_100_py
    )
  ],
  hr_any_oac_for_nnt,
  by = "chads2va_score"
)

nnt_10y_table[, treated_rate_per_100_py := untreated_rate_per_100_py * hr_adjusted]
nnt_10y_table[, rate_difference_per_100_py := treated_rate_per_100_py - untreated_rate_per_100_py]
nnt_10y_table[, measure := fifelse(rate_difference_per_100_py > 0, "NNH", "NNTB")]
nnt_10y_table[, nnt_10y := 100 / (10 * abs(rate_difference_per_100_py))]

nnt_10y_table[, treated_rate_hr_low_per_100_py := untreated_rate_per_100_py * hr_adjusted_low]
nnt_10y_table[, treated_rate_hr_high_per_100_py := untreated_rate_per_100_py * hr_adjusted_high]
nnt_10y_table[, rd_hr_low := treated_rate_hr_low_per_100_py - untreated_rate_per_100_py]
nnt_10y_table[, rd_hr_high := treated_rate_hr_high_per_100_py - untreated_rate_per_100_py]

nnt_10y_table[, nnt_10y_lower := 100 / (10 * max(abs(rd_hr_low), abs(rd_hr_high))), by = chads2va_score]
nnt_10y_table[, nnt_10y_upper := fifelse(
  rd_hr_low * rd_hr_high <= 0,
  Inf,
  100 / (10 * pmin(abs(rd_hr_low), abs(rd_hr_high)))
)]

nnt_10y_publication <- nnt_10y_table[
  ,
  .(
    chads2va_score,
    untreated_patients,
    untreated_person_years,
    untreated_stroke_events,
    untreated_rate_per_100_py,
    hr_adjusted,
    hr_adjusted_low,
    hr_adjusted_high,
    treated_rate_per_100_py,
    rate_difference_per_100_py,
    measure,
    nnt_10y,
    nnt_10y_lower,
    nnt_10y_upper
  )
]

nnt_10y_publication[, `:=`(
  untreated_person_years = round(untreated_person_years, 1),
  untreated_rate_per_100_py = round(untreated_rate_per_100_py, 2),
  hr_adjusted = round(hr_adjusted, 3),
  hr_adjusted_low = round(hr_adjusted_low, 3),
  hr_adjusted_high = round(hr_adjusted_high, 3),
  treated_rate_per_100_py = round(treated_rate_per_100_py, 2),
  rate_difference_per_100_py = round(rate_difference_per_100_py, 2),
  nnt_10y = round(nnt_10y, 1),
  nnt_10y_lower = round(nnt_10y_lower, 1),
  nnt_10y_upper = round(nnt_10y_upper, 1)
)]

print(nnt_10y_publication)

nnt_10y_html <- gt(nnt_10y_publication) |>
  tab_header(
    title = "Rate-based 10-year NNT by CHA2DS2-VA score",
    subtitle = "Crude untreated stroke rate and adjusted CHA2DS2-VA-specific hazard ratio"
  ) |>
  cols_label(
    chads2va_score = "CHA2DS2-VA score",
    untreated_patients = "Untreated patients",
    untreated_person_years = "Untreated person-years",
    untreated_stroke_events = "Untreated stroke events",
    untreated_rate_per_100_py = "Untreated rate per 100 PY",
    hr_adjusted = "Adjusted HR",
    hr_adjusted_low = "HR 95% CI low",
    hr_adjusted_high = "HR 95% CI high",
    treated_rate_per_100_py = "Estimated treated rate per 100 PY",
    rate_difference_per_100_py = "Rate difference per 100 PY",
    measure = "Measure",
    nnt_10y = "10-year NNT/NNH",
    nnt_10y_lower = "10-year N lower",
    nnt_10y_upper = "10-year N upper"
  ) |>
  tab_spanner(
    label = "Crude untreated stroke rate",
    columns = c(
      untreated_patients,
      untreated_person_years,
      untreated_stroke_events,
      untreated_rate_per_100_py
    )
  ) |>
  tab_spanner(
    label = "Adjusted hazard ratio",
    columns = c(hr_adjusted, hr_adjusted_low, hr_adjusted_high)
  ) |>
  tab_spanner(
    label = "Rate-based 10-year NNT",
    columns = c(
      treated_rate_per_100_py,
      rate_difference_per_100_py,
      measure,
      nnt_10y,
      nnt_10y_lower,
      nnt_10y_upper
    )
  ) |>
  tab_source_note(
    source_note = "10-year NNT = 100 / (10 x absolute rate difference per 100 person-years)."
  )

gtsave(nnt_10y_html, "rate_based_10y_NNT_by_CHADSVA.html")

# -----------------------------------------------------------------------------
# 12. Time-dependent VKA and DOAC exposure
# -----------------------------------------------------------------------------
# Drug class from ATC: B01AA = VKA; B01AE/B01AF = DOAC. Overlapping periods:
# DOAC takes precedence over VKA.

oac_periods_by_drug <- copy(oac_periods_dt)

oac_periods_by_drug[, drug_class := fifelse(
  substr(ATC, 1, 5) == "B01AA", "VKA",
  fifelse(substr(ATC, 1, 5) %in% c("B01AE", "B01AF"), "DOAC", NA_character_)
)]

oac_periods_by_drug <- oac_periods_by_drug[!is.na(drug_class)]

counting_base_drug <- unique(
  cox_any_oac[, .(patient_id, followup_days = max(tstop), stroke_event = max(event)), by = patient_id]
)

cox_drug_class <- tmerge(
  data1 = counting_base_drug,
  data2 = counting_base_drug,
  id = patient_id,
  tstart = 0,
  tstop = followup_days,
  event = event(followup_days, stroke_event)
)

vka_status_changes <- rbindlist(list(
  oac_periods_by_drug[drug_class == "VKA", .(patient_id, change_day = period_start_day, on_vka = 1L)],
  oac_periods_by_drug[drug_class == "VKA", .(patient_id, change_day = period_end_day, on_vka = 0L)]
))

doac_status_changes <- rbindlist(list(
  oac_periods_by_drug[drug_class == "DOAC", .(patient_id, change_day = period_start_day, on_doac = 1L)],
  oac_periods_by_drug[drug_class == "DOAC", .(patient_id, change_day = period_end_day, on_doac = 0L)]
))

vka_status_changes <- merge(
  vka_status_changes,
  counting_base_drug[, .(patient_id, followup_days)],
  by = "patient_id"
)
doac_status_changes <- merge(
  doac_status_changes,
  counting_base_drug[, .(patient_id, followup_days)],
  by = "patient_id"
)

vka_status_changes <- unique(vka_status_changes[change_day > 0 & change_day < followup_days])
doac_status_changes <- unique(doac_status_changes[change_day > 0 & change_day < followup_days])

cox_drug_class <- tmerge(
  data1 = cox_drug_class,
  data2 = vka_status_changes,
  id = patient_id,
  on_vka = tdc(change_day, on_vka)
)

cox_drug_class <- tmerge(
  data1 = cox_drug_class,
  data2 = doac_status_changes,
  id = patient_id,
  on_doac = tdc(change_day, on_doac)
)

cox_drug_class <- as.data.table(cox_drug_class)
cox_drug_class[is.na(on_vka), on_vka := 0L]
cox_drug_class[is.na(on_doac), on_doac := 0L]

cox_drug_class[, oac_drug_class := fifelse(
  on_doac == 1L, "DOAC",
  fifelse(on_vka == 1L, "VKA", "No OAC")
)]

cox_drug_class[, oac_drug_class := factor(oac_drug_class, levels = c("No OAC", "VKA", "DOAC"))]

baseline_for_drug_models <- unique(
  cox_any_oac[, .(
    patient_id,
    chads2va_baseline,
    calendar_year,
    sex,
    prior_bleeding,
    renal,
    liver,
    alcohol,
    cancer,
    dementia,
    psychiatric,
    income
  )],
  by = "patient_id"
)

cox_drug_class <- merge(cox_drug_class, baseline_for_drug_models, by = "patient_id", all.x = TRUE)
cox_drug_class <- cox_drug_class[!is.na(chads2va_baseline)]

cox_drug_class[, chads2va_baseline := factor(
  chads2va_baseline,
  levels = sort(unique(as.numeric(as.character(chads2va_baseline))))
)]

chads2va_ref_level_drug <- levels(cox_drug_class$chads2va_baseline)[1]
cox_drug_class[, chads2va_baseline := relevel(chads2va_baseline, ref = chads2va_ref_level_drug)]

# -----------------------------------------------------------------------------
# 13. Cox models — VKA and DOAC vs no OAC
# -----------------------------------------------------------------------------

cox_drug_class_adj <- coxph(
  Surv(tstart, tstop, event) ~
    oac_drug_class * chads2va_baseline +
    calendar_year +
    sex +
    prior_bleeding +
    renal +
    liver +
    alcohol +
    cancer +
    dementia +
    psychiatric +
    income,
  data = cox_drug_class
)

summary(cox_drug_class_adj)

cox_drug_class_adj_additive <- coxph(
  Surv(tstart, tstop, event) ~
    oac_drug_class +
    chads2va_baseline +
    calendar_year +
    sex +
    prior_bleeding +
    renal +
    liver +
    alcohol +
    cancer +
    dementia +
    psychiatric +
    income,
  data = cox_drug_class
)

drug_class_interaction_test <- anova(
  cox_drug_class_adj_additive,
  cox_drug_class_adj,
  test = "LRT"
)

print(drug_class_interaction_test)

hr_drug_class_adj <- extract_drug_class_hazard_ratios(
  cox_drug_class_adj,
  cox_drug_class,
  "Adjusted"
)
print(hr_drug_class_adj)

# -----------------------------------------------------------------------------
# 14. Figure — VKA and DOAC hazard ratios by CHA2DS2-VA score
# -----------------------------------------------------------------------------

plot_vka_doac <- ggplot(
  hr_drug_class_adj,
  aes(
    x = chads2va_score,
    y = hazard_ratio,
    colour = exposure,
    fill = exposure,
    group = exposure
  )
) +
  add_reference_line(max(hr_drug_class_adj$chads2va_score, na.rm = TRUE)) +
  geom_errorbar(
    aes(ymin = ci_low, ymax = ci_high),
    width = 0.12,
    linewidth = 0.7,
    position = position_drug
  ) +
  geom_line(linewidth = 0.9, position = position_drug) +
  geom_point(
    size = 3.4,
    shape = 21,
    stroke = 1.1,
    colour = "white",
    position = position_drug
  ) +
  scale_colour_manual(
    values = c("VKA" = colour_vka, "DOAC" = colour_doac),
    labels = c("VKA", "DOAC")
  ) +
  scale_fill_manual(
    values = c("VKA" = colour_vka, "DOAC" = colour_doac),
    labels = c("VKA", "DOAC")
  ) +
  scale_x_continuous(breaks = sort(unique(hr_drug_class_adj$chads2va_score))) +
  scale_y_continuous(
    limits = common_y_limits,
    breaks = common_y_breaks,
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = x_lab_chadsva, y = y_lab_hr) +
  theme_publication() +
  theme(legend.position = c(0.82, 0.18))

print(plot_vka_doac)
standard_pdf_save("HR_by_CHA2DS2_VA_VKA_DOAC.pdf", plot_vka_doac)

# -----------------------------------------------------------------------------
# 15. Formal interaction p-values
# -----------------------------------------------------------------------------

any_oac_interaction_p <- any_oac_interaction_test[2, "Pr(>|Chi|)"]
drug_class_interaction_p <- drug_class_interaction_test[2, "Pr(>|Chi|)"]

cat("\nAny OAC x CHA2DS2-VA interaction p-value:", signif(any_oac_interaction_p, 3), "\n")
cat("\nVKA/DOAC x CHA2DS2-VA interaction p-value:", signif(drug_class_interaction_p, 3), "\n")

interaction_pvalues <- data.table(
  interaction_test = c(
    "Any OAC exposure x CHA2DS2-VA",
    "OAC drug class (VKA/DOAC) x CHA2DS2-VA"
  ),
  compared_models = c(
    "OAC + CHA2DS2-VA vs OAC x CHA2DS2-VA",
    "OAC drug class + CHA2DS2-VA vs OAC drug class x CHA2DS2-VA"
  ),
  chi_square = c(
    any_oac_interaction_test[2, "Chisq"],
    drug_class_interaction_test[2, "Chisq"]
  ),
  df = c(
    any_oac_interaction_test[2, "Df"],
    drug_class_interaction_test[2, "Df"]
  ),
  p_value = c(any_oac_interaction_p, drug_class_interaction_p)
)

interaction_pvalues[, p_value_formatted := format_p(p_value)]
print(interaction_pvalues)

# -----------------------------------------------------------------------------
# 16. Unadjusted and adjusted hazard ratio table (all exposures)
# -----------------------------------------------------------------------------

cox_any_oac_unadj <- coxph(
  Surv(tstart, tstop, event) ~ on_oac * chads2va_baseline,
  data = cox_any_oac
)

cox_drug_class_unadj <- coxph(
  Surv(tstart, tstop, event) ~ oac_drug_class * chads2va_baseline,
  data = cox_drug_class
)

hr_any_oac_unadj <- extract_any_oac_hazard_ratios(cox_any_oac_unadj, cox_any_oac, "Unadjusted")
hr_any_oac_adj_table <- extract_any_oac_hazard_ratios(cox_any_oac_adj, cox_any_oac, "Adjusted")

hr_drug_unadj <- extract_drug_class_hazard_ratios(cox_drug_class_unadj, cox_drug_class, "Unadjusted")
hr_drug_adj_table <- extract_drug_class_hazard_ratios(cox_drug_class_adj, cox_drug_class, "Adjusted")

hr_combined <- rbindlist(
  list(hr_any_oac_unadj, hr_any_oac_adj_table, hr_drug_unadj, hr_drug_adj_table),
  use.names = TRUE,
  fill = TRUE
)

hr_combined[, hr_ci := sprintf(
  "%.2f (%.2f-%.2f)",
  hazard_ratio,
  ci_low,
  ci_high
)]

hr_publication_table <- dcast(
  hr_combined,
  chads2va_score + exposure ~ model,
  value.var = "hr_ci"
)

hr_publication_table[, exposure := factor(exposure, levels = c("Any OAC", "VKA", "DOAC"))]
setorder(hr_publication_table, exposure, chads2va_score)

print(hr_publication_table)

hr_publication_html <- gt(hr_publication_table) |>
  tab_header(
    title = "Hazard ratios by CHA2DS2-VA score and OAC exposure",
    subtitle = "Unadjusted and adjusted hazard ratios with 95% confidence intervals"
  ) |>
  cols_label(
    chads2va_score = "CHA2DS2-VA score",
    exposure = "Exposure",
    Unadjusted = "Unadjusted HR (95% CI)",
    Adjusted = "Adjusted HR (95% CI)"
  ) |>
  tab_source_note(
    source_note = "Reference: no OAC within each CHA2DS2-VA score stratum."
  )

gtsave(hr_publication_html, "unadjusted_adjusted_HRs_by_CHA2DS2_VA_OAC_VKA_DOAC.html")

# -----------------------------------------------------------------------------
# 17. Save numeric outputs
# -----------------------------------------------------------------------------

fwrite(hr_any_oac_adj, "HR_by_CHA2DS2_VA_any_OAC.csv")
fwrite(hr_drug_class_adj, "HR_by_CHA2DS2_VA_VKA_DOAC.csv")
fwrite(stroke_rates_table, "stroke_rates_by_chadsva_oac_status.csv")
fwrite(nnt_10y_publication, "rate_based_10y_NNT_by_CHADSVA.csv")
fwrite(interaction_pvalues, "interaction_pvalues.csv")
fwrite(hr_publication_table, "unadjusted_adjusted_HRs_by_CHA2DS2_VA_OAC_VKA_DOAC.csv")

cat("\nAnalysis complete.\n")
cat("  Figures: HR_by_CHA2DS2_VA_any_OAC.pdf, HR_by_CHA2DS2_VA_VKA_DOAC.pdf\n")
cat("  Tables: HTML and CSV files written to the working directory.\n")
