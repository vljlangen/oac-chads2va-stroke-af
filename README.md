# Does the association between oral anticoagulation and stroke risk differ by CHA₂DS₂-VA score in patients with atrial fibrillation? A nationwide cohort study

Analysis code for the manuscript above: counting-process Cox models, stratum-specific hazard ratios by CHA₂DS₂-VA score, stroke rates, rate-based 10-year NNT, and VKA/DOAC comparisons.

## Data (not included in this repository)

Place these files in the working directory under the **generic names** expected by the script (rename your local exports; internal project filenames are not used in the public code):

| File | Description |
|------|-------------|
| `cohort_baseline.sav` | One row per patient: cohort entry, follow-up dates, stroke, death, CHA₂DS₂-VA score, baseline covariates |
| `oac_prescription_periods.rds` | Oral anticoagulant purchase periods per patient (ATC code, period start/end) |

Column names in your export may differ; map them in **Section 2** of `OAC_CHA2DS2VA.R` (`cohort_name_map` and `oac_name_map`).

## Run

```r
source("OAC_CHA2DS2VA.R")
```

## Outputs

Written to the working directory:

- **Figures:** `HR_by_CHA2DS2_VA_any_OAC.pdf`, `HR_by_CHA2DS2_VA_VKA_DOAC.pdf`
- **Tables:** HTML (`gt`) and CSV (hazard ratios, stroke rates, NNT, interaction *p*-values)

## R packages

`data.table`, `survival`, `haven`, `lubridate`, `ggplot2`, `gt`; optional `showtext` and `sysfonts` for the Rosario figure font (requires a one-time internet connection).
