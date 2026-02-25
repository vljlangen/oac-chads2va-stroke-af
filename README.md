# Does the association between oral anticoagulation and stroke risk differ by CHA₂DS₂-VA score in patients with atrial fibrillation? A nationwide cohort study

Analysis code for the manuscript above: Lexis-based cohort analysis, Poisson regression, IRRs, NNT, and absolute rate differences (OAC vs no OAC, by CHA₂DS₂-VA score).

**Data:** Place `IncidentCohort.sav` in the working directory.

**Run:** `source("InteractionofOACandCHADSVA.R")`

**Outputs:** Figures (JPEG/PDF) and HTML rate tables are written to the current directory.

**R packages:** haven, dplyr, survival, knitr, gridExtra, Epi, cohorttools, ggplot2, tidyr.
