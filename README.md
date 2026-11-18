# Low-Carbohydrate and Low-Fat Dietary Patterns and Obesity-Related Cancer Risk and Survival

Analysis code for:
Zou S, et al. *Low-Carbohydrate and Low-Fat Dietary Patterns and Obesity-Related Cancer Risk and Survival: Prospective Cohort Studies.*.

The code for the revised manuscript is in the folder [`revision`](revision). Each file is named after the main figure, Extended Data figure ("FigureS*") or Supplementary Table ("tableS*") that it produces.

## 1. System requirements

| Software | Version tested | Used for |
|---|---|---|
| SAS | 9.4 (SAS Institute, Cary, NC) | Data management, Cox models, substitution, IPW and all tables |
| R | 4.5.1 | Correlation heatmaps, forest plots, mediation analysis |

- R packages (versions tested): data.table 1.18.2.1; tidyverse 2.0.0, dplyr 1.2.1, tidyr 1.3.2, stringr 1.6.0, readr 2.2.0, ggplot2 4.0.2; readxl 1.4.5; ggcorrplot 0.1.4.1; corrplot 0.95; patchwork 1.3.2; msm 1.8.2; CMAverse 0.1.0
- SAS macros from the Channing Division of Network Medicine (Brigham and Women's Hospital and Harvard Medical School): `%table1`, `%mphreg`, `%pre_pm`, `%subtype` and `%indic3`. All other macros are defined within the scripts. Brief descriptions of the Channing macros are provided in [`channing_macros.md`](channing_macros.md).
- Operating system: Linux (Channing computing cluster). The code is not platform-specific.
- No non-standard hardware is required.

## 2. Installation guide

No installation is required beyond licensed copies of SAS 9.4 and R 4.5.1. Download or clone this repository:

```
git clone https://github.com/zsyfirst/LCD_LFD_cancer.git
```

Typical install time: under 1 minute (copying the scripts only).

## 3. Demo

No demo dataset is provided. Data from the Nurses' Health Study (NHS), Nurses' Health Study II (NHSII) and Health Professionals Follow-up Study (HPFS) cannot be shared publicly because of participant confidentiality (see Data availability below). The code is shared to document every analytic step reported in the manuscript.

## 4. Instructions for use

Approved users with access to the cohort data should first update the library paths at the top of each SAS program. Run the data preparation programs first, then the analysis programs.

### Data preparation

Cohort-specific programs read the NHS, NHSII and HPFS questionnaire, diet and outcome data, apply exclusions, construct the LCD and LFD indices and their cumulative averages, and create covariates. The pooling programs then combine the three cohorts.

| Step | Files | Output |
|---|---|---|
| 1. Incidence analysis, cohort-specific | `nhs.data.inc.sas`, `nhs2.data.inc.sas`, `hpfs.data.inc.sas` | Cohort-specific datasets for obesity-related cancer incidence |
| 2. Incidence analysis, pooled | `data.inc.sas` | Pooled incidence dataset (NHS, NHSII, HPFS) |
| 3. Breast cancer survival, cohort-specific | `nhs.data.brc_postdx.sas`, `nhs2.data.brc_postdx.sas` | Post-diagnosis datasets for breast cancer survivors |
| 4. Breast cancer survival, pooled | `data_postdx_brc.sas` | Pooled breast cancer survival dataset (NHS, NHSII) |
| 5. Colorectal cancer survival, cohort-specific | `nhs.data.crc_postdx.sas`, `nhs2.data.crc_postdx.sas`, `hpfs.data.crc_postdx.sas` | Post-diagnosis datasets for colorectal cancer survivors |
| 6. Colorectal cancer survival, pooled | `data_postdx_crc.sas` | Pooled colorectal cancer survival dataset (NHS, NHSII, HPFS) |

The analysis programs below use the pooled datasets.

### Main figures

| File | Output |
|---|---|
| `Figure2_figureS2.R` | Fig. 2; Extended Data Fig. 2 (correlations among diet indices, nutrients and foods) |
| `Figure3_tableS6.sas` | Fig. 3; Supplementary Table 6 (obesity-related cancer incidence) |
| `Figure3.R` | Fig. 3 (forest plot) |
| `Figure4_table14.sas` | Fig. 4; Supplementary Table 14 (all-cause mortality, breast cancer survivors) |
| `Figure5_table15.sas` | Fig. 5; Supplementary Table 15 (all-cause mortality, colorectal cancer survivors) |

### Extended Data figures

| File | Output |
|---|---|
| `FigureS3.sas` | Extended Data Fig. 3 (restricted cubic splines) |
| `FigureS4S5.sas` | Extended Data Figs. 4–5 (cause-specific mortality, breast cancer survivors) |
| `FigureS6S7.sas` | Extended Data Figs. 6–7 (cause-specific mortality, colorectal cancer survivors) |
| `FigureS9.sas` | Extended Data Fig. 9 (associations by calendar period, 1990–2010) |
| `tableS21.*` | Extended Data Fig. 10 (dietary assessment time windows) |

### Supplementary Tables

| File | Supplementary Table |
|---|---|
| `tableS2.sas` | 2. Baseline characteristics |
| `tableS3.sas` | 3. Energy-adjusted nutrient and food intakes |
| `tableS4.sas` | 4. Characteristics at the midpoint of follow-up |
| `tableS5.sas` | 5. Case counts for the 13 component obesity-related cancers |
| `Figure3_tableS6.sas`, `tableS6_parsimony.sas` | 6. Main incidence models, including parsimonious models |
| `tableS7.sas` | 7. Cohort-specific estimates |
| `tableS8exclude_bc.sas`, `tableS8exclude_bcop.sas` | 8. Excluding the most common cancer sites |
| `tableS9_substitution_combined.sas` | 9. Isocaloric substitution (incidence and survival) |
| `tableS10*.sas` (brc, crc, endo, kidney, ovarian, panc, thyroid, othsite) | 10. Site-specific cancer incidence |
| `tableS11_subtype.sas` | 11. Heterogeneity across cancer sites |
| `tableS12_brc.sas`, `tableS12_crc.sas` | 12. Characteristics of breast and colorectal cancer survivors |
| `tableS13,R1,2,3.1.sas`, `tableS13,R1,2,3.2.sas` | 13. Comparison with SEER survival; reviewer Tables R1–R3 |
| `Figure4_table14.sas`, `Figure4_table14_parsimony.sas` | 14. Breast cancer survivors, all-cause mortality |
| `Figure5_table15.sas`, `Figure5_table15_parsimony.sas` | 15. Colorectal cancer survivors, all-cause mortality |
| `tableS16_mediation_allscores.R` | 16. Mediation through BMI |
| `tableS17_bmic9.sas`, `tableS17_withoutBMI.sas` | 17. With and without BMI adjustment |
| `tableS18.sas` | 18. Diet change and 4-year weight change |
| `tableS19.sas` | 19. Stratified analyses |
| `tableS20.1.sas` – `tableS20.3.sas` | 20. Lag analyses, continuously updated diet, no AHEI adjustment |
| `tableS21.1_*.sas` – `tableS21.3_*.sas` (brc, crc) | 21. Alternative dietary assessment windows |
| `tableS22.1.sas`, `tableS22.2.sas` | 22. Stratified by time since diagnosis |
| `tableS23.1.sas` – `tableS23.3.sas` | 23. Excluding early deaths; inverse probability of response weighting |

Supplementary Table 1 (index construction) and Supplementary Table 5 (case counts) are descriptive; the scoring algorithm is described in the Supplementary Methods.

## Data availability

Because of participant confidentiality and privacy concerns, data cannot be shared publicly. Requests to access NHS, NHSII and HPFS data are reviewed through standard controlled-access procedures; see https://www.nurseshealthstudy.org/researchers and https://hsph.harvard.edu/research/health-professionals/resources/for-external-collaborators/.

## License

This code is released under the MIT License (see `LICENSE`).

## Contact

Siyu Zou and Qi Sun, Department of Nutrition, Harvard T.H. Chan School of Public Health.
