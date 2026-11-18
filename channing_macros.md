# Channing Division macros used in this analysis

The SAS macros below were developed by the Channing Division of Network Medicine (Brigham and Women's Hospital and Harvard Medical School) and are available from the Channing Division upon request. Brief descriptions are given here; all other macros are defined within the analysis scripts.

| Macro | Purpose | Key arguments |
|---|---|---|
| `%indic3` | Creates indicator (dummy) variables for a categorical variable, with an optional indicator for missing values, and a global macro variable (`&<prefix>_`) listing the indicators for use in model statements. Used to code categorical covariates. | `vbl` (categorical variable), `reflev` (reference level), `min`/`max` or `values` (non-reference levels), `missing` (missing value), `prefix` (indicator names), `usemiss` (1 = include missing indicator; 0 = do not) |
| `%table1` | Computes age-standardized means (SD) for continuous variables and percentages for categorical variables by level of an exposure, and writes the results to a formatted table. Used for the descriptive tables (Supplementary Tables 2, 4 and 12). | `data`, `exposure` (categories coded 1, 2, …), `agegroup` (for age standardization), `varlist`, `cat` (binary variables), `poly` (multi-category variables), `noadj` (variables not age-standardized, e.g. age) |
| `%pre_pm` | Prepares a person-time dataset (cases, person-time and exposure categories by follow-up period) for stratified person-time analysis with the companion `%pm` macro, which computes incidence rates, Mantel–Haenszel rate ratios and the Breslow–Day test for trend (Breslow and Day, 1987). Used to obtain case counts and person-years. | `data`, `out`, `case` (case variables), `var` (exposure and stratification variables), `timevar` (follow-up period), `irt` (questionnaire return dates), `cutoff` (end of follow-up), `dtdth` (date of death), `dtdx` (date of diagnosis) |
| `%mphreg9` | Fits Cox proportional hazards models (PROC PHREG) for one or more exposures and models, and outputs hazard ratios with 95% confidence intervals and tests for trend. Used for the incidence and survival analyses. | exposures, covariates, strata (age, calendar period and cohort), event and time variables |
| `%subtype` | Tests whether associations of an exposure differ across disease subtypes, with an overall heterogeneity test (and optional pairwise tests) (Wang et al., 2016). Used for the test of heterogeneity across obesity-related cancer sites. | outcome subtypes, exposure(s), covariates, `pairwise` (T/F) |

## Example call

```sas
%indic3(vbl=phmsstatus, prefix=phmsstatus, reflev=1, min=2, max=4, 
                missing=0, usemiss=1,
                label1='premenopause',
                label2='pmh never use',
                label3='pmh current use',
                label4='pmh past use');
```

```sas
%table1(
   data=alldat, 
   agegroup=agegrp, 
   exposure=shift,
   varlist=ageyr bmi mets csmk cpmh thz ost null chld calc vitd prot retn alco caff,
   noadj=ageyr,
   cat=csmk cpmh thz ost null,
   rtftitle=Age-standardized characteristics of the study population of postmenopausal nurses 
            at 1988 baseline by number of years spent working in rotating night shifts,
            landscape=F,
   fn=@mets Metabolic equivalents from recreational and leisure- time activities. 
      @cpmh Postmenopausal hormone replacement therapy. 
      @chld Number of children among parous women.,
   file=testf1, 
   uselbl=F);
```

```sas
%pre_pm(data=indata, 
  out=junkdata,
  case=non_fat fatmi totmi, 
  timevar=timeper,
  irt=irt86 irt88 irt90,
  cutoff=irt92,
  dtdth=dtdth,  
  dtdtx=dtdtx,
  var=agegp bmi cholq);
```


```sas
%MPHREG9 (DATA=nhs,
	EVENT=BCCASE, TIME=BCTIME, TVAR=PERIOD,
	TIMEVAR=t84 t86 t88 t90 t92 t94 t96 t98 t00 t02 t04 t06 t08 t10 t12,
	outdat=hr_1, labels=T, id=id, agevar=agemonths,
	pyby1=&bcaaqc_,
	contvar=height,
	extrav=&bcaaqc_ &isoqc_ fambc height,
	MODEL1=&bcaaqc_ fambc ,
	MODEL2=&isoqc_ fambc height );
```

```sas
%subtype(data=cohort1, studydesign=cohort, id=id,
	exposure=alcohol, augmented=no, time=time, eventtype=cancer,
	unconstrvar=ause_p2 screen2 polyps2 cafam2
	py30ct2 py30ct3 py30ct4 py30ct5 py30ctm
	actct2 actct3 actct4 actct5 actctm
	mvit2 mvitm bmain2 bmain3 bmain4
	bmi2 bmi3 bmi4 bmi5 bmim, 
	stratavar=agemo period,
	eventtypelabel=1=high; 2=medium; 3=low,
	heterotest=heterogeneity);
```

## References

- Wang M, Spiegelman D, Kuchiba A, et al. Statistical methods for studying disease subtype heterogeneity. *Stat Med.* 2016;35(5):782–800. doi:10.1002/sim.6793.
