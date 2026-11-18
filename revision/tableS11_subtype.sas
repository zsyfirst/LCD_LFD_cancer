/*****************************************************
Purpose:       low carb, low fat, and risk of Obesity-related cancer
Study design:  prospective cohort (NHS: 1986-2022.6, NHSII: 1991-2023.6, HPFS: 1986-2024.1)
Contents:      Table S11 --- HR (Subtype Heterogeneity Evaluations)
REvision:      Disease subtype heterogeneity analysis: obesity-related cancer sites
Date:          10/2026
******************************************************/

/************************ Call in Cleaned Datasets**********************************/
%include '/lowcarb_cancer/revision/data.inc.sas';


/*
Primary heterogeneity analysis:
- Outcome: first obesity-related cancer subtype (1-13)
- Time: existing person-period time variable, tobeca
- Exposures: dietary scores per SD
- Covariates: constrained to have common associations across sites
- Baseline hazard: stratified by age, cohort, questionnaire interval
*/

%subtype(
    data         = pool,
    studydesign  = COHORT,
    id           = newid,

    augmented    = NO,

    time         = tobeca,
    eventtype    = can_obese_site,

    exposure     = lcs_sd
                   alcs_sd
                   vlcs_sd
                   ulcs_sd
                   hlcs_sd
                   lfs_sd
                   ulfs_sd
                   hlfs_sd
                   alfs_sd
                   vlfs_sd,

    constrvar    = &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam 
                   nSES,

    stratavar    = agecon cohort intervalp,

    covs         = YES,
    wald         = YES,

    eventtypelabel =
                   1=Esophageal;
                   2=Gastric;
                   3=Colorectal;
                   4=Liver;
                   5=Gallbladder;
                   6=Pancreatic;
                   7=Breast;
                   8=Endometrial;
                   9=Ovarian;
                   10=Kidney;
                   11=Thyroid;
                   12=Myeloma;
                   13=Meningioma,

    paramest     = obeca_subtype_est,
    heterotest   = obeca_subtype_hetero
);

proc print data=obeca_subtype_hetero noobs;
run;

proc print data=obeca_subtype_est noobs;
run;





