/*********************************************************************
 Purpose : Reviewer 1, Comment 2 - isocaloric substitution models
           (exploratory sensitivity analysis)

 Run once per analysis:  %let analysis = inc;  (obesity-related cancer)
                         %let analysis = brc;  (death after breast cancer)
                         %let analysis = crc;  (death after colorectal cancer)

 Design
   - Exposures are cumulative-average % of energy (post-diagnosis for brc/crc), averaged in the same way as the indices.
   - Total energy is in every model, so a contrast between two included components is an isocaloric swap that holds all other
     components constant. Only these pairwise contrasts (or contrasts with the omitted reference) are reported.
   - Model 0 (overall indices): total protein and total fat; carbohydrate (all sources) is the omitted reference.
   - Model A (overall, animal, vegetable indices): animal and plant protein, animal and plant fat; carbohydrate (all sources) is the omitted reference.
   - Model B (healthy, unhealthy indices): high-quality and refined carbohydrate, animal and plant protein, animal and plant fat;
     other carbohydrate (juice, lactose, sugar in beverages) is the omitted reference.
   - Alcohol is adjusted as a categorical covariate, as in the main models.
   - Current-diet AHEI is not adjusted for, because it shares components with the contrasts (whole grains, red meat,
     unsaturated fat). Pre-diagnosis AHEI is kept in the survival models as a confounder.
   - Contrasts are estimated with PROC PHREG ESTIMATE statements (needed for delayed entry in survival). For incidence, the
     %substitute macro is run as a cross-check of Model B.
*********************************************************************/

%let analysis = inc;                 /* inc, brc, or crc */

/* Substitution sizes (% of energy) */
%let d  = 5;                         /* carbohydrate and fat swaps   */
%let dp = 5;                         /* swaps involving protein      */
%let h  = %sysevalf(&d / 2);         /* half step for combined swaps */

/* Cumulative-average %-energy variables */
%let v_totcarb = carb_svcon;         /* total carbohydrate           */
%let v_hqcarb  = hcarb_svcon;        /* high-quality carbohydrate    */
%let v_rfcarb  = sumref_carb_svcon;  /* refined carbohydrate         */
%let v_aprot   = aprot_svcon;        /* animal protein               */
%let v_vprot   = vprot_svcon;        /* plant protein                */
%let v_afat    = afat_svcon;         /* animal fat                   */
%let v_vfat    = vfat_svcon;         /* plant fat                    */


/*------------------------------------------------------------------
  Analysis-specific data, outcome, covariates
------------------------------------------------------------------*/
%macro setup;
    %global datapath indata outdir timestmt adj;

    %let outdir = lowcarb_cancer/revision/food_items/table;

    %if &analysis = inc %then %do;
        %let datapath = lowcarb_cancer/revision/site_specific/data.sas;
        %let indata   = pool;
        %let timestmt = tobeca*can_obese(0);
    %end;
    %else %do;
        %let datapath = lowcarb_cancer/revision/&analysis./fix_postdx/data.sas;
        %let indata   = pool_timevarying;
        %let timestmt = (time_begin, time_end)*dead(0);
    %end;
%mend setup;
%setup;

%include "&datapath";

/* Main-model covariates minus total energy (modelled separately)
   and current-diet AHEI                                          */
%macro setadj;
    %if &analysis = inc %then %let adj = &alc_ &smkstatus_ &actc_ &bmib9_
            white mv asp cafh endo2c mam pmhever ocever nSES;
    %else %if &analysis = brc %then %let adj = age_dx &smkstaa_ &alcovc_post_ &postact_  &bmid9_
            white mv asp cafh pmhever ocever
            &est_status2_ &stage_ her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES;
    %else %if &analysis = crc %then %let adj = age_dx &smkstaa_ &alcovc_post_ &postact_  &bmid3cat_
            white mv asp cafh  
            &stage_cr_ &grade_cr_ &subsite_cr_ nSES;
%mend setadj;
%setadj;


/*===================================================================
  1. Components and checks
===================================================================*/
data sub_in;
    set &indata;
    totcarb = &v_totcarb;
    hqcarb  = &v_hqcarb;
    rfcarb  = &v_rfcarb;
    ocarb   = totcarb - hqcarb - rfcarb;            /* other carbohydrate */
    aprot   = &v_aprot;
    vprot   = &v_vprot;
    afat    = &v_afat;
    vfat    = &v_vfat;
    tprot   = aprot + vprot;                        /* total protein      */
    tfat    = afat + vfat;                          /* total fat          */
    /* carbohydrate + protein + fat = 100% minus alcohol             */
    pe_sum  = sum(totcarb, aprot, vprot, afat, vfat);
    %if &analysis ne inc %then %do;
        if person_time > 0;
    %end;
run;

/* Checks: ocarb should not be negative; pe_sum should be close to
   100 (minus alcohol); NMiss should be 0; each swap size (&d, &dp)
   should be within the interquartile range of the components.      */
title "&analysis: components (% of energy)";
proc means data=sub_in n nmiss mean p5 p25 median p75 p95 maxdec=1;
    var totcarb hqcarb rfcarb ocarb aprot vprot afat vfat pe_sum;
run;
title;


/*===================================================================
  1b. Model 0 - total protein and total fat; carbohydrate = reference
      Overall LCD and LFD (combined substitutions)
===================================================================*/
ods select Estimates;
ods output Estimates = subst0;
proc phreg data=sub_in;
    model &timestmt = tprot tfat &calorsvq_ &adj / ties=efron;
    strata agecon cohort intervalp;

    estimate "Overall LCD: protein + fat for carb (&d.%)"   tprot &h  tfat &h   / exp cl;
    estimate "Overall LFD: carb + protein for fat (&d.%)"   tprot &h  tfat -&d  / exp cl;
run;
ods select all;


/*===================================================================
  2. Model A - carbohydrate (all sources) is the omitted reference
     Overall, animal and vegetable indices
     "A for B" = increase A and decrease B by the same % of energy
===================================================================*/
ods select Estimates;
ods output Estimates = substA;
proc phreg data=sub_in;
    model &timestmt = aprot vprot afat vfat &calorsvq_ &adj / ties=efron;
    strata agecon cohort intervalp;

    /* Combined substitutions (one per index) */
    estimate "Vegetable LCD: plant protein + plant fat for carb (&d.%)"  vprot &h  vfat &h         / exp cl;
    estimate "Animal LCD: animal protein + animal fat for carb (&d.%)"   aprot &h  afat &h         / exp cl;
    estimate "Vegetable LFD: carb + plant protein for animal fat (&d.%)" vprot &h  afat -&d        / exp cl;
    estimate "Animal LFD: carb + animal protein for plant fat (&d.%)"    aprot &h  vfat -&d        / exp cl;

    /* Single-component substitutions */

    /* Vegetable LCD */
    estimate "Plant protein for carb (&dp.%)"       vprot  &dp            / exp cl;
    estimate "Plant fat for carb (&d.%)"            vfat   &d             / exp cl;
    /* Animal LCD */
    estimate "Animal protein for carb (&dp.%)"      aprot  &dp            / exp cl;
    estimate "Animal fat for carb (&d.%)"           afat   &d             / exp cl;
    /* Vegetable LFD */
    estimate "Carb for animal fat (&d.%)"           afat  -&d             / exp cl;
    estimate "Plant protein for animal fat (&dp.%)" vprot  &dp afat -&dp  / exp cl;
    /* Animal LFD */
    estimate "Carb for plant fat (&d.%)"            vfat  -&d             / exp cl;
    estimate "Animal protein for plant fat (&dp.%)" aprot  &dp vfat -&dp  / exp cl;

    estimate "Carb for animal protein (&d.%)"       aprot  -&d             / exp cl;
    estimate "plant fat for animal protein (&dp.%)" vfat  &dp  aprot -&dp  / exp cl;
    estimate "Carb for plant protein (&d.%)"        vprot  -&d             / exp cl;
run;
ods select all;


/*===================================================================
  3. Model B - other carbohydrate is the omitted reference
     Healthy and unhealthy indices, and source within macronutrient
===================================================================*/
ods select Estimates;
ods output Estimates = substB;
proc phreg data=sub_in;
    model &timestmt = rfcarb hqcarb aprot vprot afat vfat
                      &calorsvq_ &adj / ties=efron;
    strata agecon cohort intervalp;

    /* Combined substitutions (one per index) */
    estimate "Healthy LCD: plant protein + plant fat for refined carb (&d.%)"           vprot  &h vfat  &h rfcarb -&d / exp cl;
    estimate "Unhealthy LCD: animal protein + animal fat for high-quality carb (&d.%)"  aprot  &h afat  &h hqcarb -&d / exp cl;
    estimate "Healthy LFD: high-quality carb + plant protein for animal fat (&d.%)"     hqcarb &h vprot &h afat   -&d / exp cl;
    estimate "Unhealthy LFD: refined carb + animal protein for plant fat (&d.%)"        rfcarb &h aprot &h vfat   -&d / exp cl;

    /* Single-component substitutions */

    /* Source within macronutrient */
    estimate "Plant for animal protein (&dp.%)"               vprot  &dp aprot  -&dp / exp cl;
    estimate "Plant for animal fat (&d.%)"                    vfat   &d  afat   -&d  / exp cl;
    estimate "High-quality for refined carb (&d.%)"           hqcarb &d  rfcarb -&d  / exp cl;
    /* Healthy LCD */
    estimate "Plant protein for refined carb (&dp.%)"         vprot  &dp rfcarb -&dp / exp cl;
    estimate "Plant fat for refined carb (&d.%)"              vfat   &d  rfcarb -&d  / exp cl;
    /* Unhealthy LCD */
    estimate "Animal protein for high-quality carb (&dp.%)"   aprot  &dp hqcarb -&dp / exp cl;
    estimate "Animal fat for high-quality carb (&d.%)"        afat   &d  hqcarb -&d  / exp cl;
    /* Healthy LFD */
    estimate "High-quality carb for animal fat (&d.%)"        hqcarb &d  afat   -&d  / exp cl;
    estimate "Plant protein for animal fat (&dp.%)"           vprot  &dp afat   -&dp / exp cl;
    /* Unhealthy LFD */
    estimate "Refined carb for plant fat (&d.%)"              rfcarb &d  vfat   -&d  / exp cl;
    estimate "Animal protein for plant fat (&dp.%)"           aprot  &dp vfat   -&dp / exp cl;
    
    estimate "Animal protein for refined carb (&dp.%)"        aprot &dp  rfcarb -&dp / exp cl;
    estimate "Animal fat for refined carb (&d.%)"             afat &d    rfcarb -&d  / exp cl;
    estimate "Plant protein for high-quality carb (&dp.%)"    vprot &dp  hqcarb -&dp / exp cl;
    estimate "Plant fat for high-quality carb (&d.%)"         vfat &d    hqcarb -&d  / exp cl;

    estimate "Plant fat for animal protein (&d.%)"         vfat &d    aprot -&d  / exp cl;
    estimate "Animal fat for animal protein (&d.%)"             afat &d    aprot -&d  / exp cl;
    estimate "High-quality carb for animal protein (&dp.%)"    hqcarb &dp  aprot -&dp / exp cl;

    estimate "Plant fat for plant protein (&d.%)"         vfat &d    vprot -&d  / exp cl;
    estimate "Animal fat for plant protein (&d.%)"             afat &d    vprot -&d  / exp cl;
    estimate "High-quality carb for plant protein (&dp.%)"    hqcarb &dp  vprot -&dp / exp cl;

run;
ods select all;

/*===================================================================
  4. Add model metadata and create publication-ready results

  Model 0:
    Included nutrient components = total protein, total fat
    Omitted component            = total carbohydrate

  Model A:
    Included components = animal protein, plant protein,
                          animal fat, plant fat
    Omitted component   = total carbohydrate

  Model B:
    Included components = high-quality carb, refined carb,
                          animal protein, plant protein,
                          animal fat, plant fat
    Omitted component   = other carbohydrate

  All models additionally adjust for total energy and the
  analysis-specific covariate set.
===================================================================*/


/*------------------------------------------------------------------
  Model 0
------------------------------------------------------------------*/

data subst0_detail;
    length analysis           $4
           analysis_label     $50
           panel              $60
           model_code         $8
           model_description  $100
           included_components $300
           reference_component $60
           contrast_type      $30
           Label              $120
           hr_ci              $40
           energy_adjustment  $200
           adjustment_set     $1000;

    set subst0;

    analysis = "&analysis";

    if "&analysis" = "inc" then
        analysis_label = "Obesity-related cancer incidence";
    else if "&analysis" = "brc" then
        analysis_label = "Breast cancer survival";
    else if "&analysis" = "crc" then
        analysis_label = "Colorectal cancer survival";

    model_code = "Model 0";

    model_description =
        "Overall macronutrient substitution model";

    included_components =
        "Total protein; total fat";

    reference_component =
        "Total carbohydrate";

    panel =
        "A. Dietary-pattern substitutions";

    contrast_type =
        "Combined";

    energy_adjustment = "&calorsvq_";

    adjustment_set = symget("adj");

    hr_ci = cats(
        put(ExpEstimate, 6.2),
        " (",
        put(LowerExp, 6.2),
        ", ",
        put(UpperExp, 6.2),
        ")"
    );
run;


/*------------------------------------------------------------------
  Model A
------------------------------------------------------------------*/

data substA_detail;
    length analysis           $4
           analysis_label     $50
           panel              $60
           model_code         $8
           model_description  $100
           included_components $300
           reference_component $60
           contrast_type      $30
           Label              $120
           hr_ci              $40
           energy_adjustment  $200
           adjustment_set     $1000;

    set substA;

    analysis = "&analysis";

    if "&analysis" = "inc" then
        analysis_label = "Obesity-related cancer incidence";
    else if "&analysis" = "brc" then
        analysis_label = "Breast cancer survival";
    else if "&analysis" = "crc" then
        analysis_label = "Colorectal cancer survival";

    model_code = "Model A";

    model_description =
        "Carbohydrate not separated by quality";

    included_components =
        "Animal protein; plant protein; animal fat; plant fat";

    reference_component =
        "Total carbohydrate";

    /*
       Labels containing ":" are the combined substitutions
       corresponding directly to dietary-pattern indices.
    */

    if index(Label, ":") > 0 then do;
        panel =
            "A. Dietary-pattern substitutions";
        contrast_type =
            "Combined";
    end;

    else do;
        panel =
            "B. Individual substitutions: carbohydrate unsplit";
        contrast_type =
            "Single-component";
    end;

    energy_adjustment = "&calorsvq_";

    adjustment_set = symget("adj");

    hr_ci = cats(
        put(ExpEstimate, 6.2),
        " (",
        put(LowerExp, 6.2),
        ", ",
        put(UpperExp, 6.2),
        ")"
    );
run;


/*------------------------------------------------------------------
  Model B
------------------------------------------------------------------*/

data substB_detail;
    length analysis           $4
           analysis_label     $50
           panel              $60
           model_code         $8
           model_description  $100
           included_components $300
           reference_component $60
           contrast_type      $30
           Label              $120
           hr_ci              $40
           energy_adjustment  $200
           adjustment_set     $1000;

    set substB;

    analysis = "&analysis";

    if "&analysis" = "inc" then
        analysis_label = "Obesity-related cancer incidence";
    else if "&analysis" = "brc" then
        analysis_label = "Breast cancer survival";
    else if "&analysis" = "crc" then
        analysis_label = "Colorectal cancer survival";

    model_code = "Model B";

    model_description =
        "Carbohydrate separated into high-quality and refined sources";

    included_components =
        "High-quality carbohydrate; refined carbohydrate; "
        || "animal protein; plant protein; animal fat; plant fat";

    reference_component =
        "Other carbohydrate";

    if index(Label, ":") > 0 then do;
        panel =
            "A. Dietary-pattern substitutions";
        contrast_type =
            "Combined";
    end;

    else do;
        panel =
            "C. Individual substitutions: carbohydrate quality separated";
        contrast_type =
            "Single-component";
    end;

    energy_adjustment = "&calorsvq_";

    adjustment_set = symget("adj");

    hr_ci = cats(
        put(ExpEstimate, 6.2),
        " (",
        put(LowerExp, 6.2),
        ", ",
        put(UpperExp, 6.2),
        ")"
    );
run;


/*------------------------------------------------------------------
  Combine all results
------------------------------------------------------------------*/

data subst_&analysis;
    set subst0_detail
        substA_detail
        substB_detail;

    p_value = Probz;

    keep analysis
         analysis_label
         panel
         model_code
         model_description
         included_components
         reference_component
         contrast_type
         Label
         ExpEstimate
         LowerExp
         UpperExp
         p_value
         hr_ci
         energy_adjustment
         adjustment_set;
run;


/*------------------------------------------------------------------
  Print detailed output for QC
------------------------------------------------------------------*/

proc sort data=subst_&analysis;
    by panel model_code;
run;

title "&analysis: Detailed isocaloric substitution models";

proc print data=subst_&analysis noobs label;
    by panel model_code;

    var Label
        reference_component
        hr_ci
        p_value;

    label
        Label               = "Substitution contrast"
        reference_component = "Omitted/reference component"
        hr_ci               = "HR (95% CI)"
        p_value             = "P value";

    format p_value pvalue6.4;
run;

title;


/*------------------------------------------------------------------
  Export FULL detailed output
------------------------------------------------------------------*/

proc export
    data=subst_&analysis
    outfile="&outdir./substitution_&analysis._detailed.csv"
    dbms=csv
    replace;
run;


/*------------------------------------------------------------------
  Create compact publication table
------------------------------------------------------------------*/

data subst_&analysis._table;
    set subst_&analysis;

    keep analysis_label
         panel
         model_code
         Label
         reference_component
         hr_ci
         p_value;
run;

proc export
    data=subst_&analysis._table
    outfile="&outdir./substitution_&analysis._table.csv"
    dbms=csv
    replace;
run;