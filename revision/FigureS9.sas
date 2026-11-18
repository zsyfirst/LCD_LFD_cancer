/*****************************************************
Purpose:       low carb, low fat, and risk of Obesity-related cancer
Study design:  prospective cohort (NHS: 1986-2022.6, NHSII: 1991-2023.6, HPFS: 1986-2024.1)

Contents:      HR for Figure S9 --- HR at each interval (Obesity-related cancer incidence as outcome)
               Trend of the Associations Between Low-Fat and Low-Carbohydrate Diet Scores and Obesity-related Cancer Risk in the Pooled Dataset, 1990-2010
Revision:      add nSES for cavariates
Date:          10/2026
******************************************************/

/************************ Call in Cleaned Datasets**********************************/
%include '/lowcarb_cancer/revision/data.inc.sas';


/************************ BASELINE LCD ANALYSIS *****************************************/
proc sort data=pool; by newid intervalp; run;

/* Create baseline LCD variables (first available observation) */
data pool_baseline;
    set pool;
    by newid;
    /* Keep only first observation per person for baseline values */
    if first.newid;
    
    lcs_baseline  = lcssvcon;
    alcs_baseline = alcssvcon;
    vlcs_baseline = vlcssvcon;
    ulcs_baseline = ulcssvcon;
    hlcs_baseline = hlcssvcon;

    lfs_baseline  = lfssvcon;
    ulfs_baseline = ulfssvcon;
    hlfs_baseline = hlfssvcon;
    alfs_baseline = alfssvcon;
    vlfs_baseline = vlfssvcon;

    calor_baseline = calorsvcon;
    ahei_baseline  = aheimsvcon;
    smk_baseline   = smkstatus;
    alc_baseline   = alc;
    act_baseline   = actc;
    bmi_baseline   = bmib9;
    mam_baseline    = mam;
    endo2c_baseline = endo2c;
    endo3c_baseline = endo3c;
    screen_baseline  = screen_cancer ;
    nSES_baseline    = nSES;
    
    keep newid lcs_baseline alcs_baseline vlcs_baseline ulcs_baseline hlcs_baseline 
         lfs_baseline ulfs_baseline hlfs_baseline alfs_baseline vlfs_baseline
         calor_baseline ahei_baseline smk_baseline alc_baseline act_baseline bmi_baseline 
         mam_baseline  endo2c_baseline endo3c_baseline screen_baseline nSES_baseline ;
run;

/* Merge baseline values back to main dataset */
proc sort data=pool; by newid; run;
proc sort data=pool_baseline; by newid; run;

data pool;
    merge pool pool_baseline;
    by newid;
run;

/* Create quintiles for baseline LCD */
proc sort data=pool; by intervalp; run;

proc rank data=pool groups=5 out=pool;
    by intervalp;
    var lcs_baseline alcs_baseline vlcs_baseline ulcs_baseline hlcs_baseline 
        lfs_baseline ulfs_baseline hlfs_baseline alfs_baseline vlfs_baseline
        calor_baseline ahei_baseline;
    ranks lcs_baseq alcs_baseq vlcs_baseq ulcs_baseq hlcs_baseq 
          lfs_baseq ulfs_baseq hlfs_baseq alfs_baseq vlfs_baseq
          calor_baseq ahei_baseq;
run;

/* Convert to 1-5 scale and create indicators */
data pool;
    set pool end=_end_;
    array baseq lcs_baseq alcs_baseq vlcs_baseq ulcs_baseq hlcs_baseq 
                lfs_baseq ulfs_baseq hlfs_baseq alfs_baseq vlfs_baseq
                calor_baseq ahei_baseq;
    do over baseq;
        baseq = baseq + 1;
    end;

%indic3(vbl=lcs_baseq,  prefix=lcs_baseq,  min=2, max=5, reflev=1, missing=., usemiss=0);
%indic3(vbl=alcs_baseq, prefix=alcs_baseq, min=2, max=5, reflev=1, missing=., usemiss=0);
%indic3(vbl=vlcs_baseq, prefix=vlcs_baseq, min=2, max=5, reflev=1, missing=., usemiss=0);
%indic3(vbl=ulcs_baseq, prefix=ulcs_baseq, min=2, max=5, reflev=1, missing=., usemiss=0);
%indic3(vbl=hlcs_baseq, prefix=hlcs_baseq, min=2, max=5, reflev=1, missing=., usemiss=0);
%indic3(vbl=lfs_baseq,  prefix=lfs_baseq,  min=2, max=5, reflev=1, missing=., usemiss=0);
%indic3(vbl=ulfs_baseq, prefix=ulfs_baseq, min=2, max=5, reflev=1, missing=., usemiss=0);
%indic3(vbl=hlfs_baseq, prefix=hlfs_baseq, min=2, max=5, reflev=1, missing=., usemiss=0);
%indic3(vbl=alfs_baseq, prefix=alfs_baseq, min=2, max=5, reflev=1, missing=., usemiss=0);
%indic3(vbl=vlfs_baseq, prefix=vlfs_baseq, min=2, max=5, reflev=1, missing=., usemiss=0);
%indic3(vbl=calor_baseq, prefix=calor_baseq, min=2, max=5, reflev=1, missing=., usemiss=0);
%indic3(vbl=ahei_baseq,  prefix=ahei_baseq,  min=2, max=5, reflev=1, missing=., usemiss=0);
%indic3(vbl=smk_baseline, prefix=smk_baseline, min=2, max=3, reflev=1, missing=., usemiss=0);
%indic3(vbl=alc_baseline, prefix=alc_baseline, min=2, max=6, reflev=1, missing=., usemiss=0);
%indic3(vbl=act_baseline, prefix=act_baseline, min=2, max=5, reflev=1, missing=., usemiss=0);
%indic3(vbl=bmi_baseline, prefix=bmi_baseline, min=2, max=9, reflev=1, missing=., usemiss=1);
run;

%include 'project/LowC_ObesityCancer/pool/allcancer/macro/res.sas';
/* Run Cox models with baseline LCD */
%res(analysis=baseline.adj1, catlin=T,   
     data=pool, event=can_obese, time=tobeca, 
     timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
     id=newid, 
     tvar=intervalp,
     agevar=agecon,
     strata=agecon cohort intervalp, 
     labels=F, 
     exp1=&lcs_baseq_,   adjust1= &calor_baseq_ &ahei_baseq_ &smk_baseline_ &alc_baseline_ &act_baseline_ &bmi_baseline_ white mv asp cafh mam_baseline nSES_baseline  endo2c_baseline   pmhever ocever,
     exp2=&alcs_baseq_,  adjust2= &calor_baseq_ &ahei_baseq_ &smk_baseline_ &alc_baseline_ &act_baseline_ &bmi_baseline_ white mv asp cafh mam_baseline nSES_baseline  endo2c_baseline   pmhever ocever,
     exp3=&vlcs_baseq_,  adjust3= &calor_baseq_ &ahei_baseq_ &smk_baseline_ &alc_baseline_ &act_baseline_ &bmi_baseline_ white mv asp cafh mam_baseline nSES_baseline  endo2c_baseline   pmhever ocever,
     exp4=&ulcs_baseq_,  adjust4= &calor_baseq_ &ahei_baseq_ &smk_baseline_ &alc_baseline_ &act_baseline_ &bmi_baseline_ white mv asp cafh mam_baseline nSES_baseline  endo2c_baseline   pmhever ocever,
     exp5=&hlcs_baseq_,  adjust5= &calor_baseq_ &ahei_baseq_ &smk_baseline_ &alc_baseline_ &act_baseline_ &bmi_baseline_ white mv asp cafh mam_baseline nSES_baseline  endo2c_baseline   pmhever ocever,
     exp6=&lfs_baseq_,   adjust6= &calor_baseq_ &ahei_baseq_ &smk_baseline_ &alc_baseline_ &act_baseline_ &bmi_baseline_ white mv asp cafh mam_baseline nSES_baseline  endo2c_baseline   pmhever ocever,
     exp7=&alfs_baseq_,  adjust7= &calor_baseq_ &ahei_baseq_ &smk_baseline_ &alc_baseline_ &act_baseline_ &bmi_baseline_ white mv asp cafh mam_baseline nSES_baseline  endo2c_baseline   pmhever ocever,
     exp8=&vlfs_baseq_,  adjust8= &calor_baseq_ &ahei_baseq_ &smk_baseline_ &alc_baseline_ &act_baseline_ &bmi_baseline_ white mv asp cafh mam_baseline nSES_baseline  endo2c_baseline   pmhever ocever,
     exp9=&ulfs_baseq_,  adjust9= &calor_baseq_ &ahei_baseq_ &smk_baseline_ &alc_baseline_ &act_baseline_ &bmi_baseline_ white mv asp cafh mam_baseline nSES_baseline  endo2c_baseline   pmhever ocever,
     exp10=&hlfs_baseq_, adjust10= &calor_baseq_ &ahei_baseq_ &smk_baseline_ &alc_baseline_ &act_baseline_ &bmi_baseline_ white mv asp cafh mam_baseline nSES_baseline  endo2c_baseline   pmhever ocever);
Run;

/* Output results */
%include 'project/LowC_ObesityCancer/pool/allcancer/macro/resout.sas';
%resout(outfile=table/tableS2_baseline, dbms=tab, opt=replace);         

/************************ INTERVAL-SPECIFIC LCD ANALYSIS *****************************************/

/* identify which cohorts contribute diet data at each interval */
data pool_with_diet_flag;
    set pool;
    
    /* Flag whether diet data exists for this person at this interval */
    has_diet = 0;
    
    /* NHS1: intervals 1,3,5,7,9,11,13 (1986,1990,1994,1998,2002,2006,2010) */
    if cohort=1 and intervalp in (1,3,5,7,9,11,13) then has_diet=1;
    
    /* NHS2: intervals 4,6,8,10,12,14,16,18 (1991, 1995, 1999, 2003, 2007, 2011, 2015, 2019)*/
    if cohort=2 and intervalp in (4,6,8,10,12,14,16,18) then has_diet=1;
    
    /* HPFS: intervals 1,3,5,7,9,11,13,15,17 (1986,1990,1994,1998,2002,2006,2010, 2014, 2018.)*/
    if cohort=3 and intervalp in (1,3,5,7,9,11,13,15,17) then has_diet=1;
run;


/* Intervals where all 3 cohorts contribute (approximately same calendar years) */
*interval_analysis(3, y90);   /* ~1990: NHS1-1990, NHS2-1991, HPFS-1990 */
*interval_analysis(5, y94);   /* ~1994: NHS1-1994, NHS2-1995, HPFS-1994 */
*interval_analysis(7, y98);   /* ~1998: NHS1-1998, NHS2-1999, HPFS-1998 */
*interval_analysis(9, y02);   /* ~2002: NHS1-2002, NHS2-2003, HPFS-2002 */
*interval_analysis(11, y06);  /* ~2006: NHS1-2006, NHS2-2007, HPFS-2006 */
*interval_analysis(13, y10);  /* ~2010: NHS1-2010, NHS2-2011, HPFS-2010 */

/* Intervals with 2 cohorts */
*interval_analysis(1, y86);   /* 1986: NHS1, HPFS only */
*interval_analysis(15, y14);  /* 2014: NHS2, HPFS only */
*interval_analysis(17, y18);  /* 2018: HPFS only */


/************************ INTERVAL 3 (1990) ANALYSIS *****************************************/

/* Extract diet and covariate values at interval 3 */
proc sort data=pool_with_diet_flag; by newid intervalp; run;

data interval_values_y90;
    set pool_with_diet_flag;
    where (cohort in (1,3) and intervalp=3) or (cohort=2 and intervalp=4);
    
    /* Store diet values from this interval */
    lcs_y90 = lcssvcon;
    alcs_y90 = alcssvcon;
    vlcs_y90 = vlcssvcon;
    ulcs_y90 = ulcssvcon;
    hlcs_y90 = hlcssvcon;
    lfs_y90 = lfssvcon;
    ulfs_y90 = ulfssvcon;
    hlfs_y90 = hlfssvcon;
    alfs_y90 = alfssvcon;
    vlfs_y90 = vlfssvcon;
    
    /* Store covariate values from this interval */
    calor_y90 = calorsvcon;
    ahei_y90 = aheimsvcon;
    smk_y90 = smkstatus;
    alc_y90 = alc;
    act_y90 = actc;
    bmi_y90 = bmib9;
    mam_y90       = mam;
    endo2c_y90    = endo2c;
    endo3c_y90    = endo3c;
    screen_y90 = screen_cancer  ;
    nSES_y90 = nSES  ;

    keep newid lcs_y90 alcs_y90 vlcs_y90 ulcs_y90 hlcs_y90 
         lfs_y90 ulfs_y90 hlfs_y90 alfs_y90 vlfs_y90
         calor_y90 ahei_y90 smk_y90 alc_y90 act_y90 bmi_y90
         mam_y90 endo2c_y90 nSES_y90  endo3c_y90 screen_y90  ;
run;

/* Merge these values back to all future intervals */
data analysis_y90;
    merge pool_with_diet_flag (in=dietflag) interval_values_y90 (in=intvalue);
    by newid;
    
    /* Keep only observations at or after interval 3 */
    /* AND only people who had diet data at interval 3 */
    if dietflag and intvalue and intervalp >= 3;
run;

/* Create quintiles based on the interval-specific values */
proc rank data=analysis_y90 groups=5 out=analysis_y90;
    var lcs_y90 alcs_y90 vlcs_y90 ulcs_y90 hlcs_y90 lfs_y90 ulfs_y90 hlfs_y90 alfs_y90 vlfs_y90 calor_y90 ahei_y90;
    ranks lcs_qy90 alcs_qy90 vlcs_qy90 ulcs_qy90 hlcs_qy90 lfs_qy90 ulfs_qy90 hlfs_qy90 alfs_qy90 vlfs_qy90 calor_qy90 ahei_qy90;
run;

/* Convert to 1-5 scale and create indicator variables */
data analysis_y90;
    set analysis_y90 end=_end_;
    
    /* Convert ranks 0-4 to quintiles 1-5 */
    array qvars lcs_qy90 alcs_qy90 vlcs_qy90 ulcs_qy90 hlcs_qy90 lfs_qy90 ulfs_qy90 hlfs_qy90 alfs_qy90 vlfs_qy90 calor_qy90 ahei_qy90;
    do over qvars;
        if qvars ne . then qvars = qvars + 1;
    end;
    
    /* Create indicator variables - INSIDE the DATA step */
    %indic3(vbl=lcs_qy90,   prefix=lcs_qy90,   min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=alcs_qy90,  prefix=alcs_qy90,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=vlcs_qy90,  prefix=vlcs_qy90,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=ulcs_qy90,  prefix=ulcs_qy90,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=hlcs_qy90,  prefix=hlcs_qy90,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=lfs_qy90,   prefix=lfs_qy90,   min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=ulfs_qy90,  prefix=ulfs_qy90,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=hlfs_qy90,  prefix=hlfs_qy90,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=alfs_qy90,  prefix=alfs_qy90,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=vlfs_qy90,  prefix=vlfs_qy90,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=calor_qy90, prefix=calor_qy90, min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=ahei_qy90,  prefix=ahei_qy90,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=smk_y90,    prefix=smk_y90,    min=2, max=3, reflev=1, missing=., usemiss=0);
    %indic3(vbl=alc_y90,    prefix=alc_y90,    min=2, max=6, reflev=1, missing=., usemiss=0);
    %indic3(vbl=act_y90,    prefix=act_y90,    min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=bmi_y90,    prefix=bmi_y90,    min=2, max=9, reflev=1, missing=., usemiss=1);
run;

/* Run Cox models */
%res(analysis=interval_y90_adj0, catlin=T,   
     data=analysis_y90, event=can_obese, time=tobeca, 
     timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
     id=newid, 
     tvar=intervalp,
     agevar=agecon,
     strata=agecon cohort intervalp, 
     labels=F, 
     exp1=&lcs_qy90_,   adjust1=,
     exp2=&alcs_qy90_,  adjust2=,
     exp3=&vlcs_qy90_,  adjust3=,
     exp4=&ulcs_qy90_,  adjust4=,
     exp5=&hlcs_qy90_,  adjust5=,
     exp6=&lfs_qy90_,   adjust6=,
     exp7=&alfs_qy90_,  adjust7=,
     exp8=&vlfs_qy90_,  adjust8=,
     exp9=&ulfs_qy90_,  adjust9=,
     exp10=&hlfs_qy90_, adjust10=);
Run;

%res(analysis=interval_y90_adj1, catlin=T,   
     data=analysis_y90, event=can_obese, time=tobeca, 
     timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
     id=newid, 
     tvar=intervalp,
     agevar=agecon,
     strata=agecon cohort intervalp, 
     labels=F, 
     exp1=&lcs_qy90_,   adjust1= &calor_qy90_ &ahei_qy90_ &smk_y90_ &alc_y90_ &act_y90_ &bmi_y90_ white mv asp cafh mam_y90 endo2c_y90 nSES_y90    pmhever ocever,
     exp2=&alcs_qy90_,  adjust2= &calor_qy90_ &ahei_qy90_ &smk_y90_ &alc_y90_ &act_y90_ &bmi_y90_ white mv asp cafh mam_y90 endo2c_y90 nSES_y90    pmhever ocever,
     exp3=&vlcs_qy90_,  adjust3= &calor_qy90_ &ahei_qy90_ &smk_y90_ &alc_y90_ &act_y90_ &bmi_y90_ white mv asp cafh mam_y90 endo2c_y90 nSES_y90    pmhever ocever,
     exp4=&ulcs_qy90_,  adjust4= &calor_qy90_ &ahei_qy90_ &smk_y90_ &alc_y90_ &act_y90_ &bmi_y90_ white mv asp cafh mam_y90 endo2c_y90 nSES_y90    pmhever ocever,
     exp5=&hlcs_qy90_,  adjust5= &calor_qy90_ &ahei_qy90_ &smk_y90_ &alc_y90_ &act_y90_ &bmi_y90_ white mv asp cafh mam_y90 endo2c_y90 nSES_y90    pmhever ocever,
     exp6=&lfs_qy90_,   adjust6= &calor_qy90_ &ahei_qy90_ &smk_y90_ &alc_y90_ &act_y90_ &bmi_y90_ white mv asp cafh mam_y90 endo2c_y90 nSES_y90    pmhever ocever,
     exp7=&alfs_qy90_,  adjust7= &calor_qy90_ &ahei_qy90_ &smk_y90_ &alc_y90_ &act_y90_ &bmi_y90_ white mv asp cafh mam_y90 endo2c_y90 nSES_y90    pmhever ocever,
     exp8=&vlfs_qy90_,  adjust8= &calor_qy90_ &ahei_qy90_ &smk_y90_ &alc_y90_ &act_y90_ &bmi_y90_ white mv asp cafh mam_y90 endo2c_y90 nSES_y90    pmhever ocever,
     exp9=&ulfs_qy90_,  adjust9= &calor_qy90_ &ahei_qy90_ &smk_y90_ &alc_y90_ &act_y90_ &bmi_y90_ white mv asp cafh mam_y90 endo2c_y90 nSES_y90    pmhever ocever,
     exp10=&hlfs_qy90_, adjust10= &calor_qy90_ &ahei_qy90_ &smk_y90_ &alc_y90_ &act_y90_ &bmi_y90_ white mv asp cafh mam_y90 endo2c_y90 nSES_y90    pmhever ocever);
Run;

/* Clean up */
proc datasets library=work nolist;
    delete interval_values_y90;
quit;


/************************ INTERVAL 5 (1994) ANALYSIS *****************************************/

data interval_values_y94;
    set pool_with_diet_flag;
    where (cohort in (1,3) and intervalp=5) or (cohort=2 and intervalp=6);
    
    /* Store diet values from this interval */
    lcs_y94 = lcssvcon;
    alcs_y94 = alcssvcon;
    vlcs_y94 = vlcssvcon;
    ulcs_y94 = ulcssvcon;
    hlcs_y94 = hlcssvcon;
    lfs_y94 = lfssvcon;
    ulfs_y94 = ulfssvcon;
    hlfs_y94 = hlfssvcon;
    alfs_y94 = alfssvcon;
    vlfs_y94 = vlfssvcon;
    
    /* Store covariate values from this interval */
    calor_y94 = calorsvcon;
    ahei_y94 = aheimsvcon;
    smk_y94 = smkstatus;
    alc_y94 = alc;
    act_y94 = actc;
    bmi_y94 = bmib9;
    mam_y94       = mam;
    endo2c_y94    = endo2c;
    endo3c_y94    = endo3c;
    screen_y94 = screen_cancer  ;
    nSES_y94 = nSES;

    keep newid lcs_y94 alcs_y94 vlcs_y94 ulcs_y94 hlcs_y94 
         lfs_y94 ulfs_y94 hlfs_y94 alfs_y94 vlfs_y94
         calor_y94 ahei_y94 smk_y94 alc_y94 act_y94 bmi_y94
         mam_y94 endo2c_y94 nSES_y94  endo3c_y94 screen_y94;
run;

/* Merge these values back to all future intervals */
data analysis_y94;
    merge pool_with_diet_flag (in=dietflag) interval_values_y94 (in=intvalue);
    by newid;
    
    /* Keep only observations at or after interval 5 */
    /* AND only people who had diet data at interval 5 */
    if dietflag and intvalue and intervalp >= 5;
run;

/* Create quintiles based on the interval-specific values */
proc rank data=analysis_y94 groups=5 out=analysis_y94;
    var lcs_y94 alcs_y94 vlcs_y94 ulcs_y94 hlcs_y94 lfs_y94 ulfs_y94 hlfs_y94 alfs_y94 vlfs_y94 calor_y94 ahei_y94;
    ranks lcs_qy94 alcs_qy94 vlcs_qy94 ulcs_qy94 hlcs_qy94 lfs_qy94 ulfs_qy94 hlfs_qy94 alfs_qy94 vlfs_qy94 calor_qy94 ahei_qy94;
run;

/* Convert to 1-5 scale and create indicator variables */
data analysis_y94;
    set analysis_y94 end=_end_;
    
    /* Convert ranks 0-4 to quintiles 1-5 */
    array qvars lcs_qy94 alcs_qy94 vlcs_qy94 ulcs_qy94 hlcs_qy94 lfs_qy94 ulfs_qy94 hlfs_qy94 alfs_qy94 vlfs_qy94 calor_qy94 ahei_qy94;
    do over qvars;
        if qvars ne . then qvars = qvars + 1;
    end;
    
    /* Create indicator variables - INSIDE the DATA step */
    %indic3(vbl=lcs_qy94,   prefix=lcs_qy94,   min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=alcs_qy94,  prefix=alcs_qy94,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=vlcs_qy94,  prefix=vlcs_qy94,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=ulcs_qy94,  prefix=ulcs_qy94,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=hlcs_qy94,  prefix=hlcs_qy94,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=lfs_qy94,   prefix=lfs_qy94,   min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=ulfs_qy94,  prefix=ulfs_qy94,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=hlfs_qy94,  prefix=hlfs_qy94,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=alfs_qy94,  prefix=alfs_qy94,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=vlfs_qy94,  prefix=vlfs_qy94,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=calor_qy94, prefix=calor_qy94, min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=ahei_qy94,  prefix=ahei_qy94,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=smk_y94,    prefix=smk_y94,    min=2, max=3, reflev=1, missing=., usemiss=0);
    %indic3(vbl=alc_y94,    prefix=alc_y94,    min=2, max=6, reflev=1, missing=., usemiss=0);
    %indic3(vbl=act_y94,    prefix=act_y94,    min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=bmi_y94,    prefix=bmi_y94,    min=2, max=9, reflev=1, missing=., usemiss=1);
run;

/* Run Cox models */
%res(analysis=interval_y94_adj0, catlin=T,   
     data=analysis_y94, event=can_obese, time=tobeca, 
     timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
     id=newid, 
     tvar=intervalp,
     agevar=agecon,
     strata=agecon cohort intervalp, 
     labels=F, 
     exp1=&lcs_qy94_,   adjust1=,
     exp2=&alcs_qy94_,  adjust2=,
     exp3=&vlcs_qy94_,  adjust3=,
     exp4=&ulcs_qy94_,  adjust4=,
     exp5=&hlcs_qy94_,  adjust5=,
     exp6=&lfs_qy94_,   adjust6=,
     exp7=&alfs_qy94_,  adjust7=,
     exp8=&vlfs_qy94_,  adjust8=,
     exp9=&ulfs_qy94_,  adjust9=,
     exp10=&hlfs_qy94_, adjust10=);
Run;

%res(analysis=interval_y94_adj1, catlin=T,   
     data=analysis_y94, event=can_obese, time=tobeca, 
     timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
     id=newid, 
     tvar=intervalp,
     agevar=agecon,
     strata=agecon cohort intervalp, 
     labels=F, 
     exp1=&lcs_qy94_,   adjust1= &calor_qy94_ &ahei_qy94_ &smk_y94_ &alc_y94_ &act_y94_ &bmi_y94_ white mv asp cafh mam_y94 endo2c_y94 nSES_y94   pmhever ocever,
     exp2=&alcs_qy94_,  adjust2= &calor_qy94_ &ahei_qy94_ &smk_y94_ &alc_y94_ &act_y94_ &bmi_y94_ white mv asp cafh mam_y94 endo2c_y94 nSES_y94   pmhever ocever,
     exp3=&vlcs_qy94_,  adjust3= &calor_qy94_ &ahei_qy94_ &smk_y94_ &alc_y94_ &act_y94_ &bmi_y94_ white mv asp cafh mam_y94 endo2c_y94 nSES_y94   pmhever ocever,
     exp4=&ulcs_qy94_,  adjust4= &calor_qy94_ &ahei_qy94_ &smk_y94_ &alc_y94_ &act_y94_ &bmi_y94_ white mv asp cafh mam_y94 endo2c_y94 nSES_y94   pmhever ocever,
     exp5=&hlcs_qy94_,  adjust5= &calor_qy94_ &ahei_qy94_ &smk_y94_ &alc_y94_ &act_y94_ &bmi_y94_ white mv asp cafh mam_y94 endo2c_y94 nSES_y94   pmhever ocever,
     exp6=&lfs_qy94_,   adjust6= &calor_qy94_ &ahei_qy94_ &smk_y94_ &alc_y94_ &act_y94_ &bmi_y94_ white mv asp cafh mam_y94 endo2c_y94 nSES_y94   pmhever ocever,
     exp7=&alfs_qy94_,  adjust7= &calor_qy94_ &ahei_qy94_ &smk_y94_ &alc_y94_ &act_y94_ &bmi_y94_ white mv asp cafh mam_y94 endo2c_y94 nSES_y94   pmhever ocever,
     exp8=&vlfs_qy94_,  adjust8= &calor_qy94_ &ahei_qy94_ &smk_y94_ &alc_y94_ &act_y94_ &bmi_y94_ white mv asp cafh mam_y94 endo2c_y94 nSES_y94   pmhever ocever,
     exp9=&ulfs_qy94_,  adjust9= &calor_qy94_ &ahei_qy94_ &smk_y94_ &alc_y94_ &act_y94_ &bmi_y94_ white mv asp cafh mam_y94 endo2c_y94 nSES_y94   pmhever ocever,
     exp10=&hlfs_qy94_, adjust10= &calor_qy94_ &ahei_qy94_ &smk_y94_ &alc_y94_ &act_y94_ &bmi_y94_ white mv asp cafh mam_y94 endo2c_y94 nSES_y94   pmhever ocever);
Run;

/* Clean up */
proc datasets library=work nolist;
    delete interval_values_y94;
quit;



/************************ INTERVAL 7 (1998) ANALYSIS *****************************************/

/* Extract diet and covariate values at interval 7 */
proc sort data=pool_with_diet_flag; by newid intervalp; run;

data interval_values_y98;
    set pool_with_diet_flag;
    where (cohort in (1,3) and intervalp=7) or (cohort=2 and intervalp=8);
    
    /* Store diet values from this interval */
    lcs_y98 = lcssvcon;
    alcs_y98 = alcssvcon;
    vlcs_y98 = vlcssvcon;
    ulcs_y98 = ulcssvcon;
    hlcs_y98 = hlcssvcon;
    lfs_y98 = lfssvcon;
    ulfs_y98 = ulfssvcon;
    hlfs_y98 = hlfssvcon;
    vlfs_y98 = vlfssvcon;
    alfs_y98 = alfssvcon;
    
    /* Store covariate values from this interval */
    calor_y98 = calorsvcon;
    ahei_y98 = aheimsvcon;
    smk_y98 = smkstatus;
    alc_y98 = alc;
    act_y98 = actc;
    bmi_y98 = bmib9;
    mam_y98       = mam;
    endo2c_y98    = endo2c;
    endo3c_y98    = endo3c;
    screen_y98 = screen_cancer  ;
    nSES_y98       = nSES;

    keep newid lcs_y98 alcs_y98 vlcs_y98 ulcs_y98 hlcs_y98 
         lfs_y98 ulfs_y98 hlfs_y98 alfs_y98 vlfs_y98
         calor_y98 ahei_y98 smk_y98 alc_y98 act_y98 bmi_y98
         mam_y98 endo2c_y98 nSES_y98  endo3c_y98 screen_y98 ;
run;

/* Merge these values back to all future intervals */
data analysis_y98;
    merge pool_with_diet_flag (in=dietflag) interval_values_y98 (in=intvalue);
    by newid;
    
    /* Keep only observations at or after interval 7 */
    /* AND only people who had diet data at interval 7 */
    if dietflag and intvalue and intervalp >= 7;
run;

/* Create quintiles based on the interval-specific values */
proc rank data=analysis_y98 groups=5 out=analysis_y98;
    var lcs_y98 alcs_y98 vlcs_y98 ulcs_y98 hlcs_y98 lfs_y98 ulfs_y98 hlfs_y98 alfs_y98 vlfs_y98 calor_y98 ahei_y98;
    ranks lcs_qy98 alcs_qy98 vlcs_qy98 ulcs_qy98 hlcs_qy98 lfs_qy98 ulfs_qy98 hlfs_qy98 alfs_qy98 vlfs_qy98 calor_qy98 ahei_qy98;
run;

/* Convert to 1-5 scale and create indicator variables */
data analysis_y98;
    set analysis_y98 end=_end_;
    
    /* Convert ranks 0-4 to quintiles 1-5 */
    array qvars lcs_qy98 alcs_qy98 vlcs_qy98 ulcs_qy98 hlcs_qy98 lfs_qy98 ulfs_qy98 hlfs_qy98 alfs_qy98 vlfs_qy98 calor_qy98 ahei_qy98;
    do over qvars;
        if qvars ne . then qvars = qvars + 1;
    end;
    
    /* Create indicator variables - INSIDE the DATA step */
    %indic3(vbl=lcs_qy98,   prefix=lcs_qy98,   min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=alcs_qy98,  prefix=alcs_qy98,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=vlcs_qy98,  prefix=vlcs_qy98,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=ulcs_qy98,  prefix=ulcs_qy98,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=hlcs_qy98,  prefix=hlcs_qy98,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=lfs_qy98,   prefix=lfs_qy98,   min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=ulfs_qy98,  prefix=ulfs_qy98,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=hlfs_qy98,  prefix=hlfs_qy98,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=alfs_qy98,  prefix=alfs_qy98,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=vlfs_qy98,  prefix=vlfs_qy98,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=calor_qy98, prefix=calor_qy98, min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=ahei_qy98,  prefix=ahei_qy98,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=smk_y98,    prefix=smk_y98,    min=2, max=3, reflev=1, missing=., usemiss=0);
    %indic3(vbl=alc_y98,    prefix=alc_y98,    min=2, max=6, reflev=1, missing=., usemiss=0);
    %indic3(vbl=act_y98,    prefix=act_y98,    min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=bmi_y98,    prefix=bmi_y98,    min=2, max=9, reflev=1, missing=., usemiss=1);
run;

/* Run Cox models */
%res(analysis=interval_y98_adj0, catlin=T,   
     data=analysis_y98, event=can_obese, time=tobeca, 
     timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
     id=newid, 
     tvar=intervalp,
     agevar=agecon,
     strata=agecon cohort intervalp, 
     labels=F, 
     exp1=&lcs_qy98_,   adjust1=,
     exp2=&alcs_qy98_,  adjust2=,
     exp3=&vlcs_qy98_,  adjust3=,
     exp4=&ulcs_qy98_,  adjust4=,
     exp5=&hlcs_qy98_,  adjust5=,
     exp6=&lfs_qy98_,   adjust6=,
     exp7=&alfs_qy98_,  adjust7=,
     exp8=&vlfs_qy98_,  adjust8=,
     exp9=&ulfs_qy98_,  adjust9=,
     exp10=&hlfs_qy98_, adjust10=);
Run;

%res(analysis=interval_y98_adj1, catlin=T,   
     data=analysis_y98, event=can_obese, time=tobeca, 
     timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
     id=newid, 
     tvar=intervalp,
     agevar=agecon,
     strata=agecon cohort intervalp, 
     labels=F, 
     exp1=&lcs_qy98_,   adjust1= &calor_qy98_ &ahei_qy98_ &smk_y98_ &alc_y98_ &act_y98_ &bmi_y98_ white mv asp  cafh mam_y98 endo2c_y98 nSES_y98    pmhever ocever,
     exp2=&alcs_qy98_,  adjust2= &calor_qy98_ &ahei_qy98_ &smk_y98_ &alc_y98_ &act_y98_ &bmi_y98_ white mv asp  cafh mam_y98 endo2c_y98 nSES_y98    pmhever ocever,
     exp3=&vlcs_qy98_,  adjust3= &calor_qy98_ &ahei_qy98_ &smk_y98_ &alc_y98_ &act_y98_ &bmi_y98_ white mv asp  cafh mam_y98 endo2c_y98 nSES_y98    pmhever ocever,
     exp4=&ulcs_qy98_,  adjust4= &calor_qy98_ &ahei_qy98_ &smk_y98_ &alc_y98_ &act_y98_ &bmi_y98_ white mv asp  cafh mam_y98 endo2c_y98 nSES_y98    pmhever ocever,
     exp5=&hlcs_qy98_,  adjust5= &calor_qy98_ &ahei_qy98_ &smk_y98_ &alc_y98_ &act_y98_ &bmi_y98_ white mv asp  cafh mam_y98 endo2c_y98 nSES_y98    pmhever ocever,
     exp6=&lfs_qy98_,   adjust6= &calor_qy98_ &ahei_qy98_ &smk_y98_ &alc_y98_ &act_y98_ &bmi_y98_ white mv asp  cafh mam_y98 endo2c_y98 nSES_y98    pmhever ocever,
     exp7=&alfs_qy98_,  adjust7= &calor_qy98_ &ahei_qy98_ &smk_y98_ &alc_y98_ &act_y98_ &bmi_y98_ white mv asp  cafh mam_y98 endo2c_y98 nSES_y98    pmhever ocever,
     exp8=&vlfs_qy98_,  adjust8= &calor_qy98_ &ahei_qy98_ &smk_y98_ &alc_y98_ &act_y98_ &bmi_y98_ white mv asp  cafh mam_y98 endo2c_y98 nSES_y98    pmhever ocever,
     exp9=&ulfs_qy98_,  adjust9= &calor_qy98_ &ahei_qy98_ &smk_y98_ &alc_y98_ &act_y98_ &bmi_y98_ white mv asp  cafh mam_y98 endo2c_y98 nSES_y98    pmhever ocever,
     exp10=&hlfs_qy98_, adjust10= &calor_qy98_ &ahei_qy98_ &smk_y98_ &alc_y98_ &act_y98_ &bmi_y98_ white mv asp  cafh mam_y98 endo2c_y98 nSES_y98    pmhever ocever);
Run;

/* Clean up */
proc datasets library=work nolist;
    delete interval_values_y98;
quit;



/************************ INTERVAL 9 (2002) ANALYSIS *****************************************/

/* Extract diet and covariate values at interval 9 */
proc sort data=pool_with_diet_flag; by newid intervalp; run;

data interval_values_y02;
    set pool_with_diet_flag;
    where (cohort in (1,3) and intervalp=9) or (cohort=2 and intervalp=10);
    
    /* Store diet values from this interval */
    lcs_y02 = lcssvcon;
    alcs_y02 = alcssvcon;
    vlcs_y02 = vlcssvcon;
    ulcs_y02 = ulcssvcon;
    hlcs_y02 = hlcssvcon;
    lfs_y02 = lfssvcon;
    ulfs_y02 = ulfssvcon;
    hlfs_y02 = hlfssvcon;
    vlfs_y02 = vlfssvcon;
    alfs_y02 = alfssvcon;
    
    /* Store covariate values from this interval */
    calor_y02 = calorsvcon;
    ahei_y02 = aheimsvcon;
    smk_y02 = smkstatus;
    alc_y02 = alc;
    act_y02 = actc;
    bmi_y02 = bmib9;
    mam_y02       = mam;
    endo2c_y02    = endo2c;
    endo3c_y02    = endo3c;
    screen_y02 = screen_cancer  ;
    nSES_y02 = nSES  ;

    keep newid lcs_y02 alcs_y02 vlcs_y02 ulcs_y02 hlcs_y02 
         lfs_y02 ulfs_y02 hlfs_y02 alfs_y02 vlfs_y02
         calor_y02 ahei_y02 smk_y02 alc_y02 act_y02 bmi_y02
         mam_y02 endo2c_y02 nSES_y02  endo3c_y02 screen_y02;
run;

/* Merge these values back to all future intervals */
data analysis_y02;
    merge pool_with_diet_flag (in=dietflag) interval_values_y02 (in=intvalue);
    by newid;
    
    /* Keep only observations at or after interval 9 */
    /* AND only people who had diet data at interval 9 */
    if dietflag and intvalue and intervalp >= 9;
run;

/* Create quintiles based on the interval-specific values */
proc rank data=analysis_y02 groups=5 out=analysis_y02;
    var lcs_y02 alcs_y02 vlcs_y02 ulcs_y02 hlcs_y02 lfs_y02 ulfs_y02 hlfs_y02 alfs_y02 vlfs_y02 calor_y02 ahei_y02;
    ranks lcs_qy02 alcs_qy02 vlcs_qy02 ulcs_qy02 hlcs_qy02 lfs_qy02 ulfs_qy02 hlfs_qy02 alfs_qy02 vlfs_qy02 calor_qy02 ahei_qy02;
run;

/* Convert to 1-5 scale and create indicator variables */
data analysis_y02;
    set analysis_y02 end=_end_;
    
    /* Convert ranks 0-4 to quintiles 1-5 */
    array qvars lcs_qy02 alcs_qy02 vlcs_qy02 ulcs_qy02 hlcs_qy02 lfs_qy02 ulfs_qy02 hlfs_qy02 alfs_qy02 vlfs_qy02 calor_qy02 ahei_qy02;
    do over qvars;
        if qvars ne . then qvars = qvars + 1;
    end;
    
    /* Create indicator variables - INSIDE the DATA step */
    %indic3(vbl=lcs_qy02,   prefix=lcs_qy02,   min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=alcs_qy02,  prefix=alcs_qy02,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=vlcs_qy02,  prefix=vlcs_qy02,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=ulcs_qy02,  prefix=ulcs_qy02,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=hlcs_qy02,  prefix=hlcs_qy02,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=lfs_qy02,   prefix=lfs_qy02,   min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=ulfs_qy02,  prefix=ulfs_qy02,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=hlfs_qy02,  prefix=hlfs_qy02,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=alfs_qy02,  prefix=alfs_qy02,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=vlfs_qy02,  prefix=vlfs_qy02,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=calor_qy02, prefix=calor_qy02, min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=ahei_qy02,  prefix=ahei_qy02,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=smk_y02,    prefix=smk_y02,    min=2, max=3, reflev=1, missing=., usemiss=0);
    %indic3(vbl=alc_y02,    prefix=alc_y02,    min=2, max=6, reflev=1, missing=., usemiss=0);
    %indic3(vbl=act_y02,    prefix=act_y02,    min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=bmi_y02,    prefix=bmi_y02,    min=2, max=9, reflev=1, missing=., usemiss=1);
run;

/* Run Cox models */
%res(analysis=interval_y02_adj0, catlin=T,   
     data=analysis_y02, event=can_obese, time=tobeca, 
     timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
     id=newid, 
     tvar=intervalp,
     agevar=agecon,
     strata=agecon cohort intervalp, 
     labels=F, 
     exp1=&lcs_qy02_,   adjust1=,
     exp2=&alcs_qy02_,  adjust2=,
     exp3=&vlcs_qy02_,  adjust3=,
     exp4=&ulcs_qy02_,  adjust4=,
     exp5=&hlcs_qy02_,  adjust5=,
     exp6=&lfs_qy02_,   adjust6=,
     exp7=&alfs_qy02_,  adjust7=,
     exp8=&vlfs_qy02_,  adjust8=,
     exp9=&ulfs_qy02_,  adjust9=,
     exp10=&hlfs_qy02_, adjust10=);
Run;

%res(analysis=interval_y02_adj1, catlin=T,   
     data=analysis_y02, event=can_obese, time=tobeca, 
     timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
     id=newid, 
     tvar=intervalp,
     agevar=agecon,
     strata=agecon cohort intervalp, 
     labels=F, 
     exp1=&lcs_qy02_,   adjust1= &calor_qy02_ &ahei_qy02_ &smk_y02_ &alc_y02_ &act_y02_ &bmi_y02_ white mv asp cafh mam_y02 endo2c_y02 nSES_y02    pmhever ocever,
     exp2=&alcs_qy02_,  adjust2= &calor_qy02_ &ahei_qy02_ &smk_y02_ &alc_y02_ &act_y02_ &bmi_y02_ white mv asp cafh mam_y02 endo2c_y02 nSES_y02    pmhever ocever,
     exp3=&vlcs_qy02_,  adjust3= &calor_qy02_ &ahei_qy02_ &smk_y02_ &alc_y02_ &act_y02_ &bmi_y02_ white mv asp cafh mam_y02 endo2c_y02 nSES_y02    pmhever ocever,
     exp4=&ulcs_qy02_,  adjust4= &calor_qy02_ &ahei_qy02_ &smk_y02_ &alc_y02_ &act_y02_ &bmi_y02_ white mv asp cafh mam_y02 endo2c_y02 nSES_y02    pmhever ocever,
     exp5=&hlcs_qy02_,  adjust5= &calor_qy02_ &ahei_qy02_ &smk_y02_ &alc_y02_ &act_y02_ &bmi_y02_ white mv asp cafh mam_y02 endo2c_y02 nSES_y02    pmhever ocever,
     exp6=&lfs_qy02_,   adjust6= &calor_qy02_ &ahei_qy02_ &smk_y02_ &alc_y02_ &act_y02_ &bmi_y02_ white mv asp cafh mam_y02 endo2c_y02 nSES_y02    pmhever ocever,
     exp7=&alfs_qy02_,  adjust7= &calor_qy02_ &ahei_qy02_ &smk_y02_ &alc_y02_ &act_y02_ &bmi_y02_ white mv asp cafh mam_y02 endo2c_y02 nSES_y02    pmhever ocever,
     exp8=&vlfs_qy02_,  adjust8= &calor_qy02_ &ahei_qy02_ &smk_y02_ &alc_y02_ &act_y02_ &bmi_y02_ white mv asp cafh mam_y02 endo2c_y02 nSES_y02    pmhever ocever,
     exp9=&ulfs_qy02_,  adjust9= &calor_qy02_ &ahei_qy02_ &smk_y02_ &alc_y02_ &act_y02_ &bmi_y02_ white mv asp cafh mam_y02 endo2c_y02 nSES_y02    pmhever ocever,
     exp10=&hlfs_qy02_, adjust10= &calor_qy02_ &ahei_qy02_ &smk_y02_ &alc_y02_ &act_y02_ &bmi_y02_ white mv asp cafh mam_y02 endo2c_y02 nSES_y02    pmhever ocever);
Run;

/* Clean up */
proc datasets library=work nolist;
    delete interval_values_y02;
quit;



/************************ INTERVAL 11 (2006) ANALYSIS *****************************************/

/* Extract diet and covariate values at interval 7 */
proc sort data=pool_with_diet_flag; by newid intervalp; run;

data interval_values_y06;
    set pool_with_diet_flag;
    where (cohort in (1,3) and intervalp=11) or (cohort=2 and intervalp=12);
    
    /* Store diet values from this interval */
    lcs_y06 = lcssvcon;
    alcs_y06 = alcssvcon;
    vlcs_y06 = vlcssvcon;
    ulcs_y06 = ulcssvcon;
    hlcs_y06 = hlcssvcon;
    lfs_y06 = lfssvcon;
    ulfs_y06 = ulfssvcon;
    hlfs_y06 = hlfssvcon;
    vlfs_y06 = vlfssvcon;
    alfs_y06 = alfssvcon;
    
    /* Store covariate values from this interval */
    calor_y06 = calorsvcon;
    ahei_y06 = aheimsvcon;
    smk_y06 = smkstatus;
    alc_y06 = alc;
    act_y06 = actc;
    bmi_y06 = bmib9;
    mam_y06       = mam;
    endo2c_y06    = endo2c;
    endo3c_y06    = endo3c;
    screen_y06 = screen_cancer  ;
    nSES_y06 = nSES  ;

    keep newid lcs_y06 alcs_y06 vlcs_y06 ulcs_y06 hlcs_y06 
         lfs_y06 ulfs_y06 hlfs_y06 alfs_y06 vlfs_y06
         calor_y06 ahei_y06 smk_y06 alc_y06 act_y06 bmi_y06
         mam_y06 endo2c_y06 nSES_y06  endo3c_y06 screen_y06 ;
run;

/* Merge these values back to all future intervals */
data analysis_y06;
    merge pool_with_diet_flag (in=dietflag) interval_values_y06 (in=intvalue);
    by newid;
    
    /* Keep only observations at or after interval 11 */
    /* AND only people who had diet data at interval 11 */
    if dietflag and intvalue and intervalp >= 11;
run;

/* Create quintiles based on the interval-specific values */
proc rank data=analysis_y06 groups=5 out=analysis_y06;
    var lcs_y06 alcs_y06 vlcs_y06 ulcs_y06 hlcs_y06 lfs_y06 ulfs_y06 hlfs_y06 alfs_y06 vlfs_y06 calor_y06 ahei_y06;
    ranks lcs_qy06 alcs_qy06 vlcs_qy06 ulcs_qy06 hlcs_qy06 lfs_qy06 ulfs_qy06 hlfs_qy06 alfs_qy06 vlfs_qy06 calor_qy06 ahei_qy06;
run;

/* Convert to 1-5 scale and create indicator variables */
data analysis_y06;
    set analysis_y06 end=_end_;
    
    /* Convert ranks 0-4 to quintiles 1-5 */
    array qvars lcs_qy06 alcs_qy06 vlcs_qy06 ulcs_qy06 hlcs_qy06 lfs_qy06 ulfs_qy06 hlfs_qy06 alfs_qy06 vlfs_qy06 calor_qy06 ahei_qy06;
    do over qvars;
        if qvars ne . then qvars = qvars + 1;
    end;
    
    /* Create indicator variables - INSIDE the DATA step */
    %indic3(vbl=lcs_qy06,   prefix=lcs_qy06,   min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=alcs_qy06,  prefix=alcs_qy06,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=vlcs_qy06,  prefix=vlcs_qy06,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=ulcs_qy06,  prefix=ulcs_qy06,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=hlcs_qy06,  prefix=hlcs_qy06,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=lfs_qy06,   prefix=lfs_qy06,   min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=ulfs_qy06,  prefix=ulfs_qy06,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=hlfs_qy06,  prefix=hlfs_qy06,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=alfs_qy06,  prefix=alfs_qy06,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=vlfs_qy06,  prefix=vlfs_qy06,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=calor_qy06, prefix=calor_qy06, min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=ahei_qy06,  prefix=ahei_qy06,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=smk_y06,    prefix=smk_y06,    min=2, max=3, reflev=1, missing=., usemiss=0);
    %indic3(vbl=alc_y06,    prefix=alc_y06,    min=2, max=6, reflev=1, missing=., usemiss=0);
    %indic3(vbl=act_y06,    prefix=act_y06,    min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=bmi_y06,    prefix=bmi_y06,    min=2, max=9, reflev=1, missing=., usemiss=1);
run;

/* Run Cox models */
%res(analysis=interval_y06_adj0, catlin=T,   
     data=analysis_y06, event=can_obese, time=tobeca, 
     timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
     id=newid, 
     tvar=intervalp,
     agevar=agecon,
     strata=agecon cohort intervalp, 
     labels=F, 
     exp1=&lcs_qy06_,   adjust1=,
     exp2=&alcs_qy06_,  adjust2=,
     exp3=&vlcs_qy06_,  adjust3=,
     exp4=&ulcs_qy06_,  adjust4=,
     exp5=&hlcs_qy06_,  adjust5=,
     exp6=&lfs_qy06_,   adjust6=,
     exp7=&alfs_qy06_,  adjust7=,
     exp8=&vlfs_qy06_,  adjust8=,
     exp9=&ulfs_qy06_,  adjust9=,
     exp10=&hlfs_qy06_, adjust10=);
Run;

%res(analysis=interval_y06_adj1, catlin=T,   
     data=analysis_y06, event=can_obese, time=tobeca, 
     timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
     id=newid, 
     tvar=intervalp,
     agevar=agecon,
     strata=agecon cohort intervalp, 
     labels=F, 
     exp1=&lcs_qy06_,   adjust1= &calor_qy06_ &ahei_qy06_ &smk_y06_ &alc_y06_ &act_y06_ &bmi_y06_ white mv asp cafh mam_y06 endo2c_y06 nSES_y06    pmhever ocever,
     exp2=&alcs_qy06_,  adjust2= &calor_qy06_ &ahei_qy06_ &smk_y06_ &alc_y06_ &act_y06_ &bmi_y06_ white mv asp cafh mam_y06 endo2c_y06 nSES_y06    pmhever ocever,
     exp3=&vlcs_qy06_,  adjust3= &calor_qy06_ &ahei_qy06_ &smk_y06_ &alc_y06_ &act_y06_ &bmi_y06_ white mv asp cafh mam_y06 endo2c_y06 nSES_y06    pmhever ocever,
     exp4=&ulcs_qy06_,  adjust4= &calor_qy06_ &ahei_qy06_ &smk_y06_ &alc_y06_ &act_y06_ &bmi_y06_ white mv asp cafh mam_y06 endo2c_y06 nSES_y06    pmhever ocever,
     exp5=&hlcs_qy06_,  adjust5= &calor_qy06_ &ahei_qy06_ &smk_y06_ &alc_y06_ &act_y06_ &bmi_y06_ white mv asp cafh mam_y06 endo2c_y06 nSES_y06    pmhever ocever,
     exp6=&lfs_qy06_,   adjust6= &calor_qy06_ &ahei_qy06_ &smk_y06_ &alc_y06_ &act_y06_ &bmi_y06_ white mv asp cafh mam_y06 endo2c_y06 nSES_y06    pmhever ocever,
     exp7=&alfs_qy06_,  adjust7= &calor_qy06_ &ahei_qy06_ &smk_y06_ &alc_y06_ &act_y06_ &bmi_y06_ white mv asp cafh mam_y06 endo2c_y06 nSES_y06    pmhever ocever,
     exp8=&vlfs_qy06_,  adjust8= &calor_qy06_ &ahei_qy06_ &smk_y06_ &alc_y06_ &act_y06_ &bmi_y06_ white mv asp cafh mam_y06 endo2c_y06 nSES_y06    pmhever ocever,
     exp9=&ulfs_qy06_,  adjust9= &calor_qy06_ &ahei_qy06_ &smk_y06_ &alc_y06_ &act_y06_ &bmi_y06_ white mv asp cafh mam_y06 endo2c_y06 nSES_y06    pmhever ocever,
     exp10=&hlfs_qy06_, adjust10= &calor_qy06_ &ahei_qy06_ &smk_y06_ &alc_y06_ &act_y06_ &bmi_y06_ white mv asp cafh mam_y06 endo2c_y06 nSES_y06    pmhever ocever);
Run;

/* Clean up */
proc datasets library=work nolist;
    delete interval_values_y06;
quit;





/************************ INTERVAL 13 (2010) ANALYSIS *****************************************/

/* Extract diet and covariate values at interval 13 */
proc sort data=pool_with_diet_flag; by newid intervalp; run;

data interval_values_y10;
    set pool_with_diet_flag;
    where (cohort in (1,3) and intervalp=13) or (cohort=2 and intervalp=14);
    
    /* Store diet values from this interval */
    lcs_y10 = lcssvcon;
    alcs_y10 = alcssvcon;
    vlcs_y10 = vlcssvcon;
    ulcs_y10 = ulcssvcon;
    hlcs_y10 = hlcssvcon;
    lfs_y10 = lfssvcon;
    ulfs_y10 = ulfssvcon;
    hlfs_y10 = hlfssvcon;
    alfs_y10 = alfssvcon;
    vlfs_y10 = vlfssvcon;
    
    /* Store covariate values from this interval */
    calor_y10 = calorsvcon;
    ahei_y10 = aheimsvcon;
    smk_y10 = smkstatus;
    alc_y10 = alc;
    act_y10 = actc;
    bmi_y10 = bmib9;
    mam_y10       = mam;
    endo2c_y10    = endo2c;
    endo3c_y10    = endo3c;
    screen_y10 = screen_cancer  ;
    nSES_y10 = nSES  ;

    keep newid lcs_y10 alcs_y10 vlcs_y10 ulcs_y10 hlcs_y10 
         lfs_y10 ulfs_y10 hlfs_y10 alfs_y10 vlfs_y10
         calor_y10 ahei_y10 smk_y10 alc_y10 act_y10 bmi_y10
         mam_y10 endo2c_y10 nSES_y10  endo3c_y10 screen_y10 ;
run;

/* Merge these values back to all future intervals */
data analysis_y10;
    merge pool_with_diet_flag (in=dietflag) interval_values_y10 (in=intvalue);
    by newid;
    
    /* Keep only observations at or after interval 13 */
    /* AND only people who had diet data at interval 13 */
    if dietflag and intvalue and intervalp >= 13;
run;

/* Create quintiles based on the interval-specific values */
proc rank data=analysis_y10 groups=5 out=analysis_y10;
    var lcs_y10 alcs_y10 vlcs_y10 ulcs_y10 hlcs_y10 lfs_y10 ulfs_y10 hlfs_y10 alfs_y10 vlfs_y10 calor_y10 ahei_y10;
    ranks lcs_qy10 alcs_qy10 vlcs_qy10 ulcs_qy10 hlcs_qy10 lfs_qy10 ulfs_qy10 hlfs_qy10 alfs_qy10 vlfs_qy10 calor_qy10 ahei_qy10;
run;

/* Convert to 1-5 scale and create indicator variables */
data analysis_y10;
    set analysis_y10 end=_end_;
    
    /* Convert ranks 0-4 to quintiles 1-5 */
    array qvars lcs_qy10 alcs_qy10 vlcs_qy10 ulcs_qy10 hlcs_qy10 lfs_qy10 ulfs_qy10 hlfs_qy10 alfs_qy10 vlfs_qy10 calor_qy10 ahei_qy10;
    do over qvars;
        if qvars ne . then qvars = qvars + 1;
    end;
    
    /* Create indicator variables - INSIDE the DATA step */
    %indic3(vbl=lcs_qy10,   prefix=lcs_qy10,   min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=alcs_qy10,  prefix=alcs_qy10,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=vlcs_qy10,  prefix=vlcs_qy10,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=ulcs_qy10,  prefix=ulcs_qy10,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=hlcs_qy10,  prefix=hlcs_qy10,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=lfs_qy10,   prefix=lfs_qy10,   min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=ulfs_qy10,  prefix=ulfs_qy10,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=hlfs_qy10,  prefix=hlfs_qy10,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=alfs_qy10,  prefix=alfs_qy10,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=vlfs_qy10,  prefix=vlfs_qy10,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=calor_qy10, prefix=calor_qy10, min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=ahei_qy10,  prefix=ahei_qy10,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=smk_y10,    prefix=smk_y10,    min=2, max=3, reflev=1, missing=., usemiss=0);
    %indic3(vbl=alc_y10,    prefix=alc_y10,    min=2, max=6, reflev=1, missing=., usemiss=0);
    %indic3(vbl=act_y10,    prefix=act_y10,    min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=bmi_y10,    prefix=bmi_y10,    min=2, max=9, reflev=1, missing=., usemiss=1);
run;

/* Run Cox models */
%res(analysis=interval_y10_adj0, catlin=T,   
     data=analysis_y10, event=can_obese, time=tobeca, 
     timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
     id=newid, 
     tvar=intervalp,
     agevar=agecon,
     strata=agecon cohort intervalp, 
     labels=F, 
     exp1=&lcs_qy10_,   adjust1=,
     exp2=&alcs_qy10_,  adjust2=,
     exp3=&vlcs_qy10_,  adjust3=,
     exp4=&ulcs_qy10_,  adjust4=,
     exp5=&hlcs_qy10_,  adjust5=,
     exp6=&lfs_qy10_,   adjust6=,
     exp7=&alfs_qy10_,  adjust7=,
     exp8=&vlfs_qy10_,  adjust8=,
     exp9=&ulfs_qy10_,  adjust9=,
     exp10=&hlfs_qy10_, adjust10=);
Run;

%res(analysis=interval_y10_adj1, catlin=T,   
     data=analysis_y10, event=can_obese, time=tobeca, 
     timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
     id=newid, 
     tvar=intervalp,
     agevar=agecon,
     strata=agecon cohort intervalp, 
     labels=F, 
     exp1=&lcs_qy10_,   adjust1= &calor_qy10_ &ahei_qy10_ &smk_y10_ &alc_y10_ &act_y10_ &bmi_y10_ white mv asp cafh mam_y10 endo2c_y10 nSES_y10    pmhever ocever,
     exp2=&alcs_qy10_,  adjust2= &calor_qy10_ &ahei_qy10_ &smk_y10_ &alc_y10_ &act_y10_ &bmi_y10_ white mv asp cafh mam_y10 endo2c_y10 nSES_y10    pmhever ocever,
     exp3=&vlcs_qy10_,  adjust3= &calor_qy10_ &ahei_qy10_ &smk_y10_ &alc_y10_ &act_y10_ &bmi_y10_ white mv asp cafh mam_y10 endo2c_y10 nSES_y10    pmhever ocever,
     exp4=&ulcs_qy10_,  adjust4= &calor_qy10_ &ahei_qy10_ &smk_y10_ &alc_y10_ &act_y10_ &bmi_y10_ white mv asp cafh mam_y10 endo2c_y10 nSES_y10    pmhever ocever,
     exp5=&hlcs_qy10_,  adjust5= &calor_qy10_ &ahei_qy10_ &smk_y10_ &alc_y10_ &act_y10_ &bmi_y10_ white mv asp cafh mam_y10 endo2c_y10 nSES_y10    pmhever ocever,
     exp6=&lfs_qy10_,   adjust6= &calor_qy10_ &ahei_qy10_ &smk_y10_ &alc_y10_ &act_y10_ &bmi_y10_ white mv asp cafh mam_y10 endo2c_y10 nSES_y10    pmhever ocever,
     exp7=&alfs_qy10_,  adjust7= &calor_qy10_ &ahei_qy10_ &smk_y10_ &alc_y10_ &act_y10_ &bmi_y10_ white mv asp cafh mam_y10 endo2c_y10 nSES_y10    pmhever ocever,
     exp8=&vlfs_qy10_,  adjust8= &calor_qy10_ &ahei_qy10_ &smk_y10_ &alc_y10_ &act_y10_ &bmi_y10_ white mv asp cafh mam_y10 endo2c_y10 nSES_y10    pmhever ocever,
     exp9=&ulfs_qy10_,  adjust9= &calor_qy10_ &ahei_qy10_ &smk_y10_ &alc_y10_ &act_y10_ &bmi_y10_ white mv asp cafh mam_y10 endo2c_y10 nSES_y10    pmhever ocever,
     exp10=&hlfs_qy10_, adjust10= &calor_qy10_ &ahei_qy10_ &smk_y10_ &alc_y10_ &act_y10_ &bmi_y10_ white mv asp cafh mam_y10 endo2c_y10 nSES_y10    pmhever ocever);
Run;

/* Clean up */
proc datasets library=work nolist;
    delete interval_values_y10;
quit;






/************************ INTERVAL 1 (1986) ANALYSIS *****************************************/

/* Extract diet and covariate values at interval 1 */
proc sort data=pool_with_diet_flag; by newid intervalp; run;

data interval_values_y86;
    set pool_with_diet_flag;
    where intervalp = 1 and has_diet = 1;
    
    /* Store diet values from this interval */
    lcs_y86 = lcssvcon;
    alcs_y86 = alcssvcon;
    vlcs_y86 = vlcssvcon;
    ulcs_y86 = ulcssvcon;
    hlcs_y86 = hlcssvcon;
    lfs_y86 = lfssvcon;
    ulfs_y86 = ulfssvcon;
    hlfs_y86 = hlfssvcon;
    alfs_y86 = alfssvcon;
    vlfs_y86 = vlfssvcon;
    
    /* Store covariate values from this interval */
    calor_y86 = calorsvcon;
    ahei_y86 = aheimsvcon;
    smk_y86 = smkstatus;
    alc_y86 = alc;
    act_y86 = actc;
    bmi_y86 = bmib9;
    mam_y86       = mam;
    endo2c_y86    = endo2c;
    endo3c_y86    = endo3c;
    screen_y86 = screen_cancer  ; 
    nSES_y86 = nSES  ; 

    keep newid lcs_y86 alcs_y86 vlcs_y86 ulcs_y86 hlcs_y86 
         lfs_y86 ulfs_y86 hlfs_y86 alfs_y86 vlfs_y86
         calor_y86 ahei_y86 smk_y86 alc_y86 act_y86 bmi_y86
         mam_y86 endo2c_y86 nSES_y86  endo3c_y86 screen_y86 ;
run;

/* Merge these values back to all future intervals */
data analysis_y86;
    merge pool_with_diet_flag (in=dietflag) interval_values_y86 (in=intvalue);
    by newid;
    
    /* Keep only observations at or after interval 1 */
    /* AND only people who had diet data at interval 1 */
    if dietflag and intvalue and intervalp >= 1;
run;

/* Create quintiles based on the interval-specific values */
proc rank data=analysis_y86 groups=5 out=analysis_y86;
    var lcs_y86 alcs_y86 vlcs_y86 ulcs_y86 hlcs_y86 lfs_y86 ulfs_y86 hlfs_y86 alfs_y86 vlfs_y86 calor_y86 ahei_y86;
    ranks lcs_qy86 alcs_qy86 vlcs_qy86 ulcs_qy86 hlcs_qy86 lfs_qy86 ulfs_qy86 hlfs_qy86 alfs_qy86 vlfs_qy86 calor_qy86 ahei_qy86;
run;

/* Convert to 1-5 scale and create indicator variables */
data analysis_y86;
    set analysis_y86 end=_end_;
    
    /* Convert ranks 0-4 to quintiles 1-5 */
    array qvars lcs_qy86 alcs_qy86 vlcs_qy86 ulcs_qy86 hlcs_qy86 lfs_qy86 ulfs_qy86 hlfs_qy86 alfs_qy86 vlfs_qy86 calor_qy86 ahei_qy86;
    do over qvars;
        if qvars ne . then qvars = qvars + 1;
    end;
    
    /* Create indicator variables - INSIDE the DATA step */
    %indic3(vbl=lcs_qy86,   prefix=lcs_qy86,   min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=alcs_qy86,  prefix=alcs_qy86,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=vlcs_qy86,  prefix=vlcs_qy86,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=ulcs_qy86,  prefix=ulcs_qy86,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=hlcs_qy86,  prefix=hlcs_qy86,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=lfs_qy86,   prefix=lfs_qy86,   min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=ulfs_qy86,  prefix=ulfs_qy86,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=hlfs_qy86,  prefix=hlfs_qy86,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=alfs_qy86,  prefix=alfs_qy86,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=vlfs_qy86,  prefix=vlfs_qy86,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=calor_qy86, prefix=calor_qy86, min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=ahei_qy86,  prefix=ahei_qy86,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=smk_y86,    prefix=smk_y86,    min=2, max=3, reflev=1, missing=., usemiss=0);
    %indic3(vbl=alc_y86,    prefix=alc_y86,    min=2, max=6, reflev=1, missing=., usemiss=0);
    %indic3(vbl=act_y86,    prefix=act_y86,    min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=bmi_y86,    prefix=bmi_y86,    min=2, max=9, reflev=1, missing=., usemiss=1);
run;

/* Run Cox models */
%res(analysis=interval_y86_adj0, catlin=T,   
     data=analysis_y86, event=can_obese, time=tobeca, 
     timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
     id=newid, 
     tvar=intervalp,
     agevar=agecon,
     strata=agecon cohort intervalp, 
     labels=F, 
     exp1=&lcs_qy86_,   adjust1=,
     exp2=&alcs_qy86_,  adjust2=,
     exp3=&vlcs_qy86_,  adjust3=,
     exp4=&ulcs_qy86_,  adjust4=,
     exp5=&hlcs_qy86_,  adjust5=,
     exp6=&lfs_qy86_,   adjust6=,
     exp7=&alfs_qy86_,  adjust7=,
     exp8=&vlfs_qy86_,  adjust8=,
     exp9=&ulfs_qy86_,  adjust9=,
     exp10=&hlfs_qy86_, adjust10=);
Run;

%res(analysis=interval_y86_adj1, catlin=T,   
     data=analysis_y86, event=can_obese, time=tobeca, 
     timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
     id=newid, 
     tvar=intervalp,
     agevar=agecon,
     strata=agecon cohort intervalp, 
     labels=F, 
     exp1=&lcs_qy86_,   adjust1= &calor_qy86_ &ahei_qy86_ &smk_y86_ &alc_y86_ &act_y86_  white mv asp cafh   endo2c_y86 nSES_y86    ,
     exp2=&alcs_qy86_,  adjust2= &calor_qy86_ &ahei_qy86_ &smk_y86_ &alc_y86_ &act_y86_  white mv asp cafh   endo2c_y86 nSES_y86    ,
     exp3=&vlcs_qy86_,  adjust3= &calor_qy86_ &ahei_qy86_ &smk_y86_ &alc_y86_ &act_y86_  white mv asp cafh   endo2c_y86 nSES_y86    ,
     exp4=&ulcs_qy86_,  adjust4= &calor_qy86_ &ahei_qy86_ &smk_y86_ &alc_y86_ &act_y86_  white mv asp cafh   endo2c_y86 nSES_y86    ,
     exp5=&hlcs_qy86_,  adjust5= &calor_qy86_ &ahei_qy86_ &smk_y86_ &alc_y86_ &act_y86_  white mv asp cafh   endo2c_y86 nSES_y86    ,
     exp6=&lfs_qy86_,   adjust6= &calor_qy86_ &ahei_qy86_ &smk_y86_ &alc_y86_ &act_y86_  white mv asp cafh   endo2c_y86 nSES_y86    ,
     exp7=&alfs_qy86_,  adjust7= &calor_qy86_ &ahei_qy86_ &smk_y86_ &alc_y86_ &act_y86_  white mv asp cafh   endo2c_y86 nSES_y86    ,
     exp8=&vlfs_qy86_,  adjust8= &calor_qy86_ &ahei_qy86_ &smk_y86_ &alc_y86_ &act_y86_  white mv asp cafh   endo2c_y86 nSES_y86    ,
     exp9=&ulfs_qy86_,  adjust9= &calor_qy86_ &ahei_qy86_ &smk_y86_ &alc_y86_ &act_y86_  white mv asp cafh   endo2c_y86 nSES_y86    ,
     exp10=&hlfs_qy86_, adjust10= &calor_qy86_ &ahei_qy86_ &smk_y86_ &alc_y86_ &act_y86_  white mv asp cafh   endo2c_y86 nSES_y86    );
Run;

/* Clean up */
proc datasets library=work nolist;
    delete interval_values_y86;
quit;




/************************ INTERVAL 15 (2014) ANALYSIS *****************************************/

/* Extract diet and covariate values at interval 1 */
proc sort data=pool_with_diet_flag; by newid intervalp; run;

data interval_values_y14;
    set pool_with_diet_flag;
    where (cohort=3 and intervalp=15) or (cohort=2 and intervalp=16);
    
    /* Store diet values from this interval */
    lcs_y14 = lcssvcon;
    alcs_y14 = alcssvcon;
    vlcs_y14 = vlcssvcon;
    ulcs_y14 = ulcssvcon;
    hlcs_y14 = hlcssvcon;
    lfs_y14 = lfssvcon;
    ulfs_y14 = ulfssvcon;
    hlfs_y14 = hlfssvcon;
    alfs_y14 = alfssvcon;
    vlfs_y14 = vlfssvcon;
    
    /* Store covariate values from this interval */
    calor_y14 = calorsvcon;
    ahei_y14 = aheimsvcon;
    smk_y14 = smkstatus;
    alc_y14 = alc;
    act_y14 = actc;
    bmi_y14 = bmib9;
    mam_y14       = mam;
    endo2c_y14    = endo2c;
    endo3c_y14    = endo3c;
    screen_y14 = screen_cancer  ;
    nSES_y14 = nSES  ;

    keep newid lcs_y14 alcs_y14 vlcs_y14 ulcs_y14 hlcs_y14 
         lfs_y14 ulfs_y14 hlfs_y14 alfs_y14 vlfs_y14
         calor_y14 ahei_y14 smk_y14 alc_y14 act_y14 bmi_y14
         mam_y14 endo2c_y14 nSES_y14  endo3c_y14 screen_y14;
run;

/* Merge these values back to all future intervals */
data analysis_y14;
    merge pool_with_diet_flag (in=dietflag) interval_values_y14 (in=intvalue);
    by newid;
    
    /* Keep only observations at or after interval 15 */
    /* AND only people who had diet data at interval 15 */
    if dietflag and intvalue and intervalp >= 15;
run;

/* Create quintiles based on the interval-specific values */
proc rank data=analysis_y14 groups=5 out=analysis_y14;
    var lcs_y14 alcs_y14 vlcs_y14 ulcs_y14 hlcs_y14 lfs_y14 ulfs_y14 hlfs_y14 alfs_y14 vlfs_y14 calor_y14 ahei_y14;
    ranks lcs_qy14 alcs_qy14 vlcs_qy14 ulcs_qy14 hlcs_qy14 lfs_qy14 ulfs_qy14 hlfs_qy14 alfs_qy14 vlfs_qy14 calor_qy14 ahei_qy14;
run;

/* Convert to 1-5 scale and create indicator variables */
data analysis_y14;
    set analysis_y14 end=_end_;
    
    /* Convert ranks 0-4 to quintiles 1-5 */
    array qvars lcs_qy14 alcs_qy14 vlcs_qy14 ulcs_qy14 hlcs_qy14 lfs_qy14 ulfs_qy14 hlfs_qy14 alfs_qy14 vlfs_qy14 calor_qy14 ahei_qy14;
    do over qvars;
        if qvars ne . then qvars = qvars + 1;
    end;
    
    /* Create indicator variables - INSIDE the DATA step */
    %indic3(vbl=lcs_qy14,   prefix=lcs_qy14,   min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=alcs_qy14,  prefix=alcs_qy14,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=vlcs_qy14,  prefix=vlcs_qy14,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=ulcs_qy14,  prefix=ulcs_qy14,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=hlcs_qy14,  prefix=hlcs_qy14,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=lfs_qy14,   prefix=lfs_qy14,   min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=ulfs_qy14,  prefix=ulfs_qy14,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=hlfs_qy14,  prefix=hlfs_qy14,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=alfs_qy14,  prefix=alfs_qy14,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=vlfs_qy14,  prefix=vlfs_qy14,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=calor_qy14, prefix=calor_qy14, min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=ahei_qy14,  prefix=ahei_qy14,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=smk_y14,    prefix=smk_y14,    min=2, max=3, reflev=1, missing=., usemiss=0);
    %indic3(vbl=alc_y14,    prefix=alc_y14,    min=2, max=6, reflev=1, missing=., usemiss=0);
    %indic3(vbl=act_y14,    prefix=act_y14,    min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=bmi_y14,    prefix=bmi_y14,    min=2, max=9, reflev=1, missing=., usemiss=1);
run;

/* Run Cox models */
%res(analysis=interval_y14_adj0, catlin=T,   
     data=analysis_y14, event=can_obese, time=tobeca, 
     timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
     id=newid, 
     tvar=intervalp,
     agevar=agecon,
     strata=agecon cohort intervalp, 
     labels=F, 
     exp1=&lcs_qy14_,   adjust1=,
     exp2=&alcs_qy14_,  adjust2=,
     exp3=&vlcs_qy14_,  adjust3=,
     exp4=&ulcs_qy14_,  adjust4=,
     exp5=&hlcs_qy14_,  adjust5=,
     exp6=&lfs_qy14_,   adjust6=,
     exp7=&alfs_qy14_,  adjust7=,
     exp8=&vlfs_qy14_,  adjust8=,
     exp9=&ulfs_qy14_,  adjust9=,
     exp10=&hlfs_qy14_, adjust10=);
Run;

%res(analysis=interval_y14_adj1, catlin=T,   
     data=analysis_y14, event=can_obese, time=tobeca, 
     timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
     id=newid, 
     tvar=intervalp,
     agevar=agecon,
     strata=agecon cohort intervalp, 
     labels=F, 
     exp1=&lcs_qy14_,   adjust1= &calor_qy14_ &ahei_qy14_ &smk_y14_ &alc_y14_ &act_y14_ &bmi_y14_ white mv asp cafh mam_y14 endo2c_y14 nSES_y14    pmhever ocever,
     exp2=&alcs_qy14_,  adjust2= &calor_qy14_ &ahei_qy14_ &smk_y14_ &alc_y14_ &act_y14_ &bmi_y14_ white mv asp cafh mam_y14 endo2c_y14 nSES_y14    pmhever ocever,
     exp3=&vlcs_qy14_,  adjust3= &calor_qy14_ &ahei_qy14_ &smk_y14_ &alc_y14_ &act_y14_ &bmi_y14_ white mv asp cafh mam_y14 endo2c_y14 nSES_y14    pmhever ocever,
     exp4=&ulcs_qy14_,  adjust4= &calor_qy14_ &ahei_qy14_ &smk_y14_ &alc_y14_ &act_y14_ &bmi_y14_ white mv asp cafh mam_y14 endo2c_y14 nSES_y14    pmhever ocever,
     exp5=&hlcs_qy14_,  adjust5= &calor_qy14_ &ahei_qy14_ &smk_y14_ &alc_y14_ &act_y14_ &bmi_y14_ white mv asp cafh mam_y14 endo2c_y14 nSES_y14    pmhever ocever,
     exp6=&lfs_qy14_,   adjust6= &calor_qy14_ &ahei_qy14_ &smk_y14_ &alc_y14_ &act_y14_ &bmi_y14_ white mv asp cafh mam_y14 endo2c_y14 nSES_y14    pmhever ocever,
     exp7=&alfs_qy14_,  adjust7= &calor_qy14_ &ahei_qy14_ &smk_y14_ &alc_y14_ &act_y14_ &bmi_y14_ white mv asp cafh mam_y14 endo2c_y14 nSES_y14    pmhever ocever,
     exp8=&vlfs_qy14_,  adjust8= &calor_qy14_ &ahei_qy14_ &smk_y14_ &alc_y14_ &act_y14_ &bmi_y14_ white mv asp cafh mam_y14 endo2c_y14 nSES_y14    pmhever ocever,
     exp9=&ulfs_qy14_,  adjust9= &calor_qy14_ &ahei_qy14_ &smk_y14_ &alc_y14_ &act_y14_ &bmi_y14_ white mv asp cafh mam_y14 endo2c_y14 nSES_y14    pmhever ocever,
     exp10=&hlfs_qy14_, adjust10= &calor_qy14_ &ahei_qy14_ &smk_y14_ &alc_y14_ &act_y14_ &bmi_y14_ white mv asp cafh mam_y14 endo2c_y14 nSES_y14    pmhever ocever);
Run;

/* Clean up */
proc datasets library=work nolist;
    delete interval_values_y14;
quit;





/************************ INTERVAL 17 (2018) ANALYSIS *****************************************/

/* Extract diet and covariate values at interval 1 */
proc sort data=pool_with_diet_flag; by newid intervalp; run;

data interval_values_y18;
    set pool_with_diet_flag;
    where (cohort=3 and intervalp=17) or (cohort=2 and intervalp=18);
    
    /* Store diet values from this interval */
    lcs_y18 = lcssvcon;
    alcs_y18 = alcssvcon;
    vlcs_y18 = vlcssvcon;
    ulcs_y18 = ulcssvcon;
    hlcs_y18 = hlcssvcon;
    lfs_y18 = lfssvcon;
    ulfs_y18 = ulfssvcon;
    hlfs_y18 = hlfssvcon;
    alfs_y18 = alfssvcon;
    vlfs_y18 = vlfssvcon;
    
    /* Store covariate values from this interval */
    calor_y18 = calorsvcon;
    ahei_y18 = aheimsvcon;
    smk_y18 = smkstatus;
    alc_y18 = alc;
    act_y18 = actc;
    bmi_y18 = bmib9;
    mam_y18       = mam;
    endo2c_y18    = endo2c;
    endo3c_y18    = endo3c;
    screen_y18 = screen_cancer  ;
    nSES_y18 = nSES  ;

    keep newid lcs_y18 alcs_y18 vlcs_y18 ulcs_y18 hlcs_y18 
         lfs_y18 ulfs_y18 hlfs_y18 alfs_y18 vlfs_y18
         calor_y18 ahei_y18 smk_y18 alc_y18 act_y18 bmi_y18
         mam_y18 endo2c_y18 nSES_y18  endo3c_y18 screen_y18 ;
run;

/* Merge these values back to all future intervals */
data analysis_y18;
    merge pool_with_diet_flag (in=dietflag) interval_values_y18 (in=intvalue);
    by newid;
    
    /* Keep only observations at or after interval 17 */
    /* AND only people who had diet data at interval 17 */
    if dietflag and intvalue and intervalp >= 17;
run;

/* Create quintiles based on the interval-specific values */
proc rank data=analysis_y18 groups=5 out=analysis_y18;
    var lcs_y18 alcs_y18 vlcs_y18 ulcs_y18 hlcs_y18 lfs_y18 ulfs_y18 hlfs_y18 alfs_y18 vlfs_y18 calor_y18 ahei_y18;
    ranks lcs_qy18 alcs_qy18 vlcs_qy18 ulcs_qy18 hlcs_qy18 lfs_qy18 ulfs_qy18 hlfs_qy18 alfs_qy18 vlfs_qy18 calor_qy18 ahei_qy18;
run;

/* Convert to 1-5 scale and create indicator variables */
data analysis_y18;
    set analysis_y18 end=_end_;
    
    /* Convert ranks 0-4 to quintiles 1-5 */
    array qvars lcs_qy18 alcs_qy18 vlcs_qy18 ulcs_qy18 hlcs_qy18 lfs_qy18 ulfs_qy18 hlfs_qy18 alfs_qy18 vlfs_qy18 calor_qy18 ahei_qy18;
    do over qvars;
        if qvars ne . then qvars = qvars + 1;
    end;
    
    /* Create indicator variables - INSIDE the DATA step */
    %indic3(vbl=lcs_qy18,   prefix=lcs_qy18,   min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=alcs_qy18,  prefix=alcs_qy18,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=vlcs_qy18,  prefix=vlcs_qy18,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=ulcs_qy18,  prefix=ulcs_qy18,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=hlcs_qy18,  prefix=hlcs_qy18,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=lfs_qy18,   prefix=lfs_qy18,   min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=ulfs_qy18,  prefix=ulfs_qy18,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=hlfs_qy18,  prefix=hlfs_qy18,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=alfs_qy18,  prefix=alfs_qy18,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=vlfs_qy18,  prefix=vlfs_qy18,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=calor_qy18, prefix=calor_qy18, min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=ahei_qy18,  prefix=ahei_qy18,  min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=smk_y18,    prefix=smk_y18,    min=2, max=3, reflev=1, missing=., usemiss=0);
    %indic3(vbl=alc_y18,    prefix=alc_y18,    min=2, max=6, reflev=1, missing=., usemiss=0);
    %indic3(vbl=act_y18,    prefix=act_y18,    min=2, max=5, reflev=1, missing=., usemiss=0);
    %indic3(vbl=bmi_y18,    prefix=bmi_y18,    min=2, max=9, reflev=1, missing=., usemiss=1);
run;

/* Run Cox models */
%res(analysis=interval_y18_adj0, catlin=T,   
     data=analysis_y18, event=can_obese, time=tobeca, 
     timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
     id=newid, 
     tvar=intervalp,
     agevar=agecon,
     strata=agecon cohort intervalp, 
     labels=F, 
     exp1=&lcs_qy18_,   adjust1=,
     exp2=&alcs_qy18_,  adjust2=,
     exp3=&vlcs_qy18_,  adjust3=,
     exp4=&ulcs_qy18_,  adjust4=,
     exp5=&hlcs_qy18_,  adjust5=,
     exp6=&lfs_qy18_,   adjust6=,
     exp7=&alfs_qy18_,  adjust7=,
     exp8=&vlfs_qy18_,  adjust8=,
     exp9=&ulfs_qy18_,  adjust9=,
     exp10=&hlfs_qy18_, adjust10=);
Run;

%res(analysis=interval_y18_adj1, catlin=T,   
     data=analysis_y18, event=can_obese, time=tobeca, 
     timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
     id=newid, 
     tvar=intervalp,
     agevar=agecon,
     strata=agecon cohort intervalp, 
     labels=F, 
     exp1=&lcs_qy18_,   adjust1= &calor_qy18_ &ahei_qy18_ &smk_y18_ &alc_y18_ &act_y18_ &act_y18_  white mv asp cafh   endo2c_y18 nSES_y18  ,
     exp2=&alcs_qy18_,  adjust2= &calor_qy18_ &ahei_qy18_ &smk_y18_ &alc_y18_ &act_y18_ &act_y18_  white mv asp cafh   endo2c_y18 nSES_y18  ,
     exp3=&vlcs_qy18_,  adjust3= &calor_qy18_ &ahei_qy18_ &smk_y18_ &alc_y18_ &act_y18_ &act_y18_  white mv asp cafh   endo2c_y18 nSES_y18  ,
     exp4=&ulcs_qy18_,  adjust4= &calor_qy18_ &ahei_qy18_ &smk_y18_ &alc_y18_ &act_y18_ &act_y18_  white mv asp cafh   endo2c_y18 nSES_y18  ,
     exp5=&hlcs_qy18_,  adjust5= &calor_qy18_ &ahei_qy18_ &smk_y18_ &alc_y18_ &act_y18_ &act_y18_  white mv asp cafh   endo2c_y18 nSES_y18  ,
     exp6=&lfs_qy18_,   adjust6= &calor_qy18_ &ahei_qy18_ &smk_y18_ &alc_y18_ &act_y18_ &act_y18_  white mv asp cafh   endo2c_y18 nSES_y18  ,
     exp7=&alfs_qy18_,  adjust7= &calor_qy18_ &ahei_qy18_ &smk_y18_ &alc_y18_ &act_y18_ &act_y18_  white mv asp cafh   endo2c_y18 nSES_y18  ,
     exp8=&vlfs_qy18_,  adjust8= &calor_qy18_ &ahei_qy18_ &smk_y18_ &alc_y18_ &act_y18_ &act_y18_  white mv asp cafh   endo2c_y18 nSES_y18  ,
     exp9=&ulfs_qy18_,  adjust9= &calor_qy18_ &ahei_qy18_ &smk_y18_ &alc_y18_ &act_y18_ &act_y18_  white mv asp cafh   endo2c_y18 nSES_y18  ,
     exp10=&hlfs_qy18_, adjust10= &calor_qy18_ &ahei_qy18_ &smk_y18_ &alc_y18_ &act_y18_ &act_y18_  white mv asp cafh   endo2c_y18 nSES_y18  );
Run;

/* Clean up */
proc datasets library=work nolist;
    delete interval_values_y18;
quit;





/* Output results */
%resout(outfile=table/FigureS9_interval_specific, dbms=tab, opt=replace);