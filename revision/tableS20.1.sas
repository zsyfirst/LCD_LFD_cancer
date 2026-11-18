/*****************************************************
Purpose:       Table S20: Sensitivity analysis - lag analysis for incidence (4-year and 8-year lag)
               between assessment of dietary intake and obesity-related cancer
               
Study design:  prospective cohort (NHS, NHSII, HPFS)

******************************************************/

/************************ Call in Cleaned Datasets**********************************/
%include 'lowcarb_cancer/revision/data.inc.sas';

/************************ Data Preparation for Lag Analysis ************************/
/* 
Since the follow-up intervals are 2 years (e.g., tp01, tp02...), 
a 4-year lag corresponds to 2 intervals.
An 8-year lag corresponds to 4 intervals.

create a robust lagged dataset by merging the data with shifted intervalp values.
This ensures we are accurately capturing the exposure from exactly 4 or 8 years prior 
for the same individual (newid).
*/

* Variables to lag (exposures and covariates);
%let lag_vars = lcssvq alcssvq vlcssvq ulcssvq hlcssvq lfssvq ulfssvq hlfssvq alfssvq vlfssvq
                lcs_sd alcs_sd vlcs_sd ulcs_sd hlcs_sd lfs_sd ulfs_sd hlfs_sd alfs_sd vlfs_sd
                lcssvtrend alcssvtrend vlcssvtrend ulcssvtrend hlcssvtrend lfssvtrend ulfssvtrend hlfssvtrend alfssvtrend vlfssvtrend
                calorsvq aheimsvq smkstatus alc actc bmib9 white mv asp mifh highbp hightc cafh endo2c mam  pmhever ocever nSES ;

proc sort data=pool; by newid intervalp; run;

/* Create 4-year lag dataset (Shift intervals by +2) */
data pool_lag4(keep=newid intervalp &lag_vars);
    set pool;
    intervalp = intervalp + 2; 
run;

proc sort data=pool_lag4; by newid intervalp; run;

/* Create 8-year lag dataset (Shift intervals by +4) */
data pool_lag8(keep=newid intervalp &lag_vars);
    set pool;
    intervalp = intervalp + 4; 
run;

proc sort data=pool_lag8; by newid intervalp; run;

/* Merge into main analysis datasets */
/* 4-YEAR LAG DATASET */
data pool_4yr;
    merge pool(in=a drop=&lag_vars) pool_lag4(in=b);
    by newid intervalp;
    if a; /* keep only original follow-up periods */
    /* Remove missing lag periods */
    if nmiss(lcssvq, calorsvq)=0; 
run;

/* 8-YEAR LAG DATASET */
data pool_8yr;
    merge pool(in=a drop=&lag_vars) pool_lag8(in=b);
    by newid intervalp;
    if a; 
    /* Remove missing lag periods */
    if nmiss(lcssvq, calorsvq)=0; 
run;


/************************ HR Models (4-YEAR LAG) *****************************************/
%include 'macro/res.sas';

title "------------------ 4-YEAR LAG: CATEGORICAL ------------------";
%res(analysis=avestop.lag4.adj1,catlin=T,   
         data=pool_4yr, event=can_obese, time=tobeca, 
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         agevar=agecon,
         strata=agecon cohort intervalp, 
         labels=F, 
         exp1=&lcssvq_   ,   adjust1= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp2=&alcssvq_  ,   adjust2= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp3=&vlcssvq_  ,   adjust3= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp4=&ulcssvq_  ,   adjust4= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp5=&hlcssvq_  ,   adjust5= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES , 
         exp6=&lfssvq_   ,   adjust6= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp7=&ulfssvq_  ,   adjust7= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp8=&hlfssvq_  ,   adjust8= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp9=&alfssvq_  ,   adjust9= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp10=&vlfssvq_  ,  adjust10=&calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES );
         Run;

title "------------------ 4-YEAR LAG: CONTINUOUS ------------------";
%res(analysis=unit.lag4.adj1,catlin=T,   
         data=pool_4yr, event=can_obese, time=tobeca, 
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         agevar=agecon,
         strata=agecon cohort intervalp, 
         labels=F, 
         exp1=lcs_sd   ,   adjust1= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp2=alcs_sd  ,   adjust2= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp3=vlcs_sd  ,   adjust3= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp4=ulcs_sd  ,   adjust4= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp5=hlcs_sd  ,   adjust5= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES , 
         exp6=lfs_sd   ,   adjust6= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp7=ulfs_sd  ,   adjust7= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp8=hlfs_sd  ,   adjust8= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp9=alfs_sd  ,   adjust9= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp10=vlfs_sd  ,  adjust10=&calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ); 
         Run;

title "------------------ 4-YEAR LAG: trend ------------------";
%res(analysis=trend.lag4.adj1,catlin=T,   
         data=pool_4yr, event=can_obese, time=tobeca, 
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         agevar=agecon,
         strata=agecon cohort intervalp, 
         labels=F, 
         exp1=lcssvtrend   ,   adjust1= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp2=alcssvtrend  ,   adjust2= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp3=vlcssvtrend  ,   adjust3= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp4=ulcssvtrend  ,   adjust4= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp5=hlcssvtrend  ,   adjust5= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES , 
         exp6=lfssvtrend   ,   adjust6= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp7=ulfssvtrend  ,   adjust7= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp8=hlfssvtrend  ,   adjust8= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp9=alfssvtrend  ,   adjust9= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp10=vlfssvtrend  ,  adjust10=&calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ); 
         Run;


/************************ HR Models (8-YEAR LAG) *****************************************/
title "------------------ 8-YEAR LAG: CATEGORICAL ------------------";
%res(analysis=avestop.lag8.adj1,catlin=T,   
         data=pool_8yr, event=can_obese, time=tobeca, 
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         agevar=agecon,
         strata=agecon cohort intervalp, 
         labels=F, 
         exp1=&lcssvq_   ,   adjust1= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp2=&alcssvq_  ,   adjust2= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp3=&vlcssvq_  ,   adjust3= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp4=&ulcssvq_  ,   adjust4= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp5=&hlcssvq_  ,   adjust5= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES , 
         exp6=&lfssvq_   ,   adjust6= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp7=&ulfssvq_  ,   adjust7= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp8=&hlfssvq_  ,   adjust8= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp9=&alfssvq_  ,   adjust9= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp10=&vlfssvq_  ,  adjust10=&calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES );
         Run;

title "------------------ 8-YEAR LAG: CONTINUOUS ------------------";
%res(analysis=unit.lag8.adj1,catlin=T,   
         data=pool_8yr, event=can_obese, time=tobeca, 
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         agevar=agecon,
         strata=agecon cohort intervalp, 
         labels=F, 
         exp1=lcs_sd   ,   adjust1= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp2=alcs_sd  ,   adjust2= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp3=vlcs_sd  ,   adjust3= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp4=ulcs_sd  ,   adjust4= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp5=hlcs_sd  ,   adjust5= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES , 
         exp6=lfs_sd   ,   adjust6= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp7=ulfs_sd  ,   adjust7= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp8=hlfs_sd  ,   adjust8= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp9=alfs_sd  ,   adjust9= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp10=vlfs_sd  ,  adjust10=&calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ); 
         Run;


title "------------------ 8-YEAR LAG: trend ------------------";
%res(analysis=trend.lag8.adj1,catlin=T,   
         data=pool_8yr, event=can_obese, time=tobeca, 
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         agevar=agecon,
         strata=agecon cohort intervalp, 
         labels=F, 
         exp1=lcssvtrend   ,   adjust1= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp2=alcssvtrend  ,   adjust2= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp3=vlcssvtrend  ,   adjust3= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp4=ulcssvtrend  ,   adjust4= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp5=hlcssvtrend  ,   adjust5= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES , 
         exp6=lfssvtrend   ,   adjust6= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp7=ulfssvtrend  ,   adjust7= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp8=hlfssvtrend  ,   adjust8= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp9=alfssvtrend  ,   adjust9= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ,
         exp10=vlfssvtrend  ,  adjust10=&calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever nSES ); 
         Run;


%include 'macro/resout.sas';
/************************ HR *****************************************/

%resout(outfile=table/tableS20, dbms=tab, opt=replace);                                                                                                                             
endsas;

