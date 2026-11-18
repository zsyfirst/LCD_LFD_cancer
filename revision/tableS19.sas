/*****************************************************
Purpose:       low carb, low fat, and risk of Obesity-related cancer
Study design:  prospective cohort (NHS: 1986-2022.6, NHSII: 1991-2023.6, HPFS: 1986-2024.1)

Contents:      Table S19 --- Stratified Analysis (Obesity-related cancer as outcome)
Stratified analysis by 
        agesub - time varying age (<=65, >65)
        bmic30 - time varying BMI (<30, ≥30)
        bmic25 - time varying BMI (<25, ≥25)
        active - time-varying physical activity (<=7.5, >7.5 mets/wk)
        neverdrinker - time-varying alcohol (everdrinker, neverdrinker)
        smknever - time-varying smoking (eversmker, neversmker)
Date:          10/2026
******************************************************/


/************************ Call in Cleaned Datasets**********************************/
%include '/lowcarb_cancer/revision/data.inc.sas';

  %let noage =  &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever;
  %let nobmi =  &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever;
  %let nopa =   &calorsvq_ &aheimsvq_ &smkstatus_ &alc_        &bmib9_ white mv asp cafh endo2c mam  pmhever ocever;
  %let nodrink =&calorsvq_ &aheimsvq_ &smkstatus_       &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever;
  %let nosmk =  &calorsvq_ &aheimsvq_             &alc_ &actc_ &bmib9_ white mv asp cafh endo2c mam  pmhever ocever;
  
/************************ HR *****************************************/
%include 'macro/res.sas';

/* age stratified */
%res(analysis=agestra,catlin=T,
    data=pool, event=can_obese, time=tobeca,
    timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
    id=newid,  agevar=agecon, tvar=intervalp,
    strata=agecon cohort intervalp, 
    labels=F,
    byvar=agesub,
    exp1=&lcssvq_   ,    adjust1= &noage,
    exp2=&alcssvq_  ,    adjust2= &noage,
    exp3=&vlcssvq_  ,    adjust3= &noage,
    exp4=&ulcssvq_  ,    adjust4= &noage,
    exp5=&hlcssvq_  ,    adjust5= &noage, 
    exp6=&lfssvq_   ,    adjust6= &noage,
    exp7=&ulfssvq_  ,    adjust7= &noage,
    exp8=&hlfssvq_  ,    adjust8= &noage,
    exp9=&alfssvq_  ,    adjust9= &noage,
    exp10=&vlfssvq_  ,    adjust10= &noage);
Run;


/* BMI stratified: bmic25, bmic30 */
%res(analysis=bmistra25,catlin=T,
    data=pool, event=can_obese, time=tobeca,
    timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
    id=newid,  agevar=agecon, tvar=intervalp,
    strata=agecon cohort intervalp, 
    labels=F,
    byvar=bmic25,
    exp1=&lcssvq_   ,    adjust1= &nobmi,
    exp2=&alcssvq_  ,    adjust2= &nobmi,
    exp3=&vlcssvq_  ,    adjust3= &nobmi,
    exp4=&ulcssvq_  ,    adjust4= &nobmi,
    exp5=&hlcssvq_  ,    adjust5= &nobmi, 
    exp6=&lfssvq_   ,    adjust6= &nobmi,
    exp7=&ulfssvq_  ,    adjust7= &nobmi,
    exp8=&hlfssvq_  ,    adjust8= &nobmi,
    exp9=&alfssvq_  ,    adjust9= &nobmi,
    exp10=&vlfssvq_  ,    adjust10= &nobmi);
run;


%res(analysis=bmistra30,catlin=T,
    data=pool, event=can_obese, time=tobeca,
    timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
    id=newid,  agevar=agecon, tvar=intervalp,
    strata=agecon cohort intervalp, 
    labels=F,
    byvar=bmic30,
    exp1=&lcssvq_   ,    adjust1= &nobmi,
    exp2=&alcssvq_  ,    adjust2= &nobmi,
    exp3=&vlcssvq_  ,    adjust3= &nobmi,
    exp4=&ulcssvq_  ,    adjust4= &nobmi,
    exp5=&hlcssvq_  ,    adjust5= &nobmi, 
    exp6=&lfssvq_   ,    adjust6= &nobmi,
    exp7=&ulfssvq_  ,    adjust7= &nobmi,
    exp8=&hlfssvq_  ,    adjust8= &nobmi,
    exp9=&alfssvq_  ,    adjust9= &nobmi,
    exp10=&vlfssvq_  ,    adjust10= &nobmi);
run;


/* physical activity stratified */
%res(analysis=activestra,catlin=T,
    data=pool, event=can_obese, time=tobeca,
    timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
    id=newid,  agevar=agecon, tvar=intervalp,
    strata=agecon cohort intervalp, 
    labels=F,
    byvar=active,
    exp1=&lcssvq_   ,    adjust1= &nopa,
    exp2=&alcssvq_  ,    adjust2= &nopa,
    exp3=&vlcssvq_  ,    adjust3= &nopa,
    exp4=&ulcssvq_  ,    adjust4= &nopa,
    exp5=&hlcssvq_  ,    adjust5= &nopa, 
    exp6=&lfssvq_   ,    adjust6= &nopa,
    exp7=&ulfssvq_  ,    adjust7= &nopa,
    exp8=&hlfssvq_  ,    adjust8= &nopa,
    exp9=&alfssvq_  ,    adjust9= &nopa,
    exp10=&vlfssvq_  ,    adjust10= &nopa);
run;

/* Never drinking stratified */
%res(analysis=activestra,catlin=T,
    data=pool, event=can_obese, time=tobeca,
    timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
    id=newid,  agevar=agecon, tvar=intervalp,
    strata=agecon cohort intervalp, 
    labels=F,
    byvar=neverdrinker,
    exp1=&lcssvq_   ,    adjust1= &nodrink,
    exp2=&alcssvq_  ,    adjust2= &nodrink,
    exp3=&vlcssvq_  ,    adjust3= &nodrink,
    exp4=&ulcssvq_  ,    adjust4= &nodrink,
    exp5=&hlcssvq_  ,    adjust5= &nodrink, 
    exp6=&lfssvq_   ,    adjust6= &nodrink,
    exp7=&ulfssvq_  ,    adjust7= &nodrink,
    exp8=&hlfssvq_  ,    adjust8= &nodrink,
    exp9=&alfssvq_  ,    adjust9= &nodrink,
    exp10=&vlfssvq_  ,    adjust10= &nodrink);
run;


/* Never smoking stratified */
%res(analysis=activestra,catlin=T,
    data=pool, event=can_obese, time=tobeca,
    timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
    id=newid,  agevar=agecon, tvar=intervalp,
    strata=agecon cohort intervalp, 
    labels=F,
    byvar=smknever,
    exp1=&lcssvq_   ,    adjust1= &nosmk,
    exp2=&alcssvq_  ,    adjust2= &nosmk,
    exp3=&vlcssvq_  ,    adjust3= &nosmk,
    exp4=&ulcssvq_  ,    adjust4= &nosmk,
    exp5=&hlcssvq_  ,    adjust5= &nosmk, 
    exp6=&lfssvq_   ,    adjust6= &nosmk,
    exp7=&ulfssvq_  ,    adjust7= &nosmk,
    exp8=&hlfssvq_  ,    adjust8= &nosmk,
    exp9=&alfssvq_  ,    adjust9= &nosmk,
    exp10=&vlfssvq_  ,    adjust10= &nosmk);
run;



/**********trend***********/
/* age stratified */
%res(analysis=agestra,catlin=T,
    data=pool, event=can_obese, time=tobeca,
    timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
    id=newid,  agevar=agecon, tvar=intervalp,
    strata=agecon cohort intervalp, 
    labels=F,
    byvar=agesub,
    exp1=lcssvtrend  ,    adjust1= &noage,
    exp2=alcssvtrend ,    adjust2= &noage,
    exp3=vlcssvtrend ,    adjust3= &noage,
    exp4=ulcssvtrend ,    adjust4= &noage,
    exp5=hlcssvtrend ,    adjust5= &noage, 
    exp6=lfssvtrend  ,    adjust6= &noage,
    exp7=alfssvtrend ,    adjust7= &noage,
    exp8=vlfssvtrend ,    adjust8= &noage,
    exp9=ulfssvtrend ,    adjust9= &noage,
    exp10=hlfssvtrend ,    adjust10= &noage);
Run;


/* BMI stratified: bmic25, bmic30 */
%res(analysis=bmistra25,catlin=T,
    data=pool, event=can_obese, time=tobeca,
    timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
    id=newid,  agevar=agecon, tvar=intervalp,
    strata=agecon cohort intervalp, 
    labels=F,
    byvar=bmic25,
    exp1=lcssvtrend  ,    adjust1= &nobmi,
    exp2=alcssvtrend ,    adjust2= &nobmi,
    exp3=vlcssvtrend ,    adjust3= &nobmi,
    exp4=ulcssvtrend ,    adjust4= &nobmi,
    exp5=hlcssvtrend ,    adjust5= &nobmi, 
    exp6=lfssvtrend  ,    adjust6= &nobmi,
    exp7=alfssvtrend ,    adjust7= &nobmi,
    exp8=vlfssvtrend ,    adjust8= &nobmi,
    exp9=ulfssvtrend ,    adjust9= &nobmi,
    exp10=hlfssvtrend ,    adjust10= &nobmi);
run;


%res(analysis=bmistra30,catlin=T,
    data=pool, event=can_obese, time=tobeca,
    timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
    id=newid,  agevar=agecon, tvar=intervalp,
    strata=agecon cohort intervalp, 
    labels=F,
    byvar=bmic30,
    exp1=lcssvtrend  ,    adjust1= &nobmi,
    exp2=alcssvtrend ,    adjust2= &nobmi,
    exp3=vlcssvtrend ,    adjust3= &nobmi,
    exp4=ulcssvtrend ,    adjust4= &nobmi,
    exp5=hlcssvtrend ,    adjust5= &nobmi, 
    exp6=lfssvtrend  ,    adjust6= &nobmi,
    exp7=alfssvtrend ,    adjust7= &nobmi,
    exp8=vlfssvtrend ,    adjust8= &nobmi,
    exp9=ulfssvtrend ,    adjust9= &nobmi,
    exp10=hlfssvtrend ,    adjust10= &nobmi);
run;


/* physical activity stratified */
%res(analysis=activestra,catlin=T,
    data=pool, event=can_obese, time=tobeca,
    timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
    id=newid,  agevar=agecon, tvar=intervalp,
    strata=agecon cohort intervalp, 
    labels=F,
    byvar=active,
    exp1=lcssvtrend  ,    adjust1= &nopa,
    exp2=alcssvtrend ,    adjust2= &nopa,
    exp3=vlcssvtrend ,    adjust3= &nopa,
    exp4=ulcssvtrend ,    adjust4= &nopa,
    exp5=hlcssvtrend ,    adjust5= &nopa, 
    exp6=lfssvtrend  ,    adjust6= &nopa,
    exp7=alfssvtrend ,    adjust7= &nopa,
    exp8=vlfssvtrend ,    adjust8= &nopa,
    exp9=ulfssvtrend ,    adjust9= &nopa,
    exp10=hlfssvtrend ,    adjust10= &nopa);
run;

/* Never drinking stratified */
%res(analysis=activestra,catlin=T,
    data=pool, event=can_obese, time=tobeca,
    timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
    id=newid,  agevar=agecon, tvar=intervalp,
    strata=agecon cohort intervalp, 
    labels=F,
    byvar=neverdrinker,
    exp1=lcssvtrend  ,    adjust1= &nodrink,
    exp2=alcssvtrend ,    adjust2= &nodrink,
    exp3=vlcssvtrend ,    adjust3= &nodrink,
    exp4=ulcssvtrend ,    adjust4= &nodrink,
    exp5=hlcssvtrend ,    adjust5= &nodrink, 
    exp6=lfssvtrend  ,    adjust6= &nodrink,
    exp7=alfssvtrend ,    adjust7= &nodrink,
    exp8=vlfssvtrend ,    adjust8= &nodrink,
    exp9=ulfssvtrend ,    adjust9= &nodrink,
    exp10=hlfssvtrend ,    adjust10= &nodrink);
run;


/* Never smoking stratified */
%res(analysis=activestra,catlin=T,
    data=pool, event=can_obese, time=tobeca,
    timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
    id=newid,  agevar=agecon, tvar=intervalp,
    strata=agecon cohort intervalp, 
    labels=F,
    byvar=smknever,
    exp1=lcssvtrend  ,    adjust1= &nosmk,
    exp2=alcssvtrend ,    adjust2= &nosmk,
    exp3=vlcssvtrend ,    adjust3= &nosmk,
    exp4=ulcssvtrend ,    adjust4= &nosmk,
    exp5=hlcssvtrend ,    adjust5= &nosmk, 
    exp6=lfssvtrend  ,    adjust6= &nosmk,
    exp7=alfssvtrend ,    adjust7= &nosmk,
    exp8=vlfssvtrend ,    adjust8= &nosmk,
    exp9=ulfssvtrend ,    adjust9= &nosmk,
    exp10=hlfssvtrend ,    adjust10= &nosmk);
run;



%include 'macro/resout.sas';
/************************ HR *****************************************/

%resout(outfile=table/tableS19, dbms=tab, opt=replace);                                                                                                                             
endsas;