/*****************************************************
Purpose:       low carb, low fat, and risk of Obesity-related cancer
Study design:  prospective cohort (NHS: 1986-2022.6, NHSII: 1991-2023.6, HPFS: 1986-2024.1)
Contents:      Table S20 --- PY and HR (Obesity-related cancer as outcome, Without Adjustment for AHEI-2010)
Date:          10/2026
******************************************************/

/************************ Call in Cleaned Datasets**********************************/
%include 'lowcarb_cancer/data.inc.sas';

/************************ HR *****************************************/
%include 'macro/res.sas';

%res(analysis=avestop.adj0, catlin=T,   
         data=pool, event=can_obese, time=tobeca, 
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         agevar=agecon,
         strata=agecon cohort intervalp, 
         labels=F, 
         exp1=&lcssvq_   ,    adjust1=   ,
         exp2=&alcssvq_  ,    adjust2=   ,
         exp3=&vlcssvq_  ,    adjust3=   ,
         exp4=&ulcssvq_  ,    adjust4=   ,
         exp5=&hlcssvq_  ,    adjust5=   , 
         exp6=&lfssvq_   ,    adjust6=   ,
         exp7=&ulfssvq_  ,    adjust7=   ,
         exp8=&hlfssvq_  ,    adjust8=   ,
         exp9=&alfssvq_  ,    adjust9=   ,
         exp10=&vlfssvq_  ,    adjust10=   );
         Run;

%res(analysis=avestop.adj1,catlin=T,   
         data=pool, event=can_obese, time=tobeca, 
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         agevar=agecon,
         strata=agecon cohort intervalp, 
         labels=F, 
         exp1=&lcssvq_   ,   adjust1= &calorsvq_   &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp2=&alcssvq_  ,   adjust2= &calorsvq_   &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp3=&vlcssvq_  ,   adjust3= &calorsvq_   &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp4=&ulcssvq_  ,   adjust4= &calorsvq_   &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp5=&hlcssvq_  ,   adjust5= &calorsvq_   &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  , 
         exp6=&lfssvq_   ,   adjust6= &calorsvq_   &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp7=&ulfssvq_  ,   adjust7= &calorsvq_   &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp8=&hlfssvq_  ,   adjust8= &calorsvq_   &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp9=&alfssvq_  ,   adjust9= &calorsvq_   &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp10=&vlfssvq_  ,  adjust10=&calorsvq_   &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  );
         Run;

%res(analysis=unit.adj0, catlin=T,   
         data=pool, event=can_obese, time=tobeca, 
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         agevar=agecon,
         strata=agecon cohort intervalp, 
         labels=F, 
         exp1=lcs_sd   ,    adjust1=  ,
         exp2=alcs_sd  ,    adjust2=  ,
         exp3=vlcs_sd  ,    adjust3=  ,
         exp4=ulcs_sd  ,    adjust4=  ,
         exp5=hlcs_sd  ,    adjust5=  , 
         exp6=lfs_sd   ,    adjust6=  ,
         exp7=ulfs_sd  ,    adjust7=  ,
         exp8=hlfs_sd  ,    adjust8=  ,
         exp9=alfs_sd  ,    adjust9=  ,
         exp10=vlfs_sd  ,   adjust10=  );      
         Run;

%res(analysis=unit.adj1,catlin=T,   
         data=pool, event=can_obese, time=tobeca, 
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         agevar=agecon,
         strata=agecon cohort intervalp, 
         labels=F, 
         exp1=lcs_sd   ,   adjust1= &calorsvq_   &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp2=alcs_sd  ,   adjust2= &calorsvq_   &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp3=vlcs_sd  ,   adjust3= &calorsvq_   &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp4=ulcs_sd  ,   adjust4= &calorsvq_   &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp5=hlcs_sd  ,   adjust5= &calorsvq_   &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  , 
         exp6=lfs_sd   ,   adjust6= &calorsvq_   &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp7=ulfs_sd  ,   adjust7= &calorsvq_   &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp8=hlfs_sd  ,   adjust8= &calorsvq_   &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp9=alfs_sd  ,   adjust9= &calorsvq_   &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp10=vlfs_sd  ,  adjust10=&calorsvq_   &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ); 
         Run;

%res(analysis=trend.adj0, catlin=T,   
         data=pool, event=can_obese, time=tobeca, 
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         agevar=agecon,
         strata=agecon cohort intervalp, 
         labels=F, 
         exp1=lcssvtrend  ,     adjust1=  ,
         exp2=alcssvtrend ,     adjust2=  ,
         exp3=vlcssvtrend ,     adjust3=  ,
         exp4=ulcssvtrend ,     adjust4=  ,
         exp5=hlcssvtrend ,     adjust5=  , 
         exp6=lfssvtrend  ,     adjust6=  ,
         exp7=ulfssvtrend ,     adjust7=  ,
         exp8=hlfssvtrend ,     adjust8=  ,
         exp9=alfssvtrend ,     adjust9=  ,
         exp10=vlfssvtrend ,    adjust10=  );
         Run;

%res(analysis=trend.adj1,catlin=T,   
         data=pool, event=can_obese, time=tobeca, 
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         agevar=agecon,
         strata=agecon cohort intervalp, 
         labels=F, 
         exp1=lcssvtrend  ,    adjust1= &calorsvq_   &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp2=alcssvtrend ,    adjust2= &calorsvq_   &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp3=vlcssvtrend ,    adjust3= &calorsvq_   &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp4=ulcssvtrend ,    adjust4= &calorsvq_   &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp5=hlcssvtrend ,    adjust5= &calorsvq_   &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  , 
         exp6=lfssvtrend  ,    adjust6= &calorsvq_   &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp7=ulfssvtrend ,    adjust7= &calorsvq_   &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp8=hlfssvtrend ,    adjust8= &calorsvq_   &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp9=alfssvtrend ,    adjust9= &calorsvq_   &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp10=vlfssvtrend ,   adjust10=&calorsvq_   &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  );
         Run;

%include 'macro/resout.sas';
/************************ HR *****************************************/

%resout(outfile=table/tableS20.3, dbms=tab, opt=replace);                                                                                                                             
endsas;