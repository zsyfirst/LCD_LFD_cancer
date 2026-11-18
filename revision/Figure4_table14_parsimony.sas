/*****************************************************
Purpose:       low carb, low fat, and risk of post-dx survival outcomes
Study design:  prospective cohort (NHS: 1986-2022.6, NHSII: 1991-2023.6, HPFS: 1986-2024.1)

Contents:      Figure 4, Table S14 --- PY and HR (Outcome: All-cause mortality; Exposure: Post-dx cumulative averaged diet)
revision:      Parsimonious model
Date:          10/2026
******************************************************/


/************************ Call in Cleaned Datasets**********************************/
%include '/lowcarb_cancer/revision/data_postdx_brc.sas';


/************************ HR *****************************************/

%include '/macro/res_post.sas';

%res(analysis=avestop.adj1,catlin=T,   
         data=pool_timevarying, event=dead,  
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         tbegtend=T,
         tbeg=time_begin,
         tend=time_end,
         strata= agecon cohort intervalp, 
         labels=T, 
         exp1=&lcssvq_   ,   adjust1= age_dx &calorsvq_ &smkstaa_  &bmid9_ &stage_ &est_status2_  her2_status  &chemotx_ &hormtx_ &radiotx_ ,
         exp2=&alcssvq_  ,   adjust2= age_dx &calorsvq_ &smkstaa_  &bmid9_ &stage_ &est_status2_  her2_status  &chemotx_ &hormtx_ &radiotx_  ,
         exp3=&vlcssvq_  ,   adjust3= age_dx &calorsvq_ &smkstaa_  &bmid9_ &stage_ &est_status2_  her2_status  &chemotx_ &hormtx_ &radiotx_  ,
         exp4=&ulcssvq_  ,   adjust4= age_dx &calorsvq_ &smkstaa_  &bmid9_ &stage_ &est_status2_  her2_status  &chemotx_ &hormtx_ &radiotx_  ,
         exp5=&hlcssvq_  ,   adjust5= age_dx &calorsvq_ &smkstaa_  &bmid9_ &stage_ &est_status2_  her2_status  &chemotx_ &hormtx_ &radiotx_  , 
         exp6=&lfssvq_   ,   adjust6= age_dx &calorsvq_ &smkstaa_  &bmid9_ &stage_ &est_status2_  her2_status  &chemotx_ &hormtx_ &radiotx_  ,
         exp7=&ulfssvq_  ,   adjust7= age_dx &calorsvq_ &smkstaa_  &bmid9_ &stage_ &est_status2_  her2_status  &chemotx_ &hormtx_ &radiotx_  ,
         exp8=&hlfssvq_  ,   adjust8= age_dx &calorsvq_ &smkstaa_  &bmid9_ &stage_ &est_status2_  her2_status  &chemotx_ &hormtx_ &radiotx_  ,
         exp9=&alfssvq_  ,   adjust9= age_dx &calorsvq_ &smkstaa_  &bmid9_ &stage_ &est_status2_  her2_status  &chemotx_ &hormtx_ &radiotx_  ,
         exp10=&vlfssvq_  ,  adjust10=age_dx &calorsvq_ &smkstaa_  &bmid9_ &stage_ &est_status2_  her2_status  &chemotx_ &hormtx_ &radiotx_  );
         Run;

%res(analysis=unit.adj1,catlin=T,   
         data=pool_timevarying, event=dead,  
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         tbegtend=T,
         tbeg=time_begin,
         tend=time_end,
         strata= agecon cohort intervalp, 
         labels=T, 
         exp1=lcs_sd   ,   adjust1= age_dx &calorsvq_ &smkstaa_  &bmid9_ &stage_ &est_status2_  her2_status  &chemotx_ &hormtx_ &radiotx_  ,
         exp2=alcs_sd  ,   adjust2= age_dx &calorsvq_ &smkstaa_  &bmid9_ &stage_ &est_status2_  her2_status  &chemotx_ &hormtx_ &radiotx_  ,
         exp3=vlcs_sd  ,   adjust3= age_dx &calorsvq_ &smkstaa_  &bmid9_ &stage_ &est_status2_  her2_status  &chemotx_ &hormtx_ &radiotx_  ,
         exp4=ulcs_sd  ,   adjust4= age_dx &calorsvq_ &smkstaa_  &bmid9_ &stage_ &est_status2_  her2_status  &chemotx_ &hormtx_ &radiotx_  ,
         exp5=hlcs_sd  ,   adjust5= age_dx &calorsvq_ &smkstaa_  &bmid9_ &stage_ &est_status2_  her2_status  &chemotx_ &hormtx_ &radiotx_  , 
         exp6=lfs_sd   ,   adjust6= age_dx &calorsvq_ &smkstaa_  &bmid9_ &stage_ &est_status2_  her2_status  &chemotx_ &hormtx_ &radiotx_  ,
         exp7=ulfs_sd  ,   adjust7= age_dx &calorsvq_ &smkstaa_  &bmid9_ &stage_ &est_status2_  her2_status  &chemotx_ &hormtx_ &radiotx_  ,
         exp8=hlfs_sd  ,   adjust8= age_dx &calorsvq_ &smkstaa_  &bmid9_ &stage_ &est_status2_  her2_status  &chemotx_ &hormtx_ &radiotx_  ,
         exp9=alfs_sd  ,   adjust9= age_dx &calorsvq_ &smkstaa_  &bmid9_ &stage_ &est_status2_  her2_status  &chemotx_ &hormtx_ &radiotx_  ,
         exp10=vlfs_sd  ,  adjust10=age_dx &calorsvq_ &smkstaa_  &bmid9_ &stage_ &est_status2_  her2_status  &chemotx_ &hormtx_ &radiotx_  ); 
         Run;

%res(analysis=trend.adj1,catlin=T,   
         data=pool_timevarying, event=dead,  
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         tbegtend=T,
         tbeg=time_begin,
         tend=time_end,
         strata= agecon cohort intervalp, 
         labels=T, 
         exp1=lcssvtrend  ,    adjust1= age_dx &calorsvq_ &smkstaa_  &bmid9_ &stage_ &est_status2_  her2_status  &chemotx_ &hormtx_ &radiotx_  ,
         exp2=alcssvtrend ,    adjust2= age_dx &calorsvq_ &smkstaa_  &bmid9_ &stage_ &est_status2_  her2_status  &chemotx_ &hormtx_ &radiotx_  ,
         exp3=vlcssvtrend ,    adjust3= age_dx &calorsvq_ &smkstaa_  &bmid9_ &stage_ &est_status2_  her2_status  &chemotx_ &hormtx_ &radiotx_  ,
         exp4=ulcssvtrend ,    adjust4= age_dx &calorsvq_ &smkstaa_  &bmid9_ &stage_ &est_status2_  her2_status  &chemotx_ &hormtx_ &radiotx_  ,
         exp5=hlcssvtrend ,    adjust5= age_dx &calorsvq_ &smkstaa_  &bmid9_ &stage_ &est_status2_  her2_status  &chemotx_ &hormtx_ &radiotx_  , 
         exp6=lfssvtrend  ,    adjust6= age_dx &calorsvq_ &smkstaa_  &bmid9_ &stage_ &est_status2_  her2_status  &chemotx_ &hormtx_ &radiotx_  ,
         exp7=ulfssvtrend ,    adjust7= age_dx &calorsvq_ &smkstaa_  &bmid9_ &stage_ &est_status2_  her2_status  &chemotx_ &hormtx_ &radiotx_  ,
         exp8=hlfssvtrend ,    adjust8= age_dx &calorsvq_ &smkstaa_  &bmid9_ &stage_ &est_status2_  her2_status  &chemotx_ &hormtx_ &radiotx_  ,
         exp9=alfssvtrend ,    adjust9= age_dx &calorsvq_ &smkstaa_  &bmid9_ &stage_ &est_status2_  her2_status  &chemotx_ &hormtx_ &radiotx_  ,
         exp10=vlfssvtrend ,   adjust10=age_dx &calorsvq_ &smkstaa_  &bmid9_ &stage_ &est_status2_  her2_status  &chemotx_ &hormtx_ &radiotx_  );
         Run;


%include '/macro/resout.sas';
/************************ HR *****************************************/

%resout(outfile=table/cumLCD_dead_parsimony, dbms=tab, opt=replace);                                                                                                                             
endsas;

