/*****************************************************
Purpose:       low carb, low fat, and risk of post-dx survival outcomes
Study design:  prospective cohort (NHS: 1986-2022.6, NHSII: 1991-2023.6, HPFS: 1986-2024.1)

Contents:      Figure S4-S5 --- PY and HR (Outcome: Cause-specific Survival; Exposure: Post-dx cumulative averaged diet)
Date:          10/2026
******************************************************/


/************************ Call in Cleaned Datasets**********************************/
%include '/lowcarb_cancer/revision/data_postdx_brc.sas';


/************************ HR *****************************************/

%include '/macro/res_post.sas';

%res(analysis=deadbrc.adj0, catlin=T,   
         data=pool_timevarying, event=dead_brc,  
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         tbegtend=T,
         tbeg=time_begin,
         tend=time_end,
         strata=  cohort intervalp, 
         labels=T, 
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

%res(analysis=deadbrc.adj1,catlin=T,   
         data=pool_timevarying, event=dead_brc,  
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         tbegtend=T,
         tbeg=time_begin,
         tend=time_end,
         strata=  cohort intervalp, 
         labels=T, 
         exp1=&lcssvq_   ,   adjust1= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp2=&alcssvq_  ,   adjust2= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp3=&vlcssvq_  ,   adjust3= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp4=&ulcssvq_  ,   adjust4= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp5=&hlcssvq_  ,   adjust5= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  , 
         exp6=&lfssvq_   ,   adjust6= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp7=&ulfssvq_  ,   adjust7= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp8=&hlfssvq_  ,   adjust8= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp9=&alfssvq_  ,   adjust9= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp10=&vlfssvq_  ,  adjust10=age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  );
         Run;

%res(analysis=deadbrc_unit.adj0, catlin=T,   
         data=pool_timevarying, event=dead_brc,  
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         tbegtend=T,
         tbeg=time_begin,
         tend=time_end,
         strata=  cohort intervalp, 
         labels=T, 
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

%res(analysis=deadbrc_unit.adj1,catlin=T,   
         data=pool_timevarying, event=dead_brc,  
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         tbegtend=T,
         tbeg=time_begin,
         tend=time_end,
         strata=  cohort intervalp, 
         labels=T, 
         exp1=lcs_sd   ,   adjust1= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp2=alcs_sd  ,   adjust2= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp3=vlcs_sd  ,   adjust3= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp4=ulcs_sd  ,   adjust4= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp5=hlcs_sd  ,   adjust5= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  , 
         exp6=lfs_sd   ,   adjust6= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp7=ulfs_sd  ,   adjust7= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp8=hlfs_sd  ,   adjust8= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp9=alfs_sd  ,   adjust9= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp10=vlfs_sd  ,  adjust10=age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ); 
         Run;

%res(analysis=deadbrc_trend.adj0, catlin=T,   
         data=pool_timevarying, event=dead_brc,  
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         tbegtend=T,
         tbeg=time_begin,
         tend=time_end,
         strata=  cohort intervalp, 
         labels=T, 
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

%res(analysis=deadbrc_trend.adj1,catlin=T,   
         data=pool_timevarying, event=dead_brc,  
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         tbegtend=T,
         tbeg=time_begin,
         tend=time_end,
         strata=  cohort intervalp, 
         labels=T, 
         exp1=lcssvtrend  ,    adjust1= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp2=alcssvtrend ,    adjust2= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp3=vlcssvtrend ,    adjust3= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp4=ulcssvtrend ,    adjust4= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp5=hlcssvtrend ,    adjust5= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  , 
         exp6=lfssvtrend  ,    adjust6= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp7=ulfssvtrend ,    adjust7= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp8=hlfssvtrend ,    adjust8= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp9=alfssvtrend ,    adjust9= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp10=vlfssvtrend ,   adjust10=age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  );
         Run;



/************************ HR for CVD-Specific Death *****************************************/

%res(analysis=deadcvd.adj0, catlin=T,   
         data=pool_timevarying, event=dead_cvd,  
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         tbegtend=T,
         tbeg=time_begin,
         tend=time_end,
         strata=  cohort intervalp, 
         labels=T, 
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

%res(analysis=deadcvd.adj1,catlin=T,   
         data=pool_timevarying, event=dead_cvd,  
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         tbegtend=T,
         tbeg=time_begin,
         tend=time_end,
         strata=  cohort intervalp, 
         labels=T, 
         exp1=&lcssvq_   ,   adjust1= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp2=&alcssvq_  ,   adjust2= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp3=&vlcssvq_  ,   adjust3= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp4=&ulcssvq_  ,   adjust4= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp5=&hlcssvq_  ,   adjust5= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  , 
         exp6=&lfssvq_   ,   adjust6= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp7=&ulfssvq_  ,   adjust7= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp8=&hlfssvq_  ,   adjust8= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp9=&alfssvq_  ,   adjust9= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp10=&vlfssvq_  ,  adjust10=age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  );
         Run;

%res(analysis=deadcvd_unit.adj0, catlin=T,   
         data=pool_timevarying, event=dead_cvd,  
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         tbegtend=T,
         tbeg=time_begin,
         tend=time_end,
         strata=  cohort intervalp, 
         labels=T, 
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

%res(analysis=deadcvd_unit.adj1,catlin=T,   
         data=pool_timevarying, event=dead_cvd,  
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         tbegtend=T,
         tbeg=time_begin,
         tend=time_end,
         strata=  cohort intervalp, 
         labels=T, 
         exp1=lcs_sd   ,   adjust1= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp2=alcs_sd  ,   adjust2= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp3=vlcs_sd  ,   adjust3= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp4=ulcs_sd  ,   adjust4= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp5=hlcs_sd  ,   adjust5= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  , 
         exp6=lfs_sd   ,   adjust6= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp7=ulfs_sd  ,   adjust7= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp8=hlfs_sd  ,   adjust8= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp9=alfs_sd  ,   adjust9= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp10=vlfs_sd  ,  adjust10=age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ); 
         Run;

%res(analysis=deadcvd_trend.adj0, catlin=T,   
         data=pool_timevarying, event=dead_cvd,  
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         tbegtend=T,
         tbeg=time_begin,
         tend=time_end,
         strata=  cohort intervalp, 
         labels=T, 
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

%res(analysis=deadcvd_trend.adj1,catlin=T,   
         data=pool_timevarying, event=dead_cvd,  
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         tbegtend=T,
         tbeg=time_begin,
         tend=time_end,
         strata=  cohort intervalp, 
         labels=T, 
         exp1=lcssvtrend  ,    adjust1= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp2=alcssvtrend ,    adjust2= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp3=vlcssvtrend ,    adjust3= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp4=ulcssvtrend ,    adjust4= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp5=hlcssvtrend ,    adjust5= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  , 
         exp6=lfssvtrend  ,    adjust6= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp7=ulfssvtrend ,    adjust7= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp8=hlfssvtrend ,    adjust8= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp9=alfssvtrend ,    adjust9= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp10=vlfssvtrend ,   adjust10=age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  );
         Run;

/************************ HR for Cancer-Specific Death *****************************************/

%res(analysis=deadcan.adj0, catlin=T,   
         data=pool_timevarying, event=dead_can,  
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         tbegtend=T,
         tbeg=time_begin,
         tend=time_end,
         strata=  cohort intervalp, 
         labels=T, 
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

%res(analysis=deadcan.adj1,catlin=T,   
         data=pool_timevarying, event=dead_can,  
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         tbegtend=T,
         tbeg=time_begin,
         tend=time_end,
         strata=  cohort intervalp, 
         labels=T, 
         exp1=&lcssvq_   ,   adjust1= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp2=&alcssvq_  ,   adjust2= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp3=&vlcssvq_  ,   adjust3= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp4=&ulcssvq_  ,   adjust4= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp5=&hlcssvq_  ,   adjust5= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  , 
         exp6=&lfssvq_   ,   adjust6= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp7=&ulfssvq_  ,   adjust7= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp8=&hlfssvq_  ,   adjust8= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp9=&alfssvq_  ,   adjust9= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp10=&vlfssvq_  ,  adjust10=age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  );
         Run;

%res(analysis=deadcan_unit.adj0, catlin=T,   
         data=pool_timevarying, event=dead_can,  
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         tbegtend=T,
         tbeg=time_begin,
         tend=time_end,
         strata=  cohort intervalp, 
         labels=T, 
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

%res(analysis=deadcan_unit.adj1,catlin=T,   
         data=pool_timevarying, event=dead_can,  
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         tbegtend=T,
         tbeg=time_begin,
         tend=time_end,
         strata=  cohort intervalp, 
         labels=T, 
         exp1=lcs_sd   ,   adjust1= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp2=alcs_sd  ,   adjust2= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp3=vlcs_sd  ,   adjust3= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp4=ulcs_sd  ,   adjust4= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp5=hlcs_sd  ,   adjust5= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  , 
         exp6=lfs_sd   ,   adjust6= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp7=ulfs_sd  ,   adjust7= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp8=hlfs_sd  ,   adjust8= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp9=alfs_sd  ,   adjust9= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp10=vlfs_sd  ,  adjust10=age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ); 
         Run;

%res(analysis=deadcan_trend.adj0, catlin=T,   
         data=pool_timevarying, event=dead_can,  
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         tbegtend=T,
         tbeg=time_begin,
         tend=time_end,
         strata=  cohort intervalp, 
         labels=T, 
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

%res(analysis=deadcan_trend.adj1,catlin=T,   
         data=pool_timevarying, event=dead_can,  
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         tbegtend=T,
         tbeg=time_begin,
         tend=time_end,
         strata=  cohort intervalp, 
         labels=T, 
         exp1=lcssvtrend  ,    adjust1= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp2=alcssvtrend ,    adjust2= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp3=vlcssvtrend ,    adjust3= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp4=ulcssvtrend ,    adjust4= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp5=hlcssvtrend ,    adjust5= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  , 
         exp6=lfssvtrend  ,    adjust6= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp7=ulfssvtrend ,    adjust7= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp8=hlfssvtrend ,    adjust8= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp9=alfssvtrend ,    adjust9= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp10=vlfssvtrend ,   adjust10=age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_    &pre_aheim_fixq_ white  mvita aspa cafh   pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  );
         Run;

%include '/macro/resout.sas';
/************************ HR *****************************************/

%resout(outfile=table/cumLCD_casue_specifc_dead, dbms=tab, opt=replace);                                                                                                                             
endsas;



