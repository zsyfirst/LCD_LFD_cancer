/*****************************************************
Purpose:       Different Time Windows of Low-Carbohydrate and Low-Fat Diet Scores with All-Cause Mortality Among Breast Cancer Survivors
Study design:  prospective cohort (NHS: 1986-2022.6, NHSII: 1991-2023.6)
Contents:      Figure S10, Table S21 -- HR (Overall survival as outcome, 
                               Exposure: 1. Post-dx Single Diet)
Date:          10/2026
******************************************************/


/************************ Call in Cleaned Datasets**********************************/
%include '/lowcarb_cancer/revision/time_window/data_postdx_brc.sas';


/************************ HR *****************************************/

%include 'macro/res_post.sas';


/****************** Post-dx single diet ***********************/
%res(analysis=post_avestop.adj0, catlin=T,   
         data=pool_singlepost, event=dead, 
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         tbegtend=T,
         tbeg=time_begin,
         tend=time_end,
         strata= agecon cohort intervalp, 
         labels=T, 
         exp1=&lcs_postq_   ,    adjust1=   ,
         exp2=&alcs_postq_  ,    adjust2=   ,
         exp3=&vlcs_postq_  ,    adjust3=   ,
         exp4=&ulcs_postq_  ,    adjust4=   ,
         exp5=&hlcs_postq_  ,    adjust5=   , 
         exp6=&lfs_postq_   ,    adjust6=   ,
         exp7=&ulfs_postq_  ,    adjust7=   ,
         exp8=&hlfs_postq_  ,    adjust8=   ,
         exp9=&alfs_postq_  ,    adjust9=   ,
         exp10=&vlfs_postq_  ,    adjust10=   );
         Run;

%res(analysis=post_avestop.adj1,catlin=T,   
         data=pool_singlepost, event=dead, 
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         tbegtend=T,
         tbeg=time_begin,
         tend=time_end,
         strata= agecon cohort intervalp, 
         labels=T, 
         exp1=&lcs_postq_   ,   adjust1= age_dx &calor_postq_ &smkstaa_ &alco_postq_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa  cafh  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ /* &surg_ */ nSES ,
         exp2=&alcs_postq_  ,   adjust2= age_dx &calor_postq_ &smkstaa_ &alco_postq_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa  cafh  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ /* &surg_ */ nSES ,
         exp3=&vlcs_postq_  ,   adjust3= age_dx &calor_postq_ &smkstaa_ &alco_postq_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa  cafh  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ /* &surg_ */ nSES ,
         exp4=&ulcs_postq_  ,   adjust4= age_dx &calor_postq_ &smkstaa_ &alco_postq_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa  cafh  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ /* &surg_ */ nSES ,
         exp5=&hlcs_postq_  ,   adjust5= age_dx &calor_postq_ &smkstaa_ &alco_postq_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa  cafh  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ /* &surg_ */ nSES , 
         exp6=&lfs_postq_   ,   adjust6= age_dx &calor_postq_ &smkstaa_ &alco_postq_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa  cafh  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ /* &surg_ */ nSES ,
         exp7=&ulfs_postq_  ,   adjust7= age_dx &calor_postq_ &smkstaa_ &alco_postq_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa  cafh  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ /* &surg_ */ nSES ,
         exp8=&hlfs_postq_  ,   adjust8= age_dx &calor_postq_ &smkstaa_ &alco_postq_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa  cafh  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ /* &surg_ */ nSES ,
         exp9=&alfs_postq_  ,   adjust9= age_dx &calor_postq_ &smkstaa_ &alco_postq_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa  cafh  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ /* &surg_ */ nSES ,
         exp10=&vlfs_postq_  ,  adjust10=age_dx &calor_postq_ &smkstaa_ &alco_postq_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa  cafh  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ /* &surg_ */ nSES );
         Run;

%res(analysis=post_unit.adj0, catlin=T,   
         data=pool_singlepost, event=dead,  
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         tbegtend=T,
         tbeg=time_begin,
         tend=time_end,
         strata= agecon cohort intervalp, 
         labels=T, 
         exp1=lcs_sd_p   ,    adjust1=  ,
         exp2=alcs_sd_p  ,    adjust2=  ,
         exp3=vlcs_sd_p  ,    adjust3=  ,
         exp4=ulcs_sd_p  ,    adjust4=  ,
         exp5=hlcs_sd_p  ,    adjust5=  , 
         exp6=lfs_sd_p   ,    adjust6=  ,
         exp7=ulfs_sd_p  ,    adjust7=  ,
         exp8=hlfs_sd_p  ,    adjust8=  ,
         exp9=alfs_sd_p  ,    adjust9=  ,
         exp10=vlfs_sd_p  ,   adjust10=  );      
         Run;

%res(analysis=post_unit.adj1,catlin=T,   
         data=pool_singlepost, event=dead, 
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         tbegtend=T,
         tbeg=time_begin,
         tend=time_end,
         strata= agecon cohort intervalp, 
         labels=T, 
         exp1=lcs_sd_p   ,   adjust1= age_dx &calor_postq_ &smkstaa_ &alco_postq_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa  cafh  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ /* &surg_ */ nSES ,
         exp2=alcs_sd_p  ,   adjust2= age_dx &calor_postq_ &smkstaa_ &alco_postq_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa  cafh  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ /* &surg_ */ nSES ,
         exp3=vlcs_sd_p  ,   adjust3= age_dx &calor_postq_ &smkstaa_ &alco_postq_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa  cafh  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ /* &surg_ */ nSES ,
         exp4=ulcs_sd_p  ,   adjust4= age_dx &calor_postq_ &smkstaa_ &alco_postq_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa  cafh  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ /* &surg_ */ nSES ,
         exp5=hlcs_sd_p  ,   adjust5= age_dx &calor_postq_ &smkstaa_ &alco_postq_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa  cafh  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ /* &surg_ */ nSES , 
         exp6=lfs_sd_p   ,   adjust6= age_dx &calor_postq_ &smkstaa_ &alco_postq_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa  cafh  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ /* &surg_ */ nSES ,
         exp7=ulfs_sd_p  ,   adjust7= age_dx &calor_postq_ &smkstaa_ &alco_postq_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa  cafh  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ /* &surg_ */ nSES ,
         exp8=hlfs_sd_p  ,   adjust8= age_dx &calor_postq_ &smkstaa_ &alco_postq_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa  cafh  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ /* &surg_ */ nSES ,
         exp9=alfs_sd_p  ,   adjust9= age_dx &calor_postq_ &smkstaa_ &alco_postq_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa  cafh  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ /* &surg_ */ nSES ,
         exp10=vlfs_sd_p  ,  adjust10=age_dx &calor_postq_ &smkstaa_ &alco_postq_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa  cafh  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ /* &surg_ */ nSES ); 
         Run;

%res(analysis=post_trend.adj0, catlin=T,   
         data=pool_singlepost, event=dead, 
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         tbegtend=T,
         tbeg=time_begin,
         tend=time_end,
         strata= agecon cohort intervalp, 
         labels=T, 
         exp1=lcsposttrend  ,     adjust1=  ,
         exp2=alcsposttrend ,     adjust2=  ,
         exp3=vlcsposttrend ,     adjust3=  ,
         exp4=ulcsposttrend ,     adjust4=  ,
         exp5=hlcsposttrend ,     adjust5=  , 
         exp6=lfsposttrend  ,     adjust6=  ,
         exp7=ulfsposttrend ,     adjust7=  ,
         exp8=hlfsposttrend ,     adjust8=  ,
         exp9=alfsposttrend ,     adjust9=  ,
         exp10=vlfsposttrend ,    adjust10=  );
         Run;

%res(analysis=post_trend.adj1,catlin=T,   
         data=pool_singlepost, event=dead,  
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         tbegtend=T,
         tbeg=time_begin,
         tend=time_end,
         strata= agecon cohort intervalp, 
         labels=T, 
         exp1=lcsposttrend  ,    adjust1= age_dx &calor_postq_ &smkstaa_ &alco_postq_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa  cafh  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ /* &surg_ */ nSES  ,
         exp2=alcsposttrend ,    adjust2= age_dx &calor_postq_ &smkstaa_ &alco_postq_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa  cafh  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ /* &surg_ */ nSES  ,
         exp3=vlcsposttrend ,    adjust3= age_dx &calor_postq_ &smkstaa_ &alco_postq_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa  cafh  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ /* &surg_ */ nSES  ,
         exp4=ulcsposttrend ,    adjust4= age_dx &calor_postq_ &smkstaa_ &alco_postq_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa  cafh  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ /* &surg_ */ nSES ,
         exp5=hlcsposttrend ,    adjust5= age_dx &calor_postq_ &smkstaa_ &alco_postq_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa  cafh  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ /* &surg_ */ nSES , 
         exp6=lfsposttrend  ,    adjust6= age_dx &calor_postq_ &smkstaa_ &alco_postq_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa  cafh  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ /* &surg_ */ nSES ,
         exp7=ulfsposttrend ,    adjust7= age_dx &calor_postq_ &smkstaa_ &alco_postq_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa  cafh  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ /* &surg_ */ nSES ,
         exp8=hlfsposttrend ,    adjust8= age_dx &calor_postq_ &smkstaa_ &alco_postq_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa  cafh  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ /* &surg_ */ nSES ,
         exp9=alfsposttrend ,    adjust9= age_dx &calor_postq_ &smkstaa_ &alco_postq_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa  cafh  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ /* &surg_ */ nSES ,
         exp10=vlfsposttrend ,   adjust10=age_dx &calor_postq_ &smkstaa_ &alco_postq_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa  cafh  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ /* &surg_ */ nSES );
         Run;


%include 'macro/resout.sas';
/************************ HR *****************************************/

%resout(outfile=tableS21_post_single_dead, dbms=tab, opt=replace);                                                                                                                             
endsas;



