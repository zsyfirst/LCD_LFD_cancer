/*****************************************************
Purpose:       Different Time Windows of Low-Carbohydrate and Low-Fat Diet Scores with All-Cause Mortality Among Breast Cancer Survivors
Study design:  prospective cohort (NHS: 1984-2022.6, NHSII: 1991-2023.6)
Contents:      Figure S10, Table S21 -- HR (Overall survival as outcome, 
                               Exposure: 3. Pre-dx Diet)
Date:          10/2026
******************************************************/


/************************ Call in Cleaned Datasets**********************************/
%include 'lowcarb_cancer/revision/time_window/data_postdx_brc.sas';


/************************ HR *****************************************/

%include 'macro/res_post.sas';


/****************** pre-dx single diet ***********************/
%res(analysis=pre_avestop.adj0, catlin=T,   
         data=pool_prediag, event=dead, 
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         tbegtend=T,
         tbeg=time_begin,
         tend=time_end,
         strata=  cohort intervalp, 
         labels=T, 
         exp1=&lcs_preq_   ,    adjust1=   ,
         exp2=&alcs_preq_  ,    adjust2=   ,
         exp3=&vlcs_preq_  ,    adjust3=   ,
         exp4=&ulcs_preq_  ,    adjust4=   ,
         exp5=&hlcs_preq_  ,    adjust5=   , 
         exp6=&lfs_preq_   ,    adjust6=   ,
         exp7=&ulfs_preq_  ,    adjust7=   ,
         exp8=&hlfs_preq_  ,    adjust8=   ,
         exp9=&alfs_preq_  ,    adjust9=   ,
         exp10=&vlfs_preq_  ,    adjust10=   );
         Run;

%res(analysis=pre_avestop.adj1,catlin=T,   
         data=pool_prediag, event=dead, 
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         tbegtend=T,
         tbeg=time_begin,
         tend=time_end,
         strata=  cohort intervalp, 
         labels=T, 
         exp1=&lcs_preq_   ,   adjust1= age_dx &calor_preq_ &smkstp_ &alco_preq_ &preact_ &bmid9_   &pre_aheim_fixq_ white mv mvitp aspp cafh   pmhever ocever &est_status2_ &stage_   her2_status &chemotx_ &hormtx_ &radiotx_ &radiotx_ &surg_ nSES  ,
         exp2=&alcs_preq_  ,   adjust2= age_dx &calor_preq_ &smkstp_ &alco_preq_ &preact_ &bmid9_   &pre_aheim_fixq_ white mv mvitp aspp cafh   pmhever ocever &est_status2_ &stage_   her2_status &chemotx_ &hormtx_ &radiotx_ &radiotx_ &surg_ nSES  ,
         exp3=&vlcs_preq_  ,   adjust3= age_dx &calor_preq_ &smkstp_ &alco_preq_ &preact_ &bmid9_   &pre_aheim_fixq_ white mv mvitp aspp cafh   pmhever ocever &est_status2_ &stage_   her2_status &chemotx_ &hormtx_ &radiotx_ &radiotx_ &surg_ nSES  ,
         exp4=&ulcs_preq_  ,   adjust4= age_dx &calor_preq_ &smkstp_ &alco_preq_ &preact_ &bmid9_   &pre_aheim_fixq_ white mv mvitp aspp cafh   pmhever ocever &est_status2_ &stage_   her2_status &chemotx_ &hormtx_ &radiotx_ &radiotx_ &surg_ nSES  ,
         exp5=&hlcs_preq_  ,   adjust5= age_dx &calor_preq_ &smkstp_ &alco_preq_ &preact_ &bmid9_   &pre_aheim_fixq_ white mv mvitp aspp cafh   pmhever ocever &est_status2_ &stage_   her2_status &chemotx_ &hormtx_ &radiotx_ &radiotx_ &surg_ nSES  , 
         exp6=&lfs_preq_   ,   adjust6= age_dx &calor_preq_ &smkstp_ &alco_preq_ &preact_ &bmid9_   &pre_aheim_fixq_ white mv mvitp aspp cafh   pmhever ocever &est_status2_ &stage_   her2_status &chemotx_ &hormtx_ &radiotx_ &radiotx_ &surg_ nSES  ,
         exp7=&ulfs_preq_  ,   adjust7= age_dx &calor_preq_ &smkstp_ &alco_preq_ &preact_ &bmid9_   &pre_aheim_fixq_ white mv mvitp aspp cafh   pmhever ocever &est_status2_ &stage_   her2_status &chemotx_ &hormtx_ &radiotx_ &radiotx_ &surg_ nSES  ,
         exp8=&hlfs_preq_  ,   adjust8= age_dx &calor_preq_ &smkstp_ &alco_preq_ &preact_ &bmid9_   &pre_aheim_fixq_ white mv mvitp aspp cafh   pmhever ocever &est_status2_ &stage_   her2_status &chemotx_ &hormtx_ &radiotx_ &radiotx_ &surg_ nSES  ,
         exp9=&alfs_preq_  ,   adjust9= age_dx &calor_preq_ &smkstp_ &alco_preq_ &preact_ &bmid9_   &pre_aheim_fixq_ white mv mvitp aspp cafh   pmhever ocever &est_status2_ &stage_   her2_status &chemotx_ &hormtx_ &radiotx_ &radiotx_ &surg_ nSES  ,
         exp10=&vlfs_preq_  ,  adjust10=age_dx &calor_preq_ &smkstp_ &alco_preq_ &preact_ &bmid9_   &pre_aheim_fixq_ white mv mvitp aspp cafh   pmhever ocever &est_status2_ &stage_   her2_status &chemotx_ &hormtx_ &radiotx_ &radiotx_ &surg_ nSES  );
         Run;

%res(analysis=pre_unit.adj0, catlin=T,   
         data=pool_prediag, event=dead,  
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

%res(analysis=pre_unit.adj1,catlin=T,   
         data=pool_prediag, event=dead,  
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         tbegtend=T,
         tbeg=time_begin,
         tend=time_end,
         strata=  cohort intervalp, 
         labels=T, 
         exp1=lcs_sd_p   ,   adjust1= age_dx &calor_preq_ &smkstp_ &alco_preq_ &preact_ &bmid9_   &pre_aheim_fixq_ white mv mvitp aspp cafh   pmhever ocever &est_status2_ &stage_   her2_status &chemotx_ &hormtx_ &radiotx_ &radiotx_ &surg_ nSES  ,
         exp2=alcs_sd_p  ,   adjust2= age_dx &calor_preq_ &smkstp_ &alco_preq_ &preact_ &bmid9_   &pre_aheim_fixq_ white mv mvitp aspp cafh   pmhever ocever &est_status2_ &stage_   her2_status &chemotx_ &hormtx_ &radiotx_ &radiotx_ &surg_ nSES  ,
         exp3=vlcs_sd_p  ,   adjust3= age_dx &calor_preq_ &smkstp_ &alco_preq_ &preact_ &bmid9_   &pre_aheim_fixq_ white mv mvitp aspp cafh   pmhever ocever &est_status2_ &stage_   her2_status &chemotx_ &hormtx_ &radiotx_ &radiotx_ &surg_ nSES  ,
         exp4=ulcs_sd_p  ,   adjust4= age_dx &calor_preq_ &smkstp_ &alco_preq_ &preact_ &bmid9_   &pre_aheim_fixq_ white mv mvitp aspp cafh   pmhever ocever &est_status2_ &stage_   her2_status &chemotx_ &hormtx_ &radiotx_ &radiotx_ &surg_ nSES  ,
         exp5=hlcs_sd_p  ,   adjust5= age_dx &calor_preq_ &smkstp_ &alco_preq_ &preact_ &bmid9_   &pre_aheim_fixq_ white mv mvitp aspp cafh   pmhever ocever &est_status2_ &stage_   her2_status &chemotx_ &hormtx_ &radiotx_ &radiotx_ &surg_ nSES  , 
         exp6=lfs_sd_p   ,   adjust6= age_dx &calor_preq_ &smkstp_ &alco_preq_ &preact_ &bmid9_   &pre_aheim_fixq_ white mv mvitp aspp cafh   pmhever ocever &est_status2_ &stage_   her2_status &chemotx_ &hormtx_ &radiotx_ &radiotx_ &surg_ nSES  ,
         exp7=ulfs_sd_p  ,   adjust7= age_dx &calor_preq_ &smkstp_ &alco_preq_ &preact_ &bmid9_   &pre_aheim_fixq_ white mv mvitp aspp cafh   pmhever ocever &est_status2_ &stage_   her2_status &chemotx_ &hormtx_ &radiotx_ &radiotx_ &surg_ nSES  ,
         exp8=hlfs_sd_p  ,   adjust8= age_dx &calor_preq_ &smkstp_ &alco_preq_ &preact_ &bmid9_   &pre_aheim_fixq_ white mv mvitp aspp cafh   pmhever ocever &est_status2_ &stage_   her2_status &chemotx_ &hormtx_ &radiotx_ &radiotx_ &surg_ nSES  ,
         exp9=alfs_sd_p  ,   adjust9= age_dx &calor_preq_ &smkstp_ &alco_preq_ &preact_ &bmid9_   &pre_aheim_fixq_ white mv mvitp aspp cafh   pmhever ocever &est_status2_ &stage_   her2_status &chemotx_ &hormtx_ &radiotx_ &radiotx_ &surg_ nSES  ,
         exp10=vlfs_sd_p  ,  adjust10=age_dx &calor_preq_ &smkstp_ &alco_preq_ &preact_ &bmid9_   &pre_aheim_fixq_ white mv mvitp aspp cafh   pmhever ocever &est_status2_ &stage_   her2_status &chemotx_ &hormtx_ &radiotx_ &radiotx_ &surg_ nSES  ); 
         Run;

%res(analysis=pre_trend.adj0, catlin=T,   
         data=pool_prediag, event=dead, 
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         tbegtend=T,
         tbeg=time_begin,
         tend=time_end,
         strata=  cohort intervalp, 
         labels=T, 
         exp1=lcspretrend  ,     adjust1=  ,
         exp2=alcspretrend ,     adjust2=  ,
         exp3=vlcspretrend ,     adjust3=  ,
         exp4=ulcspretrend ,     adjust4=  ,
         exp5=hlcspretrend ,     adjust5=  , 
         exp6=lfspretrend  ,     adjust6=  ,
         exp7=ulfspretrend ,     adjust7=  ,
         exp8=hlfspretrend ,     adjust8=  ,
         exp9=alfspretrend ,     adjust9=  ,
         exp10=vlfspretrend ,    adjust10=  );
         Run;

%res(analysis=pre_trend.adj1,catlin=T,   
         data=pool_prediag, event=dead,  
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         tbegtend=T,
         tbeg=time_begin,
         tend=time_end,
         strata=  cohort intervalp, 
         labels=T, 
         exp1=lcspretrend  ,    adjust1= age_dx &calor_preq_ &smkstp_ &alco_preq_ &preact_ &bmid9_   &pre_aheim_fixq_ white mv mvitp aspp cafh   pmhever ocever &est_status2_ &stage_   her2_status &chemotx_ &hormtx_ &radiotx_ &radiotx_ &surg_ nSES   ,
         exp2=alcspretrend ,    adjust2= age_dx &calor_preq_ &smkstp_ &alco_preq_ &preact_ &bmid9_   &pre_aheim_fixq_ white mv mvitp aspp cafh   pmhever ocever &est_status2_ &stage_   her2_status &chemotx_ &hormtx_ &radiotx_ &radiotx_ &surg_ nSES   ,
         exp3=vlcspretrend ,    adjust3= age_dx &calor_preq_ &smkstp_ &alco_preq_ &preact_ &bmid9_   &pre_aheim_fixq_ white mv mvitp aspp cafh   pmhever ocever &est_status2_ &stage_   her2_status &chemotx_ &hormtx_ &radiotx_ &radiotx_ &surg_ nSES   ,
         exp4=ulcspretrend ,    adjust4= age_dx &calor_preq_ &smkstp_ &alco_preq_ &preact_ &bmid9_   &pre_aheim_fixq_ white mv mvitp aspp cafh   pmhever ocever &est_status2_ &stage_   her2_status &chemotx_ &hormtx_ &radiotx_ &radiotx_ &surg_ nSES  ,
         exp5=hlcspretrend ,    adjust5= age_dx &calor_preq_ &smkstp_ &alco_preq_ &preact_ &bmid9_   &pre_aheim_fixq_ white mv mvitp aspp cafh   pmhever ocever &est_status2_ &stage_   her2_status &chemotx_ &hormtx_ &radiotx_ &radiotx_ &surg_ nSES  , 
         exp6=lfspretrend  ,    adjust6= age_dx &calor_preq_ &smkstp_ &alco_preq_ &preact_ &bmid9_   &pre_aheim_fixq_ white mv mvitp aspp cafh   pmhever ocever &est_status2_ &stage_   her2_status &chemotx_ &hormtx_ &radiotx_ &radiotx_ &surg_ nSES  ,
         exp7=ulfspretrend ,    adjust7= age_dx &calor_preq_ &smkstp_ &alco_preq_ &preact_ &bmid9_   &pre_aheim_fixq_ white mv mvitp aspp cafh   pmhever ocever &est_status2_ &stage_   her2_status &chemotx_ &hormtx_ &radiotx_ &radiotx_ &surg_ nSES  ,
         exp8=hlfspretrend ,    adjust8= age_dx &calor_preq_ &smkstp_ &alco_preq_ &preact_ &bmid9_   &pre_aheim_fixq_ white mv mvitp aspp cafh   pmhever ocever &est_status2_ &stage_   her2_status &chemotx_ &hormtx_ &radiotx_ &radiotx_ &surg_ nSES  ,
         exp9=alfspretrend ,    adjust9= age_dx &calor_preq_ &smkstp_ &alco_preq_ &preact_ &bmid9_   &pre_aheim_fixq_ white mv mvitp aspp cafh   pmhever ocever &est_status2_ &stage_   her2_status &chemotx_ &hormtx_ &radiotx_ &radiotx_ &surg_ nSES  ,
         exp10=vlfspretrend ,   adjust10=age_dx &calor_preq_ &smkstp_ &alco_preq_ &preact_ &bmid9_   &pre_aheim_fixq_ white mv mvitp aspp cafh   pmhever ocever &est_status2_ &stage_   her2_status &chemotx_ &hormtx_ &radiotx_ &radiotx_ &surg_ nSES  );
         Run;


%include 'macro/resout.sas';
/************************ HR *****************************************/

%resout(outfile=tableS21_pre_single_dead, dbms=tab, opt=replace);                                                                                                                             
endsas;



