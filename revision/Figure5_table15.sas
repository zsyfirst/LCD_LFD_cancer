/*****************************************************
Purpose:       post-dx low carb, low fat, and mortality
Study design:  prospective cohort (NHS: 1986-2022.6, NHSII: 1991-2023.6, HPFS: 1986-2024.1)

Contents:      Figure 5, Table S15 --- PY and HR (Outcome: Overall Survival; Exposure: Post-dx cumulative averaged diet)
Revision:      Model1: age, Model2: multiple covariates, Model3: add nSES
Date:          10/2026
******************************************************/


/************************ Call in Cleaned Datasets**********************************/
%include '/lowcarb_cancer/revision/data_postdx_crc.sas';

/************************ CASE-PERSON YEARS*****************************************/
%macro casepy_postdx(data= , exp= , case= , cohortvar= cohort,
    timevar  = person_time, /* months since dx in interval */
    out      = casepy_&exp._&case 
);
    /* 1. Cohort-specific case/PY (NHS1, NHS2, HPFS) */
    proc sql;
        create table _cp_cohort as
        select 
            case 
                when &cohortvar. = 1 then "NHS1"
                when &cohortvar. = 2 then "NHS2"
                when &cohortvar. = 3 then "HPFS"
            end as cohort length=5,
            &exp.,
            sum(&case.)       as cases,
            sum(&timevar.)/12 as py
        from &data.
        where &exp. ne . and &timevar.>0
        group by &cohortvar., &exp.;
    quit;

    data _cp_cohort;
        set _cp_cohort;
        casepy = cats( put(cases, comma20.)
                      , " / "
                      , put(py, comma20. -L) );
    run;


    /* 2. pool case/PY across cohorts */
    proc sql;
        create table _cp_pool as
        select 
            "Pool"  as cohort length=5,
            &exp.,
            sum(&case.)       as cases,
            sum(&timevar.)/12 as py
        from &data.
        where &exp. ne . and &timevar.>0
        group by &exp.;
    quit;
    data _cp_pool;
        set _cp_pool;
        casepy = cats( put(cases, comma20.) , " / ", put(py, comma20. -L) );
    run;

    data &out.;
        length index $32;
        set _cp_cohort _cp_pool;
        index = "&exp";
    run;

    proc datasets noprint;
        delete _cp_cohort _cp_pool;
    run;
%mend;


%macro process_all_outcomes(data=, outcomes=, exposures=, outdir=);
    
    /* Loop through each outcome */
    %let n_outcomes = %sysfunc(countw(&outcomes));
    %do o = 1 %to &n_outcomes;
        %let outcome = %scan(&outcomes, &o);
        
        /* Loop through each exposure */
        %let n_exp = %sysfunc(countw(&exposures));
        %do e = 1 %to &n_exp;
            %let exp = %scan(&exposures, &e);
            
            /* Calculate case-PY for this exposure-outcome combination */
            %casepy_postdx(data=&data, exp=&exp, case=&outcome);
        %end;
        
        data casepy_&outcome._all;
            set 
            %do e = 1 %to &n_exp;
                %let exp = %scan(&exposures, &e);
                casepy_&exp._&outcome
            %end;
            ;
        run;
        
        proc export data=casepy_&outcome._all
            dbms=csv
            outfile="&outdir./casePY_cumLCD_&outcome..csv"
            replace;
        run;
    %end;
%mend;


%let exposures = lcssvq alcssvq vlcssvq ulcssvq hlcssvq 
                 lfssvq alfssvq vlfssvq ulfssvq hlfssvq;

%let outcomes = dead dead_crc dead_cvd dead_can dead_nocan 
                dead_brc dead_resp dead_oth;

/* Run for time-varying dataset */
%process_all_outcomes(
    data=pool_timevarying,
    outcomes=&outcomes,
    exposures=&exposures,
    outdir=lowcarb_cancer/revision/crc/table
);


/************************ HR *****************************************/
%include 'macro/res_post.sas';

%res(analysis=avestop.adj0, catlin=T,   
         data=pool_timevarying, event=dead,  
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         tbegtend=T,
         tbeg=time_begin,
         tend=time_end,
         strata= agecon cohort intervalp, 
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
         exp1=&lcssvq_   ,   adjust1= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_   ,
         exp2=&alcssvq_  ,   adjust2= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_   ,
         exp3=&vlcssvq_  ,   adjust3= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_   ,
         exp4=&ulcssvq_  ,   adjust4= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_   ,
         exp5=&hlcssvq_  ,   adjust5= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_   , 
         exp6=&lfssvq_   ,   adjust6= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_   ,
         exp7=&ulfssvq_  ,   adjust7= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_   ,
         exp8=&hlfssvq_  ,   adjust8= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_   ,
         exp9=&alfssvq_  ,   adjust9= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_   ,
         exp10=&vlfssvq_  ,  adjust10=age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_   );
         Run;

%res(analysis=avestop.adj2,catlin=T,   
         data=pool_timevarying, event=dead,  
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         tbegtend=T,
         tbeg=time_begin,
         tend=time_end,
         strata= agecon cohort intervalp, 
         labels=T, 
         exp1=&lcssvq_   ,   adjust1= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_  nSES  ,
         exp2=&alcssvq_  ,   adjust2= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_  nSES  ,
         exp3=&vlcssvq_  ,   adjust3= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_  nSES  ,
         exp4=&ulcssvq_  ,   adjust4= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_  nSES  ,
         exp5=&hlcssvq_  ,   adjust5= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_  nSES  , 
         exp6=&lfssvq_   ,   adjust6= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_  nSES  ,
         exp7=&ulfssvq_  ,   adjust7= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_  nSES  ,
         exp8=&hlfssvq_  ,   adjust8= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_  nSES  ,
         exp9=&alfssvq_  ,   adjust9= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_  nSES  ,
         exp10=&vlfssvq_  ,  adjust10=age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_  nSES  );
         Run;

%res(analysis=unit.adj0, catlin=T,   
         data=pool_timevarying, event=dead,  
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         tbegtend=T,
         tbeg=time_begin,
         tend=time_end,
         strata= agecon cohort intervalp, 
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
         exp1=lcs_sd   ,   adjust1= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_   ,
         exp2=alcs_sd  ,   adjust2= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_   ,
         exp3=vlcs_sd  ,   adjust3= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_   ,
         exp4=ulcs_sd  ,   adjust4= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_   ,
         exp5=hlcs_sd  ,   adjust5= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_   , 
         exp6=lfs_sd   ,   adjust6= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_   ,
         exp7=ulfs_sd  ,   adjust7= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_   ,
         exp8=hlfs_sd  ,   adjust8= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_   ,
         exp9=alfs_sd  ,   adjust9= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_   ,
         exp10=vlfs_sd  ,  adjust10=age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_   ); 
         Run;

%res(analysis=unit.adj2,catlin=T,   
         data=pool_timevarying, event=dead,  
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         tbegtend=T,
         tbeg=time_begin,
         tend=time_end,
         strata= agecon cohort intervalp, 
         labels=T, 
         exp1=lcs_sd   ,   adjust1= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_  nSES  ,
         exp2=alcs_sd  ,   adjust2= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_  nSES  ,
         exp3=vlcs_sd  ,   adjust3= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_  nSES  ,
         exp4=ulcs_sd  ,   adjust4= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_  nSES  ,
         exp5=hlcs_sd  ,   adjust5= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_  nSES  , 
         exp6=lfs_sd   ,   adjust6= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_  nSES  ,
         exp7=ulfs_sd  ,   adjust7= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_  nSES  ,
         exp8=hlfs_sd  ,   adjust8= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_  nSES  ,
         exp9=alfs_sd  ,   adjust9= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_  nSES  ,
         exp10=vlfs_sd  ,  adjust10=age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_  nSES  ); 
         Run;

%res(analysis=trend.adj0, catlin=T,   
         data=pool_timevarying, event=dead,  
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         tbegtend=T,
         tbeg=time_begin,
         tend=time_end,
         strata= agecon cohort intervalp, 
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
         exp1=lcssvtrend  ,    adjust1= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_   ,
         exp2=alcssvtrend ,    adjust2= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_   ,
         exp3=vlcssvtrend ,    adjust3= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_   ,
         exp4=ulcssvtrend ,    adjust4= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_   ,
         exp5=hlcssvtrend ,    adjust5= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_   , 
         exp6=lfssvtrend  ,    adjust6= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_   ,
         exp7=ulfssvtrend ,    adjust7= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_   ,
         exp8=hlfssvtrend ,    adjust8= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_   ,
         exp9=alfssvtrend ,    adjust9= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_   ,
         exp10=vlfssvtrend ,   adjust10=age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_   );
         Run;

%res(analysis=trend.adj2,catlin=T,   
         data=pool_timevarying, event=dead,  
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         tbegtend=T,
         tbeg=time_begin,
         tend=time_end,
         strata= agecon cohort intervalp, 
         labels=T, 
         exp1=lcssvtrend  ,    adjust1= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_  nSES  ,
         exp2=alcssvtrend ,    adjust2= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_  nSES  ,
         exp3=vlcssvtrend ,    adjust3= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_  nSES  ,
         exp4=ulcssvtrend ,    adjust4= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_  nSES  ,
         exp5=hlcssvtrend ,    adjust5= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_  nSES  , 
         exp6=lfssvtrend  ,    adjust6= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_  nSES  ,
         exp7=ulfssvtrend ,    adjust7= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_  nSES  ,
         exp8=hlfssvtrend ,    adjust8= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_  nSES  ,
         exp9=alfssvtrend ,    adjust9= age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_  nSES  ,
         exp10=vlfssvtrend ,   adjust10=age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_   &pre_aheim_fixq_ white mvita aspa cafh   &stage_cr_ &grade_cr_ crcduration &subsite_cr_  nSES  );
         Run;

%include 'macro/resout.sas';
/************************ HR *****************************************/

%resout(outfile=table/cumLCD_dead_postcrc, dbms=tab, opt=replace);                                                                                                                             
endsas;