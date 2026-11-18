/*****************************************************
Purpose:       Different Time Windows of Low-Carbohydrate and Low-Fat Diet Scores with All-Cause Mortality Among Breast Cancer Survivors
Study design:  prospective cohort (NHS: 1984-2022.6, NHSII: 1991-2023.6)
Contents:      Figure S10, Table S21 -- HR (Overall survival as outcome, 
                               Exposure: 1. Post-dx Diet (6 months to 4 years)
Date:          10/2026
******************************************************/


/************************ Call in Cleaned Datasets**********************************/
%include 'lowcarb_cancer/revision/data_postdx_brc.sas';



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
        
        /* QC */
        proc contents data=casepy_&outcome._all;
        run;

        proc print data=casepy_&outcome._all(obs=20);
            title "Combined case/PY results for &outcome.";
        run;


        /*----------------------------------------------------------
          Export CSV
        ----------------------------------------------------------*/
        proc export
            data=casepy_&outcome._all
            outfile="&outdir./casePY_cumLCD_&outcome..csv"
            dbms=csv
            replace;
        run;


        %put NOTE: Exported:;
        %put &outdir./casePY_cumLCD_&outcome..csv;

    %end;
%mend;


%let exposures = lcs_crcq alcs_crcq vlcs_crcq ulcs_crcq hlcs_crcq 
                 lfs_crcq alfs_crcq vlfs_crcq ulfs_crcq hlfs_crcq;

%let outcomes = dead dead_brc dead_cvd dead_can dead_nocan 
                dead_resp dead_oth;

/* Run for 6m-4y dataset */
%process_all_outcomes(
    data=pool_postcrc,
    outcomes=&outcomes,
    exposures=&exposures,
    outdir=lowcarb_cancer/revision/brc/time_window
);



/************************ HR *****************************************/

%include 'macro/res_post.sas';


/****************** Post-dx: 6 months to 4 years ***********************/
%res(analysis=avestop.adj0, catlin=T,   
         data=pool_postcrc, event=dead, 
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         tbegtend=T,
         tbeg=time_begin,
         tend=time_end,
         strata= cohort intervalp, 
         labels=T, 
         exp1=&lcs_crcq_   ,    adjust1=   ,
         exp2=&alcs_crcq_  ,    adjust2=   ,
         exp3=&vlcs_crcq_  ,    adjust3=   ,
         exp4=&ulcs_crcq_  ,    adjust4=   ,
         exp5=&hlcs_crcq_  ,    adjust5=   , 
         exp6=&lfs_crcq_   ,    adjust6=   ,
         exp7=&ulfs_crcq_  ,    adjust7=   ,
         exp8=&hlfs_crcq_  ,    adjust8=   ,
         exp9=&alfs_crcq_  ,    adjust9=   ,
         exp10=&vlfs_crcq_  ,    adjust10=   );
         Run;

%res(analysis=avestop.adj1,catlin=T,   
         data=pool_postcrc, event=dead, 
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         tbegtend=T,
         tbeg=time_begin,
         tend=time_end,
         strata= cohort intervalp, 
         labels=T, 
         exp1=&lcs_crcq_   ,   adjust1= age_dx &calor_crcq_ &smkstaa_ &alco_crcq_ &postact_ &bmid9_   &pre_aheim_fixq_  white mv mvita aspa  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES ,
         exp2=&alcs_crcq_  ,   adjust2= age_dx &calor_crcq_ &smkstaa_ &alco_crcq_ &postact_ &bmid9_   &pre_aheim_fixq_  white mv mvita aspa  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES ,
         exp3=&vlcs_crcq_  ,   adjust3= age_dx &calor_crcq_ &smkstaa_ &alco_crcq_ &postact_ &bmid9_   &pre_aheim_fixq_  white mv mvita aspa  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES ,
         exp4=&ulcs_crcq_  ,   adjust4= age_dx &calor_crcq_ &smkstaa_ &alco_crcq_ &postact_ &bmid9_   &pre_aheim_fixq_  white mv mvita aspa  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES ,
         exp5=&hlcs_crcq_  ,   adjust5= age_dx &calor_crcq_ &smkstaa_ &alco_crcq_ &postact_ &bmid9_   &pre_aheim_fixq_  white mv mvita aspa  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES , 
         exp6=&lfs_crcq_   ,   adjust6= age_dx &calor_crcq_ &smkstaa_ &alco_crcq_ &postact_ &bmid9_   &pre_aheim_fixq_  white mv mvita aspa  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES ,
         exp7=&ulfs_crcq_  ,   adjust7= age_dx &calor_crcq_ &smkstaa_ &alco_crcq_ &postact_ &bmid9_   &pre_aheim_fixq_  white mv mvita aspa  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES ,
         exp8=&hlfs_crcq_  ,   adjust8= age_dx &calor_crcq_ &smkstaa_ &alco_crcq_ &postact_ &bmid9_   &pre_aheim_fixq_  white mv mvita aspa  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES ,
         exp9=&alfs_crcq_  ,   adjust9= age_dx &calor_crcq_ &smkstaa_ &alco_crcq_ &postact_ &bmid9_   &pre_aheim_fixq_  white mv mvita aspa  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES ,
         exp10=&vlfs_crcq_  ,  adjust10=age_dx &calor_crcq_ &smkstaa_ &alco_crcq_ &postact_ &bmid9_   &pre_aheim_fixq_  white mv mvita aspa  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES );
         Run;

%res(analysis=unit.adj0, catlin=T,   
         data=pool_postcrc, event=dead,  
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         tbegtend=T,
         tbeg=time_begin,
         tend=time_end,
         strata= cohort intervalp, 
         labels=T, 
         exp1=lcs_sd_crc   ,    adjust1=  ,
         exp2=alcs_sd_crc  ,    adjust2=  ,
         exp3=vlcs_sd_crc  ,    adjust3=  ,
         exp4=ulcs_sd_crc  ,    adjust4=  ,
         exp5=hlcs_sd_crc  ,    adjust5=  , 
         exp6=lfs_sd_crc   ,    adjust6=  ,
         exp7=ulfs_sd_crc  ,    adjust7=  ,
         exp8=hlfs_sd_crc  ,    adjust8=  ,
         exp9=alfs_sd_crc  ,    adjust9=  ,
         exp10=vlfs_sd_crc  ,   adjust10=  );      
         Run;

%res(analysis=unit.adj1,catlin=T,   
         data=pool_postcrc, event=dead,  
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         tbegtend=T,
         tbeg=time_begin,
         tend=time_end,
         strata= cohort intervalp, 
         labels=T, 
         exp1=lcs_sd_crc   ,   adjust1= age_dx &calor_crcq_ &smkstaa_ &alco_crcq_ &postact_ &bmid9_   &pre_aheim_fixq_  white mv mvita aspa  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES ,
         exp2=alcs_sd_crc  ,   adjust2= age_dx &calor_crcq_ &smkstaa_ &alco_crcq_ &postact_ &bmid9_   &pre_aheim_fixq_  white mv mvita aspa  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES ,
         exp3=vlcs_sd_crc  ,   adjust3= age_dx &calor_crcq_ &smkstaa_ &alco_crcq_ &postact_ &bmid9_   &pre_aheim_fixq_  white mv mvita aspa  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES ,
         exp4=ulcs_sd_crc  ,   adjust4= age_dx &calor_crcq_ &smkstaa_ &alco_crcq_ &postact_ &bmid9_   &pre_aheim_fixq_  white mv mvita aspa  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES ,
         exp5=hlcs_sd_crc  ,   adjust5= age_dx &calor_crcq_ &smkstaa_ &alco_crcq_ &postact_ &bmid9_   &pre_aheim_fixq_  white mv mvita aspa  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES , 
         exp6=lfs_sd_crc   ,   adjust6= age_dx &calor_crcq_ &smkstaa_ &alco_crcq_ &postact_ &bmid9_   &pre_aheim_fixq_  white mv mvita aspa  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES ,
         exp7=ulfs_sd_crc  ,   adjust7= age_dx &calor_crcq_ &smkstaa_ &alco_crcq_ &postact_ &bmid9_   &pre_aheim_fixq_  white mv mvita aspa  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES ,
         exp8=hlfs_sd_crc  ,   adjust8= age_dx &calor_crcq_ &smkstaa_ &alco_crcq_ &postact_ &bmid9_   &pre_aheim_fixq_  white mv mvita aspa  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES ,
         exp9=alfs_sd_crc  ,   adjust9= age_dx &calor_crcq_ &smkstaa_ &alco_crcq_ &postact_ &bmid9_   &pre_aheim_fixq_  white mv mvita aspa  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES ,
         exp10=vlfs_sd_crc  ,  adjust10=age_dx &calor_crcq_ &smkstaa_ &alco_crcq_ &postact_ &bmid9_   &pre_aheim_fixq_  white mv mvita aspa  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES ); 
         Run;

%res(analysis=trend.adj0, catlin=T,   
         data=pool_postcrc, event=dead, 
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         tbegtend=T,
         tbeg=time_begin,
         tend=time_end,
         strata= cohort intervalp, 
         labels=T, 
         exp1=lcscrctrend  ,     adjust1=  ,
         exp2=alcscrctrend ,     adjust2=  ,
         exp3=vlcscrctrend ,     adjust3=  ,
         exp4=ulcscrctrend ,     adjust4=  ,
         exp5=hlcscrctrend ,     adjust5=  , 
         exp6=lfscrctrend  ,     adjust6=  ,
         exp7=ulfscrctrend ,     adjust7=  ,
         exp8=hlfscrctrend ,     adjust8=  ,
         exp9=alfscrctrend ,     adjust9=  ,
         exp10=vlfscrctrend ,    adjust10=  );
         Run;

%res(analysis=trend.adj1,catlin=T,   
         data=pool_postcrc, event=dead,  
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         tbegtend=T,
         tbeg=time_begin,
         tend=time_end,
         strata= cohort intervalp, 
         labels=T, 
         exp1=lcscrctrend  ,    adjust1= age_dx &calor_crcq_ &smkstaa_ &alco_crcq_ &postact_ &bmid9_   &pre_aheim_fixq_  white mv mvita aspa  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp2=alcscrctrend ,    adjust2= age_dx &calor_crcq_ &smkstaa_ &alco_crcq_ &postact_ &bmid9_   &pre_aheim_fixq_  white mv mvita aspa  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp3=vlcscrctrend ,    adjust3= age_dx &calor_crcq_ &smkstaa_ &alco_crcq_ &postact_ &bmid9_   &pre_aheim_fixq_  white mv mvita aspa  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES  ,
         exp4=ulcscrctrend ,    adjust4= age_dx &calor_crcq_ &smkstaa_ &alco_crcq_ &postact_ &bmid9_   &pre_aheim_fixq_  white mv mvita aspa  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES ,
         exp5=hlcscrctrend ,    adjust5= age_dx &calor_crcq_ &smkstaa_ &alco_crcq_ &postact_ &bmid9_   &pre_aheim_fixq_  white mv mvita aspa  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES , 
         exp6=lfscrctrend  ,    adjust6= age_dx &calor_crcq_ &smkstaa_ &alco_crcq_ &postact_ &bmid9_   &pre_aheim_fixq_  white mv mvita aspa  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES ,
         exp7=ulfscrctrend ,    adjust7= age_dx &calor_crcq_ &smkstaa_ &alco_crcq_ &postact_ &bmid9_   &pre_aheim_fixq_  white mv mvita aspa  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES ,
         exp8=hlfscrctrend ,    adjust8= age_dx &calor_crcq_ &smkstaa_ &alco_crcq_ &postact_ &bmid9_   &pre_aheim_fixq_  white mv mvita aspa  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES ,
         exp9=alfscrctrend ,    adjust9= age_dx &calor_crcq_ &smkstaa_ &alco_crcq_ &postact_ &bmid9_   &pre_aheim_fixq_  white mv mvita aspa  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES ,
         exp10=vlfscrctrend ,   adjust10=age_dx &calor_crcq_ &smkstaa_ &alco_crcq_ &postact_ &bmid9_   &pre_aheim_fixq_  white mv mvita aspa  pmhever ocever &est_status2_ &stage_  her2_status &chemotx_ &hormtx_ &radiotx_ &surg_ nSES );
         Run;



%include 'macro/resout.sas';
/************************ HR *****************************************/

%resout(outfile=tableS21_brc_LCD_dead, dbms=tab, opt=replace);                                                                                                                             
endsas;



