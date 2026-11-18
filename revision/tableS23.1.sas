/*********************************************************************
 Purpose : Reviewer comment #9 (reverse causation) - survival analyses
 Contents: Table S23
 unlagged exposure; follow-up starts 2 years after
                    the first post-dx dietary assessment (deaths and
                    person-time in that window excluded)

 Run once with site=brc and once with site=crc.
*********************************************************************/

%let site = brc;                       /* brc or crc */

/*------------------------------------------------------------------
  Site-specific paths
------------------------------------------------------------------*/
%macro setpath;
    %global datapath outdir;
    %let datapath = lowcarb_cancer/revision/data_postdx&site..sas;
    %let outdir   = lowcarb_cancer/revision/table;
%mend setpath;
%setpath;

%include "&datapath";
%include '/macro/res_post.sas';

/*------------------------------------------------------------------
  Site-specific adjustment set and time-varying covariates to lag
  (fixed covariates such as age_dx, stage, treatment, pre-dx AHEI
   are not lagged)
------------------------------------------------------------------*/
%macro setadj;
    %global adj1 tvcov;
    %if &site = brc %then %do;
        %let adj1  = age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_
                     &pre_aheim_fixq_ white mvita aspa cafh pmhever ocever
                     &est_status2_ &stage_ her2_status &chemotx_ &hormtx_
                     &radiotx_ &surg_ nSES;
        %let tvcov = calorsvq &calorsvq_    
                     bmid9 &bmid9_   nSES;
    %end;
    %else %if &site = crc %then %do;
        %let adj1  = age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid3cat_
                     &pre_aheim_fixq_ white mvita aspa cafh  
                     &stage_cr_ &grade_cr_ &subsite_cr_ nSES;
        %let tvcov = calorsvq &calorsvq_  
                     bmid3cat &bmid3cat_  nSES;
    %end;
%mend setadj;
%setadj;

/* Exposures to lag: base quartile, quartile indicators, per SD, trend */
%let explag =
    lcssvq  &lcssvq_  lcs_sd  lcssvtrend
    alcssvq &alcssvq_ alcs_sd alcssvtrend
    vlcssvq &vlcssvq_ vlcs_sd vlcssvtrend
    ulcssvq &ulcssvq_ ulcs_sd ulcssvtrend
    hlcssvq &hlcssvq_ hlcs_sd hlcssvtrend
    lfssvq  &lfssvq_  lfs_sd  lfssvtrend
    alfssvq &alfssvq_ alfs_sd alfssvtrend
    vlfssvq &vlfssvq_ vlfs_sd vlfssvtrend
    ulfssvq &ulfssvq_ ulfs_sd ulfssvtrend
    hlfssvq &hlfssvq_ hlfs_sd hlfssvtrend;

%let lagvars = &explag &tvcov;



/*===================================================================
  start follow-up 2 years after first post-dx FFQ
===================================================================*/
proc sql;
    create table pool_lm as
    select a.*, b.first_tb, b.first_tb + 24 as landmark
    from pool_timevarying as a
         inner join
         (select newid, min(time_begin) as first_tb
          from pool_timevarying
          where person_time > 0 and time_begin ne .
          group by newid) as b
         on a.newid = b.newid;
quit;

data pool_lm;
    set pool_lm;
    if time_end <= landmark then delete;          /* person-time/deaths within 2 y */
    time_begin  = max(time_begin, landmark);      /* delayed entry at landmark    */
    person_time = time_end - time_begin;
run;


/*===================================================================
  Sample sizes for each analysis (for table footnotes)
===================================================================*/
%macro nsum(data=, label=);
    proc sql;
        create table _n as
        select "&label" as analysis length=20,
               count(distinct newid) as n_people,
               sum(dead)             as deaths,
               sum(person_time)/12   as py format=comma12.
        from &data
        where person_time > 0;
    quit;
    proc append base=nsum_all data=_n force; run;
%mend nsum;

proc datasets lib=work nolist; delete nsum_all; quit;
%nsum(data=pool_timevarying, label=Main);
%nsum(data=pool_lm,          label=Landmark 2y);

title "&site: people, deaths and person-years by analysis";
proc print data=nsum_all noobs; run;
title;


/*===================================================================
  Cox models: quartiles, per SD, trend
===================================================================*/
%macro runres(data=, tag=);

%res(analysis=&tag..q.adj1, catlin=T,
     data=&data, event=dead,
     timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
     id=newid, tvar=intervalp,
     tbegtend=T, tbeg=time_begin, tend=time_end,
     strata=agecon cohort intervalp, labels=T,
     exp1=&lcssvq_,  adjust1=&adj1,   exp2=&alcssvq_, adjust2=&adj1,
     exp3=&vlcssvq_, adjust3=&adj1,   exp4=&ulcssvq_, adjust4=&adj1,
     exp5=&hlcssvq_, adjust5=&adj1,   exp6=&lfssvq_,  adjust6=&adj1,
     exp7=&ulfssvq_, adjust7=&adj1,   exp8=&hlfssvq_, adjust8=&adj1,
     exp9=&alfssvq_, adjust9=&adj1,   exp10=&vlfssvq_, adjust10=&adj1);
run;

%res(analysis=&tag..sd.adj1, catlin=T,
     data=&data, event=dead,
     timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
     id=newid, tvar=intervalp,
     tbegtend=T, tbeg=time_begin, tend=time_end,
     strata=agecon cohort intervalp, labels=T,
     exp1=lcs_sd,  adjust1=&adj1,   exp2=alcs_sd, adjust2=&adj1,
     exp3=vlcs_sd, adjust3=&adj1,   exp4=ulcs_sd, adjust4=&adj1,
     exp5=hlcs_sd, adjust5=&adj1,   exp6=lfs_sd,  adjust6=&adj1,
     exp7=ulfs_sd, adjust7=&adj1,   exp8=hlfs_sd, adjust8=&adj1,
     exp9=alfs_sd, adjust9=&adj1,   exp10=vlfs_sd, adjust10=&adj1);
run;

%res(analysis=&tag..trend.adj1, catlin=T,
     data=&data, event=dead,
     timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
     id=newid, tvar=intervalp,
     tbegtend=T, tbeg=time_begin, tend=time_end,
     strata=agecon cohort intervalp, labels=T,
     exp1=lcssvtrend,  adjust1=&adj1,   exp2=alcssvtrend, adjust2=&adj1,
     exp3=vlcssvtrend, adjust3=&adj1,   exp4=ulcssvtrend, adjust4=&adj1,
     exp5=hlcssvtrend, adjust5=&adj1,   exp6=lfssvtrend,  adjust6=&adj1,
     exp7=ulfssvtrend, adjust7=&adj1,   exp8=hlfssvtrend, adjust8=&adj1,
     exp9=alfssvtrend, adjust9=&adj1,   exp10=vlfssvtrend, adjust10=&adj1);
run;

%mend runres;

%runres(data=pool_lm,   tag=landmark2);

%include '/macro/resout.sas';
%resout(outfile=table/lag_landmark_dead_&site, dbms=tab, opt=replace);

proc export data=nsum_all dbms=csv replace
    outfile="&outdir./n_lag_landmark_dead&site..csv";
run;

endsas;
