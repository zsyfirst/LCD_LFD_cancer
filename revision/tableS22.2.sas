/*********************************************************************
 Purpose : Reviewer comment #9 (reverse causation) 
 Contents: Table S22.2
   (A) HRs stratified by time from diagnosis to the first post-dx
       dietary assessment (<=2 vs >2 years)          -- person-level
   (B) HRs stratified by follow-up time since diagnosis
       (<5, 5 to <10, >=10 years)                    -- time-varying
   (C) Tests for heterogeneity (Wald test, exposure x stratum)
   (D) Cases / person-years by stratum

 Time scale: months since diagnosis (time_begin, time_end).
 Quartiles and SDs are kept from the main analysis (not re-ranked
 within strata), so HRs are comparable across strata.
*********************************************************************/

%include 'lowcarb_cancer/revision/data_postdx_crc.sas';
%include 'macro/res_post.sas';

%let outdir = lowcarb_cancer/revision/table;
%let events = dead dead_crc dead_can dead_nocan dead_cvd dead_resp dead_oth;

%let adj1 = age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_  &bmid3cat_
            &pre_aheim_fixq_ white mvita aspa cafh
            &stage_cr_ &grade_cr_ /*crcduration*/ &subsite_cr_  nSES;

%let qexps  = lcssvq alcssvq vlcssvq ulcssvq hlcssvq
              lfssvq alfssvq vlfssvq ulfssvq hlfssvq;
%let sdexps = lcs_sd alcs_sd vlcs_sd ulcs_sd hlcs_sd
              lfs_sd alfs_sd vlfs_sd ulfs_sd hlfs_sd;


/*===================================================================
  (A) Time from diagnosis to first post-dx dietary assessment
      = start (months since dx) of the first interval with post-dx diet
===================================================================*/
proc sql;
    create table _firstdiet as
    select newid,
           min(time_begin)         as dx_to_diet,       /* months */
           calculated dx_to_diet/12 as dx_to_diet_yr
    from pool_timevarying
    where person_time > 0 and time_begin ne .
    group by newid;

    create table pool_a as
    select a.*, b.dx_to_diet,
           case when b.dx_to_diet <= 24 then 1 else 2 end as dxdiet_c
    from pool_timevarying as a
         inner join _firstdiet as b
         on a.newid = b.newid
    where a.person_time > 0;
quit;

/* Check the split before interpreting: people and deaths per group */
title "Years from diagnosis to first post-dx dietary assessment";
proc means data=_firstdiet n mean median p25 p75 min max maxdec=2;
    var dx_to_diet_yr;
run;
proc sql;
    select dxdiet_c,
           count(distinct newid)   as n_people,
           sum(dead)               as deaths,
           sum(person_time)/12     as py format=comma12.
    from pool_a
    group by dxdiet_c;
quit;
title;


/*===================================================================
  (B) Split each record at 60 and 120 months since diagnosis so the
      stratum matches follow-up time exactly (events stay in the
      final piece of each original record)
===================================================================*/
%let nev = %sysfunc(countw(&events));

data pool_b;
    set pool_a;
    array ev{&nev} &events;
    array evsave{&nev} _temporary_;
    array bound{4} _temporary_ (0 60 120 1e9);

    do k = 1 to &nev; evsave{k} = ev{k}; end;
    _tb = time_begin;
    _te = time_end;

    do c = 1 to 3;
        seg_b = max(_tb, bound{c});
        seg_e = min(_te, bound{c+1});
        if seg_e > seg_b then do;
            time_begin  = seg_b;
            time_end    = seg_e;
            person_time = seg_e - seg_b;
            tsd_c       = c;          /* 1: <5 y, 2: 5-<10 y, 3: >=10 y */
            do k = 1 to &nev;
                ev{k} = ifn(seg_e = _te, evsave{k}, 0);
            end;
            output;
        end;
    end;
    drop k c _tb _te seg_b seg_e;
run;

title "Deaths and person-years by follow-up time since diagnosis";
proc sql;
    select tsd_c,
           count(distinct newid) as n_people,
           sum(dead)             as deaths,
           sum(person_time)/12   as py format=comma12.
    from pool_b
    group by tsd_c;
quit;
title;


/*===================================================================
  Stratified HRs (quartiles, per SD, trend) via %res with byvar=
===================================================================*/
%macro strat_res(data=, byvar=, tag=);

%res(analysis=&tag..q.adj1, catlin=T,
     data=&data, event=dead,
     timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
     id=newid, tvar=intervalp,
     tbegtend=T, tbeg=time_begin, tend=time_end,
     strata=agecon cohort intervalp, labels=T, byvar=&byvar,
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
     strata=agecon cohort intervalp, labels=T, byvar=&byvar,
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
     strata=agecon cohort intervalp, labels=T, byvar=&byvar,
     exp1=lcssvtrend,  adjust1=&adj1,   exp2=alcssvtrend, adjust2=&adj1,
     exp3=vlcssvtrend, adjust3=&adj1,   exp4=ulcssvtrend, adjust4=&adj1,
     exp5=hlcssvtrend, adjust5=&adj1,   exp6=lfssvtrend,  adjust6=&adj1,
     exp7=ulfssvtrend, adjust7=&adj1,   exp8=hlfssvtrend, adjust8=&adj1,
     exp9=alfssvtrend, adjust9=&adj1,   exp10=vlfssvtrend, adjust10=&adj1);
run;

%mend strat_res;

%strat_res(data=pool_a, byvar=dxdiet_c, tag=dxdiet);   /* (A) */
%strat_res(data=pool_b, byvar=tsd_c,    tag=tsd);      /* (B) */

%include 'macro/resout.sas';
%resout(outfile=table/strat_timesincedx_dead, dbms=tab, opt=replace);


/*===================================================================
  (C) Heterogeneity: likelihood ratio test of exposure x stratum terms
      main=1 : stratum main effect included (person-level stratifier)
      main=0 : omitted, because follow-up-time strata are a function
               of the time scale and are absorbed by the baseline hazard
===================================================================*/
proc datasets lib=work nolist; delete het_all; quit;

%macro het(data=, byvar=, nlev=, main=, exps=, tag=);

    %let nx = %sysfunc(countw(&exps, %str( )));

    /* stratum indicators and all interaction terms in one pass */
    data _h;
        set &data;
        %do l = 2 %to &nlev;
            _grp&l = (&byvar = &l);
            %do e = 1 %to &nx;
                %let x = %scan(&exps, &e, %str( ));
                _i_&x._&l = &x * _grp&l;
            %end;
        %end;
    run;

    %do e = 1 %to &nx;
        %let x = %scan(&exps, &e, %str( ));
        %let grps = ;
        %let ints = ;
        %do l = 2 %to &nlev;
            %let grps = &grps _grp&l;
            %let ints = &ints _i_&x._&l;
        %end;

    /* Model without interaction */
        ods select FitStatistics;
        ods output FitStatistics = _f0;
        proc phreg data=_h;
            model (time_begin, time_end)*dead(0) =
                  &x %if &main = 1 %then &grps; &adj1 / ties=efron;
            strata agecon cohort intervalp;
        run;
    /* Model with interaction */
        ods select FitStatistics;
        ods output FitStatistics = _f1;
        proc phreg data=_h;
            model (time_begin, time_end)*dead(0) =
                  &x %if &main = 1 %then &grps; &ints &adj1 / ties=efron;
            strata agecon cohort intervalp;
        run;
        ods select all;

    /* LRT = difference in -2 log L */
        data _t;
            length stratifier exposure $20;
            merge _f0(where=(criterion = '-2 LOG L') keep=criterion withcovariates
                      rename=(withcovariates = m2ll0))
                  _f1(where=(criterion = '-2 LOG L') keep=criterion withcovariates
                      rename=(withcovariates = m2ll1));
            stratifier = "&tag";
            exposure   = "&x";
            lrt_df     = &nlev - 1;
            lrt_chisq  = m2ll0 - m2ll1;
            lrt_p      = 1 - probchi(lrt_chisq, lrt_df);
            keep stratifier exposure lrt_chisq lrt_df lrt_p;
        run;
        proc append base=het_all data=_t force; run;
   %end;   

   proc datasets lib=work nolist; delete _h _t _f0 _f1; quit;

%mend het;

%let trexps = lcssvtrend alcssvtrend vlcssvtrend ulcssvtrend hlcssvtrend
              lfssvtrend alfssvtrend vlfssvtrend ulfssvtrend hlfssvtrend;

/* Per-SD */
%het(data=pool_a, byvar=dxdiet_c, nlev=2, main=1, exps=&sdexps, tag=dxdiet);
%het(data=pool_b, byvar=tsd_c,    nlev=3, main=0, exps=&sdexps, tag=tsd);

/* Trend (median of quartile), consistent with incidence Table S13 */
%het(data=pool_a, byvar=dxdiet_c, nlev=2, main=1, exps=&trexps, tag=dxdiet);
%het(data=pool_b, byvar=tsd_c,    nlev=3, main=0, exps=&trexps, tag=tsd);

title "P for heterogeneity (exposure x stratum, likelihood ratio test)";
proc print data=het_all noobs;
    format lrt_chisq 8.2 lrt_p pvalue6.4;
run;
title;

proc export data=het_all dbms=csv replace
    outfile="&outdir./phet_timesincedx_dead.csv";
run;


/*===================================================================
  (D) Cases / person-years by stratum and quartile
===================================================================*/
%macro cp_strat(data=, byvar=, tag=, exps=);
    proc datasets lib=work nolist; delete cp_&tag; quit;
    %do e = 1 %to %sysfunc(countw(&exps, %str( )));
        %let x = %scan(&exps, &e, %str( ));
        proc sql;
            create table _cp as
            select "&x" as exposure length=20,
                   &byvar as stratum,
                   &x     as quartile,
                   sum(dead)           as cases,
                   sum(person_time)/12 as py,
                   cats(put(calculated cases, comma12.), ' / ',
                        put(calculated py, comma12.)) as casepy length=40
            from &data
            where &x ne .
            group by &byvar, &x;
        quit;
        proc append base=cp_&tag data=_cp force; run;
    %end;
    proc export data=cp_&tag dbms=csv replace
        outfile="&outdir./casePY_&tag._dead.csv";
    run;
%mend cp_strat;

%cp_strat(data=pool_a, byvar=dxdiet_c, tag=dxdiet, exps=&qexps);
%cp_strat(data=pool_b, byvar=tsd_c,    tag=tsd,    exps=&qexps);

endsas;
