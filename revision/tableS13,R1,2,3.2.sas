/*****************************************************************************************
Purpose:   Healthy-volunteer / selection-bias checks based on Reviewer 1 (comments 1 and 7)
             Part 1. Person-level dataset of ALL participants (cases and non-cases), built from BASE
             Part 2. Survival from diagnosis among ALL incident stage I-III cases  -> compare with SEER
             Part 3. Pre-diagnosis diet/lifestyle: incident cases vs cancer-free participants
             Part 4. Pre-diagnosis diet/lifestyle: post-dx FFQ responders vs non-responders
Cohort:    NHS, breast cancer (template).  See "ADAPT" notes for NHSII / HPFS / CRC.
Inputs:    nhs.data.brc_hv.sas  -> hv_nhs1, hv_long_nhs1   (needs BASE for NHS)
           nhs2.data.brc_hv.sas -> hv_nhs2, hv_long_nhs2   (needs BASE for NHSII)
Content:   Table S13, Table R1, R2, R3 in rebuttle letter
Date: 10/2026
*****************************************************************************************/

%let outdir = lowcarb_cancer/revision/table;
%let tag    = brc_pool;          /* suffix for output files */


/*========================================================================================
  PART 1. BUILD PERSON-LEVEL DATA FOR EACH COHORT, THEN STACK
========================================================================================*/

%include 'lowcarb_cancer/revision/nhs.data.brc_hv.sas';   /* hv_nhs1, hv_long_nhs1 */
%include 'lowcarb_cancer/revision/nhs2.data.brc_hv.sas';  /* hv_nhs2, hv_long_nhs2 */

/* Keep only cohort-independent variables before stacking (the raw arrays differ by cohort) */
%let hvkeep = newid id cohort incase noncase pidx pre_: next_ffq n_postffq ffq_resp in_analytic first_an responder elig_ffq alive_next resp_grp
              dtdxca dtdth surv_mo dead_any dead_site cod age_dx yeardx dxperiod agedxgp
              stage est_status her2_status chemotx radiotx hormtx surg sumstage ;


data hv_all;
    length resp_grp $30;
    set hv_nhs1(keep=&hvkeep) hv_nhs2(keep=&hvkeep);
run;

data hv_long_all;
    set hv_long_nhs1 hv_long_nhs2;
run;

/* ---- analytic person-period data (main analysis) and %indic3 macro variables ----
   Included here (not only in Part 5) so that analytic-cohort membership is read directly
   from the main analysis instead of being re-derived. */
%include 'lowcarb_cancer/revision/brc/fix_postdx/data_hv.sas';

proc sql;
    create table _an as 
    select distinct newid 
    from pool_timevarying 
    where not missing(lcssvcon)    /* has >=1 post-dx FFQ = enters main analysis */
    order by newid;
    title "Analytic cohort (main analysis) by cohort";
    select floor(newid/1000000) as cohort, count(*) as n_analytic from _an group by calculated cohort;
quit;
title;

proc sort data=hv_all; by newid; run;
data hv_all;
    merge hv_all(in=a) _an(in=b);
    by newid;
    if a;
    in_analytic_rule = in_analytic;                  /* rule re-derived in the cohort files (check only) */
    in_analytic = (b = 1 and incase = 1);            /* membership in the main analysis               */
    responder   = in_analytic;
    length resp_grp $30;
    if incase then do;
        if responder           then resp_grp = '1 In analytic cohort';
        else if not elig_ffq   then resp_grp = '0 No FFQ cycle after dx';
        else if not alive_next then resp_grp = '2 Died before next FFQ';
        else                        resp_grp = '3 Alive, no post-dx diet';
        dxera = 1 + (yeardx >= 2000);      /* 1 = before 2000, 2 = 2000 or later */
    end;
run;


/* CHECKS
   - in_analytic: totals must equal the analytic N printed above.
   - in_analytic*ffq_resp: off-diagonal cells should be ~0 after the *sv fix.
   - in_analytic*in_analytic_rule: should agree; if not, the cohort-file rule is out of date. */
proc freq data=hv_all;
    where incase = 1;
    tables cohort*resp_grp / nocol nopercent;
    tables cohort*in_analytic*ffq_resp / list missing;
    tables cohort*in_analytic*in_analytic_rule / list missing;
    title "Analytic-cohort membership among incident stage I-III breast cancer cases";
run;
title;

/*========================================================================================
  FORMATS (same as the characteristics-table program) AND THE SMD MACRO
========================================================================================*/
proc format;
    value smkstatusf  1='never smoker' 2='past smoker' 3='current smoker';
    value stagef      1='stage I' 2='stage II' 3='stage III';
    value est_statusf 1='ER-' 2='ER+' 3='Missing';
    value surgf       1='mastectomy' 2='lumpectomy' 3='unknown type' 4='no surgery' 5='Missing';
    value casef       0='Cancer-free' 1='Incident breast cancer';
    value respf 0='Non-responder' 1='Responder (in analytic cohort)';
run;

/* Standardized mean difference: (m1 - m0)/sqrt((s1^2 + s0^2)/2); unweighted (wt=one_) */
%macro smd(data=, grp=, vars=, wt=one_, out=, label=);
    %local k v n;
    %let n = %sysfunc(countw(&vars));
    %do k = 1 %to &n;
        %let v = %scan(&vars, &k);
        proc means data=&data noprint nway vardef=wgt;
            class &grp; var &v; weight &wt;
            output out=_m(keep=&grp _freq_ mean sd) mean=mean std=sd;
        run;
        data _s&k;
            length comparison $50 variable $32;
            merge _m(where=(&grp=1) rename=(mean=mean1 sd=sd1 _freq_=n1))
                  _m(where=(&grp=0) rename=(mean=mean0 sd=sd0 _freq_=n0));
            comparison = "&label";
            variable   = "&v";
            smd = (mean1 - mean0) / sqrt((sd1**2 + sd0**2)/2);
            abs_smd = abs(smd);
            keep comparison variable n1 mean1 sd1 n0 mean0 sd0 smd abs_smd;
        run;
    %end;
    data &out; set %do k=1 %to &n; _s&k %end; ; run;
    proc datasets nolist; delete _m %do k=1 %to &n; _s&k %end; ; quit;
%mend;

%let prevars = pre_ahei pre_lcs pre_alcs pre_vlcs pre_ulcs pre_hlcs pre_lfs pre_alfs pre_vlfs pre_ulfs pre_hlfs
               pre_calor pre_bmi pre_ob pre_act pre_alc pre_currsmk pre_neversmk pre_mv pre_asp pre_pmh pre_nses;



/*========================================================================================
  PART 2. SURVIVAL FROM DIAGNOSIS AMONG ALL INCIDENT STAGE I-III CASES -> compare with SEER
  SEER cause-specific survival censors deaths from other causes, so dead_site is analyzed
  with Kaplan-Meier. SEER*Explorer: breast, female; localized + regional (approx. AJCC I-III);
  matching age groups and years of diagnosis; 5- and 10-year cause-specific survival.
========================================================================================*/
%macro km_at(event=, strata=, out=);
    proc lifetest data=hv_all(where=(incase=1)) method=km timelist=60 120 reduceout outsurv=&out noprint;
        time surv_mo*&event(0);
        %if &strata ne %then strata &strata;;
    run;
    data &out; set &out; length outcome $12 strat_by $20;
        outcome = "&event"; strat_by = "&strata";
        surv_pct = 100*survival; lcl_pct = 100*sdf_lcl; ucl_pct = 100*sdf_ucl;
    run;
%mend;

proc format; value sumstagef 1='Localized' 2='Regional'; run;

%km_at(event=dead_site, strata=,         out=cs_all);
%km_at(event=dead_site, strata=stage,    out=cs_stage);
%km_at(event=dead_site, strata=dxperiod, out=cs_period);
%km_at(event=dead_site, strata=agedxgp,  out=cs_age);
%km_at(event=dead_site, strata=cohort,   out=cs_cohort);
%km_at(event=dead_site, strata=sumstage, out=cs_sumstage); /* add seer based stage  */

%km_at(event=dead_site, strata=dxera,          out=cs_era);
%km_at(event=dead_site, strata=dxera sumstage, out=cs_era_sumstage);
%km_at(event=dead_site, strata=dxera agedxgp,  out=cs_era_age);
%km_at(event=dead_site, strata=dxera stage, out=cs_era_stage);

%km_at(event=dead_any,  strata=,         out=os_all);
%km_at(event=dead_any,  strata=stage,    out=os_stage);
%km_at(event=dead_any,  strata=agedxgp,  out=os_age);
%km_at(event=dead_any,  strata=dxperiod,    out=os_period);
%km_at(event=dead_any,  strata=cohort,  out=os_cohort);
%km_at(event=dead_any, strata=sumstage, out=os_sumstage); /* add seer based stage  */


data surv_vs_seer;
    set cs_all cs_stage cs_period cs_age cs_cohort cs_sumstage
        cs_era cs_era_sumstage cs_era_age cs_era_stage
        os_all os_stage os_age os_period os_cohort os_sumstage;
    keep outcome strat_by stage sumstage dxperiod dxera agedxgp cohort timelist surv_pct lcl_pct ucl_pct;
run;
proc export data=surv_vs_seer outfile="&outdir/HV_survival_vs_SEER_&tag..csv" dbms=csv replace; run;

proc freq data=hv_all; where incase=1; tables dxera*sumstage dxera*agedxgp / norow nocol nopercent; run;


/*========================================================================================
  PART 3. PRE-DIAGNOSIS CHARACTERISTICS: INCIDENT CASES vs CANCER-FREE PARTICIPANTS
========================================================================================*/
/* case distribution of pre-diagnosis cycle, by cohort */
proc freq data=hv_long_all noprint;
    where case = 1;
    tables cohort*cycle / out=_cyc(keep=cohort cycle count) ;
run;
proc sql;
    create table _cycp as
    select cohort, cycle, count / sum(count) as p_cycle
    from _cyc group by cohort;
quit;

proc sql;
    create table _nc as
    select a.*, b.p_cycle
    from hv_long_all(where=(case=0)) as a left join _cycp as b
    on a.cohort = b.cohort and a.cycle = b.cycle;
quit;

/* weighted random draw of one cycle per person (exponential race: smallest key wins) */
data _nc;
    set _nc;
    if _n_ = 1 then call streaminit(20261014);
    if p_cycle > 0 then key = rand('exponential') / p_cycle;
    if key ne .;
run;
proc sort data=_nc; by newid key; run;
data _nc_one; set _nc; by newid; if first.newid; drop key p_cycle; run;

data hv_cmp;
    set hv_long_all(where=(case=1)) _nc_one;
    agegp = min(max(int((pre_age-30)/5), 1), 8);
    smkstatus = .;
    if pre_neversmk = 1 then smkstatus = 1;
    else if pre_currsmk = 1 then smkstatus = 3;
    else if pre_neversmk = 0 then smkstatus = 2;
    one_ = 1;
    format smkstatus smkstatusf. case casef.;
run;

/* check: cycle distribution should now be similar between groups within each cohort */
proc freq data=hv_cmp; tables cohort*cycle*case / nopercent norow nofreq; run;

%table1(data=hv_cmp,
        agegroup = agegp,
        exposure = case,
        varlist  = pre_age pre_ahei pre_lcs pre_alcs pre_vlcs pre_ulcs pre_hlcs pre_lfs pre_alfs pre_vlfs pre_ulfs pre_hlfs
                   pre_calor pre_bmi pre_ob pre_act pre_alc smkstatus pre_mv pre_asp pre_pmh pre_nses,
        noadj    = pre_age,
        cat      = pre_ob pre_mv pre_asp pre_pmh,
        poly     = smkstatus,
        rtftitle = Pre-diagnosis characteristics of participants who developed breast cancer and cancer-free participants,
        landscape= F,
        file     = %str(table/HV_cases_vs_noncases_&tag..rtf),
        dec      = 1, uselbl = T, nortf = F);

%smd(data=hv_cmp, grp=case, label=Incident breast cancer vs cancer-free,
     vars=pre_age &prevars, out=smd_case_noncase);
proc export data=smd_case_noncase outfile="&outdir/HV_SMD_cases_vs_noncases_&tag..csv" dbms=csv replace; run;


/*========================================================================================
  PART 4. RESPONDERS vs NON-RESPONDERS TO THE POST-DIAGNOSIS FFQ
  Responders     = incident cases included in the survival analysis (in_analytic = 1).
  Non-responders = incident cases NOT in the survival analysis who had >=1 FFQ cycle after
                   diagnosis and were alive at it (early death is not counted as non-response).
========================================================================================*/
proc freq data=hv_all; where incase = 1;
    tables cohort*resp_grp / out=resp_breakdown;
    title "Reasons for no post-diagnosis FFQ (report in Extended Data Fig. 1 or a Supplementary Table)";
run;
title;
proc export data=resp_breakdown outfile="&outdir/HV_response_breakdown_&tag..csv" dbms=csv replace; run;

data hv_resp;
    set hv_all;
    where incase = 1 and (responder = 1 or (elig_ffq = 1 and alive_next = 1));
    agegp  = min(max(int((age_dx-30)/5), 1), 8);
    agegp5 = min(max(int((age_dx-30)/5), 1), 10);
    stage2 = (stage = 2); stage3 = (stage = 3);
    er_pos   = (est_status = 1);  if est_status  = . then er_pos   = .;
    her2_pos = (her2_status = 1); if her2_status = . then her2_pos = .;
    chemo = (chemotx = 1); if chemotx = . then chemo = .;     /* missing stays missing */
    radio = (radiotx = 1); if radiotx = . then radio = .;
    horm  = (hormtx  = 1); if hormtx  = . then horm  = .;
    mastect = (surg = 1); if surg = . then mastect = .;
    est_status2 = est_status + 1; if est_status2 = . then est_status2 = 3;
    surg5 = surg; if surg5 = . then surg5 = 5;
    smkstatus = .;
    if pre_neversmk = 1 then smkstatus = 1;
    else if pre_currsmk = 1 then smkstatus = 3;
    else if pre_neversmk = 0 then smkstatus = 2;
    one_ = 1;
    format smkstatus smkstatusf. stage stagef. est_status2 est_statusf. surg5 surgf. responder respf.;
run;

%table1(data=hv_resp,
        agegroup = agegp,
        exposure = responder,
        varlist  = age_dx yeardx stage est_status2 her2_pos chemo radio horm surg5
                   &prevars smkstatus,
        noadj    = age_dx,
        cat      = her2_pos chemo radio horm pre_ob pre_currsmk pre_neversmk pre_mv pre_asp pre_pmh,
        poly     = stage est_status2 surg5 smkstatus,
        rtftitle = Characteristics of breast cancer cases who did and did not return a post-diagnosis FFQ,
        landscape= F,
        file     = %str(table/HV_responders_vs_nonresponders_&tag..rtf),
        dec      = 1, uselbl = T, nortf = F);

%smd(data=hv_resp, grp=responder, label=Responders vs non-responders,
     vars=age_dx yeardx stage2 stage3 er_pos her2_pos chemo radio horm mastect &prevars, out=smd_resp);
proc export data=smd_resp outfile="&outdir/HV_SMD_responders_vs_nonresponders_&tag..csv" dbms=csv replace; run;

/* Landmark analysis: mortality from the next FFQ cycle, responders vs alive non-responders */
data hv_landmark;
    set hv_resp;
    where elig_ffq = 1 and alive_next = 1;             /* common landmark = first FFQ cycle after dx */
    lm_time = surv_mo - (next_ffq - dtdxca);
    if lm_time > 0;
run;
proc phreg data=hv_landmark;
    class stage(ref=first) agedxgp dxperiod cohort / param=ref;   /* stage carries the stagef. format*/
    model lm_time*dead_any(0) = responder stage agedxgp dxperiod cohort / rl;
    title "Mortality after the next FFQ cycle: responders vs non-responders (HR for responder)";
run;
title;


