/*****************************************************************************************
  Purpose: Table S23
  IPW Cox models: post-diagnosis CUMULATIVE AVERAGE LCD/LFD scores and mortality after CRC
  Pooled NHS + NHSII + HPFS  (main model = your %res adj2)

  Selection being corrected:
    A CRC case enters the cumulative-average analysis only after returning >=1 FFQ
    after diagnosis. Weights = inverse probability of returning any post-dx FFQ,
    among cases who had the chance to (alive when an FFQ was mailed after diagnosis).

    Eligible  : elig_any = 1  (alive and under follow-up when >=1 FFQ cycle was mailed
                after dx)  -- created in each cohort file, STEP 0
    Responder : person appears in pool_timevarying (has a non-missing lcssvcon)
    Model     : Firth logistic, pre-dx / at-dx covariates only
    Weight    : stabilized, truncated at 1st/99th pct, fixed per person
    Cox       : adj2 model, WEIGHT sw, robust sandwich SE (covs(aggregate), id newid)

  ----------------------------------------------------------------------------------------
  STEP 0a  In EACH cohort file in /fix_postdx/ (nhs, nhs2, hpfs .data.crc_postdx.sas):

    (1) Comment out the period-level diet exclusion, in BOTH places
        (baseline block i=1 and follow-up block i>1):
            *%exclude(lcssvcon eq .,  skip=t);
        -> non-responders must reach pool so the response model can see them.
           The same restriction is re-applied in data_postdx.sas (STEP 0b), so the
           analytic data and the unweighted %res results do not change.

    (2) Add the eligibility flag (person level, BEFORE %beginex), keep elig_any:

        elig_any = 0;
        if dtdxca > 0 then do f = 1 to dim(ffqk);  k = ffqk{f};   -- FFQ positions
           nomd = 1014 + 24*k;        -- NHS ; NHSII 1074+24*k ; HPFS 1009+24*k
           if nomd > dtdxca and nomd < cutoff and (dtdth in (.,0) or dtdth > nomd)
              then elig_any = 1;      -- alive when an FFQ was mailed after dx
        end;

        (same loop as elig_6m4y, only the window changes: any FFQ after dx)

  STEP 0b  In data_postdx.sas, pool_timevarying step: UN-comment the restriction and
           drop zero-length rows -- both must come BEFORE proc rank:

        data pool_timevarying;
            set pool;
            if lcssvcon = . then delete;          -- re-applies the cohort-file exclusion
            if time_end <= time_begin then delete;-- zero-length final intervals
            ...
*****************************************************************************************/

%include 'lowcarb_cancer/revision/data_postdx_crc.sas';

/*========================================================================================
  1. Eligible cases (one row each, diagnosis-period row = pre-dx covariates) and response
========================================================================================*/
data _chk;
   set pool;  by newid notsorted;
   if first.newid;
   has_ffq = 0;
run;
proc sql;
   create table _resp as select distinct newid, 1 as resp from pool_timevarying;
   update _chk set has_ffq = 1 where newid in (select newid from _resp);
quit;

title 'CHECK 1: CRC cases in pool by eligibility and post-dx FFQ (any)';
proc freq data=_chk;  tables cohort*elig_any*has_ffq / list missing;  run;
/* expected: elig_any=1/has_ffq=0 = non-responders (several hundred);
             elig_any=0/has_ffq=1 should be 0 (if not, check the nominal dates)     */

proc sort data=pool out=_dxrow;  by newid descending cancase;  run;

data ipw_base;
   merge _dxrow(in=a) _resp;
   by newid;
   if a and first.newid;               /* cancase=1 row                              */
   if resp = . then resp = 0;
   if elig_any = 1 or resp = 1;        /* had the chance to return a post-dx FFQ      */
   age_dx2 = age_dx**2;
   dxyr    = dtdx_cr - 2000;
   dxyr2   = dxyr**2;

   if bmiconp = .        then bmip_c = .;
   else if bmiconp < 25  then bmip_c = 1;
   else if bmiconp < 30  then bmip_c = 2;
   else                       bmip_c = 3;

if 0 <= actconp and actconp < 7  then preact=1;
      else if actconp>=7 and actconp<15 then preact=2;
      else if actconp>=15 and actconp<25 then preact=3;
      else if actconp>=25  then preact=4;
      else  preact=2;

   if cohort = 3 then do; phmsp = 0; ocever = 0; end;   /* men */
run;

proc rank data=ipw_base out=ipw_base groups=4;
   var     nSES   calorp   pre_aheim_fix;
   ranks   nses_q calorp_q preahei_q;
run;

title 'IPW: eligible cases and responders by cohort';
proc freq data=ipw_base;  tables cohort*resp / norow nopercent;  run;

title 'IPW: missingness of response-model covariates';
proc means data=ipw_base n nmiss;
   class cohort;
   var age_dx dxyr stage_cr grade_cr subsite_cr smkstp bmiconp actconp alcp calorp
       mvitp aspp phmsp ocever white cafh nSES pre_aheim_fix;
run;

/*========================================================================================
  2. Response model and stabilized, truncated weights
========================================================================================*/
%let ipwclass = cohort stage_cr grade_cr subsite_cr preahei_q
                smkstp bmip_c preact calorp_q alcp mvitp aspp  
                white cafh nses_q;
%let ipwcov   = age_dx age_dx2 dxyr dxyr2 &ipwclass;

title 'IPW: response model (returned any post-dx FFQ)';
proc logistic data=ipw_base;
   class &ipwclass / param=ref missing;
   model resp(event='1') = &ipwcov / firth;
   output out=_ps p=p_resp;
run;

proc sql noprint;  select mean(resp) into :pbar trimmed from ipw_base;  quit;

data _ps;  set _ps;  if resp = 1;  sw = &pbar / p_resp;  run;

proc univariate data=_ps noprint;
   var sw;  output out=_pct pctlpts=1 99 pctlpre=p;
run;

data ipw_w;
   if _n_ = 1 then set _pct;
   set _ps;
   sw_raw = sw;
   sw = min(max(sw, p1), p99);
   keep newid sw sw_raw p_resp;
run;

title 'IPW: weight distribution (mean ~1)';
proc means data=ipw_w n mean std min p1 p25 median p75 p99 max;  var sw_raw sw;  run;

/* balance: eligible vs responders unweighted / weighted */
proc sort data=ipw_w;  by newid;  run;
data _bal;
   merge ipw_base(in=a) ipw_w;  by newid;  if a;
   if resp = 1 then do; grp='2 Responders weighted  '; wt=sw; output; end;
   if resp = 1 then do; grp='1 Responders unweighted'; wt=1;  output; end;
   grp='0 All eligible cases   '; wt=1; output;
run;
title 'IPW: covariate balance';
proc means data=_bal mean maxdec=3 nway;
   class grp;  weight wt;
   var age_dx dxyr stage_cr bmiconp actconp calorp aspp white nSES pre_aheim_fix;
run;

/*========================================================================================
  3. Attach weights to the time-varying analytic data
========================================================================================*/
proc sort data=pool_timevarying;  by newid;  run;
data pool_tv_ipw;
   merge pool_timevarying(in=a) ipw_w(keep=newid sw);
   by newid;
   if a;
run;

title 'IPW: analytic rows without a weight (should be 0)';
proc means data=pool_tv_ipw n nmiss;  var sw;  run;

/*========================================================================================
  4. Weighted Cox models -- same specification as %res adj2
========================================================================================*/
%let adjcov = age_dx &calorsvq_ &smkstaa_ &alcovc_post_ &postact_ &bmid9_ &pre_aheim_fixq_
              white mvita aspa cafh
              &stage_cr_ &grade_cr_ crcduration &subsite_cr_ nSES;

%macro ipwcox(exp=, keep=, label=, event=dead, wt=sw);
   proc phreg data=pool_tv_ipw covs(aggregate) nosummary;
      model (time_begin, time_end)*&event(0) = &exp &adjcov / ties=efron rl;
      strata agecon cohort intervalp;
      %if &wt ne none %then %do; weight &wt; %end;
      id newid;
      ods output ParameterEstimates=_pe;
   run;
   data _pe;
      length outcome $10 model $20 weighting $10 parameter $32;
      set _pe;
      outcome = "&event";  model = "&label";  weighting = "&wt";
      if index(upcase(parameter), upcase("&keep")) = 1;
      keep outcome model weighting parameter HazardRatio HRLowerCL HRUpperCL ProbChiSq;
   run;
   proc append base=cum_ipw_results data=_pe force;  run;
%mend;

%let scores = lcs alcs vlcs ulcs hlcs lfs ulfs hlfs alfs vlfs;

%macro runall(event=, wt=sw);
   %do i = 1 %to 10;
      %let s = %scan(&scores, &i);
      %ipwcox(exp=&&&s.svq_,     keep=&s.svq,     label=&s._Q,     event=&event, wt=&wt);
      %ipwcox(exp=&s._sd,        keep=&s._sd,     label=&s._perSD, event=&event, wt=&wt);
      %ipwcox(exp=&s.svtrend,    keep=&s.svtrend, label=&s._trend, event=&event, wt=&wt);
   %end;
%mend;

proc datasets lib=work nolist;  delete cum_ipw_results;  quit;

ods select none;
%runall(event=dead,     wt=sw);
%runall(event=dead_crc, wt=sw);
/* optional, same PHREG without weights (should match your %res adj2 closely): */
/* %runall(event=dead,     wt=none); */
/* %runall(event=dead_crc, wt=none); */
ods select all;

title 'IPW-weighted HRs: post-dx cumulative average';
proc print data=cum_ipw_results noobs;
   format HazardRatio HRLowerCL HRUpperCL 5.2 ProbChiSq pvalue6.3;
run;

proc export data=cum_ipw_results outfile='table/tableS_crc_cumLCD_ipw.csv' dbms=csv replace;
run;

endsas;
