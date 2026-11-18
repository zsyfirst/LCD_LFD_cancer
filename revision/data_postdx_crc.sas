/*****************************************************
Purpose:       post-dx low carb, low fat, and mortality
Study design:  prospective cohort (NHS: 1986-2016.6, NHSII: 1991-2019.6, HPFS: 1986-2020.1)
Endpoints:     Mortality
Exposures:     post-diagnostic cumulative lcs, alcs, vlcs, ulcs, hlcs, lfs, ulfs, hlfs (main)
               post-diagnostic simple update diet score
Inclusion:     - breast cancer / colorectal cancer (analyzed seperately)
               - ICDA codes 153 (colon cancer) and 154 (rectal cancer)
               - confirmation codes 11-19 (including cases identified from death file);
               - diagnosed between 1984 and 2020
Exclusion:     baseline cancer or cvd or diabetes, missing information on age at enrollment, missing information on post-diagnosis diet
Covariates:    age, interval, cohort, white, energy quintile, smoking, alcohol, physical activity, cancer family history, 
               hormone use, oral contraceptive use, multivitamin use, aspirin use, 
               AHEI-2010, BMI prior to diagnosis, 
               stage, duration, tumor differentiation, location.
Date:  10/2026

******************************************************/

/************************ Call in Cleaned Datasets*******************************************************/
%include 'lowcarb_cancer/revision/hpfs.data.crc_postdx.sas';
%include 'lowcarb_cancer/revision/nhs.data.crc_postdx.sas';
%include 'lowcarb_cancer/revision/nhs2.data.crc_postdx.sas';


*===================================*
*==== Read in Datasets ====*
====================================*;

data nhs1data; 
    set nhs1;           
    *time1=irt84;
    time1=irt86;
    time2=irt88;
    time3=irt90;
    time4=irt92;
    time5=irt94;
    time6=irt96;
    time7=irt98;
    time8=irt00;
    time9=irt02;
    time10=irt04;
    time11=irt06;
    time12=irt08;
    time13=irt10;
    time14=irt10+24; /* 12 */
    time15=irt10+48; /* 14 */
    time16=irt10+72; /* 16 */
    time17=irt10+96; /* 18 */
    time18=irt10+120; /* 20 */

    tp01=t86;
    tp02=t88;
    tp03=t90;
    tp04=t92;
    tp05=t94;
    tp06=t96;
    tp07=t98;
    tp08=t00;
    tp09=t02;
    tp10=t04;
    tp11=t06;
    tp12=t08;
    tp13=t10;
    tp14=0;
    tp15=0;
    tp16=0;
    tp17=0;
    tp18=0;
    
    intervalp=.;
    if tp01 eq 1 then intervalp=1;
    if tp02 eq 1 then intervalp=2;
    if tp03 eq 1 then intervalp=3;
    if tp04 eq 1 then intervalp=4;
    if tp05 eq 1 then intervalp=5;
    if tp06 eq 1 then intervalp=6;
    if tp07 eq 1 then intervalp=7;
    if tp08 eq 1 then intervalp=8;
    if tp09 eq 1 then intervalp=9;
    if tp10 eq 1 then intervalp=10;
    if tp11 eq 1 then intervalp=11;
    if tp12 eq 1 then intervalp=12;
    if tp13 eq 1 then intervalp=13;
    if tp14 eq 1 then intervalp=14;
    if tp15 eq 1 then intervalp=15;
    if tp16 eq 1 then intervalp=16;
    if tp17 eq 1 then intervalp=17;
    if tp18 eq 1 then intervalp=18;

    sex=0;
    cohort=1;
    newid=(id+1000000);                        /*IDs for pooling- NHS begins with 1*/
run;


data nhs2data; 
    set nhs2; 
    *time1=retmo91-96+12; * 86 ;
    time1=retmo91-72+12;  * 86 ;
    time2=retmo91-48+12;  * 88 ;
    time3=retmo91-24+12;  * 90 ;
    time4=retmo91+12;     * 92 ;
    time5=retmo93+12;
    time6=retmo95+12;
    time7=retmo97+12;
    time8=retmo99+12;
    time9=retmo01+12;
    time10=retmo03+12;
    time11=retmo05+12;
    time12=retmo07+12;
    time13=retmo09+12;     * 10 ;
    time14=retmo11+12;     * 12 ;
    time15=retmo13+12;     * 14 ;
    time16=retmo15+12;     * 16 ;
    time17=retmo17+12;     * 18 ;
    time18=retmo19+12;     * 20 ;

    tp01=0;
    tp02=0;
    tp03=0;
    tp04=t91;
    tp05=t93;
    tp06=t95;
    tp07=t97;
    tp08=t99;
    tp09=t01;
    tp10=t03;
    tp11=t05;
    tp12=t07;
    tp13=t09;
    tp14=t11;
    tp15=t13;
    tp16=t15;
    tp17=t17;
    tp18=t19;


    intervalp=.;
    if tp01 eq 1 then intervalp=1;
    if tp02 eq 1 then intervalp=2;
    if tp03 eq 1 then intervalp=3;
    if tp04 eq 1 then intervalp=4;
    if tp05 eq 1 then intervalp=5;
    if tp06 eq 1 then intervalp=6;
    if tp07 eq 1 then intervalp=7;
    if tp08 eq 1 then intervalp=8;
    if tp09 eq 1 then intervalp=9;
    if tp10 eq 1 then intervalp=10;
    if tp11 eq 1 then intervalp=11;
    if tp12 eq 1 then intervalp=12;
    if tp13 eq 1 then intervalp=13;
    if tp14 eq 1 then intervalp=14;
    if tp15 eq 1 then intervalp=15;
    if tp16 eq 1 then intervalp=16;
    if tp17 eq 1 then intervalp=17;
    if tp18 eq 1 then intervalp=18;

    dtdth=dtdth+12;
    dtdxchd=dtdxchd+12;
    dtdxca=dtdxca+12;

    sex=0;
    cohort=2;
    newid=(id+2000000); 

run;


data hpfsdata; 
    set hpfs; 
    *time1=rtmnyr86-24;
    time1=rtmnyr86;
    time2=rtmnyr88;
    time3=rtmnyr90;
    time4=rtmnyr92;
    time5=rtmnyr94;
    time6=rtmnyr96;
    time7=rtmnyr98;
    time8=rtmnyr00;
    time9=rtmnyr02;
    time10=rtmnyr04;
    time11=rtmnyr06;
    time12=rtmnyr08;
    time13=rtmnyr10;
    time14=rtmnyr12;
    time15=rtmnyr14;
    time16=rtmnyr16;
    time17=rtmnyr18;
    time18=rtmnyr18+24; * 20 ;

    *tp01=0;
    tp01=t86;
    tp02=t88;
    tp03=t90;
    tp04=t92;
    tp05=t94;
    tp06=t96;
    tp07=t98;
    tp08=t00;
    tp09=t02;
    tp10=t04;
    tp11=t06;
    tp12=t08;
    tp13=t10;
    tp14=t12;
    tp15=t14;
    tp16=t16;
    tp17=t18;
    tp18=0;

    intervalp=.;
    if tp01 eq 1 then intervalp=1;
    if tp02 eq 1 then intervalp=2;
    if tp03 eq 1 then intervalp=3;
    if tp04 eq 1 then intervalp=4;
    if tp05 eq 1 then intervalp=5;
    if tp06 eq 1 then intervalp=6;
    if tp07 eq 1 then intervalp=7;
    if tp08 eq 1 then intervalp=8;
    if tp09 eq 1 then intervalp=9;
    if tp10 eq 1 then intervalp=10;
    if tp11 eq 1 then intervalp=11;
    if tp12 eq 1 then intervalp=12;
    if tp13 eq 1 then intervalp=13;
    if tp14 eq 1 then intervalp=14;
    if tp15 eq 1 then intervalp=15;
    if tp16 eq 1 then intervalp=16;
    if tp17 eq 1 then intervalp=17;
    if tp18 eq 1 then intervalp=18;
   
    phmsstatus=0;
    pmhever=0; 
    ocever=0;
    breast_ca=0;

    sex=1;
    cohort=3;
    newid=(id+3000000);        /*IDs for pooling- HPFS begins with 3*/
run;


data pool;
    set nhs1data nhs2data hpfsdata;  /* Keep all variables - no exclusions yet */
run;

*proc datasets;
*    delete nhs1data nhs2data hpfsdata;
*run;


/******************************************************************************/
/* DATASET 1: Time-Varying Cumulative Average                                */
/******************************************************************************/
data pool_timevarying;
    set pool;
    if lcssvcon = . then delete;            /* Must have 6mo-4yr diet */
    if time_end <= time_begin then delete; /* zero-length final intervals*/
    age_dx2=age_dx**2;
    dtdx_cr2=dtdx_cr**2;
    *if stage_cr = . then stage_cr = 1;
run;

/* Rank quintiles within the POOLED time-varying sample */
proc sort data=pool_timevarying; by intervalp; run;

proc rank data=pool_timevarying out=pool_timevarying groups=4;
    by intervalp;
    var   pre_aheim_fix calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon ulfssvcon hlfssvcon alfssvcon vlfssvcon;
    ranks pre_aheim_fixq calorsvq aheimsvq lcssvq alcssvq vlcssvq ulcssvq hlcssvq lfssvq ulfssvq hlfssvq alfssvq vlfssvq;
    run;

/* Calculate and store standard deviation in a macro variable */
proc sql noprint; 
    select std(lcssvcon) into :std_lcs_sv from pool_timevarying; 
    select std(alcssvcon) into :std_alcs_sv from pool_timevarying;
    select std(vlcssvcon) into :std_vlcs_sv from pool_timevarying;
    select std(ulcssvcon) into :std_ulcs_sv from pool_timevarying;
    select std(hlcssvcon) into :std_hlcs_sv from pool_timevarying;
    select std(lfssvcon) into :std_lfs_sv from pool_timevarying;
    select std(ulfssvcon) into :std_ulfs_sv from pool_timevarying;
    select std(hlfssvcon) into :std_hlfs_sv from pool_timevarying;
    select std(alfssvcon) into :std_alfs_sv from pool_timevarying;
    select std(vlfssvcon) into :std_vlfs_sv from pool_timevarying;
quit;

data pool_timevarying;
    set pool_timevarying;
    array qvar pre_aheim_fixq calorsvq aheimsvq lcssvq alcssvq vlcssvq ulcssvq hlcssvq lfssvq ulfssvq hlfssvq alfssvq vlfssvq;
    do over qvar;
        qvar = qvar + 1;
        if qvar = . then qvar = 1; 
    end;


data pool_timevarying;
    set pool_timevarying end=_end_;
    
/* For SD */
lcs_sd = lcssvcon / &std_lcs_sv;
alcs_sd = alcssvcon / &std_alcs_sv;
vlcs_sd = vlcssvcon / &std_vlcs_sv;
ulcs_sd = ulcssvcon / &std_ulcs_sv;
hlcs_sd = hlcssvcon / &std_hlcs_sv;
lfs_sd = lfssvcon / &std_lfs_sv;
ulfs_sd = ulfssvcon / &std_ulfs_sv;
hlfs_sd = hlfssvcon / &std_hlfs_sv;
alfs_sd = alfssvcon / &std_alfs_sv;
vlfs_sd = vlfssvcon / &std_vlfs_sv;

/* For trend average */
if lcssvq=1 then lcssvtrend= 7;
if lcssvq=2 then lcssvtrend= 12;
if lcssvq=3 then lcssvtrend= 17;
if lcssvq=4 then lcssvtrend= 22;

if alcssvq=1 then alcssvtrend= 6;
if alcssvq=2 then alcssvtrend= 12;
if alcssvq=3 then alcssvtrend= 17;
if alcssvq=4 then alcssvtrend= 23;

if vlcssvq=1 then vlcssvtrend= 9;
if vlcssvq=2 then vlcssvtrend= 13;
if vlcssvq=3 then vlcssvtrend= 16;
if vlcssvq=4 then vlcssvtrend= 21;

if ulcssvq=1 then ulcssvtrend= 6;
if ulcssvq=2 then ulcssvtrend= 13;
if ulcssvq=3 then ulcssvtrend= 17;
if ulcssvq=4 then ulcssvtrend= 22;

if hlcssvq=1 then hlcssvtrend= 9;
if hlcssvq=2 then hlcssvtrend= 13;
if hlcssvq=3 then hlcssvtrend= 16;
if hlcssvq=4 then hlcssvtrend= 21;

if lfssvq=1 then lfssvtrend= 8;
if lfssvq=2 then lfssvtrend= 12.5;
if lfssvq=3 then lfssvtrend= 17;
if lfssvq=4 then lfssvtrend= 22;

if ulfssvq=1 then ulfssvtrend= 9;
if ulfssvq=2 then ulfssvtrend= 13;
if ulfssvq=3 then ulfssvtrend= 17;
if ulfssvq=4 then ulfssvtrend= 21;

if hlfssvq=1 then hlfssvtrend= 6;
if hlfssvq=2 then hlfssvtrend= 12;
if hlfssvq=3 then hlfssvtrend= 17;
if hlfssvq=4 then hlfssvtrend= 24.5;

if alfssvq=1 then alfssvtrend= 8;
if alfssvq=2 then alfssvtrend= 13;
if alfssvq=3 then alfssvtrend= 17;
if alfssvq=4 then alfssvtrend= 21;

if vlfssvq=1 then vlfssvtrend= 6;
if vlfssvq=2 then vlfssvtrend= 12;
if vlfssvq=3 then vlfssvtrend= 17.5;
if vlfssvq=4 then vlfssvtrend= 24;


%indic3(vbl=calorsvq,  prefix=calorsvq,     min=2, max=4, reflev=1, missing=., usemiss=0);
%indic3(vbl=pre_aheim_fixq,  prefix=pre_aheim_fixq,     min=2, max=4, reflev=1, missing=., usemiss=0);
%indic3(vbl=aheimsvq,  prefix=aheimsvq,     min=2, max=4, reflev=1, missing=., usemiss=0);
%indic3(vbl=lcssvq,    prefix=lcssvq,     min=2, max=4, reflev=1, missing=., usemiss=0);
%indic3(vbl=alcssvq,   prefix=alcssvq,     min=2, max=4, reflev=1, missing=., usemiss=0);
%indic3(vbl=vlcssvq,   prefix=vlcssvq,     min=2, max=4, reflev=1, missing=., usemiss=0);
%indic3(vbl=ulcssvq,   prefix=ulcssvq,     min=2, max=4, reflev=1, missing=., usemiss=0);
%indic3(vbl=hlcssvq,   prefix=hlcssvq,     min=2, max=4, reflev=1, missing=., usemiss=0);
%indic3(vbl=lfssvq,    prefix=lfssvq,     min=2, max=4, reflev=1, missing=., usemiss=0);
%indic3(vbl=ulfssvq,   prefix=ulfssvq,     min=2, max=4, reflev=1, missing=., usemiss=0);
%indic3(vbl=hlfssvq,   prefix=hlfssvq,     min=2, max=4, reflev=1, missing=., usemiss=0);
%indic3(vbl=alfssvq,   prefix=alfssvq,     min=2, max=4, reflev=1, missing=., usemiss=0);
%indic3(vbl=vlfssvq,   prefix=vlfssvq,     min=2, max=4, reflev=1, missing=., usemiss=0);

%indic3(vbl=agegp, prefix=agegp, min=2, max=8, reflev=1, missing=., usemiss=0);
%indic3(vbl=smkstatus, prefix=smkstatus, min=2, max=3, reflev=1, missing=., usemiss=1);
%indic3(vbl=smkstaa, prefix=smkstaa, min=2, max=3, reflev=1, missing=., usemiss=1);
%indic3(vbl=bmic3cat, prefix=bmic3cat, min=2, max=3, reflev=1, missing=., usemiss=1);
%indic3(vbl=bmib3cat, prefix=bmib3cat, min=2, max=3, reflev=1, missing=., usemiss=1);
%indic3(vbl=bmic9, prefix=bmic9, min=2, max=9, reflev=1, usemiss=1);
%indic3(vbl=bmib9, prefix=bmib9, min=2, max=9, reflev=1, usemiss=1);
%indic3(vbl=bmid9, prefix=bmid9, min=2, max=9, reflev=1, usemiss=1);

%indic3(vbl=postact, prefix=postact, min=2, max=4, reflev=1, usemiss=0);
%indic3(vbl=actc, prefix=actc, min=2, max=5, reflev=1, usemiss=0);
%indic3(vbl=alc, prefix=alc, min=2, max=6, reflev=1, missing=9, usemiss=1);
%indic3(vbl=phmsstatus, prefix=phmsstatus, min=2, max=4, reflev=1, missing=0, usemiss=1);
*if alcovc_post = . then alcovc_post = 9;
*%indic3(vbl=alcovc_post, prefix=alcovc_post, min=2, max=5, reflev=1, missing=9, usemiss=1);

%indic3(vbl=alcovc_post, prefix=alcovc_post, min=2, max=5, reflev=1, missing=., usemiss=0);
%indic3(vbl=crcduration_c, prefix=crcduration_c, min=2, max=4, reflev=1, missing=., usemiss=0);
%indic3(vbl=subsite_cr, prefix=subsite_cr, min=2, max=4, reflev=1, missing=., usemiss=0);
%indic3(vbl=stage_cr, prefix=stage_cr, min=2, max=3, reflev=1, missing=., usemiss=0);
%indic3(vbl=grade_cr, prefix=grade_cr, min=2, max=4, reflev=1, missing=., usemiss=0);
run;



******************************************************;
**      Getting Medians for  TEST FOR TREND         **;
******************************************************;


/* Getting median values according to quintiles */

%macro m(vgrp, var, data);
proc means data=&data median n nmiss;  
class &vgrp;
var &var;
run;
%mend; 

%m(lcssvq,  lcssvcon, pool_timevarying);
%m(alcssvq, alcssvcon, pool_timevarying);
%m(vlcssvq, vlcssvcon, pool_timevarying);
%m(ulcssvq, ulcssvcon, pool_timevarying);
%m(hlcssvq, hlcssvcon, pool_timevarying);
%m(lfssvq,  lfssvcon, pool_timevarying);
%m(ulfssvq, ulfssvcon, pool_timevarying);
%m(hlfssvq, hlfssvcon, pool_timevarying);
%m(alfssvq, alfssvcon, pool_timevarying);
%m(vlfssvq, vlfssvcon, pool_timevarying);


***Post-dx follow-up ;
*** >>> Median Follow up time;
proc sql;
   create table id_fup as 
   select ID, person_time, sum(person_time) as total_person_time
   from pool_timevarying
   group by ID;
quit;


**** Calculate meadian fup;
proc means data=id_fup noprint mean median p25 p75;
  var total_person_time;
  output out=fup_p
    mean=mean_pt
    median=median_pt
    p25=p25_pt
    p75=p75_pt;
run;

data fup_p_years;
  set fup_p;
  mean_years   = mean_pt/12;
  median_years = median_pt/12;
  q1_years     = p25_pt/12;
  q3_years     = p75_pt/12;
  iqr_years    = (p75_pt - p25_pt)/12;
run;

proc print data=fup_p_years noobs;
  var mean_years median_years q1_years q3_years iqr_years;
run;

