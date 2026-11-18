/*****************************************************
Purpose:       low carb, low fat, and risk of Obesity-related cancer
Study design:  prospective cohort (NHS: 1986-2022.6, NHSII: 1991-2023.6, HPFS: 1986-2024.1)
Endpoints:     Obesity-related cancer and site_specific cancer
Exposures:     lcs, alcs, vlcs, ulcs, hlcs, lfs, ulfs, hlfs (NHS: 1986-2010, NHSII: 1991-2019, HPFS: 1986-2018)
Exclusion:     baseline cancer or cvd or diabetes, exposure missing, 
               or reported implausible SFFQ data (<800 or >4200 kcal/d for men and <600 or >3500 kcal/d for women)
Covariates:    age, interval, cohort, white, energy quintile, smoking, alcohol, physical activity, cancer family history, 
               hormone use, oral contraceptive use, multivitamin use, aspirin use, 
               baseline BMI, screening

Date:  10/2026

******************************************************/

/************************ Call in Cleaned Datasets*******************************************************/
%include 'lowcarb_cancer/revision/hpfs.data.inc.sas';
%include 'lowcarb_cancer/revision/nhs.data.inc.sas';
%include 'lowcarb_cancer/revision/nhs2.data.inc.sas';

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

    sex=1;
    cohort=3;
    newid=(id+3000000);        /*IDs for pooling- HPFS begins with 3*/
run;


data pool;
    set nhs1data nhs2data hpfsdata end=_end_;
    
    keep 
    newid sex cohort fup_allca  id 

    /**Time scale**/
    tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18
    time1 time2 time3 time4 time5 time6 time7 time8 time9 time10 time11 time12 time13 time14 time15 time16 time17 time18
    
    /*interval*/
    intervalp  period
    
    /**Outcome******/
    dead tdead 
    dead_can dead_cvd dead_resp dead_oth tdead_can tdead_cvd tdead_resp tdead_oth 
    can_obese tobeca can_nonobese tnonobeca 
    cvd tcvd chd tchd str tstr 

    dtdxchd dtdxmi dtdth dtdxca dtdx_allca dtdx_obeca dtdx_nonobeca 
    can_: t: dtdx_:   
    obeca_site can_obese_site         
    /*esoph colorect stomach liver gallbladder panc
    brca endo ovca kidney thyroid myeloma menin
    lung bladder brain oral smallintest larynx nose melanoma
    diges uri myeleuk resp
    smoke r_smoke allca*/

    /**Exposure*****/
    lcsscon   alcsscon  vlcsscon  ulcsscon  hlcsscon  lfsscon  ulfsscon  hlfsscon  alfsscon  vlfsscon
    lcssvcon  alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon ulfssvcon hlfssvcon alfssvcon vlfssvcon
    
    basecarb basehcarb basesumref_carb baseprot baseaprot basevprot basetfat baseafat basevfat
    carb hcarb sumref_carb prot aprot vprot tfat afat vfat
    carb_svcon hcarb_svcon sumref_carb_svcon prot_svcon aprot_svcon vprot_svcon tfat_svcon afat_svcon vfat_svcon

    /**Covariates***/
    calorscon  aheimscon  
    calorsvcon aheimsvcon 
    agecon agesub agegp 
    white  dbfh  mifh cafh
    alcocon alc neverdrinker
    actcon  actc active
    smkstatus smknever
    bmicon bmic3cat bmib3cat bmic25 bmic30  bmic9  bmib9  
    basehighbp basehightc highbp hightc diabetes
    mv asp       
    phmsstatus pmhever ocever

    basewhgrn basergrain basefrui  baseleg basepota baseveg_noleg  basefish baserpmeats  basessb basedairy basesweets basepoultry
    whgrn_svcon rgrain_svcon frui_svcon leg_svcon pota_svcon veg_noleg_svcon fish_svcon rpmeats_svcon ssb_svcon 
    dairy_svcon swtsnochoc_svcon poul_svcon  
    mam endo2c endo3c
    nSES 
    ;
run;


proc datasets;
    delete  nhs1data nhs2data hpfsdata;
run;
* for casePY in pooldatasets;

data pool;
    set pool;
    where intervalp ne . ;
run;
proc sort data=pool; by intervalp; run;

proc rank data=pool  groups=5  out=pool;
by intervalp;
 var   calorscon aheimscon lcsscon alcsscon vlcsscon ulcsscon hlcsscon lfsscon ulfsscon hlfsscon alfsscon vlfsscon
       calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon ulfssvcon hlfssvcon alfssvcon vlfssvcon;
 ranks calorsq aheimsq lcssq alcssq vlcssq ulcssq hlcssq lfssq ulfssq hlfssq alfssq vlfssq
       calorsvq aheimsvq lcssvq alcssvq vlcssvq ulcssvq hlcssvq lfssvq ulfssvq hlfssvq alfssvq vlfssvq;
run;

/* Calculate and store standard deviation in a macro variable */
proc sql noprint; select std(lcssvcon)  into :std_lcs_sv from pool; quit;
proc sql noprint; select std(alcssvcon) into :std_alcs_sv from pool; quit;
proc sql noprint; select std(vlcssvcon) into :std_vlcs_sv from pool; quit;
proc sql noprint; select std(ulcssvcon) into :std_ulcs_sv from pool; quit;
proc sql noprint; select std(hlcssvcon) into :std_hlcs_sv from pool; quit;
proc sql noprint; select std(lfssvcon)  into :std_lfs_sv  from pool; quit;
proc sql noprint; select std(ulfssvcon) into :std_ulfs_sv from pool; quit;
proc sql noprint; select std(hlfssvcon) into :std_hlfs_sv from pool; quit;
proc sql noprint; select std(alfssvcon) into :std_alfs_sv from pool; quit;
proc sql noprint; select std(vlfssvcon) into :std_vlfs_sv from pool; quit;

proc sql noprint; select std(lcsscon)  into :std_lcs_s  from pool; quit;
proc sql noprint; select std(alcsscon) into :std_alcs_s from pool; quit;
proc sql noprint; select std(vlcsscon) into :std_vlcs_s from pool; quit;
proc sql noprint; select std(ulcsscon) into :std_ulcs_s from pool; quit;
proc sql noprint; select std(hlcsscon) into :std_hlcs_s from pool; quit;
proc sql noprint; select std(lfsscon)  into :std_lfs_s  from pool; quit;
proc sql noprint; select std(ulfsscon) into :std_ulfs_s from pool; quit;
proc sql noprint; select std(hlfsscon) into :std_hlfs_s from pool; quit;
proc sql noprint; select std(alfsscon) into :std_alfs_s from pool; quit;
proc sql noprint; select std(vlfsscon) into :std_vlfs_s from pool; quit;

data pool;
   set pool end=_end_;
   array qvar1 calorsq aheimsq lcssq alcssq vlcssq ulcssq hlcssq lfssq ulfssq hlfssq alfssq vlfssq
               calorsvq aheimsvq lcssvq alcssvq vlcssvq ulcssvq hlcssvq lfssvq ulfssvq hlfssvq alfssvq vlfssvq;
   do over qvar1;
      qvar1=qvar1+1;                
   end;

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

lcs_sd_s = lcsscon / &std_lcs_s;
alcs_sd_s = alcsscon / &std_alcs_s;
vlcs_sd_s = vlcsscon / &std_vlcs_s;
ulcs_sd_s = ulcsscon / &std_ulcs_s;
hlcs_sd_s = hlcsscon / &std_hlcs_s;
lfs_sd_s = lfsscon / &std_lfs_s;
ulfs_sd_s = ulfsscon / &std_ulfs_s;
hlfs_sd_s = hlfsscon / &std_hlfs_s;
alfs_sd_s = alfsscon / &std_alfs_s;
vlfs_sd_s = vlfsscon / &std_vlfs_s;


/* For trend average */
if lcssvq=1 then lcssvtrend= 7;
if lcssvq=2 then lcssvtrend= 12;
if lcssvq=3 then lcssvtrend= 15;
if lcssvq=4 then lcssvtrend= 18.5;
if lcssvq=5 then lcssvtrend= 23;

if alcssvq=1 then alcssvtrend= 6;
if alcssvq=2 then alcssvtrend= 11;
if alcssvq=3 then alcssvtrend= 15;
if alcssvq=4 then alcssvtrend= 19;
if alcssvq=5 then alcssvtrend= 24;

if vlcssvq=1 then vlcssvtrend= 9;
if vlcssvq=2 then vlcssvtrend= 13;
if vlcssvq=3 then vlcssvtrend= 15;
if vlcssvq=4 then vlcssvtrend= 17.5;
if vlcssvq=5 then vlcssvtrend= 21;

if ulcssvq=1 then ulcssvtrend= 7;
if ulcssvq=2 then ulcssvtrend= 12;
if ulcssvq=3 then ulcssvtrend= 15;
if ulcssvq=4 then ulcssvtrend= 18.5;
if ulcssvq=5 then ulcssvtrend= 23;

if hlcssvq=1 then hlcssvtrend= 9;
if hlcssvq=2 then hlcssvtrend= 12.5;
if hlcssvq=3 then hlcssvtrend= 15;
if hlcssvq=4 then hlcssvtrend= 17.5;
if hlcssvq=5 then hlcssvtrend= 21;

if lfssvq=1 then lfssvtrend= 8;
if lfssvq=2 then lfssvtrend= 12;
if lfssvq=3 then lfssvtrend= 15;
if lfssvq=4 then lfssvtrend= 18;
if lfssvq=5 then lfssvtrend= 22;

if ulfssvq=1 then ulfssvtrend= 9.5;
if ulfssvq=2 then ulfssvtrend= 13;
if ulfssvq=3 then ulfssvtrend= 15;
if ulfssvq=4 then ulfssvtrend= 17;
if ulfssvq=5 then ulfssvtrend= 20;

if hlfssvq=1 then hlfssvtrend= 6;
if hlfssvq=2 then hlfssvtrend= 11;
if hlfssvq=3 then hlfssvtrend= 15;
if hlfssvq=4 then hlfssvtrend= 19;
if hlfssvq=5 then hlfssvtrend= 24;

if alfssvq=1 then alfssvtrend= 9;
if alfssvq=2 then alfssvtrend= 12;
if alfssvq=3 then alfssvtrend= 15;
if alfssvq=4 then alfssvtrend= 17.5;
if alfssvq=5 then alfssvtrend= 21;

if vlfssvq=1 then vlfssvtrend= 6;
if vlfssvq=2 then vlfssvtrend= 11;
if vlfssvq=3 then vlfssvtrend= 15;
if vlfssvq=4 then vlfssvtrend= 19;
if vlfssvq=5 then vlfssvtrend= 24;

/* For trend simple updated */
if lcssq=1 then lcsstrend= 5;
if lcssq=2 then lcsstrend= 11;
if lcssq=3 then lcsstrend= 15;
if lcssq=4 then lcsstrend= 19;
if lcssq=5 then lcsstrend= 25;

if alcssq=1 then alcsstrend= 4;
if alcssq=2 then alcsstrend= 10;
if alcssq=3 then alcsstrend= 15;
if alcssq=4 then alcsstrend= 20;
if alcssq=5 then alcsstrend= 26;

if vlcssq=1 then vlcsstrend= 8;
if vlcssq=2 then vlcsstrend= 12;
if vlcssq=3 then vlcsstrend= 15;
if vlcssq=4 then vlcsstrend= 18;
if vlcssq=5 then vlcsstrend= 22;

if ulcssq=1 then ulcsstrend= 5;
if ulcssq=2 then ulcsstrend= 11;
if ulcssq=3 then ulcsstrend= 15;
if ulcssq=4 then ulcsstrend= 19;
if ulcssq=5 then ulcsstrend= 24;

if hlcssq=1 then hlcsstrend= 8;
if hlcssq=2 then hlcsstrend= 12;
if hlcssq=3 then hlcsstrend= 15;
if hlcssq=4 then hlcsstrend= 18;
if hlcssq=5 then hlcsstrend= 23;

if lfssq=1 then lfsstrend= 7;
if lfssq=2 then lfsstrend= 11;
if lfssq=3 then lfsstrend= 15;
if lfssq=4 then lfsstrend= 19;
if lfssq=5 then lfsstrend= 23;

if ulfssq=1 then ulfsstrend= 8;
if ulfssq=2 then ulfsstrend= 12;
if ulfssq=3 then ulfsstrend= 15;
if ulfssq=4 then ulfsstrend= 18;
if ulfssq=5 then ulfsstrend= 22;

if hlfssq=1 then hlfsstrend= 4;
if hlfssq=2 then hlfsstrend= 10;
if hlfssq=3 then hlfsstrend= 15;
if hlfssq=4 then hlfsstrend= 20;
if hlfssq=5 then hlfsstrend= 26;

if alfssq=1 then alfsstrend= 7;
if alfssq=2 then alfsstrend= 12;
if alfssq=3 then alfsstrend= 15;
if alfssq=4 then alfsstrend= 18;
if alfssq=5 then alfsstrend= 22;

if vlfssq=1 then vlfsstrend= 4;
if vlfssq=2 then vlfsstrend= 10;
if vlfssq=3 then vlfsstrend= 15;
if vlfssq=4 then vlfsstrend= 20;
if vlfssq=5 then vlfsstrend= 26;


%indic3(vbl=calorsq,  prefix=calorsq,     min=2, max=5, reflev=1, missing=., usemiss=0);
%indic3(vbl=aheimsq,  prefix=aheimsq,     min=2, max=5, reflev=1, missing=., usemiss=0);
%indic3(vbl=lcssq,    prefix=lcssq,     min=2, max=5, reflev=1, missing=., usemiss=0);
%indic3(vbl=alcssq,   prefix=alcssq,     min=2, max=5, reflev=1, missing=., usemiss=0);
%indic3(vbl=vlcssq,   prefix=vlcssq,     min=2, max=5, reflev=1, missing=., usemiss=0);
%indic3(vbl=ulcssq,   prefix=ulcssq,     min=2, max=5, reflev=1, missing=., usemiss=0);
%indic3(vbl=hlcssq,   prefix=hlcssq,     min=2, max=5, reflev=1, missing=., usemiss=0);
%indic3(vbl=lfssq,    prefix=lfssq,     min=2, max=5, reflev=1, missing=., usemiss=0);
%indic3(vbl=ulfssq,   prefix=ulfssq,     min=2, max=5, reflev=1, missing=., usemiss=0);
%indic3(vbl=hlfssq,   prefix=hlfssq,     min=2, max=5, reflev=1, missing=., usemiss=0);
%indic3(vbl=alfssq,   prefix=alfssq,     min=2, max=5, reflev=1, missing=., usemiss=0);
%indic3(vbl=vlfssq,   prefix=vlfssq,     min=2, max=5, reflev=1, missing=., usemiss=0);

%indic3(vbl=calorsvq,  prefix=calorsvq,     min=2, max=5, reflev=1, missing=., usemiss=0);
%indic3(vbl=aheimsvq,  prefix=aheimsvq,     min=2, max=5, reflev=1, missing=., usemiss=0);

%indic3(vbl=lcssvq,    prefix=lcssvq,     min=2, max=5, reflev=1, missing=., usemiss=0);
%indic3(vbl=alcssvq,   prefix=alcssvq,     min=2, max=5, reflev=1, missing=., usemiss=0);
%indic3(vbl=vlcssvq,   prefix=vlcssvq,     min=2, max=5, reflev=1, missing=., usemiss=0);
%indic3(vbl=ulcssvq,   prefix=ulcssvq,     min=2, max=5, reflev=1, missing=., usemiss=0);
%indic3(vbl=hlcssvq,   prefix=hlcssvq,     min=2, max=5, reflev=1, missing=., usemiss=0);
%indic3(vbl=lfssvq,    prefix=lfssvq,     min=2, max=5, reflev=1, missing=., usemiss=0);
%indic3(vbl=ulfssvq,   prefix=ulfssvq,     min=2, max=5, reflev=1, missing=., usemiss=0);
%indic3(vbl=hlfssvq,   prefix=hlfssvq,     min=2, max=5, reflev=1, missing=., usemiss=0);
%indic3(vbl=alfssvq,   prefix=alfssvq,     min=2, max=5, reflev=1, missing=., usemiss=0);
%indic3(vbl=vlfssvq,   prefix=vlfssvq,     min=2, max=5, reflev=1, missing=., usemiss=0);

%indic3(vbl= agegp,     prefix=agegp, min=2, max=8, reflev=1, missing=., usemiss=0);
%indic3(vbl= smkstatus, prefix= smkstatus, min=2, max=3, reflev=1, missing=., usemiss=0);
%indic3(vbl= bmic3cat,  prefix= bmic3cat, min=2, max=3, reflev=1, missing=., usemiss=0);
%indic3(vbl= bmib3cat,  prefix= bmib3cat, min=2, max=3, reflev=1, missing=., usemiss=0);
%indic3(vbl=bmic9, prefix=bmic9, min=2, max=9, reflev=1, usemiss=1);
%indic3(vbl=bmib9, prefix=bmib9, min=2, max=9, reflev=1, usemiss=1);

%indic3(vbl=actc, prefix=actc, min=2, max=5, reflev=1, usemiss=0,
        label1='<3 total mets/wk',
        label2='3 to < 9 mets/wk',
        label3='9 to < 18 mets/wk',
        label4='18 to < 27 mets/wk',
        label5='27+ mets/wk');

%indic3(vbl=alc, prefix=alc, min=2, max=6, reflev=1, usemiss=0,
                label1='0 g/d',
                label2='0.1-4.9 g/d',
                label3='5.0-9.9 g/d',
                label4='10.0-14.9 g/d',
                label5='15.0-29.9 g/d',
                label6='30+ g/d');

%indic3(vbl=phmsstatus, prefix=phmsstatus, min=2, max=4, reflev=1,
                missing=0, usemiss=1,
                label1='premenopause',
                label2='pmh never use',
                label3='pmh current use',
                label4='pmh past use');

if endo2c =1 or mam =1 then  screen_cancer = 1;
else screen_cancer=0;
run;


/**************************************** Descriptive statistics *****************************************/
proc contents data=pool varnum;
run;

proc freq data=pool;
    tables agegp smkstatus bmic3cat bmib3cat actc alc phmsstatus / missing;
run;


proc format;
    value obecafmt
        0  = 'No obesity-related cancer'
        1  = 'Esophageal cancer'
        2  = 'Gastric cancer'
        3  = 'Colorectal cancer'
        4  = 'Liver cancer'
        5  = 'Gallbladder cancer'
        6  = 'Pancreatic cancer'
        7  = 'Breast cancer'
        8  = 'Endometrial cancer'
        9  = 'Ovarian cancer'
        10 = 'Kidney cancer'
        11 = 'Thyroid cancer'
        12 = 'Multiple myeloma'
        13 = 'Meningioma';
run;

proc freq data=pool;
    tables can_obese_site / missing;
    format can_obese_site obecafmt.;
run;

/* Confirm that every composite event has exactly one subtype */
proc freq data=pool;
    tables can_obese*can_obese_site / missing;
run;

/* There should be no composite events with subtype=0 */
proc sql;
    select count(*) as n_problem
    from pool
    where can_obese=1 and can_obese_site=0;
quit;

/* Confirm event counts by cohort */
proc freq data=pool;
    tables cohort*can_obese_site / missing;
    where can_obese_site > 0;
    format can_obese_site obecafmt.;
run;

/* Confirm only one obesity-related cancer event per participant */
proc sql;
    create table check_multiple_events as
    select newid,
           sum(can_obese_site > 0) as n_obeca_events
    from pool
    group by newid
    having calculated n_obeca_events > 1;
quit;

proc print data=check_multiple_events(obs=20);
run;



******************************************************;
**      Getting Medians for  TEST FOR TREND         **;
******************************************************;


/* Getting median values according to quintiles */
/* pool */
%macro m(vgrp, var);
proc means data=pool median n nmiss;  
class &vgrp;
var &var;
run;
%mend; 

%m(lcssvq,  lcssvcon);
%m(alcssvq, alcssvcon);
%m(vlcssvq, vlcssvcon);
%m(ulcssvq, ulcssvcon);
%m(hlcssvq, hlcssvcon);
%m(lfssvq,  lfssvcon);
%m(ulfssvq, ulfssvcon);
%m(hlfssvq, hlfssvcon);
%m(alfssvq, alfssvcon);
%m(vlfssvq, vlfssvcon);


%m(lcssq,  lcsscon);
%m(alcssq, alcsscon);
%m(vlcssq, vlcsscon);
%m(ulcssq, ulcsscon);
%m(hlcssq, hlcsscon);
%m(lfssq,  lfsscon);
%m(ulfssq, ulfsscon);
%m(hlfssq, hlfsscon);
%m(alfssq, alfsscon);
%m(vlfssq, vlfsscon);


/****************************************************
    Calculate Median Follow-up Time
****************************************************/

proc sql;
  create table id_fup as
  select ID, sum(fup_allca) as total_fup_months
  from pool
  group by ID;
quit;

proc univariate data=id_fup noprint;
  var total_fup_months;
  output out=fup_stats q1=q1 q3=q3 median=med mean=mean;
run;

data fup_stats_years;
  set fup_stats;
  mean_years = mean/12;
  med_years  = med/12;
  q1_years   = q1/12;
  q3_years   = q3/12;
  iqr_years  = (q3-q1)/12;
run;

proc print data=fup_stats_years noobs;
  var mean_years med_years q1_years q3_years iqr_years;
run;


