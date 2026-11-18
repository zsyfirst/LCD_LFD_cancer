/*****************************************************
Purpose:       post-dx low carb, low fat, and mortality
Study design:  prospective cohort (NHS: 1986-2016.6, NHSII: 1991-2019.6, HPFS: 1986-2020.1)

Contents:      Table S12 - Characteristics of the Colorectal Cancer Survivors, Overall and By Cohort
Date:          10/2026
******************************************************/

/**** Call in Cleaned Datasets*******/

%include '/lowcarb_cancer/revision/data_postdx_crc.sas';


proc format;
value smkstatusf
1='never smker'
2='past smker'
3='curr smker'
run;

proc format;
value bmic3catf
1='bmi <24.9'
2='bmi 24.9-29.9'
3='bmi 30.0'
run;


proc format;
value alcovc_postf
1='none'
2='0.1-4.9 g/d'
3='5.0-14.9 g/d'
4='15.0-29.9 g/d'
5='30+ g/d'
run;


proc format;
value stage_crf
1='stage I'
2='stage II'
3='stage III'
run;

proc format;
value subsite_crf
1='proximal colon'
2='distal colon'
3='rectum'
4='Unspecified'
run;

proc format;
value grade_crf
1='Well Differentiated'
2='Grade 2 Moderately Diff'
3='Grade 3 Poorly Diff'
4='Unspecified'
run;

proc format;
value crcduration_cf
1='0-5 years'
2='5-10 years'
3='10-15 years'
4='>15 years'
run;

/*
proc format;
value phmsstatusf
1='pre mnp'
2='never/unknown pmhuser, post mnp'
3='curr pmh user and post mnp'
4='past pmhuser and post mnp'
run;
*/

proc sort data=pool_timevarying; by newid period; run;

data pooltable1; set pool_timevarying; 
     by newid;
     if first.newid;

    keep    cohort agecon agecond age_dx  agegp agesub white smkstatus bmicon bmi_dx bmic3cat actcon active alcosvcon alcocon neverdrinker 
            mifh cafh highbp hightc hbpd chold mv asp pmhever ocever
            alcovc_post stage_cr subsite_cr grade_cr crcduration crcduration_c 
            nSES 
            calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon ulfssvcon hlfssvcon alfssvcon vlfssvcon
            ;
    format  smkstatus smkstatusf.;
    format  bmic3cat  bmic3catf. ;
    format  alcovc_post  alcovc_postf. ;
    format  stage_cr  stage_crf. ;
    format  subsite_cr  subsite_crf. ;
    format  grade_cr  grade_crf. ;
    format  crcduration_c crcduration_cf.; 
    *format phmsstatus  phmsstatusf. ;
run;


%table1(data=pooltable1,
        agegroup =  agegp,
        exposure =  cohort,
        varlist  =  agecon age_dx agecond agesub white smkstatus smkstaa bmicon bmi_dx bmic3cat actcon active alcocon neverdrinker  alcosvcon alcovc_post 
                    mifh cafh highbp hightc hbpd chold mv asp pmhever ocever
                    crcduration crcduration_c stage_cr subsite_cr grade_cr
                    nSES
                    calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon ulfssvcon hlfssvcon alfssvcon vlfssvcon
                    ,                                                  
        noadj    =  agecon age_dx agecond,     
        cat      =  agesub white active neverdrinker mifh cafh highbp hightc hbpd chold mv asp pmhever ocever, 
        poly     =  smkstatus smkstaa bmic3cat alcovc_post crcduration_c stage_cr subsite_cr grade_cr,
        rtftitle =  Table1 Characteristics at first post-diagnosis assessment, 
        landscape=  F,                               
        file     =  %str(table/pool.table1_bycohort.rtf),
        dec      =  1,
        uselbl   =  T,
        nortf    = F);
run;

%table1(data=pooltable1,
        agegroup =  agegp,
        /*exposure =  cohort,*/
        varlist  =  agecon age_dx agecond agesub white smkstatus smkstaa bmicon bmi_dx bmic3cat actcon active alcocon neverdrinker  alcosvcon alcovc_post 
                    mifh cafh highbp hightc hbpd chold mv asp pmhever ocever
                    crcduration crcduration_c stage_cr subsite_cr grade_cr
                    nSES
                    calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon ulfssvcon hlfssvcon alfssvcon vlfssvcon
                    ,                                                
        noadj    =  agecon age_dx agecond,   
        cat      =  agesub white active neverdrinker mifh cafh highbp hightc hbpd chold mv asp pmhever ocever, 
        poly     =  smkstatus smkstaa bmic3cat alcovc_post crcduration_c stage_cr subsite_cr grade_cr,
        rtftitle =  Table1 Characteristics at first post-diagnosis assessment, 
        landscape=  F,                               
        file     =  %str(table/pool.table1_overall.rtf),
        dec      =  1,
        uselbl   =  T,
        noexp=T,
        nortf    = F);
run;