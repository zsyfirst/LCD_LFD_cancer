/*****************************************************
Purpose:       low carb, low fat, and risk of Obesity-related cancer
Study design:  prospective cohort (NHS: 1986-2022.6, NHSII: 1991-2023.6)

Contents:      Table S12 (Characteristics at first post-diagnosis assessment)
Date:          10/2026
******************************************************/

/**** Call in Cleaned Datasets*******/

%include '/lowcarb_cancer/revision/data_postdx_brc.sas';


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
4='missing'
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
value stagef
1='stage I'
2='stage II'
3='stage III'
run;

proc format;
value gradef
1='Well Differentiated'
2='Grade 2 Moderately Diff'
3='Grade 3 Poorly Diff'
4='Unspecified'
run;

proc format;
value est_statusf
1='ER-'
2='ER+'
3='Missing'
run;

proc format;
value subtypef
1='Luminal A'
2='Luminal B'
3='Her2 positive'
4='Basal-like'
5='Unclassified'
run;

proc format;
value bcduration_cf
1='0-5 years'
2='5-10 years'
3='10-15 yaers'
4='>15 years'
run;

proc format;
value surgf
1='mastectomy'
2='lumpectomy'
3='unknown type'
4='no surgery'
5='Missing'
run;

proc sort data=pool_timevarying; by newid period; run;

data pooltable1; set pool_timevarying; 
	 by newid;
     if first.newid;
    
    if radiotx = 1 then radiotx_01 = 1; else if radiotx = 2 then radiotx_01 = 0;
    if hormtx = 1 then hormtx_01 = 1; else if hormtx = 2 then hormtx_01 = 0;
    if chemotx = 1 then chemotx_01 = 1; else if chemotx = 2 then chemotx_01 = 0;
    
    first_period = period;
    first_time_since_dx = time_begin;

    keep    cohort newid first_period first_time_since_dx
    		agecon agecond agegp agesub white smkstatus 
            wtcont wtcontkg wtc4cont wtc4contkg wtc4contp 
            bmicon bmicond bmic3cat actcon active alcocon alcovc_post neverdrinker 
            mifh cafh highbp hightc hbpd chold mv asp pmhever ocever
            stage  grade subtype est_status2  radiotx_01 hormtx_01  chemotx_01 bcduration bcduration_c
            surg nSES  
            pre_aheim_fix 
            calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon ulfssvcon hlfssvcon alfssvcon vlfssvcon
            ;
    format  smkstatus smkstatusf.;
    format  bmic3cat  bmic3catf. ;
    format  alcovc_post  alcovc_postf. ;
    format  stage  stagef. ;
    format  grade  gradef. ;
    format  subtype  subtypef. ;
    format  est_status2  est_statusf. ;
    format  bcduration_c  bcduration_cf. ;
    format  surg  surgf. ; 
run;

/* CHECK the distribution of first_period and first_time_since_dx */
proc freq data=pooltable1; 
    tables first_period cohort*first_period / missing; 
run;

proc means data=pooltable1 n mean median min max;
    var first_time_since_dx bcduration;
    class cohort;
run;


%table1(data=pooltable1,
        agegroup =  agegp,
        exposure =  cohort,
        varlist  =  agecon age_dx agecond white smkstatus smkstaa wtcontkg wtc4contkg bmi_dx bmicon bmic3cat bmicond actcon active alcocon alcovc_post neverdrinker 
                    mifh cafh highbp hightc hbpd chold mv asp pmhever ocever
                    stage grade subtype est_status2 bcduration bcduration_c  radiotx_01 hormtx_01 chemotx_01 
                    surg nSES  
                    pre_aheim_fix 
                    calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon ulfssvcon hlfssvcon alfssvcon vlfssvcon
                    ,                                                  
        noadj    =  agecon age_dx agecond,     
        cat      =  white active neverdrinker mifh cafh highbp hightc  hbpd chold mv asp pmhever ocever radiotx_01 hormtx_01 chemotx_01     , 
        poly     =  smkstatus bmic3cat alcovc_post stage grade subtype est_status2 bcduration_c surg ,
        rtftitle =  Table1 characteristics at first post-diagnosis assessment, 
        landscape=  F,                               
        file     =  %str(table/table1_postdx.rtf),
        dec      =  1,
        uselbl   =  T,
        nortf    = F);
run;



%table1(data=pooltable1,
        agegroup =  agegp,
        /*exposure= cohort, */
        varlist  =  agecon age_dx agecond white smkstatus smkstaa wtcontkg wtc4contkg bmi_dx bmicon bmic3cat bmicond actcon active alcocon alcovc_post neverdrinker 
                    mifh cafh highbp hightc hbpd chold mv asp pmhever ocever
                    stage grade subtype est_status2 bcduration bcduration_c  radiotx_01 hormtx_01 chemotx_01 
                    surg nSES  
                    pre_aheim_fix
                    calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon ulfssvcon hlfssvcon alfssvcon vlfssvcon
                    ,                                                   
        noadj    =  agecon age_dx agecond ,     
        cat      =  white active neverdrinker mifh cafh highbp hightc hbpd chold mv asp pmhever ocever radiotx_01  hormtx_01 chemotx_01     , 
        poly     =  smkstatus bmic3cat alcovc_post stage grade subtype est_status2 bcduration_c surg ,
        rtftitle =  Table1 characteristics at first post-diagnosis assessment, 
        landscape=  F,                               
        file     =  %str(table/table1_overall_postdx.rtf),
        dec      =  1,
        uselbl   =  T,
        noexp=T,
        nortf    = F);
run;