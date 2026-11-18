/*****************************************************
Purpose:       low carb, low fat, and risk of Obesity-related cancer
Study design:  prospective cohort (NHS: 1986-2022.6, NHSII: 1991-2023.6, HPFS: 1986-2024.1)

Contents:      Table S4 (Covariates at the median of follow-up: 2004 for NHS and HPFS, and 2003 for NHSII [intervalp=10])
Revision:      add nSES
Date:          10/2026
******************************************************/

/**** Call in Cleaned Datasets*******/

%include 'lowcarb_cancer/revision/data.inc.sas';


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

/*
proc format;
value phmsstatusf
1='pre mnp'
2='never/unknown pmhuser, post mnp'
3='curr pmh user and post mnp'
4='past pmhuser and post mnp'
run;
*/


proc sort data=pool; by newid period; run;

data pool;
    set pool;
    if endo2c = 1 then endo2c_bin = 0;
    else if endo2c = 2 then endo2c_bin = 1;
    else endo2c_bin = .;
run;

data pooltableS2; set pool; if intervalp=10;

    keep    cohort agecon agegp agesub white smkstatus bmicon bmib3cat    
            actcon active alcocon neverdrinker 
            mifh cafh basehighbp basehightc mv asp pmhever ocever
            lcssq alcssq vlcssq ulcssq hlcssq lfssq ulfssq hlfssq alfssq vlfssq
            calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon alfssvcon vlfssvcon ulfssvcon hlfssvcon 
            basecarb basehcarb basesumref_carb baseprot baseaprot basevprot basetfat baseafat basevfat
            carb hcarb sumref_carb prot aprot vprot tfat afat vfat
            mam endo2c_bin  screen_cancer nSES ;
    format  smkstatus smkstatusf.;
    format  bmic3cat  bmic3catf. ;
    *format phmsstatus  phmsstatusf. ;
run;

%table1(data=pooltableS2,
        agegroup =  agegp,
        exposure =  cohort,
        varlist  =  agecon agesub white smkstatus bmicon bmib3cat   actcon active alcocon neverdrinker 
                    cafh mv asp pmhever ocever mam endo2c_bin  
                    calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon alfssvcon vlfssvcon ulfssvcon hlfssvcon
                    basecarb basehcarb basesumref_carb baseprot baseaprot basevprot basetfat baseafat basevfat
                    carb hcarb sumref_carb prot aprot vprot tfat afat vfat
                    mifh basehighbp basehightc screen_cancer nSES ,                                                  
        noadj    =  agecon,     
        cat      =  agesub white active neverdrinker mifh cafh basehighbp basehightc   mv asp pmhever ocever mam endo2c_bin  screen_cancer, 
        poly     =  smkstatus bmic3cat ,
        rtftitle =  tableS2 Midpoint characteristics by cohort, 
        landscape=  F,                               
        file     =  %str(table/pool.tableS2_medianfu_cohort.rtf),
        dec      =  1,
        uselbl   =  T,
        nortf    = F);
run;


%table1(data=pooltableS2,
        agegroup =  agegp,
        /*exposure =  cohort,*/
        varlist  =  agecon agesub white smkstatus bmicon bmib3cat   actcon active alcocon neverdrinker 
                    cafh mv asp pmhever ocever mam endo2c_bin  
                    calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon alfssvcon vlfssvcon ulfssvcon hlfssvcon
                    basecarb basehcarb basesumref_carb baseprot baseaprot basevprot basetfat baseafat basevfat
                    carb hcarb sumref_carb prot aprot vprot tfat afat vfat
                    mifh basehighbp basehightc screen_cancer nSES ,                                                 
        noadj    =  agecon,     
        cat      =  agesub white active neverdrinker mifh cafh basehighbp basehightc   mv asp pmhever ocever mam endo2c_bin  screen_cancer, 
        poly     =  smkstatus bmic3cat ,
        rtftitle =  tableS2 Midpoint characteristics overall, 
        landscape=  F,                               
        file     =  %str(table/tableS2_overall_medianfu.rtf),
        dec      =  1,
        uselbl   =  T,
        noexp=T,
        nortf    = F);
run;



%table1(data=pooltableS2,
        agegroup =  agegp,
        exposure =  lcssq,
        varlist  =  agecon agesub white smkstatus bmicon bmib3cat   actcon active alcocon neverdrinker 
                    cafh mv asp pmhever ocever mam endo2c_bin  
                    calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon alfssvcon vlfssvcon ulfssvcon hlfssvcon
                    basecarb basehcarb basesumref_carb baseprot baseaprot basevprot basetfat baseafat basevfat
                    carb hcarb sumref_carb prot aprot vprot tfat afat vfat
                    mifh basehighbp basehightc screen_cancer nSES ,                                                  
        noadj    =  agecon,     
        cat      =  agesub white active neverdrinker mifh cafh basehighbp basehightc   mv asp pmhever ocever mam endo2c_bin  screen_cancer, 
        poly     =  smkstatus bmic3cat ,
        rtftitle =  tableS2 Midpoint characteristics By LCD, 
        landscape=  F,                               
        file     =  %str(table/pool.tableS2_LCD_medianfu.rtf),
        dec      =  1,
        uselbl   =  T,
        nortf    = F);
run;


%table1(data=pooltableS2,
        agegroup =  agegp,
        exposure =  lfssq,
        varlist  =  agecon agesub white smkstatus bmicon bmib3cat   actcon active alcocon neverdrinker 
                    cafh mv asp pmhever ocever mam endo2c_bin  
                    calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon alfssvcon vlfssvcon ulfssvcon hlfssvcon
                    basecarb basehcarb basesumref_carb baseprot baseaprot basevprot basetfat baseafat basevfat
                    carb hcarb sumref_carb prot aprot vprot tfat afat vfat
                    mifh basehighbp basehightc screen_cancer nSES ,                                                  
        noadj    =  agecon,     
        cat      =  agesub white active neverdrinker mifh cafh basehighbp basehightc   mv asp pmhever ocever mam endo2c_bin  screen_cancer, 
        poly     =  smkstatus bmic3cat ,
        rtftitle =  tableS2 Midpoint characteristics By LFD, 
        landscape=  F,                               
        file     =  %str(table/pool.tableS2_LFD_medianfu.rtf),
        dec      =  1,
        uselbl   =  T,
        nortf    = F);
run;




%table1(data=pooltableS2,
        agegroup =  agegp,
        exposure =  hlcssq,
        varlist  =  agecon agesub white smkstatus bmicon bmib3cat   actcon active alcocon neverdrinker 
                    cafh mv asp pmhever ocever mam endo2c_bin  
                    calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon alfssvcon vlfssvcon ulfssvcon hlfssvcon
                    basecarb basehcarb basesumref_carb baseprot baseaprot basevprot basetfat baseafat basevfat
                    carb hcarb sumref_carb prot aprot vprot tfat afat vfat
                    mifh basehighbp basehightc screen_cancer nSES ,                                                  
        noadj    =  agecon,     
        cat      =  agesub white active neverdrinker mifh cafh basehighbp basehightc   mv asp pmhever ocever mam endo2c_bin  screen_cancer, 
        poly     =  smkstatus bmic3cat ,
        rtftitle =  tableS2 Midpoint characteristics By HLCD, 
        landscape=  F,                               
        file     =  %str(table/pool.tableS2_HLCD_medianfu.rtf),
        dec      =  1,
        uselbl   =  T,
        nortf    = F);
run;


%table1(data=pooltableS2,
        agegroup =  agegp,
        exposure =  hlfssq,
        varlist  =  agecon agesub white smkstatus bmicon bmib3cat   actcon active alcocon neverdrinker 
                    cafh mv asp pmhever ocever mam endo2c_bin  
                    calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon alfssvcon vlfssvcon ulfssvcon hlfssvcon
                    basecarb basehcarb basesumref_carb baseprot baseaprot basevprot basetfat baseafat basevfat
                    carb hcarb sumref_carb prot aprot vprot tfat afat vfat
                    mifh basehighbp basehightc screen_cancer nSES ,                                                  
        noadj    =  agecon,     
        cat      =  agesub white active neverdrinker mifh cafh basehighbp basehightc   mv asp pmhever ocever mam endo2c_bin  screen_cancer, 
        poly     =  smkstatus bmic3cat ,
        rtftitle =  tableS2 Midpoint characteristics By HLFD, 
        landscape=  F,                               
        file     =  %str(table/pool.tableS2_HLFD_medianfu.rtf),
        dec      =  1,
        uselbl   =  T,
        nortf    = F);
run;





%table1(data=pooltableS2,
        agegroup =  agegp,
        exposure =  ulcssq,
        varlist  =  agecon agesub white smkstatus bmicon bmib3cat   actcon active alcocon neverdrinker 
                    cafh mv asp pmhever ocever mam endo2c_bin  
                    calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon alfssvcon vlfssvcon ulfssvcon hlfssvcon
                    basecarb basehcarb basesumref_carb baseprot baseaprot basevprot basetfat baseafat basevfat
                    carb hcarb sumref_carb prot aprot vprot tfat afat vfat
                    mifh basehighbp basehightc screen_cancer nSES ,                                                  
        noadj    =  agecon,     
        cat      =  agesub white active neverdrinker mifh cafh basehighbp basehightc   mv asp pmhever ocever mam endo2c_bin  screen_cancer, 
        poly     =  smkstatus bmic3cat ,
        rtftitle =  tableS2 Midpoint characteristics By ULCD, 
        landscape=  F,                               
        file     =  %str(table/pool.tableS2_ULCD_medianfu.rtf),
        dec      =  1,
        uselbl   =  T,
        nortf    = F);
run;


%table1(data=pooltableS2,
        agegroup =  agegp,
        exposure =  ulfssq,
        varlist  =  agecon agesub white smkstatus bmicon bmib3cat   actcon active alcocon neverdrinker 
                    cafh mv asp pmhever ocever mam endo2c_bin  
                    calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon alfssvcon vlfssvcon ulfssvcon hlfssvcon
                    basecarb basehcarb basesumref_carb baseprot baseaprot basevprot basetfat baseafat basevfat
                    carb hcarb sumref_carb prot aprot vprot tfat afat vfat
                    mifh basehighbp basehightc screen_cancer nSES ,                                                  
        noadj    =  agecon,     
        cat      =  agesub white active neverdrinker mifh cafh basehighbp basehightc   mv asp pmhever ocever mam endo2c_bin  screen_cancer, 
        poly     =  smkstatus bmic3cat ,
        rtftitle =  tableS2 Midpoint characteristics By ULFD, 
        landscape=  F,                               
        file     =  %str(table/pool.tableS2_ULFD_medianfu.rtf),
        dec      =  1,
        uselbl   =  T,
        nortf    = F);
run;


%table1(data=pooltableS2,
        agegroup =  agegp,
        exposure =  alcssq,
        varlist  =  agecon agesub white smkstatus bmicon bmib3cat   actcon active alcocon neverdrinker 
                    cafh mv asp pmhever ocever mam endo2c_bin  
                    calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon alfssvcon vlfssvcon ulfssvcon hlfssvcon
                    basecarb basehcarb basesumref_carb baseprot baseaprot basevprot basetfat baseafat basevfat
                    carb hcarb sumref_carb prot aprot vprot tfat afat vfat
                    mifh basehighbp basehightc screen_cancer nSES ,                                                  
        noadj    =  agecon,     
        cat      =  agesub white active neverdrinker mifh cafh basehighbp basehightc   mv asp pmhever ocever mam endo2c_bin  screen_cancer, 
        poly     =  smkstatus bmic3cat ,
        rtftitle =  tableS2 Midpoint characteristics By ALCD, 
        landscape=  F,                               
        file     =  %str(table/pool.tableS2_ALCD_medianfu.rtf),
        dec      =  1,
        uselbl   =  T,
        nortf    = F);
run;


%table1(data=pooltableS2,
        agegroup =  agegp,
        exposure =  alfssq,
        varlist  =  agecon agesub white smkstatus bmicon bmib3cat   actcon active alcocon neverdrinker 
                    cafh mv asp pmhever ocever mam endo2c_bin  
                    calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon alfssvcon vlfssvcon ulfssvcon hlfssvcon
                    basecarb basehcarb basesumref_carb baseprot baseaprot basevprot basetfat baseafat basevfat
                    carb hcarb sumref_carb prot aprot vprot tfat afat vfat
                    mifh basehighbp basehightc screen_cancer nSES ,                                                  
        noadj    =  agecon,     
        cat      =  agesub white active neverdrinker mifh cafh basehighbp basehightc   mv asp pmhever ocever mam endo2c_bin  screen_cancer, 
        poly     =  smkstatus bmic3cat ,
        rtftitle =  tableS2 Midpoint characteristics By ALFD, 
        landscape=  F,                               
        file     =  %str(table/pool.tableS2_ALFD_medianfu.rtf),
        dec      =  1,
        uselbl   =  T,
        nortf    = F);
run;


%table1(data=pooltableS2,
        agegroup =  agegp,
        exposure =  vlcssq,
        varlist  =  agecon agesub white smkstatus bmicon bmib3cat   actcon active alcocon neverdrinker 
                    cafh mv asp pmhever ocever mam endo2c_bin  
                    calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon alfssvcon vlfssvcon ulfssvcon hlfssvcon
                    basecarb basehcarb basesumref_carb baseprot baseaprot basevprot basetfat baseafat basevfat
                    carb hcarb sumref_carb prot aprot vprot tfat afat vfat
                    mifh basehighbp basehightc screen_cancer nSES ,                                                  
        noadj    =  agecon,     
        cat      =  agesub white active neverdrinker mifh cafh basehighbp basehightc   mv asp pmhever ocever mam endo2c_bin  screen_cancer, 
        poly     =  smkstatus bmic3cat ,
        rtftitle =  tableS2 Midpoint characteristics By VLCD, 
        landscape=  F,                               
        file     =  %str(table/pool.tableS2_VLCD_medianfu.rtf),
        dec      =  1,
        uselbl   =  T,
        nortf    = F);
run;


%table1(data=pooltableS2,
        agegroup =  agegp,
        exposure =  vlfssq,
        varlist  =  agecon agesub white smkstatus bmicon bmib3cat   actcon active alcocon neverdrinker 
                    cafh mv asp pmhever ocever mam endo2c_bin  
                    calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon alfssvcon vlfssvcon ulfssvcon hlfssvcon
                    basecarb basehcarb basesumref_carb baseprot baseaprot basevprot basetfat baseafat basevfat
                    carb hcarb sumref_carb prot aprot vprot tfat afat vfat
                    mifh basehighbp basehightc screen_cancer nSES ,                                                  
        noadj    =  agecon,     
        cat      =  agesub white active neverdrinker mifh cafh basehighbp basehightc   mv asp pmhever ocever mam endo2c_bin  screen_cancer, 
        poly     =  smkstatus bmic3cat ,
        rtftitle =  tableS2 Midpoint characteristics By VLFD, 
        landscape=  F,                               
        file     =  %str(table/pool.tableS2_VLFD_medianfu.rtf),
        dec      =  1,
        uselbl   =  T,
        nortf    = F);
run;
