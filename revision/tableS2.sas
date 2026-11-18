/*****************************************************
Purpose:       low carb, low fat, and risk of Obesity-related cancer
Study design:  prospective cohort (NHS: 1986-2022.6, NHSII: 1991-2023.6, HPFS: 1986-2024.1)

Contents:      Table S2: demographic characteristics by LCD, LFD at baseline
Revision:      add nSES, by all 10 lcd, lfd scores
Date:          10/2026
******************************************************/

/**** Call in Cleaned Datasets*******/

%include '/lowcarb_cancer/revision/data.inc.sas';


proc format;
value smkstatusf
1='never smker'
2='past smker'
3='curr smker';
run;

proc format;
value bmib3catf
1='bmi <24.9'
2='bmi 24.9-29.9'
3='bmi 30.0';
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

data pooltable1; set pool; if period=1; 

    keep    cohort agecon agegp agesub white smkstatus bmicon bmib3cat  actcon active alcocon neverdrinker 
            mifh cafh basehighbp basehightc mv asp pmhever ocever
            lcssq alcssq vlcssq ulcssq hlcssq lfssq ulfssq hlfssq alfssq vlfssq
            calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon alfssvcon vlfssvcon ulfssvcon hlfssvcon 
            basecarb basehcarb basesumref_carb baseprot baseaprot basevprot basetfat baseafat basevfat
            carb hcarb sumref_carb prot aprot vprot tfat afat vfat
            mam endo2c_bin  screen_cancer  ;
    format  smkstatus smkstatusf.;
    format  bmib3cat  bmib3catf. ;
    *format phmsstatus  phmsstatusf. ;
run;

%table1(data=pooltable1,
        agegroup =  agegp,
        exposure =  lcssq,
        varlist  =  agecon agesub white smkstatus bmicon bmib3cat   actcon active alcocon neverdrinker 
                    cafh mv asp pmhever ocever mam endo2c_bin  
                    calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon alfssvcon vlfssvcon ulfssvcon hlfssvcon
                    basecarb basehcarb basesumref_carb baseprot baseaprot basevprot basetfat baseafat basevfat
                    carb hcarb sumref_carb prot aprot vprot tfat afat vfat
                    mifh basehighbp basehightc screen_cancer ,                                                  
        noadj    =  agecon,     
        cat      =  agesub white active neverdrinker mifh cafh basehighbp basehightc mv asp pmhever ocever mam endo2c_bin  screen_cancer, 
        poly     =  smkstatus bmib3cat ,
        rtftitle =  TableS1 Baseline characteristics By LCD, 
        landscape=  F,                               
        file     =  %str(table/pool.tableS1_LCD.rtf),
        dec      =  1,
        uselbl   =  T,
        nortf    = F);
run;


%table1(data=pooltable1,
        agegroup =  agegp,
        exposure =  lfssq,
        varlist  =  agecon agesub white smkstatus bmicon bmib3cat   actcon active alcocon neverdrinker 
                    cafh mv asp pmhever ocever mam endo2c_bin  
                    calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon alfssvcon vlfssvcon ulfssvcon hlfssvcon
                    basecarb basehcarb basesumref_carb baseprot baseaprot basevprot basetfat baseafat basevfat
                    carb hcarb sumref_carb prot aprot vprot tfat afat vfat
                    mifh basehighbp basehightc screen_cancer ,                                                  
        noadj    =  agecon,     
        cat      =  agesub white active neverdrinker mifh cafh basehighbp basehightc mv asp pmhever ocever mam endo2c_bin  screen_cancer, 
        poly     =  smkstatus bmib3cat ,
        rtftitle =  TableS1 Baseline characteristics By LFD, 
        landscape=  F,                               
        file     =  %str(table/pool.tableS1_LFD.rtf),
        dec      =  1,
        uselbl   =  T,
        nortf    = F);
run;


%table1(data=pooltable1,
        agegroup =  agegp,
        exposure =  hlcssq,
        varlist  =  agecon agesub white smkstatus bmicon bmib3cat   actcon active alcocon neverdrinker 
                    cafh mv asp pmhever ocever mam endo2c_bin  
                    calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon alfssvcon vlfssvcon ulfssvcon hlfssvcon
                    basecarb basehcarb basesumref_carb baseprot baseaprot basevprot basetfat baseafat basevfat
                    carb hcarb sumref_carb prot aprot vprot tfat afat vfat
                    mifh basehighbp basehightc screen_cancer ,                                                  
        noadj    =  agecon,     
        cat      =  agesub white active neverdrinker mifh cafh basehighbp basehightc mv asp pmhever ocever mam endo2c_bin  screen_cancer, 
        poly     =  smkstatus bmib3cat ,
        rtftitle =  TableS1 Baseline characteristics By HLCD, 
        landscape=  F,                               
        file     =  %str(table/pool.tableS1_HLCD.rtf),
        dec      =  1,
        uselbl   =  T,
        nortf    = F);
run;


%table1(data=pooltable1,
        agegroup =  agegp,
        exposure =  hlfssq,
        varlist  =  agecon agesub white smkstatus bmicon bmib3cat   actcon active alcocon neverdrinker 
                    cafh mv asp pmhever ocever mam endo2c_bin  
                    calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon alfssvcon vlfssvcon ulfssvcon hlfssvcon
                    basecarb basehcarb basesumref_carb baseprot baseaprot basevprot basetfat baseafat basevfat
                    carb hcarb sumref_carb prot aprot vprot tfat afat vfat
                    mifh basehighbp basehightc screen_cancer ,                                                  
        noadj    =  agecon,     
        cat      =  agesub white active neverdrinker mifh cafh basehighbp basehightc mv asp pmhever ocever mam endo2c_bin  screen_cancer, 
        poly     =  smkstatus bmib3cat ,
        rtftitle =  TableS1 Baseline characteristics By HLFD, 
        landscape=  F,                               
        file     =  %str(table/pool.tableS1_HLFD.rtf),
        dec      =  1,
        uselbl   =  T,
        nortf    = F);
run;



%table1(data=pooltable1,
        agegroup =  agegp,
        exposure =  hlcssq,
        varlist  =  agecon agesub white smkstatus bmicon bmib3cat   actcon active alcocon neverdrinker 
                    cafh mv asp pmhever ocever mam endo2c_bin  
                    calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon alfssvcon vlfssvcon ulfssvcon hlfssvcon
                    basecarb basehcarb basesumref_carb baseprot baseaprot basevprot basetfat baseafat basevfat
                    carb hcarb sumref_carb prot aprot vprot tfat afat vfat
                    mifh basehighbp basehightc screen_cancer ,                                                  
        noadj    =  agecon,     
        cat      =  agesub white active neverdrinker mifh cafh basehighbp basehightc mv asp pmhever ocever mam endo2c_bin  screen_cancer, 
        poly     =  smkstatus bmib3cat ,
        rtftitle =  TableS1 Baseline characteristics By HLCD, 
        landscape=  F,                               
        file     =  %str(table/pool.tableS1_HLCD_overall.rtf),
        dec      =  1,
        uselbl   =  T,
        noexp=T,
        nortf    = F);
run;


%table1(data=pooltable1,
        agegroup =  agegp,
        exposure =  hlfssq,
        varlist  =  agecon agesub white smkstatus bmicon bmib3cat   actcon active alcocon neverdrinker 
                    cafh mv asp pmhever ocever mam endo2c_bin  
                    calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon alfssvcon vlfssvcon ulfssvcon hlfssvcon
                    basecarb basehcarb basesumref_carb baseprot baseaprot basevprot basetfat baseafat basevfat
                    carb hcarb sumref_carb prot aprot vprot tfat afat vfat
                    mifh basehighbp basehightc screen_cancer ,                                                  
        noadj    =  agecon,     
        cat      =  agesub white active neverdrinker mifh cafh basehighbp basehightc mv asp pmhever ocever mam endo2c_bin  screen_cancer, 
        poly     =  smkstatus bmib3cat ,
        rtftitle =  TableS1 Baseline characteristics By HLFD, 
        landscape=  F,                               
        file     =  %str(table/pool.tableS1_HLFD_overall.rtf),
        dec      =  1,
        uselbl   =  T,
        noexp=T,
        nortf    = F);
run;








%table1(data=pooltable1,
        agegroup =  agegp,
        exposure =  ulcssq,
        varlist  =  agecon agesub white smkstatus bmicon bmib3cat   actcon active alcocon neverdrinker 
                    cafh mv asp pmhever ocever mam endo2c_bin  
                    calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon alfssvcon vlfssvcon ulfssvcon hlfssvcon
                    basecarb basehcarb basesumref_carb baseprot baseaprot basevprot basetfat baseafat basevfat
                    carb hcarb sumref_carb prot aprot vprot tfat afat vfat
                    mifh basehighbp basehightc screen_cancer ,                                                  
        noadj    =  agecon,     
        cat      =  agesub white active neverdrinker mifh cafh basehighbp basehightc mv asp pmhever ocever mam endo2c_bin  screen_cancer, 
        poly     =  smkstatus bmib3cat ,
        rtftitle =  TableS1 Baseline characteristics By ULCD, 
        landscape=  F,                               
        file     =  %str(table/pool.tableS1_ULCD.rtf),
        dec      =  1,
        uselbl   =  T,
        nortf    = F);
run;


%table1(data=pooltable1,
        agegroup =  agegp,
        exposure =  ulfssq,
        varlist  =  agecon agesub white smkstatus bmicon bmib3cat   actcon active alcocon neverdrinker 
                    cafh mv asp pmhever ocever mam endo2c_bin  
                    calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon alfssvcon vlfssvcon ulfssvcon hlfssvcon
                    basecarb basehcarb basesumref_carb baseprot baseaprot basevprot basetfat baseafat basevfat
                    carb hcarb sumref_carb prot aprot vprot tfat afat vfat
                    mifh basehighbp basehightc screen_cancer ,                                                  
        noadj    =  agecon,     
        cat      =  agesub white active neverdrinker mifh cafh basehighbp basehightc mv asp pmhever ocever mam endo2c_bin  screen_cancer, 
        poly     =  smkstatus bmib3cat ,
        rtftitle =  TableS1 Baseline characteristics By ULFD, 
        landscape=  F,                               
        file     =  %str(table/pool.tableS1_ULFD.rtf),
        dec      =  1,
        uselbl   =  T,
        nortf    = F);
run;



%table1(data=pooltable1,
        agegroup =  agegp,
        exposure =  ulcssq,
        varlist  =  agecon agesub white smkstatus bmicon bmib3cat   actcon active alcocon neverdrinker 
                    cafh mv asp pmhever ocever mam endo2c_bin  
                    calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon alfssvcon vlfssvcon ulfssvcon hlfssvcon
                    basecarb basehcarb basesumref_carb baseprot baseaprot basevprot basetfat baseafat basevfat
                    carb hcarb sumref_carb prot aprot vprot tfat afat vfat
                    mifh basehighbp basehightc screen_cancer ,                                                  
        noadj    =  agecon,     
        cat      =  agesub white active neverdrinker mifh cafh basehighbp basehightc mv asp pmhever ocever mam endo2c_bin  screen_cancer, 
        poly     =  smkstatus bmib3cat ,
        rtftitle =  TableS1 Baseline characteristics By ULCD, 
        landscape=  F,                               
        file     =  %str(table/pool.tableS1_ULCD_overall.rtf),
        dec      =  1,
        uselbl   =  T,
        noexp=T,
        nortf    = F);
run;


%table1(data=pooltable1,
        agegroup =  agegp,
        exposure =  ulfssq,
        varlist  =  agecon agesub white smkstatus bmicon bmib3cat   actcon active alcocon neverdrinker 
                    cafh mv asp pmhever ocever mam endo2c_bin  
                    calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon alfssvcon vlfssvcon ulfssvcon hlfssvcon
                    basecarb basehcarb basesumref_carb baseprot baseaprot basevprot basetfat baseafat basevfat
                    carb hcarb sumref_carb prot aprot vprot tfat afat vfat
                    mifh basehighbp basehightc screen_cancer ,                                                  
        noadj    =  agecon,     
        cat      =  agesub white active neverdrinker mifh cafh basehighbp basehightc mv asp pmhever ocever mam endo2c_bin  screen_cancer, 
        poly     =  smkstatus bmib3cat ,
        rtftitle =  TableS1 Baseline characteristics By ULFD, 
        landscape=  F,                               
        file     =  %str(table/pool.tableS1_ULFD_overall.rtf),
        dec      =  1,
        uselbl   =  T,
        noexp=T,
        nortf    = F);
run;


%table1(data=pooltable1,
        agegroup =  agegp,
        exposure =  alcssq,
        varlist  =  agecon agesub white smkstatus bmicon bmib3cat   actcon active alcocon neverdrinker 
                    cafh mv asp pmhever ocever mam endo2c_bin  
                    calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon alfssvcon vlfssvcon ulfssvcon hlfssvcon
                    basecarb basehcarb basesumref_carb baseprot baseaprot basevprot basetfat baseafat basevfat
                    carb hcarb sumref_carb prot aprot vprot tfat afat vfat
                    mifh basehighbp basehightc screen_cancer ,                                                  
        noadj    =  agecon,     
        cat      =  agesub white active neverdrinker mifh cafh basehighbp basehightc mv asp pmhever ocever mam endo2c_bin  screen_cancer, 
        poly     =  smkstatus bmib3cat ,
        rtftitle =  TableS1 Baseline characteristics By ALCD, 
        landscape=  F,                               
        file     =  %str(table/pool.tableS1_ALCD.rtf),
        dec      =  1,
        uselbl   =  T,
        nortf    = F);
run;


%table1(data=pooltable1,
        agegroup =  agegp,
        exposure =  alfssq,
        varlist  =  agecon agesub white smkstatus bmicon bmib3cat   actcon active alcocon neverdrinker 
                    cafh mv asp pmhever ocever mam endo2c_bin  
                    calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon alfssvcon vlfssvcon ulfssvcon hlfssvcon
                    basecarb basehcarb basesumref_carb baseprot baseaprot basevprot basetfat baseafat basevfat
                    carb hcarb sumref_carb prot aprot vprot tfat afat vfat
                    mifh basehighbp basehightc screen_cancer ,                                                  
        noadj    =  agecon,     
        cat      =  agesub white active neverdrinker mifh cafh basehighbp basehightc mv asp pmhever ocever mam endo2c_bin  screen_cancer, 
        poly     =  smkstatus bmib3cat ,
        rtftitle =  TableS1 Baseline characteristics By ALFD, 
        landscape=  F,                               
        file     =  %str(table/pool.tableS1_ALFD.rtf),
        dec      =  1,
        uselbl   =  T,
        nortf    = F);
run;



%table1(data=pooltable1,
        agegroup =  agegp,
        exposure =  vlcssq,
        varlist  =  agecon agesub white smkstatus bmicon bmib3cat   actcon active alcocon neverdrinker 
                    cafh mv asp pmhever ocever mam endo2c_bin  
                    calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon alfssvcon vlfssvcon ulfssvcon hlfssvcon
                    basecarb basehcarb basesumref_carb baseprot baseaprot basevprot basetfat baseafat basevfat
                    carb hcarb sumref_carb prot aprot vprot tfat afat vfat
                    mifh basehighbp basehightc screen_cancer ,                                                  
        noadj    =  agecon,     
        cat      =  agesub white active neverdrinker mifh cafh basehighbp basehightc mv asp pmhever ocever mam endo2c_bin  screen_cancer, 
        poly     =  smkstatus bmib3cat ,
        rtftitle =  TableS1 Baseline characteristics By VLCD, 
        landscape=  F,                               
        file     =  %str(table/pool.tableS1_VLCD.rtf),
        dec      =  1,
        uselbl   =  T,
        nortf    = F);
run;


%table1(data=pooltable1,
        agegroup =  agegp,
        exposure =  vlfssq,
        varlist  =  agecon agesub white smkstatus bmicon bmib3cat   actcon active alcocon neverdrinker 
                    cafh mv asp pmhever ocever mam endo2c_bin  
                    calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon alfssvcon vlfssvcon ulfssvcon hlfssvcon
                    basecarb basehcarb basesumref_carb baseprot baseaprot basevprot basetfat baseafat basevfat
                    carb hcarb sumref_carb prot aprot vprot tfat afat vfat
                    mifh basehighbp basehightc screen_cancer ,                                                  
        noadj    =  agecon,     
        cat      =  agesub white active neverdrinker mifh cafh basehighbp basehightc mv asp pmhever ocever mam endo2c_bin  screen_cancer, 
        poly     =  smkstatus bmib3cat ,
        rtftitle =  TableS1 Baseline characteristics By VLFD, 
        landscape=  F,                               
        file     =  %str(table/pool.tableS1_VLFD.rtf),
        dec      =  1,
        uselbl   =  T,
        nortf    = F);
run;


