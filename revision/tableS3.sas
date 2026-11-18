/*****************************************************
Purpose:       low carb, low fat, and risk of Obesity-related cancer
Study design:  prospective cohort (NHS: 1986-2022.6, NHSII: 1991-2023.6, HPFS: 1986-2024.1)

Contents:      Revision Table S3: energy-adjusted nutrient intakes across score categories
Date:          10/2026
******************************************************/

/**** Call in Cleaned Datasets*******/

%include 'lowcarb_cancer/revision/data.inc.sas';


proc sort data=pool; by newid period; run;


data pooltable1; set pool; if period=1; 

    keep    cohort agecon  agegp nSES 
            lcssq alcssq vlcssq ulcssq hlcssq lfssq ulfssq hlfssq alfssq vlfssq
            calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon alfssvcon vlfssvcon ulfssvcon hlfssvcon 
            /* energy-adjusted nutrient intakes (e.g., carbo86p, percent of total energy, already energy-adjusted)  */
            basecarb basehcarb basesumref_carb baseprot baseaprot basevprot basetfat baseafat basevfat
            carb_svcon hcarb_svcon sumref_carb_svcon prot_svcon aprot_svcon vprot_svcon tfat_svcon afat_svcon vfat_svcon

            /* food items */
            basewhgrn basergrain basefrui  baseleg basepota baseveg_noleg  basefish baserpmeats  basessb basedairy basesweets basepoultry 
            whgrn_svcon rgrain_svcon frui_svcon leg_svcon pota_svcon veg_noleg_svcon fish_svcon rpmeats_svcon ssb_svcon 
            dairy_svcon swtsnochoc_svcon poul_svcon 

            /* Individual obesity-related cancer: 13 Types, NEJM, 2016*/
            can_esoph can_stomach  can_colorect can_liver can_gallbladder can_panc  can_brca can_endo can_ovca can_kidney 
            can_thyroid can_myeloma can_menin

            /* other cancer site */
            can_oral can_nose can_larynx can_lung can_smallintest 
            can_bone can_melanoma  can_cervical  can_bladder 
            can_brain   can_myeleuk    
            ;
run;


%table1(data=pooltable1,
        agegroup =  agegp,
        exposure =  lcssq,
        varlist  =  agecon  nSES 
                    calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon alfssvcon vlfssvcon ulfssvcon hlfssvcon
                    basecarb basehcarb basesumref_carb baseprot baseaprot basevprot basetfat baseafat basevfat
                    carb_svcon hcarb_svcon sumref_carb_svcon prot_svcon aprot_svcon vprot_svcon tfat_svcon afat_svcon vfat_svcon
                     
                    /* food items */
                    basewhgrn basergrain basefrui  baseleg basepota baseveg_noleg  basefish baserpmeats  basessb basedairy basesweets basepoultry 
                    whgrn_svcon rgrain_svcon frui_svcon leg_svcon pota_svcon veg_noleg_svcon fish_svcon rpmeats_svcon ssb_svcon 
                    dairy_svcon swtsnochoc_svcon poul_svcon 

                    
                    ,                                                  
        noadj    =  agecon, 
        rtftitle =  TableS1 Baseline characteristics, 
        landscape=  F,                               
        file     =  %str(table/revision_pool.tableS1_total.rtf),
        dec      =  1,
        uselbl   =  T,
        noexp    =  T,
        nortf    = F);
run;



%table1(data=pooltable1,
        agegroup =  agegp,
        exposure =  lcssq,
        varlist  =  agecon  nSES 
                    calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon alfssvcon vlfssvcon ulfssvcon hlfssvcon
                    basecarb basehcarb basesumref_carb baseprot baseaprot basevprot basetfat baseafat basevfat
                    carb_svcon hcarb_svcon sumref_carb_svcon prot_svcon aprot_svcon vprot_svcon tfat_svcon afat_svcon vfat_svcon
                     
                    /* food items */
                    basewhgrn basergrain basefrui  baseleg basepota baseveg_noleg  basefish baserpmeats  basessb basedairy basesweets basepoultry 
                    whgrn_svcon rgrain_svcon frui_svcon leg_svcon pota_svcon veg_noleg_svcon fish_svcon rpmeats_svcon ssb_svcon 
                    dairy_svcon swtsnochoc_svcon poul_svcon 

                    
                    ,                                                  
        noadj    =  agecon, 
        rtftitle =  TableS1 Baseline characteristics By LCD, 
        landscape=  F,                               
        file     =  %str(table/revision_pool.tableS1_LCD.rtf),
        dec      =  1,
        uselbl   =  T,
        nortf    = F);
run;



%table1(data=pooltable1,
        agegroup =  agegp,
        exposure =  alcssq,
        varlist  =  agecon  nSES 
                    calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon alfssvcon vlfssvcon ulfssvcon hlfssvcon
                    basecarb basehcarb basesumref_carb baseprot baseaprot basevprot basetfat baseafat basevfat
                    carb_svcon hcarb_svcon sumref_carb_svcon prot_svcon aprot_svcon vprot_svcon tfat_svcon afat_svcon vfat_svcon
                     
                    /* food items */
                    basewhgrn basergrain basefrui  baseleg basepota baseveg_noleg  basefish baserpmeats  basessb basedairy basesweets basepoultry 
                    whgrn_svcon rgrain_svcon frui_svcon leg_svcon pota_svcon veg_noleg_svcon fish_svcon rpmeats_svcon ssb_svcon 
                    dairy_svcon swtsnochoc_svcon poul_svcon  
                    ,                                                  
        noadj    =  agecon,     
        rtftitle =  TableS1 Baseline characteristics By ALCD, 
        landscape=  F,                               
        file     =  %str(table/revision_pool.tableS1_ALCD.rtf),
        dec      =  1,
        uselbl   =  T,
        
        nortf    = F);
run;



%table1(data=pooltable1,
        agegroup =  agegp,
        exposure =  ulcssq,
        varlist  =  agecon  nSES 
                    calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon alfssvcon vlfssvcon ulfssvcon hlfssvcon
                    basecarb basehcarb basesumref_carb baseprot baseaprot basevprot basetfat baseafat basevfat
                    carb_svcon hcarb_svcon sumref_carb_svcon prot_svcon aprot_svcon vprot_svcon tfat_svcon afat_svcon vfat_svcon
                     
                    /* food items */
                    basewhgrn basergrain basefrui  baseleg basepota baseveg_noleg  basefish baserpmeats  basessb basedairy basesweets basepoultry 
                    whgrn_svcon rgrain_svcon frui_svcon leg_svcon pota_svcon veg_noleg_svcon fish_svcon rpmeats_svcon ssb_svcon 
                    dairy_svcon swtsnochoc_svcon poul_svcon 

                     
                    ,                                                  
        noadj    =  agecon,     
        
        rtftitle =  TableS1 Baseline characteristics By ULCD, 
        landscape=  F,                               
        file     =  %str(table/revision_pool.tableS1_ULCD.rtf),
        dec      =  1,
        uselbl   =  T,
        
        nortf    = F);
run;



%table1(data=pooltable1,
        agegroup =  agegp,
        exposure =  vlcssq,
        varlist  =  agecon  nSES 
                    calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon alfssvcon vlfssvcon ulfssvcon hlfssvcon
                    basecarb basehcarb basesumref_carb baseprot baseaprot basevprot basetfat baseafat basevfat
                    carb_svcon hcarb_svcon sumref_carb_svcon prot_svcon aprot_svcon vprot_svcon tfat_svcon afat_svcon vfat_svcon
                     
                    /* food items */
                    basewhgrn basergrain basefrui  baseleg basepota baseveg_noleg  basefish baserpmeats  basessb basedairy basesweets basepoultry 
                    whgrn_svcon rgrain_svcon frui_svcon leg_svcon pota_svcon veg_noleg_svcon fish_svcon rpmeats_svcon ssb_svcon 
                    dairy_svcon swtsnochoc_svcon poul_svcon 

                     
                    ,                                                  
        noadj    =  agecon,     
        
        rtftitle =  TableS1 Baseline characteristics By VLCD, 
        landscape=  F,                               
        file     =  %str(table/revision_pool.tableS1_VLCD.rtf),
        dec      =  1,
        uselbl   =  T,
        
        nortf    = F);
run;


%table1(data=pooltable1,
        agegroup =  agegp,
        exposure =  hlcssq,
        varlist  =  agecon  nSES 
                    calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon alfssvcon vlfssvcon ulfssvcon hlfssvcon
                    basecarb basehcarb basesumref_carb baseprot baseaprot basevprot basetfat baseafat basevfat
                    carb_svcon hcarb_svcon sumref_carb_svcon prot_svcon aprot_svcon vprot_svcon tfat_svcon afat_svcon vfat_svcon
                     
                    /* food items */
                    basewhgrn basergrain basefrui  baseleg basepota baseveg_noleg  basefish baserpmeats  basessb basedairy basesweets basepoultry 
                    whgrn_svcon rgrain_svcon frui_svcon leg_svcon pota_svcon veg_noleg_svcon fish_svcon rpmeats_svcon ssb_svcon 
                    dairy_svcon swtsnochoc_svcon poul_svcon 

                     
                    ,                                                  
        noadj    =  agecon,     
        
        rtftitle =  TableS1 Baseline characteristics By HLCD, 
        landscape=  F,                               
        file     =  %str(table/revision_pool.tableS1_HLCD.rtf),
        dec      =  1,
        uselbl   =  T,
        
        nortf    = F);
run;



/* By LFD */

%table1(data=pooltable1,
        agegroup =  agegp,
        exposure =  lfssq,
        varlist  =  agecon  nSES 
                    calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon alfssvcon vlfssvcon ulfssvcon hlfssvcon
                    basecarb basehcarb basesumref_carb baseprot baseaprot basevprot basetfat baseafat basevfat
                    carb_svcon hcarb_svcon sumref_carb_svcon prot_svcon aprot_svcon vprot_svcon tfat_svcon afat_svcon vfat_svcon
                     
                    /* food items */
                    basewhgrn basergrain basefrui  baseleg basepota baseveg_noleg  basefish baserpmeats  basessb basedairy basesweets basepoultry 
                    whgrn_svcon rgrain_svcon frui_svcon leg_svcon pota_svcon veg_noleg_svcon fish_svcon rpmeats_svcon ssb_svcon 
                    dairy_svcon swtsnochoc_svcon poul_svcon 

                     
                    ,                                                  
        noadj    =  agecon,     
        
        rtftitle =  TableS1 Baseline characteristics By LFD, 
        landscape=  F,                               
        file     =  %str(table/revision_pool.tableS1_LFD.rtf),
        dec      =  1,
        uselbl   =  T,
        
        nortf    = F);
run;



%table1(data=pooltable1,
        agegroup =  agegp,
        exposure =  alfssq,
        varlist  =  agecon  nSES 
                    calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon alfssvcon vlfssvcon ulfssvcon hlfssvcon
                    basecarb basehcarb basesumref_carb baseprot baseaprot basevprot basetfat baseafat basevfat
                    carb_svcon hcarb_svcon sumref_carb_svcon prot_svcon aprot_svcon vprot_svcon tfat_svcon afat_svcon vfat_svcon
                     
                    /* food items */
                    basewhgrn basergrain basefrui  baseleg basepota baseveg_noleg  basefish baserpmeats  basessb basedairy basesweets basepoultry 
                    whgrn_svcon rgrain_svcon frui_svcon leg_svcon pota_svcon veg_noleg_svcon fish_svcon rpmeats_svcon ssb_svcon 
                    dairy_svcon swtsnochoc_svcon poul_svcon 

                     
                    ,                                                  
        noadj    =  agecon,     
        
        rtftitle =  TableS1 Baseline characteristics By ALFD, 
        landscape=  F,                               
        file     =  %str(table/revision_pool.tableS1_ALFD.rtf),
        dec      =  1,
        uselbl   =  T,
        
        nortf    = F);
run;



%table1(data=pooltable1,
        agegroup =  agegp,
        exposure =  ulfssq,
        varlist  =  agecon  nSES 
                    calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon alfssvcon vlfssvcon ulfssvcon hlfssvcon
                    basecarb basehcarb basesumref_carb baseprot baseaprot basevprot basetfat baseafat basevfat
                    carb_svcon hcarb_svcon sumref_carb_svcon prot_svcon aprot_svcon vprot_svcon tfat_svcon afat_svcon vfat_svcon
                     
                    /* food items */
                    basewhgrn basergrain basefrui  baseleg basepota baseveg_noleg  basefish baserpmeats  basessb basedairy basesweets basepoultry 
                    whgrn_svcon rgrain_svcon frui_svcon leg_svcon pota_svcon veg_noleg_svcon fish_svcon rpmeats_svcon ssb_svcon 
                    dairy_svcon swtsnochoc_svcon poul_svcon 

                     
                    ,                                                  
        noadj    =  agecon,     
        
        rtftitle =  TableS1 Baseline characteristics By ULFD, 
        landscape=  F,                               
        file     =  %str(table/revision_pool.tableS1_ULFD.rtf),
        dec      =  1,
        uselbl   =  T,
        
        nortf    = F);
run;



%table1(data=pooltable1,
        agegroup =  agegp,
        exposure =  vlfssq,
        varlist  =  agecon  nSES 
                    calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon alfssvcon vlfssvcon ulfssvcon hlfssvcon
                    basecarb basehcarb basesumref_carb baseprot baseaprot basevprot basetfat baseafat basevfat
                    carb_svcon hcarb_svcon sumref_carb_svcon prot_svcon aprot_svcon vprot_svcon tfat_svcon afat_svcon vfat_svcon
                     
                    /* food items */
                    basewhgrn basergrain basefrui  baseleg basepota baseveg_noleg  basefish baserpmeats  basessb basedairy basesweets basepoultry 
                    whgrn_svcon rgrain_svcon frui_svcon leg_svcon pota_svcon veg_noleg_svcon fish_svcon rpmeats_svcon ssb_svcon 
                    dairy_svcon swtsnochoc_svcon poul_svcon 

                     
                    ,                                                  
        noadj    =  agecon,     
        
        rtftitle =  TableS1 Baseline characteristics By VLFD, 
        landscape=  F,                               
        file     =  %str(table/revision_pool.tableS1_VLFD.rtf),
        dec      =  1,
        uselbl   =  T,
        
        nortf    = F);
run;


%table1(data=pooltable1,
        agegroup =  agegp,
        exposure =  hlfssq,
        varlist  =  agecon  nSES 
                    calorsvcon aheimsvcon lcssvcon alcssvcon vlcssvcon ulcssvcon hlcssvcon lfssvcon alfssvcon vlfssvcon ulfssvcon hlfssvcon
                    basecarb basehcarb basesumref_carb baseprot baseaprot basevprot basetfat baseafat basevfat
                    carb_svcon hcarb_svcon sumref_carb_svcon prot_svcon aprot_svcon vprot_svcon tfat_svcon afat_svcon vfat_svcon
                     
                    /* food items */
                    basewhgrn basergrain basefrui  baseleg basepota baseveg_noleg  basefish baserpmeats  basessb basedairy basesweets basepoultry 
                    whgrn_svcon rgrain_svcon frui_svcon leg_svcon pota_svcon veg_noleg_svcon fish_svcon rpmeats_svcon ssb_svcon 
                    dairy_svcon swtsnochoc_svcon poul_svcon 

                     
                    ,                                                  
        noadj    =  agecon,     
        
        rtftitle =  TableS1 Baseline characteristics By HLFD, 
        landscape=  F,                               
        file     =  %str(table/revision_pool.tableS1_HLFD.rtf),
        dec      =  1,
        uselbl   =  T,
        
        nortf    = F);
run;
