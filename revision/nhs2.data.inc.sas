/*****************************************************
Purpose:       low carb, low fat, and risk of obesity-related cancer and mortality
Study design:  prospective cohort (NHS: 1986-2022.6, NHSII: 1991-2023.6, HPFS: 1986-2024.1)
Endpoints:     obesity-related cancer and mortality
Exposures:     lcs, alcs, vlcs, ulcs, hlcs, lfs, ulfs, hlfs (NHS: 1986-2010, NHSII: 1991-2019, HPFS: 1986-2018)
Exclusion:     baseline cancer or cvd or diabetes, exposure missing, 
               or reported implausible SFFQ data (<800 or >4200 kcal/d for men and <600 or >3500 kcal/d for women)
Covariates:    age, interval, cohort, white, total energy quintile, smoking, alcohol, physical activity, cancer family history, 
               hormone use, oral contraceptive use, multivitamin use, aspirin use, 
               baseline BMI, endoscopic screening, mammography 
Mediator:      time-varying BMI, 4-year weight change in kilogram
Stratifiers:   Age (<55, ≥55y), Baseline BMI (<30.0, or ≥30.0 kg/m2)
Analysis:      Survival analysis: Cox regression
               Sensitivity analysis: add SES as confounder (census tract-level neighborhood socioeconomic status)

Date:  10/2026
  
******************************************************/
   
                                                                     
*----------------1.Outcome-----------------------------;
/**** Death *****/
%deadff(file='deadff.current.nhs2');

data dead; 
    set deadff;
        dtdth=9999;
        if deadmonth > 0 then dtdth=deadmonth;
newicda=compress(nhsicda, 'E');
newicda3=substr(newicda,1,3);
if dtdth=9999 then delete;  
keep id dtdth newicda newicda3; run;

/**** Cancer *****/
/*All cancer: 1989-2023*/

data cancer;
    infile 'canc.cases' lrecl=53 recfm=d; 
    
    input
    @1      id          6.
    @7      cd      $   1.
    @10   cohort        1.
    @13     icda         $   5. 
    @20     qyr         2.
    @26     mdx1            2.
    @28     ydx1            2.
    @32     conf            2.
    @36     nhsicda         5.1
    @43     mdx         2.
    @45     ydx         2.
    @49     dtdxca          4.
    ;

    label
        id      = 'ID'
        cd      = 'Check Digit'
        cohort    = 'Cohort'
        icda        = 'Reported ICDA'
        qyr     = 'Questionnaire cycle year when reported'
        mdx1        = 'Reported month of diagnosis'
        ydx1        = 'Reported year of diagnosis'
        conf        = 'Confirmation Code'
        nhsicda     = 'NHS reviewed and confirmed ICDA code'
        mdx     = 'Confirmed Month of Diagnosis'
        ydx     = 'Confirmed Year of Diagnosis'
        dtdxca      = 'Date of dx months since 1900'
    ;
        icda2=nhsicda;

        if cohort=2;  /*include NHS2*/
        if conf in (11,12,13,14,15,16,17,18,19) then ca=1; /* ca=1: confirmed cancer ; exclude in situ ca by taking out conf code 25*/
        if ca=1 and icda2 ne . then caflag=1;

        dtdx_allca=.;
        if dtdxca>0 and ca=1 and (140<=icda2<=172 or 174<=icda2<=207 or icda2=225) then dtdx_allca=dtdxca;

    /* ---- Site groups ---- */
      /* Digestive */
      if 140 <= icda2 <= 149 and ca=1 then do; oral_ca=1;      dtdx_oral      = dtdxca; end;
      if icda2=150 and ca=1 then do; esoph_ca=1;               dtdx_esoph     = dtdxca; end;
      if icda2=151 and ca=1 then do; stomach_ca=1;             dtdx_stomach   = dtdxca; end;
      if icda2=152 and ca=1 then do; smallintest_ca=1;         dtdx_smallintest=dtdxca; end;
      if icda2 in (153,154) and ca=1 then do; colorect_ca=1;   dtdx_colorect  = dtdxca; end;
      if icda2=155 and ca=1 then do; liver_ca=1;               dtdx_liver     = dtdxca; end;
      if icda2=156 and ca=1 then do; gallbladder_ca=1;         dtdx_gallbladder=dtdxca; end;
      if icda2=157 and ca=1 then do; panc_ca=1;                dtdx_panc      = dtdxca; end;
      if 140 <= icda2 <= 159 and ca=1 then do; diges_ca=1;     dtdx_diges     = dtdxca; end;

      /* Respiratory */
      if icda2=160 and ca=1 then do; nose_ca=1;                dtdx_nose      = dtdxca; end;
      if icda2=161 and ca=1 then do; larynx_ca=1;              dtdx_larynx    = dtdxca; end;
      if icda2=162 and ca=1 then do; lung_ca=1;                dtdx_lung      = dtdxca; end;
      if 160 <= icda2 <= 165 and ca=1 then do; resp_ca=1;      dtdx_resp      = dtdxca; end;

      /* Skin / other selected */
      if icda2=170 and ca=1 then do; bone_ca=1;                dtdx_bone      = dtdxca; end;
      if icda2=172 and ca=1 then do; melanoma_ca=1;            dtdx_melanoma  = dtdxca; end;
      if icda2=174 and ca=1 then do; brca_ca=1;                dtdx_brca      = dtdxca; end;

      /* Genito-urinary Cancers */
      if icda2=180 and ca=1 then do; cervical=1;               dtdx_cervical  = dtdxca; end;
      if icda2=182 and ca=1 then do; endo_ca=1;                dtdx_endo      = dtdxca; end;
      if icda2=183 and ca=1 then do; ovca_ca=1;                dtdx_ovca      = dtdxca; end;
      if icda2=188 and ca=1 then do; bladder_ca=1;             dtdx_bladder   = dtdxca; end;
      if icda2=189 and ca=1 then do; kidney_ca=1;              dtdx_kidney    = dtdxca; end;
      if icda2 in (185,186,187,188,189) and ca=1 then do; uri_ca=1; dtdx_uri  = dtdxca; end;

      /* Other sites */
      if icda2=191 and ca=1 then do; brain_ca=1;               dtdx_brain     = dtdxca; end;
      if icda2=193 and ca=1 then do; thyroid_ca=1;             dtdx_thyroid   = dtdxca; end;

      /* Lymphatic / blood */
      if icda2=203 and ca=1 then do; myeloma_ca=1;             dtdx_myeloma   = dtdxca; end;
      if icda2=205 and ca=1 then do; myeleuk_ca=1;             dtdx_myeleuk   = dtdxca; end;

      /* Meningioma */
      if icda2=225 and ca=1 then do; menin_ca=1;               dtdx_menin     = dtdxca; end;

      /* Smoking-related cancer */
      smoke_ca = 0; dtdx_smoke_ca = .;
      if sum(of oral_ca resp_ca bladder_ca kidney_ca myeleuk_ca esoph_ca stomach_ca liver_ca colorect_ca panc_ca) >= 1 then do;
        smoke_ca = 1; dtdx_smoke_ca = dtdxca; end;
      r_smoke_ca = 0; dtdx_r_smoke_ca = .;
      if sum(oral_ca, lung_ca, esoph_ca, bladder_ca) >= 1 then do; r_smoke_ca = 1; dtdx_r_smoke_ca = dtdxca; end;


      /*obesity-related cancer: 13 Types, NEJM, 2016*/
      obese_ca=0;
      if 
         esoph_ca=1       or
         colorect_ca=1    or
         stomach_ca=1     or /*Gastric cardia cancer*/
         liver_ca=1       or
         gallbladder_ca=1 or
         panc_ca=1        or
         brca_ca=1        or     
         endo_ca=1        or
         ovca_ca=1        or
         kidney_ca=1      or
         thyroid_ca=1     or
         myeloma_ca=1     or
         menin_ca =1     
      then obese_ca=1;

     if obese_ca ne 1 and caflag=1 then nonobeca=1; /*all other cancer except obesity-related cancer*/
     run;
     proc sort data=cancer; by id dtdxca; run;


data cancer;
     set cancer end=_end_;

/*==============================================================*/
/* Obesity-related cancer; exclude brc+crc; exclude brc+crc+ovca+panc */
/*==============================================================*/

     dtdx_obeca=.;
     if obese_ca=1 and dtdxca>0 and ca=1 and icda2 in (150,151,153,154,155,156,157,174,182,183,189,193,203,225) then dtdx_obeca=dtdxca;
     dtdx_nonobeca=.;
     if nonobeca=1 and dtdxca>0 and ca=1  then dtdx_nonobeca=dtdxca;
    

     obeca_nobc=0;
     if 
         esoph_ca=1       or
         /*colorect_ca=1    or */
         stomach_ca=1     or /*Gastric cardia cancer*/
         liver_ca=1       or
         gallbladder_ca=1 or
         panc_ca=1        or
         /*brca_ca=1        or */
         endo_ca=1        or
         ovca_ca=1        or
         kidney_ca=1      or
         thyroid_ca=1     or
         myeloma_ca=1     or
         menin_ca =1     
      then obeca_nobc=1;

     dtdx_obeca_nobc=.;
     if obeca_nobc=1 and dtdxca>0 and ca=1 and icda2 in (150,151, /*153,154,*/ 155,156,157, /*174,*/ 182,183,189,193,203,225) then dtdx_obeca_nobc=dtdxca;
     
     obeca_nobcop=0;
     if 
         esoph_ca=1       or
         /*colorect_ca=1    or */
         stomach_ca=1     or /*Gastric cardia cancer*/
         liver_ca=1       or
         gallbladder_ca=1 or
         /*panc_ca=1        or */
         /*brca_ca=1        or */
         /*endo_ca=1        or */
         /*ovca_ca=1        or */
         kidney_ca=1      or
         thyroid_ca=1     or
         myeloma_ca=1     or
         menin_ca =1     
     then obeca_nobcop=1;
     dtdx_obeca_nobcop=.;
     if obeca_nobcop=1 and dtdxca>0 and ca=1 and icda2 in (150,151, /*153,154,*/ 155,156, /*157,*/ /*174,*/  /*182,183,*/ 189,193,203,225) then dtdx_obeca_nobcop=dtdxca;
     
     othsite = .; /* multiple myeloma and esophagus, stomach, gallbladder, and liver cancers*/
     if 
         esoph_ca=1       or
         stomach_ca=1     or /*Gastric cardia cancer*/
         liver_ca=1       or
         gallbladder_ca=1 or
         myeloma_ca=1     or
         menin_ca =1     
     then  othsite=1;
     dtdx_othsite = .;
     if obese_ca=1 and dtdxca>0 and ca=1 and icda2 in (150,  151, 155,156,  203 ) then dtdx_othsite=dtdxca;
    
/*==============================================================*/
/* Obesity-related cancer subtype: 13 cancer sites              */
/*==============================================================*/

obeca_site = .;

if ca=1 and dtdxca>0 then do;
    select (icda2);
        when (150)       obeca_site = 1;   /* esophagus */
        when (151)       obeca_site = 2;   /* gastric cardia: verify specificity */
        when (153,154)   obeca_site = 3;   /* colorectal */
        when (155)       obeca_site = 4;   /* liver */
        when (156)       obeca_site = 5;   /* gallbladder */
        when (157)       obeca_site = 6;   /* pancreas */
        when (174)       obeca_site = 7;   /* breast */
        when (182)       obeca_site = 8;   /* endometrium */
        when (183)       obeca_site = 9;   /* ovary */
        when (189)       obeca_site = 10;  /* kidney */
        when (193)       obeca_site = 11;  /* thyroid */
        when (203)       obeca_site = 12;  /* multiple myeloma */
        when (225)       obeca_site = 13;  /* meningioma */
        otherwise        obeca_site = .;
    end;
end;

     keep id conf dtdxca dtdx_: obeca_site 
          /*obese_ca nonobeca smoke_ca r_smoke_ca  
          esoph_ca colorect_ca stomach_ca liver_ca gallbladder_ca panc_ca kidney_ca thyroid_ca myeloma_ca menin_ca */
          ;
  run;




/**** Diabetes *****/
data diabetes;
    infile 'db8921.080423.cases' lrecl=183 recfm=d;
    
    input
    @1   id        6.
    @7    cd    $   1.
    @10   cohort        1.
    @13   icda    $    5.
    @20   qyr       2.
    @26   mdx1        2.
    @28   ydx1        2.
    @32   conf        2.
    @36   nhsicda       5.
    @43   mdx       2.
    @45   ydx       2.
    @49   dxmonth       4.
    @66   type        1.
    @67   probabil      1.
    ;

    label
        id    = 'ID'
        cd    = 'Check Digit'
        cohort    = 'Cohort'
        icda    = 'Reported ICDA'
        qyr   = 'Questionnaire cycle year when reported'
        mdx1    = 'Reported month of diagnosis'
        ydx1    = 'Reported year of diagnosis'
        conf    = 'Confirmation Code'
        nhsicda   = 'NHS reviewed and confirmed ICDA code'
        mdx   = 'Confirmed Month of Diagnosis'
        ydx   = 'Confirmed Year of Diagnosis'
        dxmonth   = 'Date of dx months since 1900'
        type    = 'Type of DM'
        probabil  = 'Ranks certainty of diabetes type status'
        ;

    dbcase=0;
    if type=2 and probabil=1 then dbcase=1;
    
    dtdxdb=.;
    if 0<dxmonth<9999 then dtdxdb=dxmonth;

    dtdxdb2=.;
    if 0<dxmonth<9999 and dbcase=1 then dtdxdb2=dxmonth;
keep id dtdxdb2;
run;



/************************************* CVD data **********************************/
/************ Read confirmed MI files, fatal & nonfatal both included ************/

/***MI***/
data mi;
%include 'mi8917.073120.cases.input';
dtdxmi=dxmonth;
keep id dtdxmi; run;

/***STROKE***/
data stroke;
%include 'str8919.031621.cases.input';
dtdxstr=dxmonth; /*total stroke*/
keep id dtdxstr; run;

/***CABG***/
data cabg;
%include 'cabg8913.060816.nodups.input';
if 0<dxmonth<9999 then dtdxcabg=dxmonth;
keep id dtdxcabg;

*Read reported cabg from %nurXX;
%nur15(keep=cabg15 cabgd15);
%nur17(keep=cabg17 cabgd17);
%nur19(keep=cabg19 cabgd19);

data update_cabg;
  merge cabg nur15 nur17 nur19;
  by id;

if dtdxcabg=. then do;
     if cabg19=1 then do;
       if cabgd19=1 then dtdxcabg=1410;
       else if cabgd19=2 then dtdxcabg=1422;
       else if cabgd19=3 then dtdxcabg=1434;
       else                   dtdxcabg=1422;
     end;

     if cabg17=1 then do;
       if cabgd17=1 then dtdxcabg=1386;
       else if cabgd17=2 then dtdxcabg=1398;
       else if cabgd17=3 then dtdxcabg=1410;
       else                   dtdxcabg=1398;
     end;

     if cabg15=1 then do;
       if cabgd15=1 then dtdxcabg=1362;
       else if cabgd15=2 then dtdxcabg=1374;
       else if cabgd15=3 then dtdxcabg=1386;
       else                   dtdxcabg=1374;
     end;
end;

keep id dtdxcabg;


/*** disease outcomes ***/    
data disease;
    merge dead cancer diabetes mi stroke update_cabg;
    by id;

    /* CVD */
    dtdxchd=.;
    if 0<dtdxmi<9999 or  0<dtdxcabg<9999 then dtdxchd=min(dtdxmi,dtdxcabg, .);
    
    dtdxcvd=.;
    if 0<dtdxchd<9999 or 0<dtdxstr<9999  then dtdxcvd=min(dtdxchd,dtdxstr,.);

    keep id dtdth dtdxca dtdx_: newicda3 dtdxdb2 dtdxcvd dtdxchd dtdxmi dtdxstr obeca_site ; run;
    proc sort; by id; run;

    PROC DATASETS;
    delete diabetes dead cancer mi stroke update_cabg cabg nur15 nur17 nur19 ; run;


*----------------2.Covariables-----------------------------;
** Endoscopy **;
*%include 'endoscopy_nhs2.sas';


%ahei2010_9119(keep=
                  ahei2010_noETOH91 ahei2010_noETOH95 ahei2010_noETOH99 ahei2010_noETOH03 ahei2010_noETOH07 
                  ahei2010_noETOH11 ahei2010_noETOH15 ahei2010_noETOH19
                  ahei91 ahei95 ahei99 ahei03 ahei07 ahei11 ahei15 ahei19);

                  ahei91=ahei2010_noETOH91;
                  ahei95=ahei2010_noETOH95;
                  ahei99=ahei2010_noETOH99;
                  ahei03=ahei2010_noETOH03;
                  ahei07=ahei2010_noETOH07;
                  ahei11=ahei2010_noETOH11;
                  ahei15=ahei2010_noETOH15;
                  ahei19=ahei2010_noETOH19;

/* Physical activity - METs/week */               
%act8917(keep=act89m  act91m  act97m  act01m  act05m  act09m act13m act17m);  

array actm{*} act89m  act91m  act97m  act01m  act05m  act09m act13m act17m;
    do i=1 to dim(actm); if actm{i}>900 then actm{i}=.; end;
    /* carry forward missing physical activity */
    do i=2 to dim(actm); if actm{i}=.  then actm{i}=actm{i-1}; end; 
run;

/* census tract-level neighborhood socioeconomic (SES) status */
%nses8917(keep=id nSES_89 nSES_91 nSES_93 nSES_95 nSES_97 nSES_99 nSES_01 nSES_03 nSES_05 nSES_07 nSES_09 nSES_11 nSES_13 nSES_15 nSES_17);
run;


/*   covariables   */
%der8919(keep= birthday height89 corr_ht   
               race race8905 white      
               retmo89  retmo91  retmo93  retmo95  retmo97  retmo99
               retmo01  retmo03  retmo05  retmo07  retmo09  retmo11  retmo13  retmo15  retmo17 retmo19 
               age89    age91    age93    age95    age97    age99
               age01    age03    age05    age07    age09    age11    age13    age15    age17    age19
               agemo89  agemo91  agemo93  agemo95  agemo97  agemo99
               agemo01  agemo03  agemo05  agemo07  agemo09  agemo11  agemo13  agemo15  agemo17  agemo19
               bmi89    bmi91    bmi93    bmi95    bmi97    bmi99
               bmi01    bmi03    bmi05    bmi07    bmi09    bmi11    bmi13    bmi15    bmi17    bmi19
               smkdr89  smkdr91  smkdr93  smkdr95  smkdr97  smkdr99
               smkdr01  smkdr03  smkdr05  smkdr07  smkdr09  smkdr11  smkdr13  smkdr15  smkdr17  smkdr19
               mnpst89  mnpst91  mnpst93  mnpst95  mnpst97  mnpst99  /* menopausal status 1.pre;\ 2.post;\ 3.dubious;\ 4.unsure;\ 5.missing*/
               mnpst01  mnpst03  mnpst05  mnpst07  mnpst09  mnpst11  mnpst13  mnpst15  mnpst17  mnpst19
               n2rh89   n2rh91   n2rh93   n2rh95   n2rh97   n2rh99  
               /* PMH status 1.never user;\2.past pmh use;\3.current pmh use;\4.past pmh use this period;\5.prior pmh use, current status unk;\9.missing*/
               n2rh01   n2rh03   n2rh05   n2rh07   n2rh09   n2rh11   n2rh13   n2rh15   n2rh17  n2rh19);

               if id=417832 then delete; *this subject died before baseline - she should not be in der8901;
               race=race8905;
               if race8905=1 then white=1; else white=0;
               run;

%nur89(keep=
       seuro89 ocauc89 hisp89 scand89 afric89 asian89 oanc89 /* Ancestry */ 
       nonwh africd asiand oancd hispd                       /* Ancestry new variables */ 
       wt89 wt1889 
       ht89 ft89 inch89 /* for height in inches*/
       mob89    /* Month Birth */  
       yob89    /* Year Birth */  

       mifh89 ammi89 afmi89 /* Family history of MI */ 
       dbfh89 mdb89 fdb89 bdb89 sdb89    /* Family history of diabetes */ 
       cafh89 ambrc89 asbrc89 bclc89 sclc89 mclc89 fclc89 smel89 fmel89 mmel89 bmel89 /* Family history of cancer */ 
       db89     /* Diabetes Not in Pregnancy */ 
       dbg89    /* Gestational diabetes */ 
       mi89     /* MI/Angina */ 
       str89    /* Stroke/TIA */ 
       can89 brcn89  ocan89  mel89 /*cancer*/
       hbp89 chol89
       mamm89 ); 

       ht89=(12*ft89)+inch89;
       /* Coding ancestry */
       if afric89=1 or asian89=1 or oanc89=1 then nonwh=1;
       else nonwh=0;
       if afric89=1 then africd=1;
       else africd=0;
       if asian89=1 then asiand=1;
       else asiand=0;
       if oanc89=1 then oancd=1;
       else oancd=0;
       if hisp89=1 then hispd=1;
       else hispd=0; 

       /* Coding family history of diabetes, MI, cancer */
        dbfh89=0;
        if mdb89=1 or fdb89=1 or bdb89=1 or sdb89=1 then dbfh89=1;
        mifh89=0;
        if ammi89 in (1,2,3,4,5) or afmi89 in (1,2,3,4,5) then mifh89=1;
        cafh89=0;
        if ambrc89=1 or asbrc89=1 or bclc89=1 or sclc89=1 or mclc89=1 or fclc89=1 or smel89=1 or fmel89=1 or mmel89=1 or bmel89=1 then cafh89=1;

        can89=0; if brcn89=1 or ocan89=1 or mel89=1 then can89=1;
        run;

%nur91 (keep=wt91 hbp91 db91 chol91 mi91 ang91 str91 
             oc91 can91 brcn91 ocan91 mel91 mamm91 batch91 sigm91 );
        can91=0; if brcn91=1 or ocan91=1 then can91=1;
        run;


%serv91 (keep=
/* whole grains */  cer_s91  ckcer_s91 oat_s91 dkbr_s91 brice_s91 otgrn_s91 bran_s91 oatbr_s91 wgerm_s91 popc_s91 tort_s91
/* refined grains */  whbr_s91 engl_s91 muff_s91 wrice_s91 pasta_s91 pcake_s91  
/* fruits */ rais_s91 prune_s91 ban_s91 cant_s91 appl_s91 avo_s91 peach_s91 oran_s91 grfr_s91 straw_s91 blueb_s91  
/* vegetables */ tom_s91 toj_s91 tosau_s91 chsau_s91 brocc_s91 cabb_s91 cauli_s91 bruss_s91 rcar_s91 ccar_s91 mixv_s91 yam_s91 osqua_s91 eggpl_s91 kale_s91 cspin_s91 rspin_s91 ilett_s91 rlett_s91 cel_s91 beet_s91 oniog_s91 oniov_s91 grpep_s91 garl_s91 corn_s91
/* legumes */ sbean_s91 tofu_s91 bean_s91 peas_s91  
/* tea-coffee */ tea_s91 dcaf_s91 coff_s91 

/* potatoes */  fries_s91 pot_s91 pchip_s91
/* dairy */ skim_s91 whole_s91 sherb_s91 icecr_s91 yog_s91 cream_s91 sour_s91 cotch_s91 crmch_s91 otch_s91
/* eggs */ egg_s91
/* fish */ ctuna_s91 dkfsh_s91 ofish_s91 shrim_s91
/* fruit juices */ aj_s91 oj_s91 grj_s91 othj_s91  

/* ssb */  cola_s91 cnoc_s91 otsug_s91 punch_s91   
/* sweets */ choco_s91 cdyw_s91 cdywo_s91 cokh_s91 cokr_s91 brwn_s91 donut_s91 cakh_s91 cakr_s91 pieh_s91 pier_s91 srolr_s91 srolh_s91 jam_s91
/* diet bev */ lccaf_s91 lcnoc_s91 lcoth_s91
/* meat */ chwi_s91 chwo_s91 turk_s91 bacon_s91 dog_s91 procm_s91 hamb_s91 bmix_s91 pmain_s91 bmain_s91 livb_s91 livc_s91
/* rpmeats_s91 ssb_s91 swtsvg_s91 swtsnochoc_s91 */
rwine_s91 wwine_s91
);  /* continuous variable, servings/day */

%n91_dt(keep=  cerbr91d );

data serv91;
merge serv91 n91_dt;
by id;

/*** derive daily servings of whole grain, refined grains infomation provided by David Jacobs*****/
if cerbr91d in (0, 4, 11, 12, 13, 15, 16, 17, 19, 20, 21, 22, 23, 25, 29, 30, 31, 33, 35, 37, 38, 42, 43, 45, 46, 47, 59, 60, 
61, 62, 63, 64, 65, 66, 68, 73, 74, 75, 76, 79,80, 81, 82, 83, 84, 85, 86, 88, 95, 96, 97, 98, 101, 109,
111, 112, 113, 114, 115, 119, 123, 124, 127, 128, 129, 131, 132, 133, 134, 146, 147, 148, 149, 150, 151, 152, 154, 155, 157,  158,
130, 164, 165, 166, 167, 170, 173, 176,  177, 182, 184, 185, 186, 188/*,190,197,198,997*/ )

then rcer_s91=cer_s91;
      
else wcer_s91=cer_s91;  

if rcer_s91=. then rcer_s91=0;      
if wcer_s91=. then wcer_s91=0;   
/* take sum of caffeinated and decaffeinated tea & coffee */
tea_s91  = SUM(tea_s91);
coff_s91  = sum(coff_s91, dcaf_s91);
/** whole grains **/
whgrn_s91 =sum(wcer_s91, dkbr_s91, popc_s91, oat_s91, ckcer_s91, wgerm_s91, brice_s91, oatbr_s91, bran_s91, otgrn_s91, popc_s91);
/** refined grains **/
rgrain_s91 =sum(rcer_s91, whbr_s91, engl_s91, muff_s91, wrice_s91, pasta_s91, pcake_s91);
/* fruit and vegetable*/
frui_s91 =sum(rais_s91,  ban_s91, cant_s91,avo_s91, appl_s91,peach_s91, straw_s91,blueb_s91, prune_s91, oran_s91, grfr_s91);
vege_s91 =sum(tom_s91, toj_s91, tosau_s91, chsau_s91, tofu_s91,sbean_s91,brocc_s91, cabb_s91,cauli_s91,bruss_s91, rcar_s91,ccar_s91, corn_s91,
                  peas_s91, mixv_s91, bean_s91,osqua_s91,eggpl_s91,  yam_s91,cspin_s91,rspin_s91,kale_s91,ilett_s91,  
                  rlett_s91,  cel_s91, beet_s91,oniog_s91,oniov_s91);
fv_s91 =sum(frui_s91, vege_s91);
/* potatoes */
pota_s91=sum(fries_s91, pot_s91, pchip_s91);
/*legumes */
leg_s91=sum(sbean_s91, tofu_s91, bean_s91, peas_s91);
/* fruit juices */
fj_s91=sum(aj_s91, oj_s91, grj_s91, othj_s91  );
/* vegetables minus potatoes and legumes */
veg_noleg_s91=vege_s91-leg_s91;
/* dairy and fish */
dairy_s91=sum(skim_s91, whole_s91, sherb_s91, icecr_s91, yog_s91, cream_s91, sour_s91 ,cotch_s91 ,crmch_s91, otch_s91);
fish_s91=sum(ctuna_s91, dkfsh_s91, ofish_s91, shrim_s91);
/* poultry */
poul_s91=sum(chwi_s91, chwo_s91, turk_s91);
/* red and processed meat */
rpmeats_s91= sum(bacon_s91, dog_s91, procm_s91, hamb_s91, bmix_s91, pmain_s91, bmain_s91);
/* ssb*/
ssb_s91  = sum( cola_s91, cnoc_s91, otsug_s91, punch_s91);

/*** SWEETS & DESERTS  ***/
swtsvg_s91  = sum(choco_s91, cdyw_s91, cdywo_s91, cokh_s91, cokr_s91, brwn_s91, donut_s91, cakh_s91, cakr_s91, pieh_s91, pier_s91, srolr_s91, srolh_s91, jam_s91);
swtsnochoc_s91  = sum( cdyw_s91, cdywo_s91, cokh_s91, cokr_s91, brwn_s91, donut_s91, cakh_s91, cakr_s91, pieh_s91, pier_s91, srolr_s91, srolh_s91, jam_s91);
run;


%nur93 (keep=wt93 db93 hbp93 chol93 mi93 ang93 str93 
             oc93 brcn93  ocan93 mel93 can93  
             mmi93  fmi93  mifh93 mamm93 batch93 sigm93 );
        can93=0; if brcn93=1 or ocan93=1 or mel93=1 then can93=1;
        mifh93=0; if mmi93=1 or fmi93=1 then mifh93=1;
        run;

%nur95 (keep=wt95 hbp95 db95 chol95 mi95 ang95 str95 cabg95 
             oc95 brcn95  ocan95 mel95 can95 mamm95 batch95 sigm95 );  
        can95=0; if brcn95=1 or ocan95=1 or mel95=1 then can95=1;
        run;

%serv95 (keep=
/* whole grains */  cer_s95 ckcer_s95 oat_s95 dkbr_s95 brice_s95 otgrn_s95 bran_s95 oatbr_s95 wgerm_s95 popc_s95 tort_s95
/* fruits */ rais_s95 prune_s95 ban_s95 cant_s95 avo_s95 appl_s95 peach_s95 oran_s95 grfr_s95 straw_s95 blueb_s95  
/* vegetables */ tom_s95 toj_s95 tosau_s95 salsa_s95 brocc_s95 cabb_s95 cauli_s95 bruss_s95 rcar_s95 ccar_s95 mixv_s95 yam_s95 osqua_s95 eggpl_s95 kale_s95 cspin_s95 rspin_s95 ilett_s95 rlett_s95 cel_s95 grpep_s95 oniog_s95 oniov_s95 garl_s95 corn_s95
/* legumes */ sbean_s95 tofu_s95 bean_s95 peas_s95  
/* tea-coffee */ tea_s95 dtea_s95 dcaf_s95 coff_s95 

/* potatoes */  fries_s95 pot_s95 pchip_s95
/* dairy */ skim_s95 m1or2_s95 whole_s95 sherb_s95 icecr_s95 flyog_s95 plyog_s95 cream_s95 cotch_s95 crmch_s95 otch_s95
/* eggs */ egg_s95
/* fish */ ctuna_s95 dkfsh_s95 bfsh_s95 ofish_s95 shrim_s95
/* fruit juices */ aj_s95 oj_s95 grj_s95 othj_s95  

/* ssb */  cola_s95 cnoc_s95 otsug_s95 punch_s95   
/* sweets */ choco_s95 cdyw_s95 cdywo_s95 cokh_s95 cokr_s95 brwn_s95 donut_s95 cakh_s95 cakr_s95 pieh_s95 pier_s95 srolr_s95 srolh_s95 jam_s95
/* meat */ chwi_s95 chwo_s95 bacon_s95 bfdog_s95 ctdog_s95 procm_s95 hamb_s95 hambl_s95 bmix_s95 pmain_s95 bmain_s95 livb_s95 livc_s95
/* rpmeats_s95 ssb_s95 swtsvg_s95 swtsnochoc_s95 */
rwine_s95 wwine_s95
);  

%n95_dt(keep=   cerbr95d );

data serv95;
merge serv95 n95_dt;
by id;

if cerbr95d in (0 , 4, 11, 12, 13, 15, 16, 17, 19, 20, 21, 22, 23, 25, 29, 30, 31, 33, 35, 37, 38, 42, 43, 45, 46,  47, 59, 60, 
61, 62, 63, 64, 65, 66, 68, 73, 74, 75, 76, 79, 80, 81, 82, 83, 84, 85, 86, 88, 95, 96, 97, 98, 101, 109,
111, 112, 113, 114, 115, 119, 123, 124, 127, 128, 129, 130, 131, 132, 133, 134, 146, 147, 148, 149, 150, 151, 152, 153,
154, 155, 157, 158, 164, 165, 166, 167,  170, 173,176, 177, 182, 184, 185, 186, 188/*,190,197,198,997*/)

then rcer_s95=cer_s95;      
else wcer_s95=cer_s95;    
if rcer_s95=. then rcer_s95=0;   
if wcer_s95=. then wcer_s95=0;   

/* take sum of caffeinated and decaffeinated tea & coffee */
tea_s95  = SUM(tea_s95, dtea_s95);
coff_s95  = sum(coff_s95, dcaf_s95);
/** whole grains **/
whgrn_s95 =sum(wcer_s95, dkbr_s95, popc_s95, oat_s95, ckcer_s95, wgerm_s95, brice_s95, oatbr_s95, bran_s95, otgrn_s95);
/* fruit and vegetable*/
frui_s95 =sum(rais_s95,  ban_s95, cant_s95, avo_s95, appl_s95,peach_s95, straw_s95,blueb_s95, prune_s95, oran_s95, grfr_s95);
vege_s95 =sum(tom_s95  ,   toj_s95  ,   tosau_s95,   salsa_s95,   tofu_s95 ,   sbean_s95,   brocc_s95,   
               cabb_s95 ,   cauli_s95,   bruss_s95,   rcar_s95 ,   ccar_s95 ,   corn_s95 ,   peas_s95 ,   
               mixv_s95 ,   bean_s95 ,   osqua_s95,   eggpl_s95,   yam_s95  ,   cspin_s95,   rspin_s95,   
               kale_s95 ,   ilett_s95,   rlett_s95,   cel_s95  ,   grpep_s95,   oniog_s95,   oniov_s95 );
fv_s95 =sum(frui_s95, vege_s95);
/* potatoes */
pota_s95=sum(fries_s95, pot_s95, pchip_s95);
/*legumes */
leg_s95=sum(sbean_s95, tofu_s95, bean_s95, peas_s95);
/* fruit juices */
fj_s95=sum(aj_s95, oj_s95, grj_s95, othj_s95  );
/* vegetables minus potatoes and legumes */
veg_noleg_s95=vege_s95-leg_s95;
/* dairy and fish */
dairy_s95=sum(skim_s95, m1or2_s95, whole_s95, sherb_s95, icecr_s95, flyog_s95, plyog_s95, cream_s95, cotch_s95, crmch_s95, otch_s95);
fish_s95=sum(ctuna_s95, dkfsh_s95, bfsh_s95, ofish_s95, shrim_s95);
/* poultry */
poul_s95=sum(chwi_s95, chwo_s95, ctdog_s95);
/* red and processed meat */
rpmeats_s95= sum(bacon_s95, bfdog_s95, procm_s95, hamb_s95, hambl_s95, bmix_s95, pmain_s95, bmain_s95);
/* ssb*/
ssb_s95  = sum( cola_s95, cnoc_s95, otsug_s95, punch_s95);
/*** SWEETS & DESERTS  ***/
swtsvg_s95  = sum(choco_s95, cdyw_s95, cdywo_s95, cokh_s95, cokr_s95, brwn_s95, donut_s95, cakh_s95, cakr_s95, pieh_s95, pier_s95, srolr_s95, srolh_s95, jam_s95);
swtsnochoc_s95  = sum(cdyw_s95, cdywo_s95, cokh_s95, cokr_s95, brwn_s95, donut_s95, cakh_s95, cakr_s95, pieh_s95, pier_s95, srolr_s95, srolh_s95, jam_s95);
run;


%nur97 (keep=hbp97 db97 chol97 mi97 ang97 str97 cabg97 
             dbfh97 pdb97  sbdb97  
             mifh97 mmi97 fmi97 sbmi97 strfh97 mstr97 fstr97 sbstr97
             oc97 brcn97  ocan97 mel97 can97 
             cafh97 mbrcn97 sbrc197 sbrc297 sov97 mov97 pmel97 smel97 sclc197 sclc297 pclc97
             mamno97 mamsc97 mamsy97 mampt97 batch97 sigmx97 sigsc97 sigsy97 sigpt97 );

            if mmi97=1  or fmi97=1  or sbmi97=1  then mifh97=1; else mifh97=0;
            if mstr97=1 or fstr97=1 or sbstr97=1 then strfh97=1; else  strfh97=0;
            if pdb97=1 or sbdb97=1 then dbfh97=1; else dbfh97=0;
            if mbrcn97=1 or sbrc197=1 or sbrc297=1 or sov97=1 or mov97=1 or pmel97=1 or smel97=1 or sclc197=1 or sclc297=1 or pclc97=1 then cafh97=1; else cafh97=0;
            can97=0; if brcn97=1 or ocan97=1 or mel97=1 then can97=1;
            run;


%nur99 (keep=wt99 hbp99 db99 chol99 mi99 ang99 str99 cabg99 
             oc99 brcn99 ocan99 mel99 can99
             mamno99 mamsc99 mamsy99 mampt99 q99 sigmx99 sigsc99 sigsy99 sigpt99 );
             can99=0; if brcn99=1 or ocan99=1 or mel99=1 then can99=1;
             run;
%serv99 (keep=
/* whole grains */  cer_s99 ckcer_s99 oat_s99 dkbr_s99 brice_s99 otgrn_s99 bran_s99 oatbr_s99 wgerm_s99 popc_s99 tort_s99  
/* fruits */ rais_s99 prune_s99 ban_s99 cant_s99 avo_s99 appl_s99 peach_s99 oran_s99 grfr_s99 straw_s99 blueb_s99  
/* vegetables */ tom_s99 toj_s99 tosau_s99 salsa_s99 brocc_s99 cabb_s99 cauli_s99 bruss_s99 rcar_s99 ccar_s99 mixv_s99 yam_s99 osqua_s99 eggpl_s99 kale_s99 cspin_s99 rspin_s99 ilett_s99 rlett_s99 cel_s99 grpep_s99 oniog_s99 oniov_s99 garl_s99 corn_s99
/* legumes */ sbean_s99 tofu_s99 bean_s99 peas_s99 soy_s99  

/* potatoes */  fries_s99 pot_s99 pchip_s99
/* dairy */ skim_s99 m2_s99 whole_s99 sherb_s99 icecr_s99 flyog_s99 plyog_s99 cream_s99 cotch_s99 crmch_s99 otch_s99
/* eggs */ egg_s99 eggwh_s99
/* fish */ ctuna_s99 dkfsh_s99 bfsh_s99 ofish_s99 shrim_s99
/* fruit juices */ aj_s99 oj_s99 ojca_s99 prunj_s99 othj_s99  

/* tea-coffee */ tea_s99 dtea_s99 dcaf_s99 coff_s99 
/* ssb */  cola_s99 otsnc_s99 otsug_s99 punch_s99   
/* sweets */ choco_s99 cdyw_s99 cdywo_s99 coknf_s99 cokh_s99 cokr_s99 brwn_s99 donut_s99 cakr_s99 cakh_s99 pieh_s99 srolf_s99 srolr_s99 srolh_s99 jam_s99
/* meat */ chksa_s99 chwi_s99 chwo_s99 bacon_s99 bfdog_s99 ctdog_s99 pmsan_s99 procm_s99 hamb_s99 hambl_s99 bmix_s99 pmain_s99 bmain_s99 livb_s99 livc_s99
/* rpmeats_s99 ssb_s99 swtsvg_s99 swtsnochoc_s99  */
rwine_s99 wwine_s99
);  


%n99_dt(keep=   cerbr99d );

data serv99;
merge serv99 n99_dt;
by id;

if cerbr99d in (0, 4, 11, 12, 13, 15, 16, 17, 19, 20, 21, 22, 23, 25, 29, 30, 31, 33, 35, 37, 38, 42, 43, 45, 46, 47, 59, 60, 
61, 62, 63, 64, 65, 66, 68, 73, 74, 75, 76, 79,80, 81, 82, 83, 84, 85, 86, 88, 95, 96, 97, 98, 101, 109,
111, 112, 113, 114, 115, 119, 123, 124, 127, 128, 129, 131, 132, 133, 134, 146, 147, 148, 149, 150, 151, 152, 154, 155, 157,  158,
130, 164, 165, 166, 167, 170, 173, 176,  177, 182, 184, 185, 186, 188/*,190,197,198,997*/ )

then rcer_s99=cer_s99;      
else wcer_s99=cer_s99;      
if rcer_s99=. then rcer_s99=0;   
if wcer_s99=. then wcer_s99=0;   

/* take sum of caffeinated and decaffeinated tea & coffee */
tea_s99  = SUM(tea_s99, dtea_s99);
coff_s99  = sum(coff_s99, dcaf_s99);
/** whole grains **/
whgrn_s99 =sum(wcer_s99, dkbr_s99, popc_s99, oat_s99, ckcer_s99, wgerm_s99, brice_s99, oatbr_s99, bran_s99, otgrn_s99);

/* fruit and vegetable*/
frui_s99 =sum(rais_s99,  ban_s99, cant_s99, avo_s99, appl_s99,peach_s99, straw_s99,blueb_s99, prune_s99, oran_s99, grfr_s99);
vege_s99 =sum(tom_s99   ,  toj_s99   ,  tosau_s99 ,  salsa_s99 ,  tofu_s99  ,  sbean_s99 , soy_s99, brocc_s99 ,  cabb_s99  ,  
                 cauli_s99 ,  bruss_s99 ,  rcar_s99  ,  ccar_s99  ,  corn_s99,  peas_s99  ,  mixv_s99  ,  bean_s99  ,  
                 yam_s99  ,  osqua_s99 ,  eggpl_s99 ,  kale_s99  ,  cspin_s99 ,  rspin_s99 ,  ilett_s99 ,  rlett_s99 ,  
                 cel_s99   ,  grpep_s99 ,  oniog_s99 ,  oniov_s99 );
fv_s99 =sum(frui_s99, vege_s99);

/* potatoes */
pota_s99=sum(fries_s99, pot_s99, pchip_s99);
/*legumes */
leg_s99=sum(sbean_s99, tofu_s99, bean_s99, peas_s99, soy_s99);
/* fruit juices */
fj_s99=sum(aj_s99, oj_s99, ojca_s99, prunj_s99, othj_s99 );
/* vegetables minus potatoes and legumes */
veg_noleg_s99=vege_s99-leg_s99;
/* dairy and fish */
dairy_s99=sum(skim_s99, m2_s99, whole_s99, sherb_s99, icecr_s99, flyog_s99, plyog_s99, cream_s99, cotch_s99, crmch_s99, otch_s99);
fish_s99=sum(ctuna_s99, dkfsh_s99, bfsh_s99, ofish_s99, shrim_s99);
egg_s99=sum(egg_s99, eggwh_s99);
/* poultry */
poul_s99=sum(chksa_s99, chwi_s99, chwo_s99, ctdog_s99 );
/* red and processed meat */
rpmeats_s99 = sum(bacon_s99,  bfdog_s99,  pmsan_s99,  procm_s99,  hamb_s99,  hambl_s99,  bmix_s99,  pmain_s99,  bmain_s99);
/* ssb*/
ssb_s99   = sum( cola_s99,  otsnc_s99,  otsug_s99,  punch_s99);
/*** SWEETS & DESERTS  ***/
swtsvg_s99   = sum(choco_s99,  cdyw_s99,  cdywo_s99,  coknf_s99,  cokh_s99,  cokr_s99,  brwn_s99,  donut_s99,  cakr_s99,  cakh_s99,  pieh_s99,  srolf_s99,  srolr_s99,  srolh_s99,  jam_s99);
swtsnochoc_s99   = sum(cdyw_s99,  cdywo_s99,  coknf_s99,  cokh_s99,  cokr_s99,  brwn_s99,  donut_s99,  cakr_s99,  cakh_s99,  pieh_s99,  srolf_s99,  srolr_s99,  srolh_s99,  jam_s99);
run;

%nur01 (keep=wt01 db01 hbp01 chol01 mi01 ang01 str01 cabg01 
             dbfh01 mdb01 fdb01 sdb01   
             mifh01 mmi01 fmi01 smi01  strfh01 mstr01 fstr01 sstr01
             oc01  brcn01  ocan01 mel01 can01
             cafh01 mov01 sov01 mbrcn01 sbrc101 sbrc201 pclc01 sclc101 sclc201 fpro01 bpro01 mut01 smut01 ppan01 span01 smel01 pmel01
             mamno01 mamsc01 mamsy01 mampt01 
             q01         sigm01  vbld01 diarr01 fecal01 hxclc01 rout01 bariu01 abdpn01 q23ap01 );
             if mmi01=1 or fmi01=1 or smi01=1 then mifh01=1; else mifh01=0;
             if mstr01=1 or fstr01=1 or sstr01=1 then strfh01=1; else strfh01=0;
             if fdb01=1 or mdb01=1 or sdb01=1 then dbfh01=1; else dbfh01=0;
             if mov01=1 or sov01=1 or mbrcn01=1 or sbrc101=1 or sbrc201=1 or pclc01=1 or sclc101=1 or sclc201=1 or 
                fpro01=1 or bpro01=1 or mut01=1 or smut01=1 or ppan01=1 or span01=1 or smel01=1 or pmel01=1 then cafh01=1; else cafh01=0; 
             can01=0; if brcn01=1 or ocan01=1 or mel01=1 then can01=1;
             run;

%nur03(keep= wt03 db03 hbp03 chol03 mi03 ang03 str03 cabg03 brcn03 ocan03 mel03 can03
             cerbr03d
             q03 colsc03  sigm03  vbld03 diarr03 fecal03 hxclc03  rout03 bariu03 abdpn03 q12ap03 virt03 cpolp03 );  
             can03=0; if brcn03=1 or ocan03=1 or mel03=1 then can03=1;
             run;

%serv03 (keep=
/* whole grains */  cer_s03 ckcer_s03 oat_s03 dkbr_s03 ryebr_s03 brice_s03 oatbr_s03 bran_s03 wgerm_s03 popc_s03 ffpop_s03 tort_s03   
/* fruits */ rais_s03 grape_s03 prune_s03 ban_s03 cant_s03 appl_s03 apsau_s03 peach_s03 oran_s03 grfr_s03 straw_s03 blueb_s03  
/* legumes */ sbean_s03 tofu_s03 bean_s03 peas_s03 soy_s03 

/* potatoes */  fries_s03 pot_s03 pchip_s03
/* dairy */ skim_s03 m1or2_s03 whole_s03 sherb_s03 icecr_s03 flyog_s03 plyog_s03 cream_s03 cotch_s03 crmch_s03 otch_s03
/* eggs */ egg_s03 eggwh_s03
/* fish */ ctuna_s03 dkfsh_s03 bfsh_s03 ofish_s03 shrim_s03
/* fruit juices */ aj_s03 oj_s03 ojca_s03 prunj_s03 othj_s03  

/* vegetables */ tom_s03 toj_s03 tosau_s03 salsa_s03 brocc_s03 cauli_s03 cabb_s03 bruss_s03 rcar_s03 ccar_s03 mixv_s03 yam_s03 osqua_s03 eggpl_s03 kale_s03 cspin_s03 rspin_s03 ilett_s03 rlett_s03 cel_s03 grpep_s03 oniog_s03 oniov_s03 garl_s03 corn_s03
/* tea-coffee */ tea_s03 dtea_s03 dcaf_s03 coff_s03 
/* ssb */  cola_s03 otsug_s03 punch_s03   
/* sweets */ choc_s03 cdyw_s03 cdywo_s03 coknf_s03 cokr_s03 cokh_s03 brwn_s03 donut_s03 cakr_s03 cakh_s03 pieh_s03 srolf_s03 srolr_s03 srolh_s03 jam_s03
/* meat */ chksa_s03 chwi_s03 chwo_s03 bacon_s03 bfdog_s03 ctdog_s03 pmsan_s03 procm_s03 hamb_s03 hambl_s03 bmix_s03 pmain_s03 bmain_s03 livb_s03 livc_s03
rwine_s03 wwine_s03
/* choco_s03 rpmeats_s03 ssb_s03 swtsvg_s03 swtsnochoc_s03 */
);  

data serv03;
  merge serv03 nur03(keep=id cerbr03d);
  by id;

if cerbr03d in (0, 4, 11, 12, 13, 15, 16, 17, 19, 20, 21, 22, 23, 25, 29, 30, 31, 33, 35, 37, 38, 42, 43, 45, 46, 47, 59, 60, 
61, 62, 63, 64, 65, 66, 68, 73, 74, 75, 76, 79,80, 81, 82, 83, 84, 85, 86, 88, 95, 96, 97, 98, 101, 109,
111, 112, 113, 114, 115, 119, 123, 124, 127, 128, 129, 131, 132, 133, 134, 146, 147, 148, 149, 150, 151, 152, 154, 155, 157,  158,
130, 164, 165, 166, 167, 170, 173, 176,  177, 182, 184, 185, 186, 188/*,190,197,198,997*/ )

then rcer_s03=cer_s03;
else wcer_s03=cer_s03;
if rcer_s03=. then rcer_s03=0;      
if wcer_s03=. then wcer_s03=0;

choco_s03=choc_s03;  /* after 2003 the name is different */

/* take sum of caffeinated and decaffeinated tea & coffee */
tea_s03  = SUM(tea_s03, dtea_s03);
coff_s03  = sum(coff_s03, dcaf_s03);
/** whole grains **/
whgrn_s03 =sum(wcer_s03, dkbr_s03, popc_s03, oat_s03, ckcer_s03, wgerm_s03, brice_s03, oatbr_s03, bran_s03);
/* fruit and vegetable*/
frui_s03 =sum(rais_s03 ,  prune_s03,            ban_s03  ,  cant_s03 ,   /*avo_s03,*/         appl_s03 ,  
                          oran_s03 ,                      grfr_s03 ,            straw_s03,  blueb_s03,  peach_s03  );
vege_s03 =sum(tom_s03,  toj_s03,   tosau_s03,  salsa_s03,   tofu_s03, sbean_s03,  bean_s03,    peas_s03,  soy_s03,
                brocc_s03,cauli_s03,    cabb_s03,  bruss_s03,   rcar_s03,  ccar_s03,  corn_s03,    mixv_s03,  
                 yam_s03,osqua_s03,   eggpl_s03,   kale_s03,  cspin_s03, rspin_s03, ilett_s03,   rlett_s03,  
                  cel_s03,grpep_s03,   oniog_s03,  oniov_s03 );
fv_s03 =sum(frui_s03, vege_s03);
/* potatoes */
pota_s03=sum(fries_s03, pot_s03, pchip_s03);
/*legumes */
leg_s03=sum(sbean_s03, tofu_s03, bean_s03, peas_s03, soy_s03);
/* fruit juices */
fj_s03=sum(aj_s03, oj_s03, ojca_s03, prunj_s03, othj_s03  );
/* vegetables minus potatoes and legumes */
veg_noleg_s03=vege_s03-leg_s03;
/* dairy and fish */
dairy_s03=sum(skim_s03, m1or2_s03, whole_s03, sherb_s03, icecr_s03, flyog_s03, plyog_s03, cream_s03, cotch_s03, crmch_s03, otch_s03);
fish_s03=sum(ctuna_s03, dkfsh_s03, bfsh_s03, ofish_s03, shrim_s03);
egg_s03=sum(egg_s03, eggwh_s03);
/* poultry */
poul_s03=sum(chksa_s03, chwi_s03, chwo_s03, ctdog_s03 );
/* red and processed meat */
rpmeats_s03 = sum(bacon_s03, bfdog_s03, pmsan_s03, procm_s03, hamb_s03, hambl_s03, bmix_s03, pmain_s03, bmain_s03);
/* ssb*/
ssb_s03   = sum(cola_s03, otsug_s03, punch_s03);
/*** SWEETS & DESERTS  ***/
swtsvg_s03   = sum(choco_s03, cdyw_s03, cdywo_s03, coknf_s03, cokh_s03, cokr_s03, brwn_s03, donut_s03, cakr_s03, cakh_s03, pieh_s03, srolf_s03, srolr_s03, srolh_s03, jam_s03);
swtsnochoc_s03   = sum(cdyw_s03, cdywo_s03, coknf_s03, cokh_s03, cokr_s03, brwn_s03, donut_s03, cakr_s03, cakh_s03, pieh_s03, srolf_s03, srolr_s03, srolh_s03, jam_s03);
run;

%nur05 (keep=wt05 db05 hbp05 chol05 mi05 ang05 str05 cabg05 
             dbfh05 mdb05 fdb05 sdb05  
             mifh05 mdchd05 fdchd05 mdstr05 fdstr05 strfh05
             oc05 brcn05 ocan05 mel05 can05 
             cafh05 mov05 sov05 mbrcn05 sbrc105 sbrc205 sclc105 pclc05 smel05 pmel05 sov05 mov05 
             mamno05 mamsc05 mamsy05 mampt05 
             q05 colsc05 sigsc05 vbld05 diarr05 fecal05 hxclc05 rout05 bariu05 abdpn05 q13ap05 virt05 cpolp05 virtc05
             signo05 sig9405 sig9505 sig9705 sig9905 sig0105 sig0305 sigpt05
             colno05 col9405 col9505 col9705 col9905 col0105 col0305 colpt05);

             if mdb05=1 or fdb05=1 or sdb05=1 then dbfh05=1; else dbfh05=0;
             if mdchd05=1 or fdchd05=1 then mifh05=1; else mifh05=0;
             if mdstr05=1 or fdstr05=1 then strfh05=1; else strfh05=0;
             can05=0; if brcn05=1 or ocan05=1 or mel05=1 then can05=1;
             if mov05=1 or sov05=1 or mbrcn05=1 or sbrc105=1 or sbrc205=1 or sclc105=1 or
                pclc05=1 or smel05=1 or pmel05=1 or sov05=1 or mov05=1 then cafh05=1; else cafh05=0;
             run;

%nur07(keep= db07 hbp07 chol07 mi07 ang07 str07 cabg07 
             brcn07 ocan07 mel07 can07
             cerbr07d mamno07 mamsc07 mamsy07 mampt07 
             q07 colsc07 sigsc07 vbld07 diarr07 fecal07 hxclc07 rout07 bariu07 abdpn07 q14apt07 virt07 cpolp07 virtc07 );
             can07=0; if brcn07=1 or ocan07=1 or mel07=1 then can07=1;
             run;
%serv07(keep=  
mchoc_s07 /*milk chocolate*/
dchoc_s07 /*dark chocolate*/
cdyw_s07  /*Candy bar*/
cdywo_s07/*Candy wo Chocolate*/
coknf_s07/*Cookies, Fat Free*/
cokr_s07/*Cookies, Commercial*/
cokh_s07/*Cookies, Home Baked*/
/*brwn_s07 Brownies   merged to cokh_s07 in NHS2*/
donut_s07/*Doughnuts*/
cake_s07/*Cake*/
piehr_s07/*Pie, Ready-Made*/
jam_s07 /*Jams, Jellies*/
pbut_s07 /*Peanut buter*/
ffpop_s07/*Popcorn: Fat Free*/
popc_s07/*Popcorn: Regular*/
srolf_s07/*Sweetroll Fat Free*/
srolr_s07/*Sweetroll, Coffee Cake-Ready Made*/
srolh_s07/*Sweetroll, Coffee Cake-Home Made*/
pretz_s07/*Pretzels*/
brbar_s07/*Breakfast Bars*/
enbar_s07/*Energy Bars*/
lcbar_s07/*Low Carb Bars*/
pnut_s07/*Peanuts*/
wnut_s07/*Walnuts*/
onut_s07/*Other Nuts*/

pchip_s07 /*chips*/
icecr_s07 /*Ice cream*/
plyog_s07 /*DAIRY - Yogurt-plain, lo-carb*/
flyog_s07 /*Yogurt*/

rais_s07 prune_s07 ban_s07 cant_s07 avo_s07 appl_s07 oran_s07 grfr_s07 straw_s07 blueb_s07 peach_s07
tom_s07 toj_s07 tosau_s07 salsa_s07 tofu_s07 sbean_s07 bean_s07 peas_s07 soy_s07 cabb_s07 cauli_s07 bruss_s07 rcar_s07 ccar_s07 kale_s07  osqua_s07 
yam_s07 cspin_s07 rspin_s07 rlett_s07 ilett_s07 corn_s07 mixv_s07 eggpl_s07 cel_s07 grpep_s07 oniog_s07 oniov_s07 brocc_s07
cer_s07 dkbr_s07 popc_s07 oat_s07 ckcer_s07 brice_s07 oatbr_s07 bran_s07
bmix_s07 pmain_s07 bmain_s07  pmsan_s07 procm_s07 bacon_s07 bfdog_s07 hamb_s07 hambl_s07
cola_s07 otsug_s07 punch_s07
tea_s07 dtea_s07 coff_s07 dcaf_s07 rwine_s07 wwine_s07
/* fv_s07 fruit_s07 veg_s07 whgrn_s07 cerbr_s07 rpmeats_s07 ssb_s07 wcer_s07  swtsvg_s07 */
/* potatoes */  fries_s07 pot_s07 pchip_s07
/* dairy */ skim_s07 m1or2_s07 whole_s07 sherb_s07 icecr_s07 flyog_s07 plyog_s07 cream_s07 cotch_s07 crmch_s07 otch_s07
/* eggs */ egg_s07 eggom_s07
/* fish */ ctuna_s07 dkfsh_s07 bfsh_s07 ofish_s07 shrim_s07
/* fruit juices */ aj_s07 oj_s07 ojca_s07 prunj_s07 othj_s07  
chksa_s07 chwi_s07 chwo_s07 ctdog_s07
);

data serv07;
merge serv07 nur07(keep=id cerbr07d);
by id;

if cerbr07d in (0, 4, 11, 12, 13, 15, 16, 17, 19, 20, 21, 22, 23, 25, 29, 30, 31, 33, 35, 37, 38, 42, 43, 45, 46, 47, 59, 60, 
61, 62, 63, 64, 65, 66, 68, 73, 74, 75, 76, 79,80, 81, 82, 83, 84, 85, 86, 88, 95, 96, 97, 98, 101, 109,
111, 112, 113, 114, 115, 119, 123, 124, 127, 128, 129, 131, 132, 133, 134, 146, 147, 148, 149, 150, 151, 152, 154, 155, 157,  158,
130, 164, 165, 166, 167, 170, 173, 176,  177, 182, 184, 185, 186, 188 )

then rcer_s07=cer_s07;      
else wcer_s07=cer_s07;  
if rcer_s07=. then rcer_s07=0;      
if wcer_s07=. then wcer_s07=0;     

choco_s07=SUM(mchoc_s07, dchoc_s07);
/* take sum of caffeinated and decaffeinated tea & coffee */
tea_s07  = SUM(tea_s07, dtea_s07);
coffee_s07  = sum(coff_s07, dcaf_s07);
/* whole grain*/
whgrn_s07  = SUM(wcer_s07, dkbr_s07, popc_s07, oat_s07, ckcer_s07, brice_s07, oatbr_s07, bran_s07);
/* fruit and vegetables*/
fruit_s07 = SUM(rais_s07, prune_s07, ban_s07, cant_s07, avo_s07, appl_s07, oran_s07, grfr_s07, straw_s07, blueb_s07, peach_s07);  
vege_s07   = SUM(tom_s07,  toj_s07,   tosau_s07,  salsa_s07,   tofu_s07, sbean_s07,  bean_s07,    peas_s07,  soy_s07,
                brocc_s07,cauli_s07,    cabb_s07,  bruss_s07,   rcar_s07,  ccar_s07,  corn_s07,   mixv_s07,  
                 yam_s07,osqua_s07,   eggpl_s07,   kale_s07,  cspin_s07, rspin_s07, ilett_s07,   rlett_s07,  
                  cel_s07,grpep_s07,   oniog_s07,  oniov_s07    );
fv_s07    = SUM(fruit_s07, vege_s07);

/* potatoes */
pota_s07=sum(fries_s07, pot_s07, pchip_s07);
/*legumes */
leg_s07=sum(sbean_s07, tofu_s07, bean_s07, peas_s07, soy_s07);
/* fruit juices */
fj_s07=sum(aj_s07, oj_s07, ojca_s07, prunj_s07, othj_s07  );
/* vegetables minus potatoes and legumes */
veg_noleg_s07=vege_s07-leg_s07;
/* dairy and fish */
dairy_s07=sum(skim_s07, m1or2_s07, whole_s07, sherb_s07, icecr_s07, flyog_s07, plyog_s07, cream_s07, cotch_s07, crmch_s07, otch_s07);
fish_s07=sum(ctuna_s07, dkfsh_s07, bfsh_s07, ofish_s07, shrim_s07);
egg_s07=sum(egg_s07, eggom_s07);
/* poultry */
poul_s07=sum(chksa_s07, chwi_s07, chwo_s07, ctdog_s07 );
/* Total meat intake - servings added up*/
rpmeats_s07= sum(bmix_s07, pmain_s07, bmain_s07,  pmsan_s07, procm_s07, bacon_s07, bfdog_s07, hamb_s07, hambl_s07);
/* ssb*/
ssb_s07    = SUM(cola_s07, otsug_s07, punch_s07);
/*** SWEETS & DESERTS  ***/
swtsvg_s07= sum(mchoc_s07, dchoc_s07, cdyw_s07, cdywo_s07, coknf_s07, cokh_s07, cokr_s07, donut_s07, cake_s07, piehr_s07, srolf_s07, srolr_s07, srolh_s07, icecr_s07);
swtsnochoc_s07 = sum(cdyw_s07 ,  cdywo_s07,  coknf_s07,  cokr_s07 ,  cokh_s07 ,   donut_s07,  cake_s07 ,  piehr_s07 ,  srolf_s07, srolr_s07, srolh_s07 ,icecr_s07 );
run;

%nur09 (keep=wt09 db09 hbp09 chol09 mi09 str09 ang09 cabg09 
             dbfh09 pdb09 sdb09 mdchd09 fdchd09 mifh09 mdstr09 fdstr09 strfh09
             oc09 brcn09  ocan09 mel09 can09 
             cafh09 mov09 sov09 mbrcn09 sbrcn09 pclc09 sclc109 sclc209 pmel09 smel09
             mamno09 mamsc09 mamsy09 mampt09 
             q09 colsc09 sigsc09 vbld09 diarr09 fecal09 hxclc09 rout09 bariu09 abdpn09 q14ap09  virt09 cpolp09 virtc09 );

             if pdb09=1 or sdb09=1 then dbfh09=1; else dbfh09=0;
             if mdchd09=1 or fdchd09=1 then mifh09=1; else mifh09=0;
             if mdstr09=1 or fdstr09=1 then strfh09=1; else strfh09=0;
             if mov09=1 or sov09=1 or mbrcn09=1 or sbrcn09=1 or pclc09=1 or sclc109=1 or
                sclc209=1 or pmel09=1 or smel09=1 then cafh09=1; else cafh09=0;
             can09=0; if brcn09=1 or ocan09=1 or mel09=1 then can09=1;
             run;

%nur11(keep= wt11 db11 hbp11 chol11 
             mi11 str11 cabg11 ang11
             brcn11 ut11 ov11 colc11 mel11 ocan11 can11
             cerbr11d mamno11 mamsc11 mamsy11 mampt11 
             q11 colsc11 sigsc11 vbld11 diarr11 fecal11 hxclc11 rout11 bariu11 abdpn11 colrsnpt11 virt11 cpolp11 virtc11 );
             can11=0; if brcn11=1 or ut11=1 or ov11=1 or colc11=1 or mel11=1 or ocan11=1 then can11=1;
             run;

%serv11(keep= 
mchoc_s11 /*milk chocolate*/
dchoc_s11 /*dark chocolate*/

cdyw_s11  /*Candy bar*/
cdywo_s11/*Candy wo Chocolate*/
coknf_s11/*Cookies, Fat Free*/
cokr_s11/*Cookies, Commercial*/
cokh_s11/*Cookies, Home Baked*/
donut_s11/*Doughnuts*/
cake_s11 /*Cake, Homemade*/
piehr_s11/*Pie, Ready-Made*/
jam_s11 /*Jams, Jellies*/
pbut_s11 /*Peanut buter*/
ffpop_s11/*Popcorn: Fat Free*/
popc_s11/*Popcorn: Regular*/
sroll_s11/*Sweetroll Fat Free*/
brbar_s11/*Breakfast Bars*/
enbar_s11/*Energy Bars*/
hpbars_s11/*Protein Bars*/
pretz_s11/*Pretzels*/
slimf_s11/*Slimfast*/
ensur_s11/*Ensure*/
mdrfr_s11 /*Mixed Dried Fruit*/
pnut_s11/*Peanuts*/
wnut_s11/*Walnuts*/
onut_s11/*Other Nuts*/

pchip_s11 /*chips*/
icecr_s11 /*Ice cream*/
plyog_s11 /*DAIRY - Yogurt-plain, lo-carb*/
flyog_s11 /*Yogurt*/


rais_s11  prune_s11  ban_s11  cant_s11  avo_s11  appl_s11  oran_s11  grfr_s11  straw_s11  blueb_s11  peach_s11   
brocc_s11  cabb_s11  cauli_s11  bruss_s11  kale_s11 tom_s11  toj_s11  tosau_s11  ketch_s11
rcar_s11  ccar_s11  osqua_s11  yam_s11 cspin_s11  rspin_s11  rlett_s11  ilett_s11
tofu_s11  sbean_s11  peas_s11  bean_s11  soy_s11 corn_s11  mixv_s11  eggpl_s11  cel_s11  grpep_s11  oniog_s11  oniov_s11
cer_s11   dkbr_s11  popc_s11  oat_s11  ckcer_s11  brice_s11  oatbr_s11   wgerm_s11
bmix_s11  pmain_s11  bmain_s11   pmsan_s11  procm_s11  bacon_s11  bfdog_s11  hamb_s11  hambl_s11
cola_s11  otsug_s11  punch_s11
tea_s11 dtea_s11 coff_s11 dcaf_s11  rwine_s11 wwine_s11

/* potatoes */  fries_s11 pot_s11 pchip_s11
/* dairy */ skim_s11 m1or2_s11 whole_s11 sherb_s11 icecr_s11 flyog_s11 plyog_s11 cream_s11 cotch_s11 crmch_s11 otch_s11
/* eggs */ egg_s11 eggom_s11
/* fish */ ctuna_s11 dkfsh_s11 bfsh_s11 ofish_s11 shrim_s11
/* fruit juices */ aj_s11 oj_s11 ojca_s11 prunj_s11 othj_s11  

chksa_s11 chwi_s11 chwo_s11 ctdog_s11
/* fv_s11 fruit_s11 veg_s11 whgrn_s11 rpmeats_s11 ssb_s11 cerbr_s11 wcer_s11 swtsvg_s11 */
);

data serv11;
merge serv11 nur11(keep=id cerbr11d);
by id;

if cerbr11d in (0, 4, 11, 12, 13, 15, 16, 17, 19, 20, 21, 22, 23, 25, 29, 30, 31, 33, 35, 37, 38, 42, 43, 45, 46, 47, 59, 60, 
61, 62, 63, 64, 65, 66, 68, 73, 74, 75, 76, 79,80, 81, 82, 83, 84, 85, 86, 88, 95, 96, 97, 98, 101, 109,
111, 112, 113, 114, 115, 119, 123, 124, 127, 128, 129, 131, 132, 133, 134, 146, 147, 148, 149, 150, 151, 152, 154, 155, 157,  158,
130, 164, 165, 166, 167, 170, 173, 176,  177, 182, 184, 185, 186, 188 )
  then rcer_s11=cer_s11;      
  else wcer_s11=cer_s11;
  if rcer_s11=. then rcer_s11=0;
  if wcer_s11=. then wcer_s11=0;  

choco_s11=SUM(mchoc_s11, dchoc_s11);

/* take sum of caffeinated and decaffeinated tea & coffee */
tea_s11  = SUM(tea_s11, dtea_s11);
coffee_s11  = sum(coff_s11, dcaf_s11);
/* whole grain*/
whgrn_s11  = SUM(wcer_s11, dkbr_s11, popc_s11, oat_s11, ckcer_s11,  brice_s11, oatbr_s11,  wgerm_s11);
/* fruit and vegetables*/
fruit_s11  = SUM(rais_s11, prune_s11, ban_s11, cant_s11, avo_s11, appl_s11, oran_s11, grfr_s11, straw_s11, blueb_s11, peach_s11);
vege_s11    = SUM(brocc_s11, cabb_s11,   cauli_s11,  bruss_s11, kale_s11,   tom_s11,   toj_s11,   tosau_s11, ketch_s11,
                  rcar_s11,  ccar_s11,   osqua_s11,  yam_s11,  cspin_s11,  rspin_s11, rlett_s11, ilett_s11,
                  tofu_s11,  sbean_s11,  peas_s11,   bean_s11,  soy_s11,    corn_s11,  mixv_s11,  eggpl_s11, cel_s11,  grpep_s11,  oniog_s11,  oniov_s11);
fv_s11     = SUM(fruit_s11, vege_s11);

/* potatoes */
pota_s11=sum(fries_s11, pot_s11, pchip_s11);
/*legumes */
leg_s11=sum(sbean_s11, tofu_s11, bean_s11, peas_s11, soy_s11);
/* fruit juices */
fj_s11=sum(aj_s11, oj_s11, ojca_s11, prunj_s11, othj_s11  );
/* vegetables minus potatoes and legumes */
veg_noleg_s11=vege_s11-leg_s11;
/* dairy and fish */
dairy_s11=sum(skim_s11, m1or2_s11, whole_s11, sherb_s11, icecr_s11, flyog_s11, plyog_s11, cream_s11, cotch_s11, crmch_s11, otch_s11);
fish_s11=sum(ctuna_s11, dkfsh_s11, bfsh_s11, ofish_s11, shrim_s11);
egg_s11=sum(egg_s11, eggom_s11);
/* poultry */
poul_s11=sum(chksa_s11, chwi_s11, chwo_s11, ctdog_s11 );
/* Total meat intake - servings added up*/
rpmeats_s11=sum(bmix_s11, pmain_s11, bmain_s11,  pmsan_s11, procm_s11, bacon_s11, bfdog_s11, hamb_s11, hambl_s11);
/* ssb*/
ssb_s11    = SUM(cola_s11, otsug_s11, punch_s11);
/*** SWEETS & DESERTS  ***/
swtsvg_s11= sum(mchoc_s11, dchoc_s11, cdyw_s11, cdywo_s11, coknf_s11, cokh_s11, cokr_s11, donut_s11, cake_s11, piehr_s11, sroll_s11, icecr_s11);
swtsnochoc_s11 = sum(cdyw_s11 ,  cdywo_s11,  coknf_s11,  cokr_s11 ,  cokh_s11 ,   donut_s11,  cake_s11 ,  piehr_s11 ,  sroll_s11 ,icecr_s11 );
run;

%nur13 (keep=wt13 db13 hbp13 chol13 mi13 str13 cabg13 ang13 
             dbfh13 hxdb13 mdchd13 fdchd13 mifh13 mdstr13 fdstr13 strfh13
             brcn13 ut13 ov13 colc13 mel13 ocan13 can13 
             cafh13 hxov13 hxbc13 hxmel13 hxcol13
             mamno13 mamsc13 mamsy13 mampt13 
             q13 colsc13 sigsc13 vbld13 diarr13 fecal13 hxclc13 rout13 bariu13 abdpn13 colrsnpt13 virt13 cpolp13 virtc13 );

            if hxdb13=2 then dbfh13=1;else dbfh13=0;
            if mdchd13=1 or fdchd13=1 then mifh13=1; else mifh13=0;
            if mdstr13=1 or fdstr13=1 then strfh13=1; else strfh13=0;
            if hxov13=2 or hxov13=2 or hxbc13=2 or hxmel13=2 or hxcol13=2 then cafh13=1;else cafh13=0;
            can13=0; if brcn13=1 or ut13=1 or ov13=1 or colc13=1 or mel13=1 or ocan13=1 then can13=1;
            run;


%nur15(keep=wt15 db15 hbp15  chol15  mi15  ang15 str15 cabg15  
            brcn15  ocan15 ut15 ov15 mel15 can15 
            cerbr15d mamno15 mamsc15 mamsy15 mampt15 
            q15 colsc15 sigsc15 vbld15 diarr15 fecal15 hxclc15 rout15 fecdna15 abdpn15 colrsnpt15 virt15 cpolp15 virtc15 );
            can15=0; if brcn15=1 or ocan15=1 or ut15=1 or ov15=1 or mel15=1 then can15=1;
            run;


%serv15 (keep= 
mchoc_s15 /*milk chocolate*/
dchoc_s15 /*dark chocolate*/

cdyw_s15 /*Candy bar*/
cdywo_s15/*Candy wo Chocolate*/
coknf_s15/*Cookies, Fat Free*/
cokr_s15/*Cookies, Commercial*/
cokh_s15/*Cookies, Home Baked*/
donut_s15/*Doughnuts*/
cake_s15 /*Cake, Homemade*/
piehr_s15/*Pie, Ready-Made*/
jam_s15  /*Jams, Jellies*/
pbut_s15 /*Peanut buter*/
ffpop_s15/*Popcorn: Fat Free*/
popc_s15/*Popcorn: Regular*/
sroll_s15/*Sweetroll Fat Free*/
brbar_s15/*Breakfast Bars*/
powrb_s15/*Energy Bars*/
pretz_s15/*Pretzels*/
slimf_s15/*Slimfast*/
ensur_s15/*Ensure*/
mdrfr_s15 /*Mixed Dried Fruit*/
pnut_s15/*Peanuts*/
wnut_s15/*Walnuts*/
onut_s15/*Other Nuts*/

pchip_s15 /*chips*/
icecr_s15 /*Ice cream*/
plyog_s15 /*DAIRY - Yogurt-plain, lo-carb*/
flyog_s15 /*Yogurt*/


rais_s15  prune_s15  ban_s15  cant_s15  avo_s15  appl_s15  oran_s15  grfr_s15  straw_s15  blueb_s15  peach_s15   
brocc_s15  cabb_s15  cauli_s15  bruss_s15  kale_s15 tom_s15  toj_s15  tosau_s15  ketch_s15
rcar_s15  ccar_s15  osqua_s15  yam_s15 cspin_s15  rspin_s15  rlett_s15  ilett_s15
tofu_s15  sbean_s15  peas_s15  bean_s15  soy_s15 corn_s15  mixv_s15  eggpl_s15  cel_s15  grpep_s15  oniog_s15  oniov_s15
cer_s15   dkbr_s15  popc_s15  oat_s15  ckcer_s15  brice_s15  oatbr_s15   wgerm_s15
bmix_s15  pmain_s15  bmain_s15   pmsan_s15  procm_s15  bacon_s15  bfdog_s15  hamb_s15  hambl_s15
cola_s15  otsug_s15  punch_s15
tea_s15 dtea_s15 coff_s15 dcaf_s15  rwine_s15 wwine_s15
/* potatoes */  fries_s15 pot_s15 pchip_s15
/* dairy */ skim_s15 m1or2_s15 whole_s15 sherb_s15 icecr_s15 flyog_s15 plyog_s15 cream_s15 cotch_s15 crmch_s15 otch_s15
/* eggs */ egg_s15 eggom_s15
/* fish */ ctuna_s15 dkfsh_s15 bfsh_s15 ofish_s15 shrim_s15
/* fruit juices */ aj_s15 oj_s15 ojca_s15  othj_s15  

chksa_s15 chwi_s15 chwo_s15 ctdog_s15
/* fv_s15 fruit_s15 veg_s15 whgrn_s15 rpmeats_s15 ssb_s15 cerbr_s15 wcer_s15 swtsvg_s15 */
);

data serv15;
merge serv15 nur15(keep=id cerbr15d);
by id;

if cerbr15d in (0, 4, 11, 12, 13, 15, 16, 17, 19, 20, 21, 22, 23, 25, 29, 30, 31, 33, 35, 37, 38, 42, 43, 45, 46, 47, 59, 60, 
61, 62, 63, 64, 65, 66, 68, 73, 74, 75, 76, 79,80, 81, 82, 83, 84, 85, 86, 88, 95, 96, 97, 98, 101, 109,
111, 112, 113, 114, 115, 119, 123, 124, 127, 128, 129, 131, 132, 133, 134, 146, 147, 148, 149, 150, 151, 152, 154, 155, 157,  158,
130, 164, 165, 166, 167, 170, 173, 176,  177, 182, 184, 185, 186, 188 )
  then rcer_s15=cer_s15;      
  else wcer_s15=cer_s15;  
  if rcer_s15=. then rcer_s15=0;
  if wcer_s15=. then wcer_s15=0;  

choco_s15=SUM(mchoc_s15, dchoc_s15);

/* take sum of caffeinated and decaffeinated tea & coffee */
tea_s15  = SUM(tea_s15, dtea_s15);
coffee_s15  = sum(coff_s15, dcaf_s15);
/* whole grain*/
whgrn_s15  = SUM(wcer_s15, dkbr_s15, popc_s15, oat_s15, ckcer_s15,  brice_s15, oatbr_s15,  wgerm_s15);
/* fruit and vegetables*/
fruit_s15  = SUM(rais_s15, prune_s15, ban_s15, cant_s15, avo_s15, appl_s15, oran_s15, grfr_s15, straw_s15, blueb_s15, peach_s15);
vege_s15    = SUM(brocc_s15, cabb_s15,   cauli_s15,  bruss_s15, kale_s15,   tom_s15,   toj_s15,   tosau_s15, ketch_s15,
                  rcar_s15,  ccar_s15,   osqua_s15,  yam_s15,  cspin_s15,  rspin_s15, rlett_s15, ilett_s15,
                  tofu_s15,  sbean_s15,  peas_s15,   bean_s15,  soy_s15,    corn_s15,  mixv_s15,  eggpl_s15, cel_s15,  grpep_s15,  oniog_s15,  oniov_s15);
fv_s15     = SUM(fruit_s15, vege_s15);
/* potatoes */
pota_s15=sum(fries_s15, pot_s15, pchip_s15);
/*legumes */
leg_s15=sum(sbean_s15, tofu_s15, bean_s15, peas_s15, soy_s15);
/* fruit juices */
fj_s15=sum(aj_s15, oj_s15, ojca_s15, othj_s15  );
/* vegetables minus potatoes and legumes */
veg_noleg_s15=vege_s15-leg_s15;
/* dairy and fish */
dairy_s15=sum(skim_s15, m1or2_s15, whole_s15, sherb_s15, icecr_s15, flyog_s15, plyog_s15, cream_s15, cotch_s15, crmch_s15, otch_s15);
fish_s15=sum(ctuna_s15, dkfsh_s15, bfsh_s15, ofish_s15, shrim_s15);
egg_s15=sum(egg_s15, eggom_s15);
/* poultry */
poul_s15=sum(chksa_s15, chwi_s15, chwo_s15, ctdog_s15 );
/* Total meat intake - servings added up*/
rpmeats_s15=sum(bmix_s15, pmain_s15, bmain_s15,  pmsan_s15, procm_s15, bacon_s15, bfdog_s15, hamb_s15, hambl_s15);
/* ssb*/
ssb_s15    = SUM(cola_s15, otsug_s15, punch_s15);
/*** SWEETS & DESERTS  ***/
swtsvg_s15= sum(mchoc_s15, dchoc_s15, cdyw_s15, cdywo_s15, coknf_s15, cokh_s15, cokr_s15, donut_s15, cake_s15, piehr_s15, sroll_s15, icecr_s15);
swtsnochoc_s15 = sum(cdyw_s15 ,  cdywo_s15,  coknf_s15,  cokr_s15 ,  cokh_s15 ,   donut_s15,  cake_s15 ,  piehr_s15 ,  sroll_s15 ,icecr_s15 );
run;

%nur17(keep=wt17 db17 hbp17  chol17  
            mi17  ang17  str17 cabg17 
            brcn17  ocan17 ut17 ov17 mel17 can17
            mamno17 mamsc17 mamsy17 mampt17 
            q17 colsc17 sigsc17 vbld17 diarr17 fecal17 hxclc17 rout17   abdpn17 colrsnpt17 virt17 cpolp17 virtc17 fecdna17 );
            can17=0; if brcn17=1 or ocan17=1 or ut17=1 or ov17=1 or mel17=1 then can17=1;
            run;

%nur19(keep=wt19 db19 hbp19  chol19  
            mi19  ang19  str19 cabg19 
            brcn19  ocan19 /*ut19*/ ov19 mel19 can19 cerbr19d
            mamno19 mamsc19 mamsy19 mampt19 
            q19 colsc19 sigsc19 vbld19 diarr19 fecal19 hxclc19 rout19 bariu19  abdpn19 colrsnpt19 virt19 cpolp19 virtc19 fecdna19 );
            can19=0; if brcn19=1 or ocan19=1 /*or ut19=1*/ or ov19=1 or mel19=1 then can19=1;
            run;


%serv19 (keep= 
mchoc_s19 /*milk chocolate*/
dchoc_s19 /*dark chocolate*/

cdyw_s19 /*Candy bar*/
cdywo_s19/*Candy wo Chocolate*/
/*coknf_s19*/ /*Cookies, Fat Free*/
cokr_s19/*Cookies, Commercial*/
cokh_s19/*Cookies, Home Baked*/
donut_s19/*Doughnuts*/
cake_s19 /*Cake, Homemade*/
piehr_s19/*Pie, Ready-Made*/
jam_s19  /*Jams, Jellies*/
pbut_s19 /*Peanut buter*/
ffpop_s19/*Popcorn: Fat Free*/
popc_s19/*Popcorn: Regular*/
sroll_s19/*Sweetroll Fat Free*/
brbar_s19/*Breakfast Bars*/
powrb_s19/*Energy Bars*/
pretz_s19/*Pretzels*/
slimf_s19/*Slimfast*/
ensur_s19/*Ensure*/
mdrfr_s19 /*Mixed Dried Fruit*/
pnut_s19/*Peanuts*/
wnut_s19/*Walnuts*/
onut_s19/*Other Nuts*/

pchip_s19 /*chips*/
icecr_s19 /*Ice cream*/
plyog_s19 /*DAIRY - Yogurt-plain, lo-carb*/
flyog_s19 /*Yogurt*/


rais_s19  prune_s19  ban_s19  cant_s19  avo_s19  appl_s19  oran_s19  grfr_s19  straw_s19  blueb_s19  peach_s19   
brocc_s19  cabb_s19  cauli_s19  bruss_s19  kale_s19 tom_s19  toj_s19  tosau_s19  ketch_s19
rcar_s19  ccar_s19  osqua_s19  yam_s19 cspin_s19  rspin_s19  rlett_s19  ilett_s19
tofu_s19  sbean_s19  peas_s19  bean_s19  soy_s19 corn_s19  mixv_s19  eggpl_s19    grpep_s19  oniog_s19  oniov_s19
cer_s19   dkbr_s19  popc_s19  oat_s19  ckcer_s19  brice_s19  oatbr_s19   wgerm_s19
bmix_s19  pmain_s19  bmain_s19   pmsan_s19  procm_s19  bacon_s19  bfdog_s19  hamb_s19  hambl_s19
cola_s19  otsug_s19  punch_s19
tea_s19 dtea_s19 coff_s19 dcaf_s19  rwine_s19 wwine_s19
/* potatoes */  fries_s19 pot_s19 pchip_s19
/* dairy */ skim_s19 m1or2_s19 whole_s19 sherb_s19 icecr_s19 flyog_s19 plyog_s19 cream_s19 cotch_s19 crmch_s19 otch_s19
/* eggs */ egg_s19 eggom_s19
/* fish */ ctuna_s19 dkfsh_s19 bfsh_s19 ofish_s19 shrim_s19
/* fruit juices */ aj_s19 oj_s19 ojca_s19  othj_s19  

chksa_s19 chwi_s19 chwo_s19 ctdog_s19
/* fv_s19 fruit_s19 veg_s19 whgrn_s19 rpmeats_s19 ssb_s19 cerbr_s19 wcer_s19 swtsvg_s19 */
);

data serv19;
merge serv19 nur19(keep=id cerbr19d);
by id;

if cerbr19d in (0, 4, 11, 12, 13, 19, 16, 17, 19, 20, 21, 22, 23, 25, 29, 30, 31, 33, 35, 37, 38, 42, 43, 45, 46, 47, 59, 60, 
61, 62, 63, 64, 65, 66, 68, 73, 74, 75, 76, 79,80, 81, 82, 83, 84, 85, 86, 88, 95, 96, 97, 98, 101, 109,
111, 112, 113, 114, 119, 119, 123, 124, 127, 128, 129, 131, 132, 133, 134, 146, 147, 148, 149, 190, 191, 192, 194, 195, 197,  198,
130, 164, 165, 166, 167, 170, 173, 176,  177, 182, 184, 185, 186, 188 )
  then rcer_s19=cer_s19;      
  else wcer_s19=cer_s19;  
  if rcer_s19=. then rcer_s19=0;
  if wcer_s19=. then wcer_s19=0;  

choco_s19=SUM(mchoc_s19, dchoc_s19);

/* take sum of caffeinated and decaffeinated tea & coffee */
tea_s19  = SUM(tea_s19, dtea_s19);
coffee_s19  = sum(coff_s19, dcaf_s19);
/* whole grain*/
whgrn_s19  = SUM(wcer_s19, dkbr_s19, popc_s19, oat_s19, ckcer_s19,  brice_s19, oatbr_s19,  wgerm_s19);
/* fruit and vegetables*/
fruit_s19  = SUM(rais_s19, prune_s19, ban_s19, cant_s19, avo_s19, appl_s19, oran_s19, grfr_s19, straw_s19, blueb_s19, peach_s19);
vege_s19    = SUM(brocc_s19, cabb_s19,   cauli_s19,  bruss_s19, kale_s19,   tom_s19,   toj_s19,   tosau_s19, ketch_s19,
                  rcar_s19,  ccar_s19,   osqua_s19,  yam_s19,  cspin_s19,  rspin_s19, rlett_s19, ilett_s19,
                  tofu_s19,  sbean_s19,  peas_s19,   bean_s19,  soy_s19,    corn_s19,  mixv_s19,  eggpl_s19,  grpep_s19,  oniog_s19,  oniov_s19);
fv_s19     = SUM(fruit_s19, vege_s19);
/* potatoes */
pota_s19=sum(fries_s19, pot_s19, pchip_s19);
/*legumes */
leg_s19=sum(sbean_s19, tofu_s19, bean_s19, peas_s19, soy_s19);
/* fruit juices */
fj_s19=sum(aj_s19, oj_s19, ojca_s19, othj_s19  );
/* vegetables minus potatoes and legumes */
veg_noleg_s19=vege_s19-leg_s19;
/* dairy and fish */
dairy_s19=sum(skim_s19, m1or2_s19, whole_s19, sherb_s19, icecr_s19, flyog_s19, plyog_s19, cream_s19, cotch_s19, crmch_s19, otch_s19);
fish_s19=sum(ctuna_s19, dkfsh_s19, bfsh_s19, ofish_s19, shrim_s19);
egg_s19=sum(egg_s19, eggom_s19);
/* poultry */
poul_s19=sum(chksa_s19, chwi_s19, chwo_s19, ctdog_s19 );
/* Total meat intake - servings added up*/
rpmeats_s19=sum(bmix_s19, pmain_s19, bmain_s19,  pmsan_s19, procm_s19, bacon_s19, bfdog_s19, hamb_s19, hambl_s19);
/* ssb*/
ssb_s19    = SUM(cola_s19, otsug_s19, punch_s19);
/*** SWEETS & DESERTS  ***/
swtsvg_s19= sum(mchoc_s19, dchoc_s19, cdyw_s19, cdywo_s19, cokh_s19, cokr_s19, donut_s19, cake_s19, piehr_s19, sroll_s19, icecr_s19);
swtsnochoc_s19 = sum(cdyw_s19 ,  cdywo_s19, cokr_s19 ,  cokh_s19 ,   donut_s19,  cake_s19 ,  piehr_s19 ,  sroll_s19 ,icecr_s19 );
run;


/* oc8911 latest file date 2017/4/28; Access on 1/2/2023 */
%oc8911(keep=nocu89 nocu91 nocu93 nocu95 nocu97 nocu99
             nocu01 nocu03 nocu05 nocu07 nocu09 nocu11);  /* OC no longer asked in questionnaires after 2011 */
/* multivitamin use status 1-current, 0-non, 9-unknown */
%supp8915 (keep=    mvitu89 mvitu91 mvitu93 mvitu95 mvitu97 
                    mvitu99 mvitu01 mvitu03 mvitu05 mvitu07 
                    mvitu09 mvitu11 mvitu13 mvitu15);

/**********READING IN  DATA ON ASPIRIN USE AND USE OF NON-STEROIDAL-INFLAMMATORY DRUGS (NSAID)**********/
%meds8919(keep=
aspu89  aspu91  aspu93  aspu95  aspu97  aspu99  aspu01  aspu03  aspu05  aspu07  aspu09  aspu11  aspu13 aspu15  aspu17 aspu19
nsaiu89 nsaiu91 nsaiu93 nsaiu95 nsaiu97 nsaiu99 nsaiu01 nsaiu03 nsaiu05 nsaiu07 nsaiu09 nsaiu11 nsaiu13 nsaiu15 nsaiu17 nsaiu19);


/*** Nutrient Intake from 1991-2015, in grams ***/
%n91_nts(keep=id alco91n calor91n na91n  nblnk91n
  carbo91n
  pot_carb91n/*Carbohydrate from potato*/
  leg_carb91n/*Carbohdyrate from legumes*/
  frtcb_njs91n/*Fruit carbohydrate minus juice sugar*/
  vegcb_npl91n/*Vegetable Carbohydrate minus potato and legume carbohydrate*/
  sumref_carb91n/*Sum of carbohydrate from potato, added sugar in foods, and refined grain, foods only*/
  sum_whgcarb91n/*Sum Carbohydrate From Wholegrain*/
  sumoth_crb_njs91n/*Sum of carbohydrate from fruit minus juice sugar, plus carb from vegetable, legume, and lact*/
  mill_grn91n/*Carbohydrate from milled whole grain*/
  int_grn91n/*Carbohydrate from intact whole grain*/
  ref_grn91n/*Carbohydrate from refined grain*/
  adds_cer91n/********Added sugar minus cereal fruit for Carb Quality Variables only********/
  adds_fds91n/********Added Sugar From Foods for Carb Quality Variables only********/

  tfat91n
  afat91n/*Animal Fat*/
  vfat91n/*Vegetable Fat*/
  dfat91n/*Dairy Fat*/
  sat91n/*Total Saturated Fat*/
  mon91n/*Total Monounsaturated Fat*/
  mufa_pl91n/*Plant Monounsaturated Fat*/
  mufa_a91n/*Animal Monounsaturated Fat*/
  poly91n/*Total Polyunsaturated Fat*/
  omega91n/*Omega 3 (f205+f226) gm, no alpha 18:3*/
  trnss91n/*Total Trans Fat*/

  prot91n
  aprot91n/*Animal Protein*/
  vprot91n/*Vegetable Protein*/
  dprot91n/*Dairy Protein*/
  );

data n91_nts;
set n91_nts;

array gram_cp carbo91n pot_carb91n leg_carb91n frtcb_njs91n vegcb_npl91n sumref_carb91n sum_whgcarb91n 
sumoth_crb_njs91n mill_grn91n int_grn91n ref_grn91n prot91n aprot91n vprot91n dprot91n adds_cer91n adds_fds91n;
array perc_cp carbo91p pot_carb91p leg_carb91p frtcb_njs91p vegcb_npl91p sumref_carb91p sum_whgcarb91p 
sumoth_crb_njs91p mill_grn91p int_grn91p ref_grn91p prot91p aprot91p vprot91p dprot91p adds_cer91p adds_fds91p;

array gram_f  tfat91n afat91n vfat91n dfat91n sat91n mon91n mufa_pl91n mufa_a91n poly91n omega91n trnss91n;
array perc_f  tfat91p afat91p vfat91p dfat91p sat91p mon91p mufa_pl91p mufa_a91p poly91p omega91p trnss91p;

do over gram_cp;
perc_cp=(4*gram_cp)/calor91n*100;
end;

do over gram_f;
perc_f=(9*gram_f)/calor91n*100;
end;

hfat91p=sum(mufa_pl91p,poly91p);      /*healthy fat as % of energy*/
ufat91p=sum(sat91p,mufa_a91p,trnss91p);   /*unhealthy fat as % of energy*/

hcarbo91p=sum(leg_carb91p,frtcb_njs91p,vegcb_npl91p,sum_whgcarb91p);   /*healthy carb as % of energy*/
run;

proc rank data=n91_nts groups=11 out=n91_nts;       * assign 11 groups based on these values (high prot=10, low prot = 0);
var    prot91p  aprot91p  vprot91p  tfat91p  afat91p  vfat91p  sum_whgcarb91p  sumoth_crb_njs91p  hfat91p  ufat91p;
ranks  prot91pd aprot91pd vprot91pd tfat91pd afat91pd vfat91pd sum_whgcarb91pd sumoth_crb_njs91pd hfat91pd ufat91pd;
run;

proc rank data=n91_nts descending groups=11 out=n91_nts;    * high carb = 0, low carb = 10;
var    carbo91p  sumref_carb91p  adds_cer91p  adds_fds91p  hcarbo91p;
ranks  carbo91pd sumref_carb91pd adds_cer91pd adds_fds91pd hcarbo91pd;
run;


data n91_nts;
set n91_nts;            
/* low carb score: more prot/fat (higher group ranking), less carb (higher group ranking)*/
lcs91=sum(carbo91pd,prot91pd,tfat91pd);       * low carb score = sum of total carb, prot, fat;
alcs91=sum(carbo91pd,aprot91pd,afat91pd);     * animal-based low carb score = sum or carb, animal protein, animal fat;
vlcs91=sum(carbo91pd,vprot91pd,vfat91pd);     * vege(plant-based) low carb score = sum of carb, vege protein, vege fat;

* healthy low carb score: more plant-based fat/protein (higher group ranking), less unhealthy carb (higher group ranking);
* unhealthy low carb score: more animal-based fat/protein (higher group ranking), less healthy carb (higher group ranking);
hlcs91=sum(sumref_carb91pd,vprot91pd,vfat91pd);   * healthy(high-quality) low carb score = sum of a reversed rank of unhealthful carb, vege protein, vege fat;
ulcs91=sum(hcarbo91pd     ,aprot91pd,afat91pd);   * unhealthy(low-quality) low carb score = sum of a reversed rank of healthy carb, animal protein, animal fat;

* carb quality: more healthy carb (whole grain, fruits=>higher group ranking), less unhealthy carb (higher group ranking);
* sugar quality ;

cqs91=sum(sum_whgcarb91pd,sumoth_crb_njs91pd,sumref_carb91pd);    
sug_cqs91=sum(adds_cer91pd,adds_fds91pd);
run;

proc rank data=n91_nts groups=11 out=n91_nts;
var    carbo91p  sumref_carb91p  hcarbo91p   prot91p  aprot91p  vprot91p   ;
ranks  carbo91pdd sumref_carb91pdd hcarbo91pdd  prot91pdd aprot91pdd vprot91pdd  ;
run;

proc rank data=n91_nts descending groups=11 out=n91_nts ;
var    tfat91p  afat91p  vfat91p  ;
ranks  tfat91pdd afat91pdd vfat91pdd ;
run;

data n91_nts;
set n91_nts;

lfs91=sum(carbo91pdd,prot91pdd,tfat91pdd);
alfs91=sum(carbo91pdd,aprot91pdd,vfat91pdd); /* updated: from tfat91pdd to vfat91pdd */
vlfs91=sum(carbo91pdd,vprot91pdd,afat91pdd); /* updated: from tfat91pdd to afat91pdd */

hlfs91=sum(hcarbo91pdd, vprot91pdd, afat91pdd);
ulfs91=sum(sumref_carb91pdd,aprot91pdd,vfat91pdd);
run;



%n95_nts(keep=id alco95n calor95n na95n
  carbo95n
  pot_carb95n/*Carbohydrate from potato*/
  leg_carb95n/*Carbohdyrate from legumes*/
  frtcb_njs95n/*Fruit carbohydrate minus juice sugar*/
  vegcb_npl95n/*Vegetable Carbohydrate minus potato and legume carbohydrate*/
  sumref_carb95n/*Sum of carbohydrate from potato, added sugar in foods, and refined grain, foods only*/
  sum_whgcarb95n/*Sum Carbohydrate From Wholegrain*/
  sumoth_crb_njs95n/*Sum of carbohydrate from fruit minus juice sugar, plus carb from vegetable, legume, and lact*/
  mill_grn95n/*Carbohydrate from milled whole grain*/
  int_grn95n/*Carbohydrate from intact whole grain*/
  ref_grn95n/*Carbohydrate from refined grain*/
  adds_cer95n/********Added sugar minus cereal fruit for Carb Quality Variables only********/
  adds_fds95n/********Added Sugar From Foods for Carb Quality Variables only********/

  tfat95n
  afat95n/*Animal Fat*/
  vfat95n/*Vegetable Fat*/
  dfat95n/*Dairy Fat*/
  sat95n/*Total Saturated Fat*/
  mon95n/*Total Monounsaturated Fat*/
  mufa_pl95n/*Plant Monounsaturated Fat*/
  mufa_a95n/*Animal Monounsaturated Fat*/
  poly95n/*Total Polyunsaturated Fat*/
  omega95n/*Omega 3 (f205+f226) gm, no alpha 18:3*/
  trnss95n/*Total Trans Fat*/

  prot95n
  aprot95n/*Animal Protein*/
  vprot95n/*Vegetable Protein*/
  dprot95n/*Dairy Protein*/
  );

data n95_nts;
set n95_nts;

array gram_cp carbo95n pot_carb95n leg_carb95n frtcb_njs95n vegcb_npl95n sumref_carb95n sum_whgcarb95n sumoth_crb_njs95n mill_grn95n int_grn95n ref_grn95n prot95n aprot95n vprot95n dprot95n adds_cer95n adds_fds95n;
array perc_cp carbo95p pot_carb95p leg_carb95p frtcb_njs95p vegcb_npl95p sumref_carb95p sum_whgcarb95p sumoth_crb_njs95p mill_grn95p int_grn95p ref_grn95p prot95p aprot95p vprot95p dprot95p adds_cer95p adds_fds95p;

array gram_f  tfat95n afat95n vfat95n dfat95n sat95n mon95n mufa_pl95n mufa_a95n poly95n omega95n trnss95n;
array perc_f  tfat95p afat95p vfat95p dfat95p sat95p mon95p mufa_pl95p mufa_a95p poly95p omega95p trnss95p;

do over gram_cp;
perc_cp=(4*gram_cp)/calor95n*100;
end;

do over gram_f;
perc_f=(9*gram_f)/calor95n*100;
end;

hfat95p=sum(mufa_pl95p,poly95p);
ufat95p=sum(sat95p,mufa_a95p,trnss95p);

hcarbo95p=sum(leg_carb95p,frtcb_njs95p,vegcb_npl95p,sum_whgcarb95p);

run;

proc rank data=n95_nts groups=11 out=n95_nts;
var    prot95p  aprot95p  vprot95p  tfat95p  afat95p  vfat95p  sum_whgcarb95p  sumoth_crb_njs95p  hfat95p  ufat95p;
ranks  prot95pd aprot95pd vprot95pd tfat95pd afat95pd vfat95pd sum_whgcarb95pd sumoth_crb_njs95pd hfat95pd ufat95pd;
run;

proc rank data=n95_nts descending groups=11 out=n95_nts;
var    carbo95p  sumref_carb95p  adds_cer95p  adds_fds95p  hcarbo95p;
ranks  carbo95pd sumref_carb95pd adds_cer95pd adds_fds95pd hcarbo95pd;
run;

data n95_nts;
set n95_nts;

lcs95=sum(carbo95pd,prot95pd,tfat95pd);
alcs95=sum(carbo95pd,aprot95pd,afat95pd);
vlcs95=sum(carbo95pd,vprot95pd,vfat95pd);

hlcs95=sum(sumref_carb95pd,vprot95pd,vfat95pd);
ulcs95=sum(hcarbo95pd     ,aprot95pd,afat95pd);

cqs95=sum(sum_whgcarb95pd,sumoth_crb_njs95pd,sumref_carb95pd);
sug_cqs95=sum(adds_cer95pd,adds_fds95pd);
run;

proc rank data=n95_nts groups=11 out=n95_nts;
var    carbo95p  sumref_carb95p  hcarbo95p   prot95p  aprot95p  vprot95p   ;
ranks  carbo95pdd sumref_carb95pdd hcarbo95pdd  prot95pdd aprot95pdd vprot95pdd  ;
run;

proc rank data=n95_nts descending groups=11 out=n95_nts ;
var    tfat95p  afat95p  vfat95p  ;
ranks  tfat95pdd afat95pdd vfat95pdd ;
run;

data n95_nts;
set n95_nts;

lfs95=sum(carbo95pdd,prot95pdd,tfat95pdd);
alfs95=sum(carbo95pdd,aprot95pdd,vfat95pdd);
vlfs95=sum(carbo95pdd,vprot95pdd,afat95pdd);

hlfs95=sum(hcarbo95pdd, vprot95pdd, afat95pdd);
ulfs95=sum(sumref_carb95pdd,aprot95pdd,vfat95pdd);
run;

%n99_nts(keep=id alco99n calor99n na99n
  carbo99n
  pot_carb99n/*Carbohydrate from potato*/
  leg_carb99n/*Carbohdyrate from legumes*/
  frtcb_njs99n/*Fruit carbohydrate minus juice sugar*/
  vegcb_npl99n/*Vegetable Carbohydrate minus potato and legume carbohydrate*/
  sumref_carb99n/*Sum of carbohydrate from potato, added sugar in foods, and refined grain, foods only*/
  sum_whgcarb99n/*Sum Carbohydrate From Wholegrain*/
  sumoth_crb_njs99n/*Sum of carbohydrate from fruit minus juice sugar, plus carb from vegetable, legume, and lact*/
  mill_grn99n/*Carbohydrate from milled whole grain*/
  int_grn99n/*Carbohydrate from intact whole grain*/
  ref_grn99n/*Carbohydrate from refined grain*/
  adds_cer99n/********Added sugar minus cereal fruit for Carb Quality Variables only********/
  adds_fds99n/********Added Sugar From Foods for Carb Quality Variables only********/

  tfat99n
  afat99n/*Animal Fat*/
  vfat99n/*Vegetable Fat*/
  dfat99n/*Dairy Fat*/
  sat99n/*Total Saturated Fat*/
  mon99n/*Total Monounsaturated Fat*/
  mufa_pl99n/*Plant Monounsaturated Fat*/
  mufa_a99n/*Animal Monounsaturated Fat*/
  poly99n/*Total Polyunsaturated Fat*/
  omega99n/*Omega 3 (f205+f226) gm, no alpha 18:3*/
  trn0099n/*Total Trans Fat*/

  prot99n
  aprot99n/*Animal Protein*/
  vprot99n/*Vegetable Protein*/
  dprot99n/*Dairy Protein*/
  );

data n99_nts;
set n99_nts;

array gram_cp carbo99n pot_carb99n leg_carb99n frtcb_njs99n vegcb_npl99n sumref_carb99n sum_whgcarb99n sumoth_crb_njs99n mill_grn99n int_grn99n ref_grn99n prot99n aprot99n vprot99n dprot99n adds_cer99n adds_fds99n;
array perc_cp carbo99p pot_carb99p leg_carb99p frtcb_njs99p vegcb_npl99p sumref_carb99p sum_whgcarb99p sumoth_crb_njs99p mill_grn99p int_grn99p ref_grn99p prot99p aprot99p vprot99p dprot99p adds_cer99p adds_fds99p;

array gram_f  tfat99n afat99n vfat99n dfat99n sat99n mon99n mufa_pl99n mufa_a99n poly99n omega99n trn0099n;
array perc_f  tfat99p afat99p vfat99p dfat99p sat99p mon99p mufa_pl99p mufa_a99p poly99p omega99p trnss99p;

do over gram_cp;
perc_cp=(4*gram_cp)/calor99n*100;
end;

do over gram_f;
perc_f=(9*gram_f)/calor99n*100;
end;

hfat99p=sum(mufa_pl99p,poly99p);
ufat99p=sum(sat99p,mufa_a99p,trnss99p);

hcarbo99p=sum(leg_carb99p,frtcb_njs99p,vegcb_npl99p,sum_whgcarb99p);
run;

proc rank data=n99_nts groups=11 out=n99_nts;
var    prot99p  aprot99p  vprot99p  tfat99p  afat99p  vfat99p  sum_whgcarb99p  sumoth_crb_njs99p  hfat99p  ufat99p;
ranks  prot99pd aprot99pd vprot99pd tfat99pd afat99pd vfat99pd sum_whgcarb99pd sumoth_crb_njs99pd hfat99pd ufat99pd;
run;

proc rank data=n99_nts descending groups=11 out=n99_nts;
var    carbo99p  sumref_carb99p  adds_cer99p  adds_fds99p  hcarbo99p;
ranks  carbo99pd sumref_carb99pd adds_cer99pd adds_fds99pd hcarbo99pd;
run;

data n99_nts;
set n99_nts;

lcs99=sum(carbo99pd,prot99pd,tfat99pd);
alcs99=sum(carbo99pd,aprot99pd,afat99pd);
vlcs99=sum(carbo99pd,vprot99pd,vfat99pd);

hlcs99=sum(sumref_carb99pd,vprot99pd,vfat99pd);
ulcs99=sum(hcarbo99pd     ,aprot99pd,afat99pd);

cqs99=sum(sum_whgcarb99pd,sumoth_crb_njs99pd,sumref_carb99pd);
sug_cqs99=sum(adds_cer99pd,adds_fds99pd);
run;

proc rank data=n99_nts groups=11 out=n99_nts;
var    carbo99p  sumref_carb99p  hcarbo99p   prot99p  aprot99p  vprot99p   ;
ranks  carbo99pdd sumref_carb99pdd hcarbo99pdd  prot99pdd aprot99pdd vprot99pdd  ;
run;

proc rank data=n99_nts descending groups=11 out=n99_nts ;
var    tfat99p  afat99p  vfat99p  ;
ranks  tfat99pdd afat99pdd vfat99pdd ;
run;

data n99_nts;
set n99_nts;

lfs99=sum(carbo99pdd,prot99pdd,tfat99pdd);
alfs99=sum(carbo99pdd,aprot99pdd,vfat99pdd);
vlfs99=sum(carbo99pdd,vprot99pdd,afat99pdd);

hlfs99=sum(hcarbo99pdd, vprot99pdd, afat99pdd);
ulfs99=sum(sumref_carb99pdd,aprot99pdd,vfat99pdd);
run;

%n03_nts(keep=id alco03n calor03n na03n
  carbo03n
  pot_carb03n/*Carbohydrate from potato*/
  leg_carb03n/*Carbohdyrate from legumes*/
  frtcb_njs03n/*Fruit carbohydrate minus juice sugar*/
  vegcb_npl03n/*Vegetable Carbohydrate minus potato and legume carbohydrate*/
  sumref_carb03n/*Sum of carbohydrate from potato, added sugar in foods, and refined grain, foods only*/
  sum_whgcarb03n/*Sum Carbohydrate From Wholegrain*/
  sumoth_crb_njs03n/*Sum of carbohydrate from fruit minus juice sugar, plus carb from vegetable, legume, and lact*/
  mill_grn03n/*Carbohydrate from milled whole grain*/
  int_grn03n/*Carbohydrate from intact whole grain*/
  ref_grn03n/*Carbohydrate from refined grain*/
  adds_cer03n/********Added sugar minus cereal fruit for Carb Quality Variables only********/
  adds_fds03n/********Added Sugar From Foods for Carb Quality Variables only********/

  tfat03n
  afat03n/*Animal Fat*/
  vfat03n/*Vegetable Fat*/
  dfat03n/*Dairy Fat*/
  sat03n/*Total Saturated Fat*/
  mon03n/*Total Monounsaturated Fat*/
  mufa_pl03n/*Plant Monounsaturated Fat*/
  mufa_a03n/*Animal Monounsaturated Fat*/
  poly03n/*Total Polyunsaturated Fat*/
  omega03n/*Omega 3 (f205+f226) gm, no alpha 18:3*/
  trn0203n/*Total Trans Fat*/

  prot03n
  aprot03n/*Animal Protein*/
  vprot03n/*Vegetable Protein*/
  dprot03n/*Dairy Protein*/
  );

data n03_nts;
set n03_nts;

array gram_cp carbo03n pot_carb03n leg_carb03n frtcb_njs03n vegcb_npl03n sumref_carb03n sum_whgcarb03n sumoth_crb_njs03n mill_grn03n int_grn03n ref_grn03n prot03n aprot03n vprot03n dprot03n adds_cer03n adds_fds03n;
array perc_cp carbo03p pot_carb03p leg_carb03p frtcb_njs03p vegcb_npl03p sumref_carb03p sum_whgcarb03p sumoth_crb_njs03p mill_grn03p int_grn03p ref_grn03p prot03p aprot03p vprot03p dprot03p adds_cer03p adds_fds03p;

array gram_f  tfat03n afat03n vfat03n dfat03n sat03n mon03n mufa_pl03n mufa_a03n poly03n omega03n trn0203n;
array perc_f  tfat03p afat03p vfat03p dfat03p sat03p mon03p mufa_pl03p mufa_a03p poly03p omega03p trnss03p;

do over gram_cp;
perc_cp=(4*gram_cp)/calor03n*100;
end;

do over gram_f;
perc_f=(9*gram_f)/calor03n*100;
end;

hfat03p=sum(mufa_pl03p,poly03p);
ufat03p=sum(sat03p,mufa_a03p,trnss03p);

hcarbo03p=sum(leg_carb03p,frtcb_njs03p,vegcb_npl03p,sum_whgcarb03p);
run;

proc rank data=n03_nts groups=11 out=n03_nts;
var    prot03p  aprot03p  vprot03p  tfat03p  afat03p  vfat03p  sum_whgcarb03p  sumoth_crb_njs03p  hfat03p  ufat03p;
ranks  prot03pd aprot03pd vprot03pd tfat03pd afat03pd vfat03pd sum_whgcarb03pd sumoth_crb_njs03pd hfat03pd ufat03pd;
run;

proc rank data=n03_nts descending groups=11 out=n03_nts;
var    carbo03p  sumref_carb03p  adds_cer03p  adds_fds03p  hcarbo03p;
ranks  carbo03pd sumref_carb03pd adds_cer03pd adds_fds03pd hcarbo03pd;
run;


data n03_nts;
set n03_nts;

lcs03=sum(carbo03pd,prot03pd,tfat03pd);
alcs03=sum(carbo03pd,aprot03pd,afat03pd);
vlcs03=sum(carbo03pd,vprot03pd,vfat03pd);

hlcs03=sum(sumref_carb03pd,vprot03pd,vfat03pd);
ulcs03=sum(hcarbo03pd     ,aprot03pd,afat03pd);

cqs03=sum(sum_whgcarb03pd,sumoth_crb_njs03pd,sumref_carb03pd);
sug_cqs03=sum(adds_cer03pd,adds_fds03pd);
run;

proc rank data=n03_nts groups=11 out=n03_nts;
var    carbo03p  sumref_carb03p  hcarbo03p   prot03p  aprot03p  vprot03p   ;
ranks  carbo03pdd sumref_carb03pdd hcarbo03pdd  prot03pdd aprot03pdd vprot03pdd  ;
run;

proc rank data=n03_nts descending groups=11 out=n03_nts ;
var    tfat03p  afat03p  vfat03p  ;
ranks  tfat03pdd afat03pdd vfat03pdd ;
run;

data n03_nts;
set n03_nts;

lfs03=sum(carbo03pdd,prot03pdd,tfat03pdd);
alfs03=sum(carbo03pdd,aprot03pdd,vfat03pdd);
vlfs03=sum(carbo03pdd,vprot03pdd,afat03pdd);

hlfs03=sum(hcarbo03pdd, vprot03pdd, afat03pdd);
ulfs03=sum(sumref_carb03pdd,aprot03pdd,vfat03pdd);
run;

%n07_nts(keep=id alco07n calor07n na07n
  carbo07n
  pot_carb07n/*Carbohydrate from potato*/
  leg_carb07n/*Carbohdyrate from legumes*/
  frtcb_njs07n/*Fruit carbohydrate minus juice sugar*/
  vegcb_npl07n/*Vegetable Carbohydrate minus potato and legume carbohydrate*/
  sumref_carb07n/*Sum of carbohydrate from potato, added sugar in foods, and refined grain, foods only*/
  sum_whgcarb07n/*Sum Carbohydrate From Wholegrain*/
  sumoth_crb_njs07n/*Sum of carbohydrate from fruit minus juice sugar, plus carb from vegetable, legume, and lact*/
  mill_grn07n/*Carbohydrate from milled whole grain*/
  int_grn07n/*Carbohydrate from intact whole grain*/
  ref_grn07n/*Carbohydrate from refined grain*/
  adds_cer07n/********Added sugar minus cereal fruit for Carb Quality Variables only********/
  adds_fds07n/********Added Sugar From Foods for Carb Quality Variables only********/

  tfat07n
  afat07n/*Animal Fat*/
  vfat07n/*Vegetable Fat*/
  dfat07n/*Dairy Fat*/
  sat07n/*Total Saturated Fat*/
  mon07n/*Total Monounsaturated Fat*/
  mufa_pl07n/*Plant Monounsaturated Fat*/
  mufa_a07n/*Animal Monounsaturated Fat*/
  poly07n/*Total Polyunsaturated Fat*/
  omega07n/*Omega 3 (f205+f226) gm, no alpha 18:3*/
  trn0707n/*Total Trans Fat*/

  prot07n
  aprot07n/*Animal Protein*/
  vprot07n/*Vegetable Protein*/
  dprot07n/*Dairy Protein*/
  );

data n07_nts;
set n07_nts;

array gram_cp carbo07n pot_carb07n leg_carb07n frtcb_njs07n vegcb_npl07n sumref_carb07n sum_whgcarb07n sumoth_crb_njs07n mill_grn07n int_grn07n ref_grn07n prot07n aprot07n vprot07n dprot07n adds_cer07n adds_fds07n;
array perc_cp carbo07p pot_carb07p leg_carb07p frtcb_njs07p vegcb_npl07p sumref_carb07p sum_whgcarb07p sumoth_crb_njs07p mill_grn07p int_grn07p ref_grn07p prot07p aprot07p vprot07p dprot07p adds_cer07p adds_fds07p;

array gram_f  tfat07n afat07n vfat07n dfat07n sat07n mon07n mufa_pl07n mufa_a07n poly07n omega07n trn0707n;
array perc_f  tfat07p afat07p vfat07p dfat07p sat07p mon07p mufa_pl07p mufa_a07p poly07p omega07p trnss07p;

do over gram_cp;
perc_cp=(4*gram_cp)/calor07n*100;
end;

do over gram_f;
perc_f=(9*gram_f)/calor07n*100;
end;

hfat07p=sum(mufa_pl07p,poly07p);
ufat07p=sum(sat07p,mufa_a07p,trnss07p);

hcarbo07p=sum(leg_carb07p,frtcb_njs07p,vegcb_npl07p,sum_whgcarb07p);
run;

proc rank data=n07_nts groups=11 out=n07_nts;
var    prot07p  aprot07p  vprot07p  tfat07p  afat07p  vfat07p  sum_whgcarb07p  sumoth_crb_njs07p  hfat07p  ufat07p;
ranks  prot07pd aprot07pd vprot07pd tfat07pd afat07pd vfat07pd sum_whgcarb07pd sumoth_crb_njs07pd hfat07pd ufat07pd;
run;

proc rank data=n07_nts descending groups=11 out=n07_nts;
var    carbo07p  sumref_carb07p  adds_cer07p  adds_fds07p  hcarbo07p;
ranks  carbo07pd sumref_carb07pd adds_cer07pd adds_fds07pd hcarbo07pd;
run;


data n07_nts;
set n07_nts;

lcs07=sum(carbo07pd,prot07pd,tfat07pd);
alcs07=sum(carbo07pd,aprot07pd,afat07pd);
vlcs07=sum(carbo07pd,vprot07pd,vfat07pd);

hlcs07=sum(sumref_carb07pd,vprot07pd,vfat07pd);
ulcs07=sum(hcarbo07pd     ,aprot07pd,afat07pd);

cqs07=sum(sum_whgcarb07pd,sumoth_crb_njs07pd,sumref_carb07pd);
sug_cqs07=sum(adds_cer07pd,adds_fds07pd);
run;

proc rank data=n07_nts groups=11 out=n07_nts;
var    carbo07p  sumref_carb07p  hcarbo07p   prot07p  aprot07p  vprot07p   ;
ranks  carbo07pdd sumref_carb07pdd hcarbo07pdd  prot07pdd aprot07pdd vprot07pdd  ;
run;

proc rank data=n07_nts descending groups=11 out=n07_nts ;
var    tfat07p  afat07p  vfat07p  ;
ranks  tfat07pdd afat07pdd vfat07pdd ;
run;

data n07_nts;
set n07_nts;

lfs07=sum(carbo07pdd,prot07pdd,tfat07pdd);
alfs07=sum(carbo07pdd,aprot07pdd,vfat07pdd);
vlfs07=sum(carbo07pdd,vprot07pdd,afat07pdd);

hlfs07=sum(hcarbo07pdd, vprot07pdd, afat07pdd);
ulfs07=sum(sumref_carb07pdd,aprot07pdd,vfat07pdd);
run;

%n11_nts(keep=id alco11n calor11n na11n
  carbo11n
  pot_carb11n/*Carbohydrate from potato*/
  leg_carb11n/*Carbohdyrate from legumes*/
  frtcb_njs11n/*Fruit carbohydrate minus juice sugar*/
  vegcb_npl11n/*Vegetable Carbohydrate minus potato and legume carbohydrate*/
  sumref_carb11n/*Sum of carbohydrate from potato, added sugar in foods, and refined grain, foods only*/
  sum_whgcarb11n/*Sum Carbohydrate From Wholegrain*/
  sumoth_crb_njs11n/*Sum of carbohydrate from fruit minus juice sugar, plus carb from vegetable, legume, and lact*/
  mill_grn11n/*Carbohydrate from milled whole grain*/
  int_grn11n/*Carbohydrate from intact whole grain*/
  ref_grn11n/*Carbohydrate from refined grain*/
  adds_cer11n/********Added sugar minus cereal fruit for Carb Quality Variables only********/
  adds_fds11n/********Added Sugar From Foods for Carb Quality Variables only********/

  tfat11n
  afat11n/*Animal Fat*/
  vfat11n/*Vegetable Fat*/
  dfat11n/*Dairy Fat*/
  sat11n/*Total Saturated Fat*/
  mon11n/*Total Monounsaturated Fat*/
  mufa_pl11n/*Plant Monounsaturated Fat*/
  mufa_a11n/*Animal Monounsaturated Fat*/
  poly11n/*Total Polyunsaturated Fat*/
  omega11n/*Omega 3 (f205+f226) gm, no alpha 18:3*/
  trn1111n/*Total Trans Fat*/

  prot11n
  aprot11n/*Animal Protein*/
  vprot11n/*Vegetable Protein*/
  dprot11n/*Dairy Protein*/
  );

data n11_nts;
set n11_nts;

array gram_cp carbo11n pot_carb11n leg_carb11n frtcb_njs11n vegcb_npl11n sumref_carb11n sum_whgcarb11n sumoth_crb_njs11n mill_grn11n int_grn11n ref_grn11n prot11n aprot11n vprot11n dprot11n adds_cer11n adds_fds11n;
array perc_cp carbo11p pot_carb11p leg_carb11p frtcb_njs11p vegcb_npl11p sumref_carb11p sum_whgcarb11p sumoth_crb_njs11p mill_grn11p int_grn11p ref_grn11p prot11p aprot11p vprot11p dprot11p adds_cer11p adds_fds11p;

array gram_f  tfat11n afat11n vfat11n dfat11n sat11n mon11n mufa_pl11n mufa_a11n poly11n omega11n trn1111n;
array perc_f  tfat11p afat11p vfat11p dfat11p sat11p mon11p mufa_pl11p mufa_a11p poly11p omega11p trnss11p;

do over gram_cp;
perc_cp=(4*gram_cp)/calor11n*100;
end;

do over gram_f;
perc_f=(9*gram_f)/calor11n*100;
end;

hfat11p=sum(mufa_pl11p,poly11p);
ufat11p=sum(sat11p,mufa_a11p,trnss11p);

hcarbo11p=sum(leg_carb11p,frtcb_njs11p,vegcb_npl11p,sum_whgcarb11p);
run;

proc rank data=n11_nts groups=11 out=n11_nts;
var    prot11p  aprot11p  vprot11p  tfat11p  afat11p  vfat11p  sum_whgcarb11p  sumoth_crb_njs11p  hfat11p  ufat11p;
ranks  prot11pd aprot11pd vprot11pd tfat11pd afat11pd vfat11pd sum_whgcarb11pd sumoth_crb_njs11pd hfat11pd ufat11pd;
run;

proc rank data=n11_nts descending groups=11 out=n11_nts;
var    carbo11p  sumref_carb11p  adds_cer11p  adds_fds11p  hcarbo11p;
ranks  carbo11pd sumref_carb11pd adds_cer11pd adds_fds11pd hcarbo11pd;
run;


data n11_nts;
set n11_nts;

lcs11=sum(carbo11pd,prot11pd,tfat11pd);
alcs11=sum(carbo11pd,aprot11pd,afat11pd);
vlcs11=sum(carbo11pd,vprot11pd,vfat11pd);

hlcs11=sum(sumref_carb11pd,vprot11pd,vfat11pd);
ulcs11=sum(hcarbo11pd     ,aprot11pd,afat11pd);

cqs11=sum(sum_whgcarb11pd,sumoth_crb_njs11pd,sumref_carb11pd);
sug_cqs11=sum(adds_cer11pd,adds_fds11pd);
run;

proc rank data=n11_nts groups=11 out=n11_nts;
var    carbo11p  sumref_carb11p  hcarbo11p   prot11p  aprot11p  vprot11p   ;
ranks  carbo11pdd sumref_carb11pdd hcarbo11pdd  prot11pdd aprot11pdd vprot11pdd  ;
run;

proc rank data=n11_nts descending groups=11 out=n11_nts ;
var    tfat11p  afat11p  vfat11p  ;
ranks  tfat11pdd afat11pdd vfat11pdd ;
run;

data n11_nts;
set n11_nts;

lfs11=sum(carbo11pdd,prot11pdd,tfat11pdd);
alfs11=sum(carbo11pdd,aprot11pdd,vfat11pdd);
vlfs11=sum(carbo11pdd,vprot11pdd,afat11pdd);
hlfs11=sum(hcarbo11pdd, vprot11pdd, afat11pdd);
ulfs11=sum(sumref_carb11pdd,aprot11pdd,vfat11pdd);
run;

%n15_nts(keep=id alco15n calor15n na15n
  carbo15n
  pot_carb15n/*Carbohydrate from potato*/
  leg_carb15n/*Carbohdyrate from legumes*/
  frtcb_njs15n/*Fruit carbohydrate minus juice sugar*/
  vegcb_npl15n/*Vegetable Carbohydrate minus potato and legume carbohydrate*/
  sumref_carb15n/*Sum of carbohydrate from potato, added sugar in foods, and refined grain, foods only*/
  sum_whgcarb15n/*Sum Carbohydrate From Wholegrain*/
  sumoth_crb_njs15n/*Sum of carbohydrate from fruit minus juice sugar, plus carb from vegetable, legume, and lact*/
  mill_grn15n/*Carbohydrate from milled whole grain*/
  int_grn15n/*Carbohydrate from intact whole grain*/
  ref_grn15n/*Carbohydrate from refined grain*/
  adds_cer15n/********Added sugar minus cereal fruit for Carb Quality Variables only********/
  adds_fds15n/********Added Sugar From Foods for Carb Quality Variables only********/

  tfat15n
  afat15n/*Animal Fat*/
  vfat15n/*Vegetable Fat*/
  dfat15n/*Dairy Fat*/
  sat15n/*Total Saturated Fat*/
  mon15n/*Total Monounsaturated Fat*/
  mufa_pl15n/*Plant Monounsaturated Fat*/
  mufa_a15n/*Animal Monounsaturated Fat*/
  poly15n/*Total Polyunsaturated Fat*/
  omega15n/*Omega 3 (f205+f226) gm, no alpha 18:3*/
  trn1115n/*Total Trans Fat*/

  prot15n
  aprot15n/*Animal Protein*/
  vprot15n/*Vegetable Protein*/
  dprot15n/*Dairy Protein*/
  );

data n15_nts;
set n15_nts;

array gram_cp carbo15n pot_carb15n leg_carb15n frtcb_njs15n vegcb_npl15n sumref_carb15n sum_whgcarb15n sumoth_crb_njs15n mill_grn15n int_grn15n ref_grn15n prot15n aprot15n vprot15n dprot15n adds_cer15n adds_fds15n;
array perc_cp carbo15p pot_carb15p leg_carb15p frtcb_njs15p vegcb_npl15p sumref_carb15p sum_whgcarb15p sumoth_crb_njs15p mill_grn15p int_grn15p ref_grn15p prot15p aprot15p vprot15p dprot15p adds_cer15p adds_fds15p;

array gram_f  tfat15n afat15n vfat15n dfat15n sat15n mon15n mufa_pl15n mufa_a15n poly15n omega15n trn1115n;
array perc_f  tfat15p afat15p vfat15p dfat15p sat15p mon15p mufa_pl15p mufa_a15p poly15p omega15p trnss15p;

do over gram_cp;
perc_cp=(4*gram_cp)/calor15n*100;
end;

do over gram_f;
perc_f=(9*gram_f)/calor15n*100;
end;

hfat15p=sum(mufa_pl15p,poly15p);
ufat15p=sum(sat15p,mufa_a15p,trnss15p);

hcarbo15p=sum(leg_carb15p,frtcb_njs15p,vegcb_npl15p,sum_whgcarb15p);
run;

proc rank data=n15_nts groups=11 out=n15_nts;
var    prot15p  aprot15p  vprot15p  tfat15p  afat15p  vfat15p  sum_whgcarb15p  sumoth_crb_njs15p  hfat15p  ufat15p;
ranks  prot15pd aprot15pd vprot15pd tfat15pd afat15pd vfat15pd sum_whgcarb15pd sumoth_crb_njs15pd hfat15pd ufat15pd;
run;

proc rank data=n15_nts descending groups=11 out=n15_nts;
var    carbo15p  sumref_carb15p  adds_cer15p  adds_fds15p  hcarbo15p;
ranks  carbo15pd sumref_carb15pd adds_cer15pd adds_fds15pd hcarbo15pd;
run;

data n15_nts;
set n15_nts;

lcs15=sum(carbo15pd,prot15pd,tfat15pd);
alcs15=sum(carbo15pd,aprot15pd,afat15pd);
vlcs15=sum(carbo15pd,vprot15pd,vfat15pd);

hlcs15=sum(sumref_carb15pd,vprot15pd,vfat15pd);
ulcs15=sum(hcarbo15pd     ,aprot15pd,afat15pd);

cqs15=sum(sum_whgcarb15pd,sumoth_crb_njs15pd,sumref_carb15pd);
sug_cqs15=sum(adds_cer15pd,adds_fds15pd);
run;

proc rank data=n15_nts groups=11 out=n15_nts;
var    carbo15p  sumref_carb15p  hcarbo15p   prot15p  aprot15p  vprot15p   ;
ranks  carbo15pdd sumref_carb15pdd hcarbo15pdd  prot15pdd aprot15pdd vprot15pdd  ;
run;

proc rank data=n15_nts descending groups=11 out=n15_nts ;
var    tfat15p  afat15p  vfat15p  ;
ranks  tfat15pdd afat15pdd vfat15pdd ;
run;

data n15_nts;
set n15_nts;

lfs15=sum(carbo15pdd,prot15pdd,tfat15pdd);
alfs15=sum(carbo15pdd,aprot15pdd,vfat15pdd);
vlfs15=sum(carbo15pdd,vprot15pdd,afat15pdd);
hlfs15=sum(hcarbo15pdd, vprot15pdd, afat15pdd);
ulfs15=sum(sumref_carb15pdd,aprot15pdd,vfat15pdd);
run;


%n19_nts(keep=id alco19n calor19n na19n
  carbo19n
  pot_carb19n/*Carbohydrate from potato*/
  leg_carb19n/*Carbohdyrate from legumes*/
  frtcb_njs19n/*Fruit carbohydrate minus juice sugar*/
  vegcb_npl19n/*Vegetable Carbohydrate minus potato and legume carbohydrate*/
  sumref_carb19n/*Sum of carbohydrate from potato, added sugar in foods, and refined grain, foods only*/
  sum_whgcarb19n/*Sum Carbohydrate From Wholegrain*/
  sumoth_crb_njs19n/*Sum of carbohydrate from fruit minus juice sugar, plus carb from vegetable, legume, and lact*/
  mill_grn19n/*Carbohydrate from milled whole grain*/
  int_grn19n/*Carbohydrate from intact whole grain*/
  ref_grn19n/*Carbohydrate from refined grain*/
  adds_cer19n/********Added sugar minus cereal fruit for Carb Quality Variables only********/
  adds_fds19n/********Added Sugar From Foods for Carb Quality Variables only********/

  tfat19n
  afat19n/*Animal Fat*/
  vfat19n/*Vegetable Fat*/
  dfat19n/*Dairy Fat*/
  sat19n/*Total Saturated Fat*/
  mon19n/*Total Monounsaturated Fat*/
  mufa_pl19n/*Plant Monounsaturated Fat*/
  mufa_a19n/*Animal Monounsaturated Fat*/
  poly19n/*Total Polyunsaturated Fat*/
  omega19n/*Omega 3 (f205+f226) gm, no alpha 18:3*/
  /*trn1119n*/ /*Total Trans Fat*/

  prot19n
  aprot19n/*Animal Protein*/
  vprot19n/*Vegetable Protein*/
  dprot19n/*Dairy Protein*/
  );

data n19_nts;
set n19_nts;

array gram_cp carbo19n pot_carb19n leg_carb19n frtcb_njs19n vegcb_npl19n sumref_carb19n sum_whgcarb19n sumoth_crb_njs19n mill_grn19n int_grn19n ref_grn19n prot19n aprot19n vprot19n dprot19n adds_cer19n adds_fds19n;
array perc_cp carbo19p pot_carb19p leg_carb19p frtcb_njs19p vegcb_npl19p sumref_carb19p sum_whgcarb19p sumoth_crb_njs19p mill_grn19p int_grn19p ref_grn19p prot19p aprot19p vprot19p dprot19p adds_cer19p adds_fds19p;

array gram_f  tfat19n afat19n vfat19n dfat19n sat19n mon19n mufa_pl19n mufa_a19n poly19n omega19n /*trn1119n*/;
array perc_f  tfat19p afat19p vfat19p dfat19p sat19p mon19p mufa_pl19p mufa_a19p poly19p omega19p /*trnss19p*/;

do over gram_cp;
perc_cp=(4*gram_cp)/calor19n*100;
end;

do over gram_f;
perc_f=(9*gram_f)/calor19n*100;
end;

hfat19p=sum(mufa_pl19p,poly19p);
ufat19p=sum(sat19p,mufa_a19p/*,trnss19p*/);

hcarbo19p=sum(leg_carb19p,frtcb_njs19p,vegcb_npl19p,sum_whgcarb19p);
run;

proc rank data=n19_nts groups=11 out=n19_nts;
var    prot19p  aprot19p  vprot19p  tfat19p  afat19p  vfat19p  sum_whgcarb19p  sumoth_crb_njs19p  hfat19p  ufat19p;
ranks  prot19pd aprot19pd vprot19pd tfat19pd afat19pd vfat19pd sum_whgcarb19pd sumoth_crb_njs19pd hfat19pd ufat19pd;
run;

proc rank data=n19_nts descending groups=11 out=n19_nts;
var    carbo19p  sumref_carb19p  adds_cer19p  adds_fds19p  hcarbo19p;
ranks  carbo19pd sumref_carb19pd adds_cer19pd adds_fds19pd hcarbo19pd;
run;

data n19_nts;
set n19_nts;

lcs19=sum(carbo19pd,prot19pd,tfat19pd);
alcs19=sum(carbo19pd,aprot19pd,afat19pd);
vlcs19=sum(carbo19pd,vprot19pd,vfat19pd);

hlcs19=sum(sumref_carb19pd,vprot19pd,vfat19pd);
ulcs19=sum(hcarbo19pd     ,aprot19pd,afat19pd);

cqs19=sum(sum_whgcarb19pd,sumoth_crb_njs19pd,sumref_carb19pd);
sug_cqs19=sum(adds_cer19pd,adds_fds19pd);
run;

proc rank data=n19_nts groups=11 out=n19_nts;
var    carbo19p  sumref_carb19p  hcarbo19p   prot19p  aprot19p  vprot19p   ;
ranks  carbo19pdd sumref_carb19pdd hcarbo19pdd  prot19pdd aprot19pdd vprot19pdd  ;
run;

proc rank data=n19_nts descending groups=11 out=n19_nts ;
var    tfat19p  afat19p  vfat19p  ;
ranks  tfat19pdd afat19pdd vfat19pdd ;
run;

data n19_nts;
set n19_nts;

lfs19=sum(carbo19pdd,prot19pdd,tfat19pdd);
alfs19=sum(carbo19pdd,aprot19pdd,vfat19pdd);
vlfs19=sum(carbo19pdd,vprot19pdd,afat19pdd);
hlfs19=sum(hcarbo19pdd, vprot19pdd, afat19pdd);
ulfs19=sum(sumref_carb19pdd,aprot19pdd,vfat19pdd);
run;


*----------------3.Merge-----------------------------;
 data base;         
 merge der8919(in=mst) disease  ahei2010_9119 act8917 supp8915 meds8919 oc8911  
       nur89 nur91 nur93 nur95 nur97 nur99  nur01 nur03 nur05 nur07 nur09 nur11 nur13 nur15 nur17 nur19
       serv91 serv95 serv99 serv03 serv07 serv11 serv15 serv19
       n91_nts(in=ffq) n95_nts n99_nts n03_nts n07_nts n11_nts n15_nts n19_nts 
       nses8917 end=_end_;
       by id;    
       exrec=1; if first.id and mst then exrec=0;
       baseffq=1; if  ffq then baseffq=0; 

/* Family history of diabetes, MI, stroke, and cancer */
    if dbfh89=1 or dbfh97=1  or dbfh01=1  or dbfh05=1 or dbfh09=1 or dbfh13=1  then dbfh=1;else dbfh=0; /* last dbfh 13*/
    if mifh89=1 or mifh93=1 or mifh97=1 or mifh01=1 or mifh05=1  or mifh09=1 or mifh13=1 then mifh=1; else mifh=0;
    if strfh97=1 or strfh01=1 or strfh05=1 or strfh09=1 or strfh13=1 then strfh=1; else strfh=0;
    if cafh89=1 or cafh97=1 or cafh01=1 or cafh05=1 or cafh09=1 or cafh13=1 then cafh=1; else cafh=0;  

/************** cancer: any confirmed cancer (from registry and self-reported) *****************/
  array canc_prev {*} can89  can91  can93  can95  can97  can99  can01  can03  can05  can07  can09  can11 can13  can15 can17 can19;
  array canc_new  {*} cann89 cann91 cann93 cann95 cann97 cann99 cann01 cann03 cann05 cann07 cann09 cann11 cann13 cann15 cann17 cann19;

   /* set first cancer */
   if dtdx_allca > 0 then first_dtdx_ca = dtdx_allca; /* any confirmed cancer */
   else if dtdxca > 0 then first_dtdx_ca = dtdxca;
   else first_dtdx_ca = .;

    array qqyear {*} retmo89 retmo91 retmo93 retmo95   retmo97   retmo99    retmo01   retmo03   retmo05   retmo07   retmo09   
                     retmo11 retmo13 retmo15 retmo17 retmo19; 

    do i = 1 to dim(qqyear);
        if qqyear{i} > 0 and first_dtdx_ca > 0 and first_dtdx_ca <= qqyear{i} then canc_new{i} = 1;
        else if canc_prev{i} = 1 then canc_new{i} = 1;  /* keep self-reported if earlier */
    end;

  /* cumulative cancer status: once cancer, always cancer */
    do i = 2 to dim(canc_new);
        if canc_new{i-1} = 1 then canc_new{i} = 1;
        if canc_prev{i-1} = 1 then canc_prev{i} = 1;
    end;

  /* Replace old canc array with new registry-based cancer status */
    do i = 1 to dim(canc_new);
        canc_prev{i} = canc_new{i};
        if canc_prev{i} ne 1 then canc_prev{i} = 0;
    end;


/* chronic disease status (flag for stopping diet updates) */
    /*array canc_prev {*} can89  can91  can93  can95  can97  can99  can01  can03  can05  can07  can09  can11 can13  can15 can17 cabg19;*/
    array cabg    {*} xx     xx     xx     cabg95 cabg97 cabg99 cabg01 cabg03 cabg05 cabg07 cabg09 cabg11 cabg13 cabg15 cabg17 cabg19;
    array mi      {*} mi89   mi91   mi93   mi95   mi97   mi99   mi01   mi03   mi05   mi07   mi09   mi11 mi13 mi15 mi17 mi19;
    array strk    {*} str89  str91  str93  str95  str97  str99  str01  str03  str05  str07  str09  str11 str13  str15  str17 str19;
    array hbp     {*} hbp89  hbp91  hbp93  hbp95  hbp97  hbp99  hbp01  hbp03  hbp05  hbp07  hbp09  hbp11 hbp13  hbp15 hbp17 hbp19;
    array chol    {*} chol89 chol91 chol93 chol95 chol97 chol99 chol01 chol03 chol05 chol07 chol09 chol11 chol13 chol15 chol17 chol19;
    array diab    {*} db89   db91   db93   db95   db97   db99   db01   db03   db05   db07   db09   db11 db13   db15 db17 db19;
    array chrodis {*} chrodis89   chrodis91   chrodis93   chrodis95   chrodis97   chrodis99   chrodis01   chrodis03   
                      chrodis05   chrodis07   chrodis09   chrodis11 chrodis13   chrodis15  chrodis17 chrodis19;     
   
do j=2 to DIM(hbp);
      *First, make intermediate outcomes cumulative varibles: i.e. once cancer, always cancer;  
      if hbp{j-1}=1    then hbp{j}=1;
      if chol{j-1}=1   then chol{j}=1;
      if diab{j-1}=1   then diab{j}=1;
      if mi{j-1}=1     then mi{j}=1;
      if cabg{j-1}=1   then cabg{j}=1;
      if strk{j-1}=1   then strk{j}=1;
      if canc_prev{j-1}=1   then canc_prev{j}=1;
    end;

    do i=1 to DIM(hbp);
        if canc_prev {i} ne 1    then canc_prev {i} = 0;
        if mi   {i} ne 1    then mi   {i} = 0;
        if hbp  {i} ne 1    then hbp  {i} = 0;
        if chol {i} ne 1    then chol {i} = 0;
        if diab {i} ne 1    then diab {i} = 0;
        if strk {i} ne 1    then strk {i} = 0;
        if cabg {i} ne 1    then cabg {i} = 0;
        if mi(i) =1 or strk(i)=1 or cabg(i)=1 or dtdxcvd >0 or diab(i)=1 or canc_prev(i)=1  then chrodis(i)=1; else chrodis(i)=0; 
    end;

    do i = 2 to dim(chrodis); if chrodis{i-1} = 1 then chrodis{i} = 1; 
    end; /* Make chrodis cumulative */

*** nSES;
if nSES_89 = . then nSES_89m = -0.426; else nSES_89m = nSES_89;
if nSES_91 = . then nSES_91m = nSES_89m; else nSES_91m = nSES_91;

/***************** Diet scores: carry forward + STOP updating after cancer/CVD *******************/
/* updated scores */ 
array dis   {*} chrodis91 chrodis95 chrodis99 chrodis03 chrodis07 chrodis11 chrodis15 chrodis19;
array calor {*} calor91n calor95n calor99n calor03n calor07n calor11n calor15n calor19n; 
array aheim {*} ahei91 ahei95 ahei99 ahei03 ahei07 ahei11 ahei15 ahei19;
array lcs   {*} lcs91 lcs95 lcs99 lcs03 lcs07 lcs11 lcs15 lcs19;
array alcs  {*} alcs91 alcs95 alcs99 alcs03 alcs07 alcs11 alcs15 alcs19;
array vlcs  {*} vlcs91 vlcs95 vlcs99 vlcs03 vlcs07 vlcs11 vlcs15 vlcs19;
array ulcs  {*} ulcs91 ulcs95 ulcs99 ulcs03 ulcs07 ulcs11 ulcs15 ulcs19;
array hlcs  {*} hlcs91 hlcs95 hlcs99 hlcs03 hlcs07 hlcs11 hlcs15 hlcs19;
array lfs   {*} lfs91 lfs95 lfs99 lfs03 lfs07 lfs11 lfs15 lfs19;
array ulfs  {*} ulfs91 ulfs95 ulfs99 ulfs03 ulfs07 ulfs11 ulfs15 ulfs19;
array hlfs  {*} hlfs91 hlfs95 hlfs99 hlfs03 hlfs07 hlfs11 hlfs15 hlfs19;
array alfs {*} alfs91 alfs95 alfs99 alfs03 alfs07 alfs11 alfs15 alfs19;
array vlfs {*} vlfs91 vlfs95 vlfs99 vlfs03 vlfs07 vlfs11 vlfs15 vlfs19;

array    carbp  {*}   carbo91p         carbo95p       carbo99p         carbo03p       carbo07p         carbo11p       carbo15p         carbo19p   ;
array    hcarbp  {*}   hcarbo91p        hcarbo95p      hcarbo99p        hcarbo03p      hcarbo07p        hcarbo11p      hcarbo15p        hcarbo19p  ;
array sumref_carbp {*}  sumref_carb91p   sumref_carb95p sumref_carb99p   sumref_carb03p sumref_carb07p   sumref_carb11p  sumref_carb15p   sumref_carb19p      ;
array     protp  {*}   prot91p          prot95p        prot99p          prot03p        prot07p          prot11p        prot15p          prot19p;
array    aprotp  {*}   aprot91p         aprot95p       aprot99p         aprot03p       aprot07p         aprot11p       aprot15p         aprot19p ;
array    vprotp  {*}   vprot91p         vprot95p       vprot99p         vprot03p       vprot07p         vprot11p       vprot15p         vprot19p ;
array    tfatp   {*}   tfat91p          tfat95p        tfat99p          tfat03p        tfat07p          tfat11p        tfat15p          tfat19p;
array    afatp   {*}   afat91p          afat95p        afat99p          afat03p        afat07p          afat11p        afat15p          afat19p;
array    vfatp   {*}   vfat91p          vfat95p        vfat99p          vfat03p        vfat07p          vfat11p        vfat15p          vfat19p;

/* carry forward for updated  scores */
do i=2 to dim(calor);
if calor{i}=.   then calor{i}=calor{i-1};
if aheim{i}=.   then aheim{i}=aheim{i-1};
if lcs{i}=.     then lcs{i}=lcs{i-1};
if alcs{i}=.    then alcs{i}=alcs{i-1};
if vlcs{i}=.    then vlcs{i}=vlcs{i-1};
if ulcs{i}=.    then ulcs{i}=ulcs{i-1};
if hlcs{i}=.    then hlcs{i}=hlcs{i-1};
if lfs{i}=.     then lfs{i}=lfs{i-1};
if ulfs{i}=.    then ulfs{i}=ulfs{i-1};
if hlfs{i}=.    then hlfs{i}=hlfs{i-1};
if alfs{i}=.   then alfs{i}=alfs{i-1};
if vlfs{i}=.   then vlfs{i}=vlfs{i-1};

if carbp{i}=.     then carbp{i}=carbp{i-1};
if hcarbp{i}=.    then hcarbp{i}=hcarbp{i-1};
if sumref_carbp{i}=.    then sumref_carbp{i}=sumref_carbp{i-1};
if protp{i}=.   then protp{i}=protp{i-1};
if aprotp{i}=.   then aprotp{i}=aprotp{i-1};
if vprotp{i}=.     then vprotp{i}=vprotp{i-1};
if tfatp{i}=.    then tfatp{i}=tfatp{i-1};
if afatp{i}=.    then afatp{i}=afatp{i-1};
if vfatp{i}=.   then vfatp{i}=vfatp{i-1};
end;

/* stop updated scores */
array calors {*} calor91ns calor95ns calor99ns calor03ns calor07ns calor11ns calor15ns calor19ns; 
array aheims {*} ahei91s ahei95s ahei99s ahei03s ahei07s ahei11s ahei15s ahei19s;
array lcss {*} lcs91s lcs95s lcs99s lcs03s lcs07s lcs11s lcs15s lcs19s;
array alcss {*} alcs91s alcs95s alcs99s alcs03s alcs07s alcs11s alcs15s alcs19s;
array vlcss {*} vlcs91s vlcs95s vlcs99s vlcs03s vlcs07s vlcs11s vlcs15s vlcs19s;
array ulcss {*} ulcs91s ulcs95s ulcs99s ulcs03s ulcs07s ulcs11s ulcs15s ulcs19s;
array hlcss {*} hlcs91s hlcs95s hlcs99s hlcs03s hlcs07s hlcs11s hlcs15s hlcs19s;
array lfss {*} lfs91s lfs95s lfs99s lfs03s lfs07s lfs11s lfs15s lfs19s;
array ulfss {*} ulfs91s ulfs95s ulfs99s ulfs03s ulfs07s ulfs11s ulfs15s ulfs19s;
array hlfss {*} hlfs91s hlfs95s hlfs99s hlfs03s hlfs07s hlfs11s hlfs15s hlfs19s;
array alfss {*} alfs91s alfs95s alfs99s alfs03s alfs07s alfs11s alfs15s alfs19s;
array vlfss {*} vlfs91s vlfs95s vlfs99s vlfs03s vlfs07s vlfs11s vlfs15s vlfs19s;

array    carbs  {*}   carbo91s         carbo95s       carbo99s         carbo03s       carbo07s         carbo11s       carbo15s         carbo19s   ;
array    hcarbs  {*}   hcarbo91s        hcarbo95s      hcarbo99s        hcarbo03s      hcarbo07s        hcarbo11s      hcarbo15s        hcarbo19s  ;
array sumref_carbs {*}  sumref_carb91s   sumref_carb95s sumref_carb99s   sumref_carb03s sumref_carb07s   sumref_carb11s  sumref_carb15s   sumref_carb19s      ;
array     prots  {*}   prot91s          prot95s        prot99s          prot03s        prot07s          prot11s        prot15s          prot19s; 
array    aprots  {*}   aprot91s         aprot95s       aprot99s         aprot03s       aprot07s         aprot11s       aprot15s         aprot19s ;
array    vprots  {*}   vprot91s         vprot95s       vprot99s         vprot03s       vprot07s         vprot11s       vprot15s         vprot19s ;
array    tfats   {*}   tfat91s          tfat95s        tfat99s          tfat03s        tfat07s          tfat11s        tfat15s          tfat19s; 
array    afats   {*}   afat91s          afat95s        afat99s          afat03s        afat07s          afat11s        afat15s          afat19s; 
array    vfats   {*}   vfat91s          vfat95s        vfat99s          vfat03s        vfat07s          vfat11s        vfat15s          vfat19s; 


/* stop updating scores after cancer or CVD */
do i=1 to dim(calor);
  calors{i}=calor{i};  
  aheims{i}=aheim{i};
  lcss{i}=lcs{i};  
  alcss{i}=alcs{i};
  vlcss{i}=vlcs{i};
  ulcss{i}=ulcs{i};
  hlcss{i}=hlcs{i};
  lfss{i}=lfs{i};
  ulfss{i}=ulfs{i};
  hlfss{i}=hlfs{i};
  alfss{i}=alfs{i};
  vlfss{i}=vlfs{i};
end;

do i=2 to dim(calor);
if  dis{i}=1    then calors{i}=calors{i-1};
if  dis{i}=1    then aheims{i}=aheims{i-1};
if  dis{i}=1    then lcss{i}=lcss{i-1};
if  dis{i}=1    then alcss{i}=alcss{i-1};
if  dis{i}=1    then vlcss{i}=vlcss{i-1};
if  dis{i}=1    then ulcss{i}=ulcss{i-1};
if  dis{i}=1    then hlcss{i}=hlcss{i-1};
if  dis{i}=1    then lfss{i}=lfss{i-1};
if  dis{i}=1    then ulfss{i}=ulfss{i-1};
if  dis{i}=1    then hlfss{i}=hlfss{i-1};
if  dis{i}=1    then alfss{i}=alfss{i-1};
if  dis{i}=1    then vlfss{i}=vlfss{i-1};
end;

/*calculate cumulative average of updated score**/
array nut{8,12}
calor91ns ahei91s lcs91s alcs91s vlcs91s ulcs91s hlcs91s lfs91s ulfs91s hlfs91s alfs91s vlfs91s
calor95ns ahei95s lcs95s alcs95s vlcs95s ulcs95s hlcs95s lfs95s ulfs95s hlfs95s alfs95s vlfs95s
calor99ns ahei99s lcs99s alcs99s vlcs99s ulcs99s hlcs99s lfs99s ulfs99s hlfs99s alfs99s vlfs99s
calor03ns ahei03s lcs03s alcs03s vlcs03s ulcs03s hlcs03s lfs03s ulfs03s hlfs03s alfs03s vlfs03s
calor07ns ahei07s lcs07s alcs07s vlcs07s ulcs07s hlcs07s lfs07s ulfs07s hlfs07s alfs07s vlfs07s
calor11ns ahei11s lcs11s alcs11s vlcs11s ulcs11s hlcs11s lfs11s ulfs11s hlfs11s alfs11s vlfs11s
calor15ns ahei15s lcs15s alcs15s vlcs15s ulcs15s hlcs15s lfs15s ulfs15s hlfs15s alfs15s vlfs15s
calor19ns ahei19s lcs19s alcs19s vlcs19s ulcs19s hlcs19s lfs19s ulfs19s hlfs19s alfs19s vlfs19s
;
array nutv{8,12}
calor91nsv ahei91sv lcs91sv alcs91sv vlcs91sv ulcs91sv hlcs91sv lfs91sv ulfs91sv hlfs91sv alfs91sv vlfs91sv
calor95nsv ahei95sv lcs95sv alcs95sv vlcs95sv ulcs95sv hlcs95sv lfs95sv ulfs95sv hlfs95sv alfs95sv vlfs95sv
calor99nsv ahei99sv lcs99sv alcs99sv vlcs99sv ulcs99sv hlcs99sv lfs99sv ulfs99sv hlfs99sv alfs99sv vlfs99sv
calor03nsv ahei03sv lcs03sv alcs03sv vlcs03sv ulcs03sv hlcs03sv lfs03sv ulfs03sv hlfs03sv alfs03sv vlfs03sv
calor07nsv ahei07sv lcs07sv alcs07sv vlcs07sv ulcs07sv hlcs07sv lfs07sv ulfs07sv hlfs07sv alfs07sv vlfs07sv
calor11nsv ahei11sv lcs11sv alcs11sv vlcs11sv ulcs11sv hlcs11sv lfs11sv ulfs11sv hlfs11sv alfs11sv vlfs11sv
calor15nsv ahei15sv lcs15sv alcs15sv vlcs15sv ulcs15sv hlcs15sv lfs15sv ulfs15sv hlfs15sv alfs15sv vlfs15sv
calor19nsv ahei19sv lcs19sv alcs19sv vlcs19sv ulcs19sv hlcs19sv lfs19sv ulfs19sv hlfs19sv alfs19sv vlfs19sv
;

/* calculate cumulative average intake */
do j=1 to 12;
  nutv{1,j}=nut{1,j};
 do i=2 to 8;
  sumvar=0;
   n=0;
     do k=1 to i;
      if (nut{k,j} ne .) then do;
      n=n+1; sumvar=sumvar+nut{k,j};
      end;
     end;
   if n=0 then nutv{i,j}=nut{1,j};
   else nutv{i,j}=sumvar/n;
  end;
end;


/* stop updating scores after cancer or CVD */
do i=1 to dim(carbp);
  carbs{i}=carbp{i};  
  hcarbs{i}=hcarbp{i};
  sumref_carbs{i}=sumref_carbp{i};
  prots{i}=protp{i};
  aprots{i}=aprotp{i};
  vprots{i}=vprotp{i};
  tfats{i}=tfatp{i};
  afats{i}=afatp{i};
  vfats{i}=vfatp{i};
end;

do i=2 to dim(carbp);
if  dis{i}=1    then carbs{i}=carbs{i-1};
if  dis{i}=1    then hcarbs{i}=hcarbs{i-1};
if  dis{i}=1    then sumref_carbs{i}=sumref_carbs{i-1};
if  dis{i}=1    then prots{i}=prots{i-1};
if  dis{i}=1    then aprots{i}=aprots{i-1};
if  dis{i}=1    then vprots{i}=vprots{i-1};
if  dis{i}=1    then tfats{i}=tfats{i-1};
if  dis{i}=1    then afats{i}=afats{i-1};
if  dis{i}=1    then vfats{i}=vfats{i-1};
end;

/*calculate cumulative average of updated score**/
array macron{8,9}
carbo91ns hcarbo91s  sumref_carb91s prot91s aprot91s vprot91s tfat91s  afat91s vfat91s
carbo95ns hcarbo95s  sumref_carb95s prot95s aprot95s vprot95s tfat95s  afat95s vfat95s
carbo99ns hcarbo99s  sumref_carb99s prot99s aprot99s vprot99s tfat99s  afat99s vfat99s
carbo03ns hcarbo03s  sumref_carb03s prot03s aprot03s vprot03s tfat03s  afat03s vfat03s
carbo07ns hcarbo07s  sumref_carb07s prot07s aprot07s vprot07s tfat07s  afat07s vfat07s
carbo11ns hcarbo11s  sumref_carb11s prot11s aprot11s vprot11s tfat11s  afat11s vfat11s
carbo15ns hcarbo15s  sumref_carb15s prot15s aprot15s vprot15s tfat15s  afat15s vfat15s
carbo19ns hcarbo19s  sumref_carb19s prot19s aprot19s vprot19s tfat19s  afat19s vfat19s
;
array macronv{8,9}
carbo91nsv hcarbo91sv  sumref_carb91sv prot91sv aprot91sv vprot91sv tfat91sv  afat91sv vfat91sv
carbo95nsv hcarbo95sv  sumref_carb95sv prot95sv aprot95sv vprot95sv tfat95sv  afat95sv vfat95sv
carbo99nsv hcarbo99sv  sumref_carb99sv prot99sv aprot99sv vprot99sv tfat99sv  afat99sv vfat99sv
carbo03nsv hcarbo03sv  sumref_carb03sv prot03sv aprot03sv vprot03sv tfat03sv  afat03sv vfat03sv
carbo07nsv hcarbo07sv  sumref_carb07sv prot07sv aprot07sv vprot07sv tfat07sv  afat07sv vfat07sv
carbo11nsv hcarbo11sv  sumref_carb11sv prot11sv aprot11sv vprot11sv tfat11sv  afat11sv vfat11sv
carbo15nsv hcarbo15sv  sumref_carb15sv prot15sv aprot15sv vprot15sv tfat15sv  afat15sv vfat15sv
carbo19nsv hcarbo19sv  sumref_carb19sv prot19sv aprot19sv vprot19sv tfat19sv  afat19sv vfat19sv
;


/* calculate cumulative average intake */
do j=1 to 9;
  macronv{1,j}=macron{1,j};
 do i=2 to 8;
  sumvar=0;
   n=0;
     do k=1 to i;
      if (macron{k,j} ne .) then do;
      n=n+1; sumvar=sumvar+macron{k,j};
      end;
     end;
   if n=0 then macronv{i,j}=macron{1,j};
   else macronv{i,j}=sumvar/n;
  end;
end;


array    whgrn_o  {*}   whgrn_s91       whgrn_s95      whgrn_s99       whgrn_s03       whgrn_s07      whgrn_s11    whgrn_s15   whgrn_s19;
array    rgrain_o {*}   rgrain_s91      rgrain_s95     rgrain_s99      rgrain_s03      rgrain_s07     rgrain_s11   rgrain_s15  rgrain_s19;
array    frui_o   {*}   frui_s91    frui_s95    frui_s99 frui_s03 frui_s07 frui_s11 frui_s15 frui_s19 ;
array     leg_o   {*}   leg_s91        leg_s95         leg_s99        leg_s03        leg_s07       leg_s11     leg_s15  leg_s19 ;
array    pota_o   {*}   pota_s91       pota_s95        pota_s99       pota_s03       pota_s07      pota_s11    pota_s15 pota_s19 ;
array  veg_noleg_o{*}   veg_noleg_s91   veg_noleg_s95  veg_noleg_s99       veg_noleg_s03       veg_noleg_s07      veg_noleg_s11    veg_noleg_s15 veg_noleg_s19 ;
array    fish_o   {*}   fish_s91        fish_s95       fish_s99        fish_s03        fish_s07       fish_s11     fish_s15  fish_s19 ;
array   rpmeats_o {*}   rpmeats_s91     rpmeats_s95    rpmeats_s99        rpmeats_s03        rpmeats_s07       rpmeats_s11     rpmeats_s15 rpmeats_s19 ;
array    ssb_o    {*}   ssb_s91         ssb_s95        ssb_s99        ssb_s03        ssb_s07       ssb_s11     ssb_s15 ssb_s19 ;
array    dairy_o   {*}   dairy_s91        dairy_s95       dairy_s99        dairy_s03        dairy_s07       dairy_s11     dairy_s15 dairy_s19 ;
array swtsnochoc_o {*}   swtsnochoc_s91     swtsnochoc_s95    swtsnochoc_s99        swtsnochoc_s03        swtsnochoc_s07       swtsnochoc_s11     swtsnochoc_s15  swtsnochoc_s19 ;
array    poul_o    {*}   poul_s91         poul_s95        poul_s99        poul_s03        poul_s07       poul_s11     poul_s15  poul_s19 ;

array    whgrn_os  {*}   whgrn_s91s       whgrn_s95s      whgrn_s99s       whgrn_s03s       whgrn_s07s      whgrn_s11s    whgrn_s15s whgrn_s19s ;
array    rgrain_os {*}   rgrain_s91s      rgrain_s95s     rgrain_s99s      rgrain_s03s      rgrain_s07s     rgrain_s11s   rgrain_s15s rgrain_s19s ;
array    frui_os   {*}   frui_s91s    frui_s95s    frui_s99s frui_s03s frui_s07s frui_s11s frui_s15s frui_s19s ;
array     leg_os   {*}   leg_s91s        leg_s95s         leg_s99s        leg_s03s        leg_s07s       leg_s11s     leg_s15s  leg_s19s ;
array    pota_os   {*}   pota_s91s       pota_s95s        pota_s99s       pota_s03s       pota_s07s      pota_s11s    pota_s15s pota_s19s ;
array  veg_noleg_os{*}   veg_noleg_s91s   veg_noleg_s95s  veg_noleg_s99s       veg_noleg_s03s       veg_noleg_s07s      veg_noleg_s11s    veg_noleg_s15s veg_noleg_s19s ;
array    fish_os   {*}   fish_s91s        fish_s95s       fish_s99s        fish_s03s        fish_s07s       fish_s11s     fish_s15s fish_s19s ;
array   rpmeats_os {*}   rpmeats_s91s     rpmeats_s95s    rpmeats_s99s        rpmeats_s03s        rpmeats_s07s       rpmeats_s11s     rpmeats_s15s rpmeats_s19s ;
array    ssb_os    {*}   ssb_s91s         ssb_s95s        ssb_s99s        ssb_s03s        ssb_s07s       ssb_s11s     ssb_s15s ssb_s19s ;
array    dairy_os   {*}   dairy_s91s        dairy_s95s       dairy_s99s        dairy_s03s        dairy_s07s       dairy_s11s     dairy_s15s dairy_s19s ;
array swtsnochoc_os {*}   swtsnochoc_s91s     swtsnochoc_s95s    swtsnochoc_s99s        swtsnochoc_s03s        swtsnochoc_s07s       swtsnochoc_s11s     swtsnochoc_s15s  swtsnochoc_s19s ;
array    poul_os    {*}   poul_s91s         poul_s95s        poul_s99s        poul_s03s        poul_s07s       poul_s11s     poul_s15s poul_s19s ;

do i=2 to dim(whgrn_o);
  whgrn_os{i}=whgrn_o{i};
  rgrain_os{i}=rgrain_o{i};
  frui_os{i}=frui_o{i};
  leg_os{i}=leg_o{i};
  pota_os{i}=pota_o{i};
  veg_noleg_os{i}=veg_noleg_o{i};
  fish_os{i} = fish_o{i};
  rpmeats_os{i}=rpmeats_o{i};
  ssb_os{i}=ssb_o{i};
  dairy_os{i}=dairy_o{i};
  swtsnochoc_os{i}=swtsnochoc_o{i};
  poul_os{i}=poul_o{i};
end;


do i=2 to dim(whgrn_o);
if  dis{i}=1    then whgrn_os{i}=whgrn_os{i-1};
if  dis{i}=1    then rgrain_os{i}=rgrain_os{i-1};
if  dis{i}=1    then frui_os{i}=frui_os{i-1};
if  dis{i}=1    then leg_os{i}=leg_os{i-1};
if  dis{i}=1    then pota_os{i}=pota_os{i-1};
if  dis{i}=1    then veg_noleg_os{i}=veg_noleg_os{i-1};
if  dis{i}=1    then fish_os{i}=fish_os{i-1};
if  dis{i}=1    then rpmeats_os{i}=rpmeats_os{i-1};
if  dis{i}=1    then ssb_os{i}=ssb_os{i-1};
if  dis{i}=1    then dairy_os{i}=dairy_os{i-1};
if  dis{i}=1    then swtsnochoc_os{i}=swtsnochoc_os{i-1};
if  dis{i}=1    then poul_os{i}=poul_os{i-1};
end;


/* calculate cumulative average intake */
array foodlist{8,12}
whgrn_s91s  rgrain_s91s  frui_s91s  leg_s91s pota_s91s  veg_noleg_s91s  fish_s91s  rpmeats_s91s  ssb_s91s  dairy_s91s swtsnochoc_s91s  poul_s91s
whgrn_s95s  rgrain_s95s  frui_s95s  leg_s95s pota_s95s  veg_noleg_s95s  fish_s95s  rpmeats_s95s  ssb_s95s  dairy_s95s swtsnochoc_s95s  poul_s95s
whgrn_s99s  rgrain_s99s  frui_s99s  leg_s99s pota_s99s  veg_noleg_s99s  fish_s99s  rpmeats_s99s  ssb_s99s  dairy_s99s swtsnochoc_s99s  poul_s99s
whgrn_s03s  rgrain_s03s  frui_s03s  leg_s03s pota_s03s  veg_noleg_s03s  fish_s03s  rpmeats_s03s  ssb_s03s  dairy_s03s swtsnochoc_s03s  poul_s03s
whgrn_s07s  rgrain_s07s  frui_s07s  leg_s07s pota_s07s  veg_noleg_s07s  fish_s07s  rpmeats_s07s  ssb_s07s  dairy_s07s swtsnochoc_s07s  poul_s07s
whgrn_s11s  rgrain_s11s  frui_s11s  leg_s11s pota_s11s  veg_noleg_s11s  fish_s11s  rpmeats_s11s  ssb_s11s  dairy_s11s swtsnochoc_s11s  poul_s11s
whgrn_s15s  rgrain_s15s  frui_s15s  leg_s15s pota_s15s  veg_noleg_s15s  fish_s15s  rpmeats_s15s  ssb_s15s  dairy_s15s swtsnochoc_s15s  poul_s15s
whgrn_s19s  rgrain_s19s  frui_s19s  leg_s19s pota_s19s  veg_noleg_s19s  fish_s19s  rpmeats_s19s  ssb_s19s  dairy_s19s swtsnochoc_s19s  poul_s19s
;
array foodlistv{8,12}
whgrn_s91sv  rgrain_s91sv  frui_s91sv  leg_s91sv pota_s91sv  veg_noleg_s91sv  fish_s91sv  rpmeats_s91sv  ssb_s91sv  dairy_s91sv swtsnochoc_s91sv  poul_s91sv
whgrn_s95sv  rgrain_s95sv  frui_s95sv  leg_s95sv pota_s95sv  veg_noleg_s95sv  fish_s95sv  rpmeats_s95sv  ssb_s95sv  dairy_s95sv swtsnochoc_s95sv  poul_s95sv
whgrn_s99sv  rgrain_s99sv  frui_s99sv  leg_s99sv pota_s99sv  veg_noleg_s99sv  fish_s99sv  rpmeats_s99sv  ssb_s99sv  dairy_s99sv swtsnochoc_s99sv  poul_s99sv
whgrn_s03sv  rgrain_s03sv  frui_s03sv  leg_s03sv pota_s03sv  veg_noleg_s03sv  fish_s03sv  rpmeats_s03sv  ssb_s03sv  dairy_s03sv swtsnochoc_s03sv  poul_s03sv
whgrn_s07sv  rgrain_s07sv  frui_s07sv  leg_s07sv pota_s07sv  veg_noleg_s07sv  fish_s07sv  rpmeats_s07sv  ssb_s07sv  dairy_s07sv swtsnochoc_s07sv  poul_s07sv
whgrn_s11sv  rgrain_s11sv  frui_s11sv  leg_s11sv pota_s11sv  veg_noleg_s11sv  fish_s11sv  rpmeats_s11sv  ssb_s11sv  dairy_s11sv swtsnochoc_s11sv  poul_s11sv
whgrn_s15sv  rgrain_s15sv  frui_s15sv  leg_s15sv pota_s15sv  veg_noleg_s15sv  fish_s15sv  rpmeats_s15sv  ssb_s15sv  dairy_s15sv swtsnochoc_s15sv  poul_s15sv
whgrn_s19sv  rgrain_s19sv  frui_s19sv  leg_s19sv pota_s19sv  veg_noleg_s19sv  fish_s19sv  rpmeats_s19sv  ssb_s19sv  dairy_s19sv swtsnochoc_s19sv  poul_s19sv
;


do j=1 to 12;
  foodlistv{1,j}=foodlist{1,j};
 do i=2 to 8;
  sumvar=0;
   n=0;
     do k=1 to i;
      if (foodlist{k,j} ne .) then do;
      n=n+1; sumvar=sumvar+foodlist{k,j};
      end;
     end;
   if n=0 then foodlistv{i,j}=foodlist{1,j};
   else foodlistv{i,j}=sumvar/n;
  end;
end;

* MAMMOGRAPHY USE *;
* use one i+1 mem memography status *;
**********************;
/*Mammo in last 2 yrs -   
  1 = no
  2 = yes screening
  3 = yes symptoms
  4 = pass through */
if mamm91 in (2,3) then mam89=1; if mamm91=1 then mam89=0; 
if mamm93 in (2,3) then mam91=1; if mamm93=1 then mam91=0;
if mamm95 in (2,3) then mam93=1; if mamm95=1 then mam93=0;
if mamsc97=1 or mamsy97=1 then mam95=1; if mamno97=1 then mam95=0;
if mamsc99=1 or mamsy99=1 then mam97=1; if mamno99=1 then mam97=0;
if mamsc01=1 or mamsy01=1 then mam99=1; if mamno01=1 then mam99=0;
mam01=mam99; /*not collected in 03*/
if mamsc05=1 or mamsy05=1 then mam03=1; if mamno05=1 then mam03=0;
if mamsc07=1 or mamsy07=1 then mam05=1; if mamno07=1 then mam05=0;
if mamsc09=1 or mamsy09=1 then mam07=1; if mamno09=1 then mam07=0;
if mamsc11=1 or mamsy11=1 then mam09=1; if mamno11=1 then mam09=0;
if mamsc13=1 or mamsy13=1 then mam11=1; if mamno13=1 then mam11=0;
if mamsc15=1 or mamsy15=1 then mam13=1; if mamno15=1 then mam13=0;
if mamsc17=1 or mamsy17=1 then mam15=1; if mamno17=1 then mam15=0;
if mamsc19=1 or mamsy19=1 then mam17=1; if mamno19=1 then mam17=0;


/***********/
/*endoscopy*/
/***********/
if sigm91=1 then endo91=1; /*endo no*/
else if sigm91 in (2,3) then endo91=2; /*endo yes*/
if sigm93=1 then endo93=1; /*endo no*/
else if sigm93 in (2,3) then endo93=2; /*endo yes*/
if sigm95=1 then endo95=1; /*endo no*/
else if sigm95 in (2,3) then endo95=2; /*endo yes*/
if sigsc97=1 or sigsy97=1 then endo97=2; /*endo yes*/
else if sigmx97=1 then endo97=1;
if sigsc99=1 or sigsy99=1 then endo99=2; /*endo yes*/
else if sigmx99=1 then endo99=1;
if sigm01=1 then endo01=2;
else if sigm01=2 then endo01=1;
if colsc03=2 or sigm03=2 then endo03=2;
else if colsc03=1 or sigm03=1 then endo03=1;
if virtc05=2 or colsc05=2 or sigsc05=2 or sig0305=1 or col0305=1 then endo05=2;
else if virtc05=1 or colsc05=1 or sigsc05=1 then endo05=1;
if virtc07=2 or colsc07=2 or sigsc07=2 then endo07=2;
else if virtc07=1 or colsc07=1 or sigsc07=1 then endo07=1;
if virtc09=2 or colsc09=2 or sigsc09=2 then endo09=2;
else if virtc09=1 or colsc09=1 or sigsc09=1 then endo09=1;
if virtc11=2 or colsc11=2 or sigsc11=2 then endo11=2;
else if virtc11=1 or colsc11=1 or sigsc11=1 then endo11=1;
if virtc13=2 or colsc13=2 or sigsc13=2 then endo13=2;
else if virtc13=1 or colsc13=1 or sigsc13=1 then endo13=1;
if virtc15=2 or colsc15=2 or sigsc15=2 then endo15=2;
else if virtc15=1 or colsc15=1 or sigsc15=1 then endo15=1;
if virtc17=2 or colsc17=2 or sigsc17=2 then endo17=2;
else if virtc17=1 or colsc17=1 or sigsc17=1 then endo17=1;
if virtc19=2 or colsc19=2 or sigsc19=2 then endo19=2;
else if virtc19=1 or colsc19=1 or sigsc19=1 then endo19=1;

*-----------*
 Use 2005 Qx
*-----------*;
/*For only missing variables of the past Qx, we use 2005 info*/
/*We do not impute endo[i] if person already endo[i-1]=yes or endo[i+1]=yes*/
/*never*/
if (endo91=2 or endo93=2 or endo95=2 or
    endo97=2 or endo99=2 or endo01=2 or endo03=2 or endo05=2)
   then atleast=1; /*yes edo at least once in some year. They cannot be "Never endo"*/
array arendo[*] endo91 endo93 endo95 endo97 endo99 endo01 endo03 endo05;
if (atleast ne 1 and signo05=1 and colno05=1) then do; /*never in 2005 qx*/
   do i=1 to dim(arendo);
      if arendo[i]=. then arendo[i]=1; /*no*/
   end;
end;

/*we fill from the latest information*/
/*2003/Later endo can fill endo05*/
if (endo03 ne 2 and endo05=. and (sig0305=1 or col0305=1))
   then endo05=2; /*yes*/
else if (endo05=. and sig0305 ne 1 and col0305 ne 1 and q05 in (1,2))
   then endo05=1; /*no*/

/*2001-2003 endo can fill endo03*/
if (endo01 ne 2 and endo05 ne 2 and endo03=. and (sig0105=1 or col0105=1))
   then endo03=2; /*yes*/
else if (endo03=. and sig0105 ne 1 and col0105 ne 1 and q05 in (1,2))
   then endo03=1; /*no*/

/*1999-2001 endo can fill endo01*/
if (endo99 ne 2 and endo03 ne 2 and endo01=. and (sig9905=1 or col9905=1))
   then endo01=2; /*yes*/
else if (endo01=. and sig9905 ne 1 and col9905 ne 1 and q05 in (1,2))
   then endo01=1; /*no*/

/*1997-1999 endo can fill endo99*/
if (endo97 ne 2 and endo01 ne 2 and endo99=. and (sig9705=1 or col9705=1))
   then endo99=2; /*yes*/
else if (endo99=. and sig9705 ne 1 and col9705 ne 1 and q05 in (1,2))
   then endo99=1; /*no*/

/*1995-1997 endo can fill endo97*/
if (endo95 ne 2 and endo99 ne 2 and endo97=. and (sig9505=1 or col9505=1))
   then endo97=2; /*yes*/
else if (endo97=. and sig9505 ne 1 and col9505 ne 1 and q05 in (1,2))
   then endo97=1; /*no*/

/*1994/Earlier, we just do not know exact year. but can fill no-endo*/
if (endo95=. and sig9405 ne 1 and col9405 ne 1 and q05 in (1,2)) then endo95=1; /*no*/
else if (endo93=. and sig9405 ne 1 and col9405 ne 1 and q05 in (1,2)) then endo93=1; /*no*/
else if (endo91=. and sig9405 ne 1 and col9405 ne 1 and q05 in (1,2)) then endo91=1; /*no*/

*****fill in if missing*******;
if endo91=. and .<batch91<800 then endo91=1; /*one answered long qx but missing, we assume no endo*/
if endo93=. and .<batch93<800 then endo93=1; /*one answered long qx but missing, we assume no endo*/
if endo95=. and .<batch95<800 then endo95=1; /*one answered long qx but missing, we assume no endo*/
if endo97=. and .<batch97<800 then endo97=1; /*one answered long qx but missing, we assume no endo*/
if endo99=. and q99 in (1,2) then endo99=1;
if endo01=. and q01 in (1,3) then endo01=1;
if endo03=. and q03 in (1,2) then endo03=1;
if endo05=. and q05 in (1,2) then endo05=1;
if endo07=. and q07 in (1,2,3,4,5) then endo07=1;
if endo09=. and q09 in (1,2) then endo09=1;
if endo11=. and q11 in (1,2,3,4,5) then endo11=1;
if endo13=. and q13 in (1,2,3) then endo13=1;
if endo15=. and q15 in (1,2,3,4) then endo15=1;
if endo19=. and q15 in (1,2,3) then endo19=1;

*--------------------------------------*
 reason for endoscopy, make dichotomous
*--------------------------------------*;
   /*1991-*/
   if sigm91=3 then reasoned91=1; /*symptom. NOTE, Precise reasons were not asked in 1991 Qx*/
   if sigm91=2 then reasoned91=2; /*screening*/
   /*1993-*/
   if sigm93=3 then reasoned93=1; /*symptom. NOTE, Precise reasons were not asked in 1993 Qx*/
   if sigm93=2 then reasoned93=2; /*screening*/
   /*1995-*/
   if sigm95=3 then reasoned95=1; /*symptom. NOTE, Precise reasons were not asked in 1995 Qx*/
   if sigm95=2 then reasoned95=2; /*screening*/
   /*1997-*/
   if sigsy97=1 then reasoned97=1; /*symptom. NOTE, Precise reasons were not asked in 1997 Qx*/
   if sigsc97=1 then reasoned97=2; /*screening*/
   /*1999-*/
   if sigsy99=1 then reasoned99=1; /*symptom. NOTE, Precise reasons were not asked in 1999 Qx*/
   if sigsc99=1 then reasoned99=2; /*screening*/
   /*2001*/
   if (vbld01=1 or fecal01=1 or abdpn01=1 or diarr01=1 or bariu01=1) then reasoned01=1; /*symptom*/
   else if hxclc01=1 or rout01=1 then reasoned01=2;                 /*screening or family hx*/
   /*2003*/
   if (vbld03=1 or fecal03=1 or abdpn03=1 or diarr03=1 or bariu03=1 or cpolp03=1 or virt03=1) then reasoned03=1; /*symptom*/
   else if (hxclc03=1 or rout03=1) then reasoned03=2;         /*screening or family hx*/
   /*2005*/
   if (vbld05=1 or fecal05=1 or abdpn05=1 or diarr05=1 or bariu05=1 or cpolp05=1 or virt05=1) then reasoned05=1; /*symptom*/
   else if (hxclc05=1 or rout05=1) then reasoned05=2; /*screening or family hx*/
   /*2007*/
   if (vbld07=1 or fecal07=1 or abdpn07=1 or diarr07=1 or bariu07=1 or cpolp07=1 or virt07=1) then reasoned07=1; /*symptom*/
   else if (hxclc07=1 or rout07=1) then reasoned07=2; /*screening or family hx*/
   /*2009*/
   if (vbld09=1 or fecal09=1 or abdpn09=1 or diarr09=1 or bariu09=1 or cpolp09=1 or virt09=1) then reasoned09=1; /*symptom*/
   else if (hxclc09=1 or rout09=1) then reasoned09=2; /*screening or family hx*/
   /*2011*/
   if (vbld11=1 or fecal11=1 or abdpn11=1 or diarr11=1 or bariu11=1 or cpolp11=1 or virt11=1) then reasoned11=1; /*symptom*/
   else if (hxclc11=1 or rout11=1) then reasoned11=2; /*screening or family hx*/
   /*2013*/
   if (vbld13=1 or fecal13=1 or abdpn13=1 or diarr13=1 or bariu13=1 or cpolp13=1  or virt13=1) then reasoned13=1; /*symptom*/
   else if (hxclc13=1 or rout13=1) then reasoned13=2; /*screening or family hx*/
   /*2015*/
   if (vbld15=1 or fecal15=1 or abdpn15=1 or diarr15=1 or fecdna15=1 or cpolp15=1  or virt15=1) then reasoned15=1; /*symptom*/
   else if (hxclc15=1 or rout15=1) then reasoned15=2; /*screening or family hx*/
   /*2017*/
   if (vbld17=1 or fecal17=1 or abdpn17=1 or diarr17=1 or fecdna17=1 or cpolp17=1  or virt17=1) then reasoned17=1; /*symptom*/
   else if (hxclc17=1 or rout17=1) then reasoned17=2; /*screening or family hx*/
   /*2019*/
   if (vbld19=1 or fecal19=1 or abdpn19=1 or diarr19=1 or bariu19=1 or cpolp19=1  or virt19=1) then reasoned19=1; /*symptom*/
   else if (hxclc19=1 or rout19=1) then reasoned19=2; /*screening or family hx*/

array endoxx [*]  endo91 endo93 endo95 endo97 endo99 endo01 endo03 endo05 endo07 endo09 endo11 endo13 endo15 endo17 endo19 ; /*1=no endo, 2=yes endo*/
array reasonxx[*] reasoned91 reasoned93 reasoned95 reasoned97 reasoned99 reasoned01 reasoned03 reasoned05 reasoned07 reasoned09 reasoned11 reasoned13 reasoned15 reasoned17 reasoned19 ;/*At the time of qx 1=symptom,2=screening*/


/***1. time varying, updated every 2 years, no endo, endo for symptoms, screening endo***/
array endo3c [*] endo3c91 endo3c93 endo3c95 endo3c97 endo3c99 endo3c01 endo3c03 endo3c05 endo3c07 endo3c09 endo3c11 endo3c13 endo3c15 endo3c17 endo3c19 ; /*1=no endo, 2=endo for symptoms, 3=screening endo*/

do i=1 to dim(endo3c);
   if endoxx[i]=1 then endo3c[i]=1; *no endo;
   else if endoxx[i]=2 and reasonxx[i]=1 then endo3c[i]=2; *endo for symptoms;
   else if endoxx[i]=2 and reasonxx[i]=2 then endo3c[i]=3; *endo for screening;
   else endo3c[i]=.; *missing;
end;

/***2. ever screening endo, eendo3c***/
array eendo3c [*] eendo3c91 eendo3c93 eendo3c95 eendo3c97 eendo3c99 eendo3c01 eendo3c03 eendo3c05 eendo3c07 eendo3c09 eendo3c11 eendo3c13 eendo3c15 eendo3c17 eendo3c19;
do i=1 to dim(eendo3c);
  if endo3c[i]=3 then do;
    do j=i to dim(eendo3c);
      eendo3c[j]=3;
    end;
  end;
end;

do i=1 to dim(eendo3c);
  if eendo3c[i]=. and endo3c[i]=2 then do;
    do j=i to dim(eendo3c);
      if eendo3c[j]=. then eendo3c[j]=2;
    end;
  end;
end;

do i=1 to dim(eendo3c);
  if eendo3c[i]=. and endo3c[i]=1 then do;
    do j=i to dim(eendo3c);
      if eendo3c[j]=. and endo3c[j]=1 then eendo3c[j]=1;
    end;
  end;
end;

do i=2 to dim(eendo3c);*recode endo status after missing;
  if endo3c[i-1]=. and eendo3c[i]=1 then do;
    do j=i to dim(eendo3c);
      if eendo3c[j]=1 then eendo3c[j]=.;
    end;
  end;
end;

/***3. ever screening endo： 1=no screening endo, 2=screening endo  ***/
array eendo2c [*] eendo2c91    eendo2c93    eendo2c95    eendo2c97    eendo2c99    eendo2c01    eendo2c03    eendo2c05    eendo2c07    eendo2c09   eendo2c11 eendo2c13 eendo2c15 eendo2c17 eendo2c19 ;
do i=1 to dim(eendo2c);
  if endo3c[i]=3 then do;
    do j=i to dim(eendo2c);
      eendo2c[j]=2; *screening endo;
    end;
  end;
end;

do i=1 to dim(eendo2c);
  if eendo2c[i]=. and (endo3c[i]=2 or endo3c[i]=1) then do;
    eendo2c[i]=1; *no screening endo;
  end;
end;

do i=2 to dim(eendo2c); *recode endo status after missing;
  if endo3c[i-1]=. and eendo2c[i]=1 then do;
    do j=i to dim(eendo2c);
      if eendo2c[j]=1 then eendo2c[j]=.;
    end;
  end;
end;
run;

PROC DATASETS;
delete    disease  ahei2010_9119  der8919  act8917 supp8915  meds8919  oc8911
          nur89 nur91 nur93 nur95 nur97 nur99  nur01 nur03 nur05 nur07 nur09 nur11 nur13 nur15 nur17 nur19
          n91_nts n95_nts n99_nts n03_nts n07_nts n11_nts n15_nts n19_nts 
          ;
RUN;

data nhs2;
     set base end=_end_;

/*Siyu, revise the cutoff: NHS1 2020.6, NHS2 2021.6, HPFS 2020.1 */
cutoff=1458;    /*** Set cutoff at 2021.6 ***/
array   irt     {16}    retmo91   retmo93   retmo95   retmo97   retmo99    retmo01   retmo03   retmo05   retmo07   retmo09   
                        retmo11 retmo13 retmo15  retmo17 retmo19 cutoff;
array   tvar    {15}    t91   t93  t95  t97  t99  t01  t03  t05  t07  t09  t11  t13  t15 t17  t19 ;
array   bmi     {15}    bmi91 bmi93 bmi95 bmi97 bmi99 bmi01 bmi03 bmi05 bmi07 bmi09 bmi11 bmi13 bmi15  bmi17 bmi19 ;
array   age     {15}    age91   age93   age95   age97   age99   age01   age03   age05   age07   age09   age11   age13   age15  age17   age19;
array   phms    {15}    mnpst91 mnpst93 mnpst95 mnpst97 mnpst99 mnpst01 mnpst03 mnpst05 mnpst07 mnpst09 mnpst11 mnpst13 mnpst15  mnpst17 mnpst19;
array   hor     {15}    n2rh91    n2rh93    n2rh95    n2rh97    n2rh99     n2rh01    n2rh03    n2rh05    n2rh07    n2rh09    n2rh11   n2rh13   n2rh15 n2rh17   n2rh19 ;
array   smok    {15}    smkdr91   smkdr93   smkdr95   smkdr97   smkdr99    smkdr01   smkdr03   smkdr05   smkdr07   smkdr09   smkdr11   smkdr13  smkdr15 smkdr17  smkdr19;
array   alcon   {15}    alco91n   alco91n   alco95n     alco95n   alco99n   alco99n    alco03n   alco03n   alco07n   alco07n   alco11n   alco11n alco15n alco15n alco19n;
array nSESArr   {15}    nSES_91m nSES_93 nSES_95 nSES_97 nSES_99 nSES_01 nSES_03 nSES_05 nSES_07 nSES_09 nSES_11 nSES_13 nSES_15 nSES_17 nSES_17;

array   act     {15}    act91m    act91m    act91m      act97m    act97m     act01m    act01m    act05m    act05m    act09m    act09m    act13m   act13m  act17m act17m;
array   mvit    {15}    mvitu91  mvitu93   mvitu95    mvitu97   mvitu99  mvitu01    mvitu03   mvitu05   mvitu07   mvitu09  mvitu11   mvitu13  mvitu15  mvitu15 mvitu15;
array   aspi    {15}    aspu91   aspu93    aspu95     aspu97    aspu99   aspu01     aspu03    aspu05    aspu07    aspu09   aspu11    aspu13   aspu15 aspu17 aspu19;
array   ocuse   {15}    nocu91   nocu93    nocu95     nocu97    nocu99   nocu01     nocu03    nocu05    nocu07    nocu09   nocu11    nocu11   nocu11 nocu11 nocu11;
array   hbp     {15}    hbp91  hbp93  hbp95  hbp97  hbp99  hbp01  hbp03  hbp05  hbp07  hbp09  hbp11 hbp13  hbp15 hbp17  hbp19;
array   chol    {15}    chol91 chol93 chol95 chol97 chol99 chol01 chol03 chol05 chol07 chol09 chol11 chol13 chol15 chol17 chol19;
array   diab    {15}    db91   db93   db95   db97   db99   db01   db03   db05   db07   db09   db11 db13  db15 db17  db19;
array   canc    {15}    can91  can93  can95  can97  can99  can01  can03  can05  can07  can09  can11 can13  can15 can17 can19;

array    carbp  {*}   carbo91p  carbo91p carbo95p  carbo95p carbo99p  carbo99p carbo03p  carbo03p carbo07p  carbo07p   carbo11p  carbo11p carbo15p  carbo15p  carbo19p   ;
array    hcarbp  {*}   hcarbo91p  hcarbo91p hcarbo95p  hcarbo95p hcarbo99p  hcarbo99p hcarbo03p  hcarbo03p hcarbo07p  hcarbo07p   hcarbo11p  hcarbo11p hcarbo15p  hcarbo15p  hcarbo19p   ;
array sumref_carbp {*}  sumref_carb91p  sumref_carb91p sumref_carb95p  sumref_carb95p sumref_carb99p  sumref_carb99p sumref_carb03p  sumref_carb03p sumref_carb07p  sumref_carb07p   sumref_carb11p  sumref_carb11p sumref_carb15p  sumref_carb15p  sumref_carb19p   ;
array     protp  {*}   prot91p  prot91p prot95p  prot95p prot99p  prot99p prot03p  prot03p prot07p  prot07p   prot11p  prot11p prot15p  prot15p  prot19p   ;
array    aprotp  {*}   aprot91p  aprot91p aprot95p  aprot95p aprot99p  aprot99p aprot03p  aprot03p aprot07p  aprot07p   aprot11p  aprot11p aprot15p  aprot15p  aprot19p   ;
array    vprotp  {*}   vprot91p  vprot91p vprot95p  vprot95p vprot99p  vprot99p vprot03p  vprot03p vprot07p  vprot07p   vprot11p  vprot11p vprot15p  vprot15p  vprot19p   ;
array    tfatp   {*}   tfat91p  tfat91p tfat95p  tfat95p tfat99p  tfat99p tfat03p  tfat03p tfat07p  tfat07p   tfat11p  tfat11p tfat15p  tfat15p  tfat19p   ;
array    afatp   {*}   afat91p  afat91p afat95p  afat95p afat99p  afat99p afat03p  afat03p afat07p  afat07p   afat11p  afat11p afat15p  afat15p  afat19p   ;
array    vfatp   {*}   vfat91p  vfat91p vfat95p  vfat95p vfat99p  vfat99p vfat03p  vfat03p vfat07p  vfat07p   vfat11p  vfat11p vfat15p  vfat15p  vfat19p   ;

 array eendo3c   {*} eendo3c91 eendo3c93 eendo3c95 eendo3c97 eendo3c99 eendo3c01 eendo3c03 eendo3c05 eendo3c07 eendo3c09 eendo3c11 eendo3c13 eendo3c15 eendo3c17 eendo3c19 ;
 array eendo2c   {*} eendo2c91    eendo2c93    eendo2c95    eendo2c97    eendo2c99    eendo2c01    eendo2c03    eendo2c05    eendo2c07    eendo2c09   eendo2c11 eendo2c13 eendo2c15 eendo2c17 eendo2c19 ;
 array mama      {*} mam91 mam93 mam95 mam97 mam99 mam01 mam03 mam05 mam07 mam09 mam11 mam13 mam15 mam17 mam19 ;

do j=2 to DIM(bmi);
      if age{j}=.      then age{j}=age{j-1}+2;
      if bmi{j}=.      then bmi{j}=bmi{j-1};
      if phms{j}=.      then phms{j}=phms{j-1}; 
      if hor{j}=.      then hor{j}=hor{j-1}; 
      if smok{j}=.     then smok{j}=smok{j-1}; /* smoking status range 0-15 0.missing 1.never 2-8.past 9-15.curr*/
      if alcon{j}=.    then alcon{j}=alcon{j-1};
      if act{j}=.      then act{j}=act{j-1};
      if mvit{j}=.     then mvit{j}=mvit{j-1};
      if aspi{j}=.     then aspi{j}=aspi{j-1};
      if ocuse{j}=.     then ocuse{j}=ocuse{j-1};
      *if marriage{j}=. then marriage{j}=marriage{j-1};
      *if alone{j}=.    then alone{j}=alone{j-1};
      *if antihp{j}=.   then antihp{j}=antihp{j-1};
      *if antitc{j}=.   then antitc{j}=antitc{j-1};
      if canc{j}=.   then canc{j}=canc{j-1};
      if nSESArr{j}=.   then  nSESArr{j}= nSESArr{j-1};
end;

/* -single s-stop updated v-average sv-stop updated average*/
/* -s */
array calors {15} calor91ns calor91ns calor95ns calor95ns calor99ns calor99ns calor03ns calor03ns calor07ns calor07ns calor11ns calor11ns calor15ns calor15ns calor19ns;
array aheims {15} ahei91s ahei91s ahei95s ahei95s ahei99s ahei99s ahei03s ahei03s ahei07s ahei07s ahei11s ahei11s ahei15s ahei15s ahei19s; 
array lcss {15} lcs91s lcs91s lcs95s lcs95s lcs99s lcs99s lcs03s lcs03s lcs07s lcs07s lcs11s lcs11s lcs15s lcs15s lcs19s;
array alcss {15} alcs91s alcs91s alcs95s alcs95s alcs99s alcs99s alcs03s alcs03s alcs07s alcs07s alcs11s alcs11s alcs15s alcs15s alcs19s;
array vlcss {15} vlcs91s vlcs91s vlcs95s vlcs95s vlcs99s vlcs99s vlcs03s vlcs03s vlcs07s vlcs07s vlcs11s vlcs11s vlcs15s vlcs15s vlcs19s;
array ulcss {15} ulcs91s ulcs91s ulcs95s ulcs95s ulcs99s ulcs99s ulcs03s ulcs03s ulcs07s ulcs07s ulcs11s ulcs11s ulcs15s ulcs15s ulcs19s; 
array hlcss {15} hlcs91s hlcs91s hlcs95s hlcs95s hlcs99s hlcs99s hlcs03s hlcs03s hlcs07s hlcs07s hlcs11s hlcs11s hlcs15s hlcs15s hlcs19s;
array lfss {15} lfs91s lfs91s lfs95s lfs95s lfs99s lfs99s lfs03s lfs03s lfs07s lfs07s lfs11s lfs11s lfs15s lfs15s lfs19s;
array ulfss {15} ulfs91s ulfs91s ulfs95s ulfs95s ulfs99s ulfs99s ulfs03s ulfs03s ulfs07s ulfs07s ulfs11s ulfs11s ulfs15s ulfs15s ulfs19s;
array hlfss {15} hlfs91s hlfs91s hlfs95s hlfs95s hlfs99s hlfs99s hlfs03s hlfs03s hlfs07s hlfs07s hlfs11s hlfs11s hlfs15s hlfs15s hlfs19s;
array alfss {15} alfs91s alfs91s alfs95s alfs95s alfs99s alfs99s alfs03s alfs03s alfs07s alfs07s alfs11s alfs11s alfs15s alfs15s alfs19s;
array vlfss {15} vlfs91s vlfs91s vlfs95s vlfs95s vlfs99s vlfs99s vlfs03s vlfs03s vlfs07s vlfs07s vlfs11s vlfs11s vlfs15s vlfs15s vlfs19s;

/* -sv */
array calorsv {15} calor91nsv calor91nsv calor95nsv calor95nsv calor99nsv calor99nsv calor03nsv calor03nsv calor07nsv calor07nsv calor11nsv calor11nsv calor15nsv calor15nsv calor19nsv;
array aheimsv {15} ahei91sv ahei91sv ahei95sv ahei95sv ahei99sv ahei99sv ahei03sv ahei03sv ahei07sv ahei07sv ahei11sv ahei11sv ahei15sv ahei15sv ahei19sv; 
array lcssv {15} lcs91sv lcs91sv lcs95sv lcs95sv lcs99sv lcs99sv lcs03sv lcs03sv lcs07sv lcs07sv lcs11sv lcs11sv lcs15sv lcs15sv lcs19sv;
array alcssv {15} alcs91sv alcs91sv alcs95sv alcs95sv alcs99sv alcs99sv alcs03sv alcs03sv alcs07sv alcs07sv alcs11sv alcs11sv alcs15sv alcs15sv alcs19sv;
array vlcssv {15} vlcs91sv vlcs91sv vlcs95sv vlcs95sv vlcs99sv vlcs99sv vlcs03sv vlcs03sv vlcs07sv vlcs07sv vlcs11sv vlcs11sv vlcs15sv vlcs15sv vlcs19sv;
array ulcssv {15} ulcs91sv ulcs91sv ulcs95sv ulcs95sv ulcs99sv ulcs99sv ulcs03sv ulcs03sv ulcs07sv ulcs07sv ulcs11sv ulcs11sv ulcs15sv ulcs15sv ulcs19sv;
array hlcssv {15} hlcs91sv hlcs91sv hlcs95sv hlcs95sv hlcs99sv hlcs99sv hlcs03sv hlcs03sv hlcs07sv hlcs07sv hlcs11sv hlcs11sv hlcs15sv hlcs15sv hlcs19sv;
array lfssv {15} lfs91sv lfs91sv lfs95sv lfs95sv lfs99sv lfs99sv lfs03sv lfs03sv lfs07sv lfs07sv lfs11sv lfs11sv lfs15sv lfs15sv lfs19sv;
array ulfssv {15} ulfs91sv ulfs91sv ulfs95sv ulfs95sv ulfs99sv ulfs99sv ulfs03sv ulfs03sv ulfs07sv ulfs07sv ulfs11sv ulfs11sv ulfs15sv ulfs15sv ulfs19sv;
array hlfssv {15} hlfs91sv hlfs91sv hlfs95sv hlfs95sv hlfs99sv hlfs99sv hlfs03sv hlfs03sv hlfs07sv hlfs07sv hlfs11sv hlfs11sv hlfs15sv hlfs15sv hlfs19sv;
array alfssv {15} alfs91sv alfs91sv alfs95sv alfs95sv alfs99sv alfs99sv alfs03sv alfs03sv alfs07sv alfs07sv alfs11sv alfs11sv alfs15sv alfs15sv alfs19sv;
array vlfssv {15} vlfs91sv vlfs91sv vlfs95sv vlfs95sv vlfs99sv vlfs99sv vlfs03sv vlfs03sv vlfs07sv vlfs07sv vlfs11sv vlfs11sv vlfs15sv vlfs15sv vlfs19sv;


array whgrn_sv {15}  whgrn_s91sv whgrn_s91sv whgrn_s95sv whgrn_s95sv whgrn_s99sv whgrn_s99sv whgrn_s03sv whgrn_s03sv whgrn_s07sv whgrn_s07sv whgrn_s11sv whgrn_s11sv whgrn_s15sv whgrn_s15sv whgrn_s19sv;
array rgrain_sv{15}  rgrain_s91sv rgrain_s91sv rgrain_s95sv rgrain_s95sv rgrain_s99sv rgrain_s99sv rgrain_s03sv rgrain_s03sv rgrain_s07sv rgrain_s07sv rgrain_s11sv rgrain_s11sv rgrain_s15sv rgrain_s15sv rgrain_s19sv;   
array frui_sv  {15}  frui_s91sv    frui_s91sv    frui_s95sv    frui_s95sv    frui_s99sv    frui_s99sv    frui_s03sv    frui_s03sv    frui_s07sv    frui_s07sv    frui_s11sv    frui_s11sv    frui_s15sv frui_s15sv frui_s19sv;
array leg_sv   {15}  leg_s91sv   leg_s91sv   leg_s95sv   leg_s95sv   leg_s99sv   leg_s99sv   leg_s03sv   leg_s03sv   leg_s07sv   leg_s07sv   leg_s11sv   leg_s11sv   leg_s15sv leg_s15sv leg_s19sv;
array pota_sv  {15}  pota_s91sv   pota_s91sv   pota_s95sv   pota_s95sv   pota_s99sv   pota_s99sv   pota_s03sv   pota_s03sv   pota_s07sv   pota_s07sv   pota_s11sv   pota_s11sv   pota_s15sv pota_s15sv pota_s19sv;
array veg_noleg_sv {15}  veg_noleg_s91sv   veg_noleg_s91sv   veg_noleg_s95sv   veg_noleg_s95sv   veg_noleg_s99sv   veg_noleg_s99sv   veg_noleg_s03sv   veg_noleg_s03sv   veg_noleg_s07sv   veg_noleg_s07sv   veg_noleg_s11sv   veg_noleg_s11sv   veg_noleg_s15sv veg_noleg_s15sv veg_noleg_s19sv;
array fish_sv  {15}  fish_s91sv   fish_s91sv   fish_s95sv   fish_s95sv   fish_s99sv   fish_s99sv   fish_s03sv   fish_s03sv   fish_s07sv   fish_s07sv   fish_s11sv   fish_s11sv   fish_s15sv fish_s15sv fish_s19sv;
array rpmeats_sv   {15}  rpmeats_s91sv    rpmeats_s91sv    rpmeats_s95sv    rpmeats_s95sv    rpmeats_s99sv    rpmeats_s99sv    rpmeats_s03sv    rpmeats_s03sv    rpmeats_s07sv    rpmeats_s07sv    rpmeats_s11sv    rpmeats_s11sv    rpmeats_s15sv rpmeats_s15sv rpmeats_s19sv;
array ssb_sv  {15}  ssb_s91sv   ssb_s91sv   ssb_s95sv   ssb_s95sv   ssb_s99sv   ssb_s99sv   ssb_s03sv   ssb_s03sv   ssb_s07sv   ssb_s07sv   ssb_s11sv   ssb_s11sv   ssb_s15sv ssb_s15sv ssb_s19sv;
array dairy_sv{15}  dairy_s91sv   dairy_s91sv   dairy_s95sv   dairy_s95sv   dairy_s99sv   dairy_s99sv   dairy_s03sv   dairy_s03sv   dairy_s07sv   dairy_s07sv   dairy_s11sv   dairy_s11sv   dairy_s15sv dairy_s15sv dairy_s19sv;
array swtsnochoc_sv{15} swtsnochoc_s91sv swtsnochoc_s91sv swtsnochoc_s95sv swtsnochoc_s95sv swtsnochoc_s99sv swtsnochoc_s99sv swtsnochoc_s03sv swtsnochoc_s03sv swtsnochoc_s07sv swtsnochoc_s07sv swtsnochoc_s11sv swtsnochoc_s11sv swtsnochoc_s15sv swtsnochoc_s15sv swtsnochoc_s19sv ;
array poul_sv {15} poul_s91sv poul_s91sv poul_s95sv poul_s95sv poul_s99sv poul_s99sv poul_s03sv poul_s03sv poul_s07sv poul_s07sv poul_s11sv poul_s11sv poul_s15sv poul_s15sv poul_s19sv ;


array carbsv {15} carbo91nsv carbo91nsv carbo95nsv carbo95nsv carbo99nsv carbo99nsv carbo03nsv carbo03nsv carbo07nsv carbo07nsv carbo11nsv carbo11nsv carbo15nsv carbo15nsv carbo19nsv;
array hcarbsv {15} hcarbo91sv hcarbo91sv hcarbo95sv hcarbo95sv hcarbo99sv hcarbo99sv hcarbo03sv hcarbo03sv hcarbo07sv hcarbo07sv hcarbo11sv hcarbo11sv hcarbo15sv hcarbo15sv hcarbo19sv; 
array sumref_carbsv {15} sumref_carb91sv sumref_carb91sv sumref_carb95sv sumref_carb95sv sumref_carb99sv sumref_carb99sv sumref_carb03sv sumref_carb03sv sumref_carb07sv sumref_carb07sv sumref_carb11sv sumref_carb11sv sumref_carb15sv sumref_carb15sv sumref_carb19sv;
array protsv {15} prot91sv prot91sv prot95sv prot95sv prot99sv prot99sv prot03sv prot03sv prot07sv prot07sv prot11sv prot11sv prot15sv prot15sv prot19sv;
array aprotsv {15} aprot91sv aprot91sv aprot95sv aprot95sv aprot99sv aprot99sv aprot03sv aprot03sv aprot07sv aprot07sv aprot11sv aprot11sv aprot15sv aprot15sv aprot19sv;
array vprotsv {15} vprot91sv vprot91sv vprot95sv vprot95sv vprot99sv vprot99sv vprot03sv vprot03sv vprot07sv vprot07sv vprot11sv vprot11sv vprot15sv vprot15sv vprot19sv;
array tfatsv {15} tfat91sv tfat91sv tfat95sv tfat95sv tfat99sv tfat99sv tfat03sv tfat03sv tfat07sv tfat07sv tfat11sv tfat11sv tfat15sv tfat15sv tfat19sv;
array afatsv {15} afat91sv afat91sv afat95sv afat95sv afat99sv afat99sv afat03sv afat03sv afat07sv afat07sv afat11sv afat11sv afat15sv afat15sv afat19sv;
array vfatsv {15} vfat91sv vfat91sv vfat95sv vfat95sv vfat99sv vfat99sv vfat03sv vfat03sv vfat07sv vfat07sv vfat11sv vfat11sv vfat15sv vfat15sv vfat19sv;


/*************************************************************************
***** If an irt date is before June of that qq year or after or equal ****
***** to the next qq year it is incorrect and should be defaulted to  ****
***** June of that qq year.    Make time period indicator tvar=0.     ****
*************************************************************************/

/* From 1991: create ltf using the original irt values */
ltf=.;
do i=1 to dim(irt)-1; if irt{i}>0 then ltf=irt{i}; end;
if ltf=irt{dim(irt)-1} then ltf=.;
/* Correct irt values */
do i=1 to dim(irt)-1;
    if (irt{i} <  (1074+24*i)) or (irt{i} >= (1098+24*i))  then irt{i}=1074+24*i;
end;

/* INCIDENT CANCER MACRO */
%macro make_incident_cancers(sitelist);
   %local k site dvar tvar evar;
   %do k = 1 %to %sysfunc(countw(&sitelist));
      %let site = %scan(&sitelist, &k);
      %let dvar = dtdx_&site; 
      %let evar = can_&site;
      %let tvar = t&site;

      /* initialize event/time for interval */
      &evar = 0; 
      &tvar = irt{i+1} - irt{i};

      if &dvar > 0 then do;
         if irt{i} < &dvar <= irt{i+1} then do;
            &evar = 1;
            &tvar = &dvar - irt{i};
         end;
      end;

      if irt{i} < dtdth < irt{i+1} then &tvar = min(&tvar, dtdth - irt{i});
   %end;
%mend make_incident_cancers;


%beginex();
/*************************************************************/
/*                DO-LOOP OVER TIME PERIODS                  */

   
    do i=1 to DIM(irt)-1;     
      
    period=i;
      do j=1 to DIM(tvar);
      tvar{j}=0;
      end;
    tvar{i}=1; 


/****** OUTCOME ******/
    /*** death ***/
    dead = 0; tdead = irt{i+1} - irt{i};
    if irt{i} < dtdth <= irt{i+1} then do; dead = 1; tdead = dtdth - irt{i}; end;
    /* Cause-specific death */
    dead_can   = 0; tdead_can=irt{i+1}-irt{i}; 
    if irt{i}< dtdth <=irt{i+1} and (140<=newicda3<=208) then do; dead_can=1; tdead_can=dtdth-irt{i}; end;

    dead_cvd   = 0; tdead_cvd=irt{i+1}-irt{i};
    if irt{i}< dtdth <=irt{i+1} and (390<=newicda3<=458) then do; dead_cvd=1; tdead_cvd=dtdth-irt{i}; end;

    dead_resp=0;  tdead_resp=irt{i+1}-irt{i};        
    if irt{i}<dtdth<=irt{i+1} and (460<=newicda3<520) then do; dead_resp=1; tdead_resp=dtdth-irt{i}; end;

    dead_oth=0;tdead_oth=irt{i+1}-irt{i}; 
    if irt{i}<dtdth<=irt{i+1} and (dead_can=0 and dead_cvd=0 and dead_resp=0 ) then do; dead_oth=1; tdead_oth=dtdth-irt{i}; end;

    /*** CVD, CHD ***/
    cvd = 0; tcvd = irt{i+1} - irt{i};
    if irt{i} < dtdxcvd <= irt{i+1} then do; cvd = 1; tcvd = dtdxcvd - irt{i}; end;
    if dtdth > irt{i} and dtdth < irt{i+1} then tcvd = min(tcvd, dtdth - irt{i});

    chd = 0; tchd = irt{i+1} - irt{i};
    if irt{i} < dtdxchd <= irt{i+1} then do; chd = 1; tchd = dtdxchd - irt{i}; end;
    if dtdth > irt{i} and dtdth < irt{i+1} then tchd = min(tchd, dtdth - irt{i});

    str = 0; tstr = irt{i+1} - irt{i};
    if irt{i} < dtdxstr <= irt{i+1} then do; str = 1; tstr = dtdxstr - irt{i}; end;
    if dtdth > irt{i} and dtdth < irt{i+1} then tstr = min(tstr, dtdth - irt{i});
    
    /*** obesity-related cancer ***/
    can_obese = 0; tobeca = irt{i+1} - irt{i};
    if irt{i} < dtdx_obeca <= irt{i+1} then do; can_obese = 1; tobeca = dtdx_obeca - irt{i}; end;
    if dtdth > irt{i} and dtdth < irt{i+1} then tobeca = min(tobeca, dtdth - irt{i});

    can_nonobese = 0; tnonobeca = irt{i+1} - irt{i};
    if irt{i} < dtdx_nonobeca <= irt{i+1} then do; can_nonobese = 1; tnonobeca = dtdx_nonobeca - irt{i}; end;
    if dtdth > irt{i} and dtdth < irt{i+1} then tnonobeca = min(tnonobeca, dtdth - irt{i});

    /*********************************************************************/
    /* Categorical subtype outcome for %subtype macro                   */
    /* 0 = no obesity-related cancer event in this interval             */
    /* 1-13 = first obesity-related cancer subtype                      */
    /*********************************************************************/
    can_obese_site = 0;
    if can_obese = 1 then do;
        can_obese_site = obeca_site;
    end;

    /*** Site-specific cancer ***/
    %make_incident_cancers(
    /* Individual sites */
    oral esoph stomach smallintest colorect liver gallbladder panc
    nose larynx lung 
    bone melanoma brca cervical endo ovca bladder kidney
    brain thyroid myeloma myeleuk menin

    obeca_nobcop  obeca_nobc othsite

    /* Grouped sites */
    diges resp uri                            
    smoke_ca r_smoke_ca                            
    allca /*obeca nonobeca*/                     
);

**** create follow-up time (cancer-free to dx);
end_cal = min(irt{i+1}, cutoff, ltf);
fup_allca = 0;
if end_cal ne . and end_cal > irt{i} then do;

      /* time to event within interval */
      fup_allca = end_cal - irt{i};
      /* Censor at death if occurs before end_cal */
      if 0 < dtdth <= end_cal then fup_allca = min(fup_allca, dtdth - irt{i});
      /* Censor at incident cancer diagnosis (ALL cancer example) */
      if 0 < dtdx_allca <= end_cal then fup_allca = min(fup_allca, dtdx_allca - irt{i});
      fup_allca = max(0, fup_allca);
    end;
   /* Ensure no contribution after diagnosis date (even if later intervals exist) */
    if 0 < dtdx_allca and irt{i} >= dtdx_allca then fup_allca = 0;


/******EXPOSURE******/
/* -s */
calorscon=calors{i};
aheimscon=aheims{i};
lcsscon=lcss{i};
alcsscon=alcss{i};
vlcsscon=vlcss{i};
ulcsscon=ulcss{i};
hlcsscon=hlcss{i};
lfsscon=lfss{i};
ulfsscon=ulfss{i};
hlfsscon=hlfss{i};
alfsscon=alfss{i};
vlfsscon=vlfss{i};

carb=carbp {i};
hcarb=hcarbp {i};
sumref_carb=sumref_carbp {i};
prot=protp {i};
aprot=aprotp {i};
vprot=vprotp {i};
tfat=tfatp {i};
afat=afatp {i};
vfat=vfatp {i};
basecarb=carbp {1};
basehcarb=hcarbp {1};
basesumref_carb=sumref_carbp {1};
baseprot=protp {1};
baseaprot=aprotp {1};
basevprot=vprotp {1};
basetfat=tfatp {1};
baseafat=afatp {1};
basevfat=vfatp {1};
    carb_svcon=carbsv {i};
    hcarb_svcon=hcarbsv {i};
    sumref_carb_svcon=sumref_carbsv {i};
    prot_svcon=protsv {i};
    aprot_svcon=aprotsv {i};
    vprot_svcon=vprotsv {i};
    tfat_svcon=tfatsv {i};
    afat_svcon=afatsv {i};
    vfat_svcon=vfatsv {i};

whgrn_svcon=whgrn_sv {i};
rgrain_svcon=rgrain_sv {i};
frui_svcon=frui_sv {i};
leg_svcon=leg_sv {i};
pota_svcon=pota_sv {i};
veg_noleg_svcon=veg_noleg_sv {i};
fish_svcon=fish_sv {i};
rpmeats_svcon=rpmeats_sv {i};
ssb_svcon=ssb_sv {i};
dairy_svcon=dairy_sv {i};
swtsnochoc_svcon=swtsnochoc_sv {i};
poul_svcon=poul_sv {i};

/* -sv */
calorsvcon=calorsv{i};
aheimsvcon=aheimsv{i};
lcssvcon=lcssv{i};
alcssvcon=alcssv{i};
vlcssvcon=vlcssv{i};
ulcssvcon=ulcssv{i};
hlcssvcon=hlcssv{i};
lfssvcon=lfssv{i};
ulfssvcon=ulfssv{i};
hlfssvcon=hlfssv{i};
alfssvcon =alfssv{i};
vlfssvcon =vlfssv{i};

/* AGE GROUP */
     agecon=age{i};
     agegp=int((age{i}-20)/5);     **** Define the agegp in the i-th period;
     if agegp<=0 then agegp=1;
     if agegp>6 then agegp=6;

     %indic3(vbl=agegp, prefix=agegp, reflev=1, min=2, max=6, missing=., usemiss=0,
                label1='age<30',
                label2='age30-34',
                label3='age35-39',
                label4='age40-44',
                label5='age45-50',
                label6='age > 50');
     if agecon>65 then agesub=1; else agesub=0;

/* menopausal status and PMH use */

        if      phms{i}=1         then phmsstatus=1; /*pre menopause*/
        else if hor{i}=1          then phmsstatus=2; /*never use*/
        else if hor{i}=3          then phmsstatus=3; /*current use*/
        else if hor{i} in (2,4)   then phmsstatus=4; /*past use*/
        else                           phmsstatus=1;

        %indic3(vbl=phmsstatus, prefix=phmsstatus, min=2, max=4, reflev=1,
                missing=9, usemiss=0,
                label1='premenopause',
                label2='pmh never use',
                label3='pmh current use',
                label4='pmh past use');

        if phmsstatus=1 then postpau=0; else postpau=1;
        if phmsstatus in (3,4) then pmhever=1; else pmhever=0;

/***smoking status***/
        if smok{i}=1 then smkk=1;  *never;
        else if smok{i} in(2,3,4,5,6,7,8) then smkk=2; *past smoker;
        else if smok{i} in (9,10) then smkk=3; *current, 1-14 cigs ;
        else if smok{i}=11 then smkk=4; *current, 15-24 cigs ;
        else if smok{i} in(12,13,14) then smkk=5; *current, 25+ cigs ;
        else smkk=1; 

     %indic3(vbl=smkk, prefix=smkk, min=2,max=5, reflev=1, missing=., usemiss=0,
              label1='never smker',
              label2='past smker',
              label3='curr 1-14 smker',
              label4='curr 15+ smker',
              label5='curr 25+ smker');

      if smkk=2 then smkstatus=2;
      else if smkk>=3 then smkstatus=3;
      else smkstatus=1;
      %indic3(vbl=smkstatus, prefix=smkstatus, min=2, max=3, reflev=1, missing=., usemiss=0);

      if smkstatus=1  then smknever=1; else smknever=0;

 /***category of exercise***/
       actcon=act{i};
       
       if 0 <= actcon and actcon < 3 then actc=1;
       else if actcon>=3 and actcon<9 then actc=2;
       else if actcon>=9 and actcon<18 then actc=3;
       else if actcon>=18 and actcon<27 then actc=4;
       else if actcon>=27 then actc=5;
       else  actc=3;

%indic3(vbl=actc, prefix=actc, min=2, max=5, reflev=1, missing=., usemiss=0,
        label1='<3 total mets/wk',
        label2='3 to < 9 mets/wk',
        label3='9 to < 18 mets/wk',
        label4='18 to < 27 mets/wk',
        label5='27+ mets/wk');

if actcon>7.5 then active=1; else active=0;

/****** Indicator Alcohol consumption ******/
       alcocon=alcon{i}; /***continuous variable, g/day***/

            alc=1;
            if alcocon=0.0                          then  alc=1;
            else if alcocon>0.0  and  alcocon<5.0   then  alc=2;
            else if alcocon>=5.0 and  alcocon<10.0  then  alc=3;
            else if alcocon>=10.0 and alcocon<15.0  then  alc=4;
            else if alcocon>=15.0 and alcocon<30.0  then  alc=5;
            else if alcocon>=30.0                   then  alc=6;

       %indic3(vbl=alc, prefix=alc, min=2, max=6, reflev=1, usemiss=0,
                label1='0 g/d',
                label2='0.1-4.9 g/d',
                label3='5.0-9.9 g/d',
                label4='10.0-14.9 g/d',
                label5='15.0-29.9 g/d',
                label6='30+ g/d');

        if alc in (2,3,4,5,6) then neverdrinker=0; else neverdrinker=1;

/*** Indicator for BMI ***/
         bmicon=bmi{i};

         if      25<= bmicon and bmicon < 30 then bmic3cat=2;
         else if 30<= bmicon     then bmic3cat=3;
         else                         bmic3cat=1;
         %indic3(vbl= bmic3cat, prefix= bmic3cat, min=2, max=3, reflev=1, missing=., usemiss=0);

         if bmicon>=25 then bmic25=1;  else  bmic25=0;
         if bmicon>=30 then bmic30=1;  else  bmic30=0;

        if       0< bmicon<21    then bmic9=1;
        else if 21<=bmicon<23    then bmic9=2;
        else if 23<=bmicon<25    then bmic9=3;
        else if 25<=bmicon<27    then bmic9=4;
        else if 27<=bmicon<29    then bmic9=5;
        else if 29<=bmicon<31    then bmic9=6;
        else if 31<=bmicon<33    then bmic9=7;
        else if 33<=bmicon<35    then bmic9=8;
        else if     bmicon>=35   then bmic9=9;
        else bmic9=.;
        %indic3(vbl=bmic9, prefix=bmic9, min=2, max=9, reflev=1, usemiss=1);


         bmibcon=bmi{1};

         if      25<= bmibcon <30 then bmib3cat=2;
         else if 30<= bmibcon     then bmib3cat=3;
         else                          bmib3cat=1;
         %indic3(vbl= bmib3cat, prefix= bmib3cat, min=2, max=3, reflev=1, missing=., usemiss=0);

        if       0< bmibcon<21    then bmib9=1;
        else if 21<=bmibcon<23    then bmib9=2;
        else if 23<=bmibcon<25    then bmib9=3;
        else if 25<=bmibcon<27    then bmib9=4;
        else if 27<=bmibcon<29    then bmib9=5;
        else if 29<=bmibcon<31    then bmib9=6;
        else if 31<=bmibcon<33    then bmib9=7;
        else if 33<=bmibcon<35    then bmib9=8;
        else if     bmibcon>=35   then bmib9=9;
        else bmib9=.;

        %indic3(vbl=bmib9, prefix=bmib9, min=2, max=9, reflev=1, usemiss=1);


/*** Indicator for OC use ******/
        if ocuse{i}=1 then ocu=1;
        else if ocuse{i}=2 then ocu=2;
        else if ocuse{i}=3 then ocu=3;
        else                 ocu=1;

%indic3(vbl=ocu, prefix=ocu, min=2, max=3, reflev=1,
                missing=9, usemiss=0,
                label1='never OC',
                label2='past OC',
                label3='current OC');

        ocever=ocuse{i};
        if ocever=2 or ocever=3 then ocever=1;else ocever=0;

/*** aspirin ***/
            select(aspi{i});
                when (1)     asp=1;
                otherwise    asp=0;
            end;

/*** multiple vitamin ***/
             select(mvit{i});
                when (1) mv=1;
                otherwise  mv=0;
            end;

/***family history ***/
       if mifh ne 1 then mifh=0;
       if dbfh ne 1 then dbfh=0;
       if cafh ne 1 then cafh=0;
       
/****** Indicator for History of antihypertensive user ******
            select(antihp{i});
                when (1)     antihbp=1;
                otherwise    antihbp=0;
            end; */   
            if hbp{i}=1 then highbp=1; else highbp=0;
            if hbp{1}=1 then basehighbp=1; else basehighbp=0;

/****** Indicator for TC lowering ******
            select(antitc{i});
                when (1)     antihtc=1;
                otherwise    antihtc=0;
            end; */
            if chol{i}=1 then hightc=1; else hightc=0;            
            if chol{1}=1 then basehightc=1; else basehightc=0;

/****** Indicator for History of diabetes *******/
          if diab{i}=1 then diabetes=1; else diabetes=0;

/****** Indicator for History of cancer ******/
          if canc{i}=1 then cancer=1; else cancer=0;

/****** endoscopy for screening purpose ******/
    endo2c=eendo2c{i}; /* ever screening endo, 1=no screening endo, 2=screening endo */
    endo3c=eendo3c{i}; /* ever endo3c, 1=no endo, 2=endo for symptoms, 3=screening endo */
    if endo2c = . then endo2c = 1;
/****** Mammogram in last 2 yrs *******/
   mam=mama{i};                     
   if mam = . then mam = 0;

/** nSES ***/
nSES=nSESArr{i};

 /* baseline food*/
basewhgrn=whgrn_s91;
basergrain=rgrain_s91;
basefrui=frui_s91;
baseleg=leg_s91;
basepota=pota_s91;
baseveg_noleg=veg_noleg_s91;
basefish=fish_s91;
baserpmeats=rpmeats_s91;
basessb=ssb_s91;
basedairy=dairy_s91;
basesweets=swtsnochoc_s91;
basepoultry=poul_s91;

/****************  BASELINE EXCLUSIONS ********************************************/

if i=1 then do;
   %exclude(exrec eq 1);                 *multiple records and not in master file;
   %exclude(baseffq eq 1);               * no baseline ffq;
   %exclude(lcs91 eq .);
   %exclude(alcs91 eq .);
   %exclude(vlcs91 eq .);
   %exclude(ulcs91 eq .);
   %exclude(hlcs91 eq .);
   %exclude(lfs91 eq .);
   %exclude(ulfs91 eq .);
   %exclude(hlfs91 eq .);

   %exclude(ltf eq irt{1});

   %exclude(0 lt dtdth   le irt{i});
   %exclude(0 lt dtdx_allca  le irt{i}); 
   %exclude(can91   eq 1);

   %exclude(0 lt dtdxcvd le irt{i}); 
   %exclude(mi91   eq 1);
   %exclude(str91   eq 1);

   %exclude(0 lt dtdxdb2 le irt{i});
   %exclude(db91    eq 1); 

   %exclude(birthday eq .);
   %exclude(age91 eq .);

   %exclude(calor91n eq .);
   %exclude(calor91n < 500);
   %exclude(calor91n > 3500);
   %exclude(ahei91 eq .);

   %output();
  end;

/* EXCLUSIONS DURING FOLLOW-UP */
   if i>1 then do;

   %exclude(irt{i-1} lt dtdth le irt{i} );       *censor death; 
   %exclude(irt{i-1} lt dtdx_obeca le irt{i} );  *censor obesity-related cancer: stop follow-up after obesity-related cancer;
   %output();
  end;
end;        
/* END OF DO-LOOP OVER TIME PERIODs */
%endex();

proc sort data=nhs2; by period; run;

proc rank data=nhs2 out=nhs2 group=5;
by period;
var 
calorscon  calorsvcon 
aheimscon  aheimsvcon
lcsscon    lcssvcon
alcsscon   alcssvcon
vlcsscon   vlcssvcon
ulcsscon   ulcssvcon
hlcsscon   hlcssvcon
lfsscon    lfssvcon
ulfsscon   ulfssvcon
hlfsscon   hlfssvcon
alfsscon   alfssvcon
vlfsscon   vlfssvcon
;

ranks
calorsq  calorsvq
aheimsq  aheimsvq
lcssq    lcssvq
alcssq   alcssvq
vlcssq   vlcssvq
ulcssq   ulcssvq
hlcssq   hlcssvq
lfssq    lfssvq
ulfssq   ulfssvq
hlfssq   hlfssvq
alfssq   alfssvq
vlfssq   vlfssvq
;
run;

data nhs2;
set  nhs2 end=_end_;

array qvar 
calorsq  calorsvq
aheimsq  aheimsvq
lcssq    lcssvq
alcssq   alcssvq
vlcssq   vlcssvq
ulcssq   ulcssvq
hlcssq   hlcssvq
lfssq    lfssvq
ulfssq   ulfssvq
hlfssq   hlfssvq
alfssq   alfssvq
vlfssq   vlfssvq
;

do over qvar;
  qvar=qvar+1;
  if qvar=. then qvar=1; 
end;

keep id period 
     /**Time scale**/
     retmo91   retmo93   retmo95   retmo97   retmo99    retmo01   retmo03   retmo05   retmo07   retmo09   
     retmo11 retmo13 retmo15 retmo17 retmo19 cutoff
     t91   t93  t95  t97  t99  t01  t03  t05  t07  t09  t11  t13  t15 t17  t19

     dtdxchd dtdxmi dtdth dtdxca  dtdx_allca dtdx_obeca  dtdx_nonobeca fup_allca
     /** Outcome **/
     dead tdead
     dead_can dead_cvd dead_resp dead_oth tdead_can tdead_cvd tdead_resp tdead_oth 
     can_obese tobeca can_nonobese tnonobeca 
     cvd tcvd chd tchd str tstr
     /** Site-specific cancer **/
     can_: t: dtdx_:    obeca_site can_obese_site     

     /** Exposures (stop-updated, stop-updated cumulative) **/
     lcsscon    lcssvcon   alcsscon   alcssvcon  vlcsscon   vlcssvcon  ulcsscon   ulcssvcon  hlcsscon   hlcssvcon
     lfsscon    lfssvcon   ulfsscon   ulfssvcon  hlfsscon   hlfssvcon  alfsscon   alfssvcon  vlfsscon   vlfssvcon

     calorsq  calorsvq aheimsq  aheimsvq
     lcssq    lcssvq alcssq   alcssvq vlcssq   vlcssvq ulcssq   ulcssvq hlcssq   hlcssvq
     lfssq    lfssvq ulfssq   ulfssvq hlfssq   hlfssvq alfssq   alfssvq vlfssq   vlfssvq

     basecarb basehcarb basesumref_carb baseprot baseaprot basevprot basetfat baseafat basevfat
     carb hcarb sumref_carb prot aprot vprot tfat afat vfat
     carb_svcon hcarb_svcon sumref_carb_svcon prot_svcon aprot_svcon vprot_svcon tfat_svcon afat_svcon vfat_svcon

     /* === Covariates === */
     calorscon  calorsvcon  aheimscon  aheimsvcon
     white dbfh mifh cafh basehighbp basehightc highbp hightc diabetes mv asp cancer
     agecon agegp  agesub  bmicon alcocon actcon
   
     bmic3cat &bmic3cat_ bmib3cat &bmib3cat_ bmic25 bmic30 bmic9 &bmic9_ bmib9 &bmib9_ 
     smkstatus &smkstatus_ smknever 
     alc  &alc_ neverdrinker  actc &actc_ active phmsstatus &phmsstatus_ postpau pmhever ocever
     basewhgrn basergrain basefrui  baseleg basepota baseveg_noleg  basefish baserpmeats  basessb basedairy basesweets basepoultry 
     whgrn_svcon rgrain_svcon frui_svcon leg_svcon pota_svcon veg_noleg_svcon fish_svcon rpmeats_svcon ssb_svcon 
     dairy_svcon swtsnochoc_svcon poul_svcon 
     mam endo2c endo3c 
     nSES
; 
run;

/*
proc means data=nhs2 mean std median n nmiss min max;
var 
calorscon  calorsvcon  aheimscon  aheimsvcon agecon bmicon alcocon actcon 

lcsscon    lcssvcon   alcsscon   alcssvcon  vlcsscon   vlcssvcon  ulcsscon   ulcssvcon  hlcsscon   hlcssvcon
lfsscon    lfssvcon   ulfsscon   ulfssvcon  hlfsscon   hlfssvcon;
run;

proc freq data=nhs1;
tables 
agegp white bmib3cat smkstatus alc actc mifh basehighbp basehightc mv asp phmsstatus ocever
;
run;
*/