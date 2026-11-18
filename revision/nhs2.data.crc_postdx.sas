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
               stage, tumor differentiation, location.
revision:      add SES to covariates 
Date:  10/2026

******************************************************/

/**** Cancer *****/
/*All cancer: 1989-2023*/

data allcancer;
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
      colorect_ca=0; if icda2 in (153,154) and ca=1 then do; colorect_ca=1;   dtdx_colorect  = dtdxca; end;
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
     proc sort data=allcancer; by id dtdxca; run;


     data allcancer;
     set allcancer end=_end_;

     dtdx_obeca=.;
     if obese_ca=1 and dtdxca>0 and ca=1 and icda2 in (150,151,153,154,155,156,157,174,182,183,189,193,203,225) then dtdx_obeca=dtdxca;
     dtdx_nonobeca=.;
     if nonobeca=1 and dtdxca>0 and ca=1  then dtdx_nonobeca=dtdxca;
  
     keep id conf dtdxca colorect_ca brca_ca dtdx_colorect dtdx_brca 
          /*obese_ca nonobeca smoke_ca r_smoke_ca  
          esoph_ca colorect_ca stomach_ca liver_ca gallbladder_ca panc_ca kidney_ca thyroid_ca myeloma_ca menin_ca */
          ;
  run;



/************** Colorectal Cancer **************/

data crc; *colc.cases.dd;
  %include 'colc8921.021422.cases.input'; /* use this old version */
run;

/******************************
id              1-6   id number for NHS2 
icda            7-9      icd code  
nhsicda        36-40   NHS reviewed and confirmed ICDA code
dxmonth        49-52   Date of dx months since 1900

tstage     129 Extent of disease through bowel wall
        $label 1.T1 invades submucosa;
        2.T2 invades muscularis;
        3.T3 invades thru muscularis to subserosa;
        4.T4 invades other organs, structures, perforation;
        5.Tis carcinoma in situ, invasion of lamina propria;
        6.Tx unknown
nstage          i1              132             Number of lymph nodes involved
        $label 1.N0 - no regional lymph node metastases;
        2.N1 - 1-3 regional lymph node metastases;
        3.N2 - 4+ regional lymph node metastases;
        4.positive lymph nodes: unknown number;
        5.Nx unknown status
mets          i1          135     Distant metastases
        $label 1.M0-no distant metastases;
        2.M1-distant metastases

histgr        i1              123       Histologic Grade
        $label 1.well differentiated;
        2.grade II;
        3.grade III poor;
        4.unspecified

site    30     1= ascending (right)
               2= transverse
               3= descending, splenic flexure (left)
               4= sigmoid
               5= rectum
               6= cecum
               7= site unspecified
               8= rectosigmoid              
*************************************/

data crc_nhs2;
  set crc;
  dtdxca = dxmonth;
  if conf in (11,12,13,14,15,16,17,18,19)  then cr=1; else cr=0;
  conf1119=0; if 11<=conf<=19 then conf1119=1;
  
  *if tstage=5 and nstage=1 and mets=1 then stage_cr=1;  /* Tis, N0, M0 */
  *else if tstage in (1,2) and nstage=1 and mets=1 then stage_cr=1;  /* T1-T2, N0, M0 */
  *else if tstage in (3,4) and nstage=1 and mets=1 then stage_cr=2;  /* T3-T4, N0, M0 */
  *else if nstage in (2,3,4) then stage_cr=3;  /* Any T, N1-N2 or unknown positive */
  *else if mets=2 then stage_cr=4;  /* Any T, any N, M1 */
  *else stage_cr=.;  /* Unknown/unspecified */

  /* AJCC stage: M first, then N, then T (unknown M treated as M0) -- same rule as NHS and HPFS */
  if      mets = 2                       then stage_cr = 4;   /* any T, any N, M1          */
  else if nstage in (2,3,4)              then stage_cr = 3;   /* any T, N1-N2 / N+ unknown */
  else if nstage = 1 and tstage in (3,4) then stage_cr = 2;   /* T3-T4, N0                 */
  else if nstage = 1 and tstage in (1,2) then stage_cr = 1;   /* T1-T2, N0                 */
  else if nstage = 1 and tstage = 5      then stage_cr = 1;   /* Tis, N0                   */
  else                                        stage_cr = .;   /* Nx, or N0 with Tx         */

  /* SEER summary stage (needs only node status) */
  if      stage_cr = 4       then sumstage = 3;   /* distant (excluded) */
  else if nstage = 1         then sumstage = 1;   /* localized          */
  else if nstage in (2,3,4)  then sumstage = 2;   /* regional           */
  else                            sumstage = .;


  /* GRADE: Create grade_cr from histgr */
  if histgr=1 then grade_cr=1;      /* Well differentiated */
  else if histgr=2 then grade_cr=2; /* Grade II (moderately differentiated) */
  else if histgr=3 then grade_cr=3; /* Grade III (poorly differentiated) */
  else grade_cr=4;                  /* Unspecified/unknown */
  
  /* SUBSITE: Create subsite_cr from site variable */
    if      site in (1,2,6) then subsite_cr=1; /* Proximal */
    else if site in (3,4)   then subsite_cr=2; /* Distal */
    else if site in (5,8)   then subsite_cr=3; /* Rectum */
    else subsite_cr=4;                        /* Unknown */
  
  dtdx_cr = int(dxmonth/12 + 1900);
  
keep id cr conf1119 stage_cr sumstage tstage nstage mets grade_cr subsite_cr dtdxca dtdx_cr;
run;
proc sort data=crc_nhs2; by id dtdxca; run;
proc sort data=crc_nhs2 nodupkey; by id; run;

proc means data=crc_nhs2; var _numeric_; run;
proc freq data=crc_nhs2; tables cr*(stage_cr grade_cr subsite_cr) / missing; run;
proc freq data=crc_nhs2;tables dtdx_cr;run; /*dtdx_cr 1989 to 2023;*/

proc freq data=crc_nhs2; where conf1119 = 1;
  tables nstage*mets*stage_cr / list missing;      /* how many N+ M1 move from III to IV */
  tables cr*(stage_cr sumstage grade_cr subsite_cr) / missing;
run;

                                                                             
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
    merge dead diabetes mi stroke update_cabg;
    by id;

    /* CVD */
    dtdxchd=.;
    if 0<dtdxmi<9999 or  0<dtdxcabg<9999 then dtdxchd=min(dtdxmi,dtdxcabg, .);
    
    dtdxcvd=.;
    if 0<dtdxchd<9999 or 0<dtdxstr<9999  then dtdxcvd=min(dtdxchd,dtdxstr,.);

    keep id dtdth newicda3 dtdxdb2 dtdxcvd dtdxchd dtdxmi dtdxstr; run;
    proc sort; by id; run;

    PROC DATASETS;
    delete diabetes dead mi stroke update_cabg cabg nur15 nur17 nur19; run;


*----------------2.Covariables-----------------------------;
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
       hbp89 chol89); 

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
             oc91 can91 brcn91 ocan91 mel91);
        can91=0; if brcn91=1 or ocan91=1 then can91=1;
        run;

%nur93 (keep=wt93 db93 hbp93 chol93 mi93 ang93 str93 
             oc93 brcn93  ocan93 mel93 can93  
             mmi93  fmi93  mifh93);
        can93=0; if brcn93=1 or ocan93=1 or mel93=1 then can93=1;
        mifh93=0; if mmi93=1 or fmi93=1 then mifh93=1;
        run;

%nur95 (keep=wt95 hbp95 db95 chol95 mi95 ang95 str95 cabg95 
             oc95 brcn95  ocan95 mel95 can95);  
        can95=0; if brcn95=1 or ocan95=1 or mel95=1 then can95=1;
        run;


%nur97 (keep=hbp97 db97 chol97 mi97 ang97 str97 cabg97 
             dbfh97 pdb97  sbdb97  
             mifh97 mmi97 fmi97 sbmi97 strfh97 mstr97 fstr97 sbstr97
             oc97 brcn97  ocan97 mel97 can97 
             cafh97 mbrcn97 sbrc197 sbrc297 sov97 mov97 pmel97 smel97 sclc197 sclc297 pclc97
             );

            if mmi97=1  or fmi97=1  or sbmi97=1  then mifh97=1; else mifh97=0;
            if mstr97=1 or fstr97=1 or sbstr97=1 then strfh97=1; else  strfh97=0;
            if pdb97=1 or sbdb97=1 then dbfh97=1; else dbfh97=0;
            if mbrcn97=1 or sbrc197=1 or sbrc297=1 or sov97=1 or mov97=1 or pmel97=1 or smel97=1 or sclc197=1 or sclc297=1 or pclc97=1 then cafh97=1; else cafh97=0;
            can97=0; if brcn97=1 or ocan97=1 or mel97=1 then can97=1;
            run;


%nur99 (keep=wt99 hbp99 db99 chol99 mi99 ang99 str99 cabg99 
             oc99 brcn99 ocan99 mel99 can99);
             can99=0; if brcn99=1 or ocan99=1 or mel99=1 then can99=1;
             run;

%nur01 (keep=wt01 db01 hbp01 chol01 mi01 ang01 str01 cabg01 
             dbfh01 mdb01 fdb01 sdb01   
             mifh01 mmi01 fmi01 smi01  strfh01 mstr01 fstr01 sstr01
             oc01  brcn01  ocan01 mel01 can01
             cafh01 mov01 sov01 mbrcn01 sbrc101 sbrc201 pclc01 sclc101 sclc201 fpro01 bpro01 mut01 smut01 ppan01 span01 smel01 pmel01);
             if mmi01=1 or fmi01=1 or smi01=1 then mifh01=1; else mifh01=0;
             if mstr01=1 or fstr01=1 or sstr01=1 then strfh01=1; else strfh01=0;
             if fdb01=1 or mdb01=1 or sdb01=1 then dbfh01=1; else dbfh01=0;
             if mov01=1 or sov01=1 or mbrcn01=1 or sbrc101=1 or sbrc201=1 or pclc01=1 or sclc101=1 or sclc201=1 or 
                fpro01=1 or bpro01=1 or mut01=1 or smut01=1 or ppan01=1 or span01=1 or smel01=1 or pmel01=1 then cafh01=1; else cafh01=0; 
             can01=0; if brcn01=1 or ocan01=1 or mel01=1 then can01=1;
             run;

%nur03(keep= wt03 db03 hbp03 chol03 mi03 ang03 str03 cabg03 brcn03 ocan03 mel03 can03
             cerbr03d);  
             can03=0; if brcn03=1 or ocan03=1 or mel03=1 then can03=1;
             run;

%nur05 (keep=wt05 db05 hbp05 chol05 mi05 ang05 str05 cabg05 
             dbfh05 mdb05 fdb05 sdb05  
             mifh05 mdchd05 fdchd05 mdstr05 fdstr05 strfh05
             oc05 brcn05 ocan05 mel05 can05 
             cafh05 mov05 sov05 mbrcn05 sbrc105 sbrc205 sclc105 pclc05 smel05 pmel05 sov05 mov05 );

             if mdb05=1 or fdb05=1 or sdb05=1 then dbfh05=1; else dbfh05=0;
             if mdchd05=1 or fdchd05=1 then mifh05=1; else mifh05=0;
             if mdstr05=1 or fdstr05=1 then strfh05=1; else strfh05=0;
             can05=0; if brcn05=1 or ocan05=1 or mel05=1 then can05=1;
             if mov05=1 or sov05=1 or mbrcn05=1 or sbrc105=1 or sbrc205=1 or sclc105=1 or
                pclc05=1 or smel05=1 or pmel05=1 or sov05=1 or mov05=1 then cafh05=1; else cafh05=0;
             run;

%nur07(keep= db07 hbp07 chol07 mi07 ang07 str07 cabg07 
             brcn07 ocan07 mel07 can07
             cerbr07d);
             can07=0; if brcn07=1 or ocan07=1 or mel07=1 then can07=1;
             run;

%nur09 (keep=wt09 db09 hbp09 chol09 mi09 str09 ang09 cabg09 
             dbfh09 pdb09 sdb09 mdchd09 fdchd09 mifh09 mdstr09 fdstr09 strfh09
             oc09 brcn09  ocan09 mel09 can09 
             cafh09 mov09 sov09 mbrcn09 sbrcn09 pclc09 sclc109 sclc209 pmel09 smel09);

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
             cerbr11d);
             can11=0; if brcn11=1 or ut11=1 or ov11=1 or colc11=1 or mel11=1 or ocan11=1 then can11=1;
             run;

%nur13 (keep=wt13 db13 hbp13 chol13 mi13 str13 cabg13 ang13 
             dbfh13 hxdb13 mdchd13 fdchd13 mifh13 mdstr13 fdstr13 strfh13
             brcn13 ut13 ov13 colc13 mel13 ocan13 can13 
             cafh13 hxov13 hxbc13 hxmel13 hxcol13);

            if hxdb13=2 then dbfh13=1;else dbfh13=0;
            if mdchd13=1 or fdchd13=1 then mifh13=1; else mifh13=0;
            if mdstr13=1 or fdstr13=1 then strfh13=1; else strfh13=0;
            if hxov13=2 or hxov13=2 or hxbc13=2 or hxmel13=2 or hxcol13=2 then cafh13=1;else cafh13=0;
            can13=0; if brcn13=1 or ut13=1 or ov13=1 or colc13=1 or mel13=1 or ocan13=1 then can13=1;
            run;


%nur15(keep=wt15 db15 hbp15  chol15  mi15  ang15 str15 cabg15  
            brcn15  ocan15 ut15 ov15 mel15 can15 
            cerbr15d);
            can15=0; if brcn15=1 or ocan15=1 or ut15=1 or ov15=1 or mel15=1 then can15=1;
            run;


%nur17(keep=wt17 db17 hbp17  chol17  
            mi17  ang17  str17 cabg17 
            brcn17  ocan17 ut17 ov17 mel17 can17);
            can17=0; if brcn17=1 or ocan17=1 or ut17=1 or ov17=1 or mel17=1 then can17=1;
            run;

%nur19(keep=wt19 db19 hbp19  chol19  
            mi19  ang19  str19 cabg19 
            brcn19  ocan19 /*ut19*/ ov19 mel19 can19);
            can19=0; if brcn19=1 or ocan19=1 /*or ut19=1*/ or ov19=1 or mel19=1 then can19=1;
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
alfs91=sum(carbo91pdd,aprot91pdd,vfat91pdd);
vlfs91=sum(carbo91pdd,vprot91pdd,afat91pdd);

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
 merge der8919(in=mst) disease crc_nhs2(in=crc) ahei2010_9119   act8917 supp8915 meds8919 oc8911
       nur89 nur91 nur93 nur95 nur97 nur99  nur01 nur03 nur05 nur07 nur09 nur11 nur13 nur15 nur17 nur19
       n91_nts n95_nts n99_nts n03_nts n07_nts n11_nts n15_nts n19_nts 
       nses8917 end=_end_;
       by id;    
       exrec=1; if first.id and mst then exrec=0; 
       crc_cases=1; if crc then crc_cases =0;

    /*pathology info for CRC cases*/
      %indic3(vbl=subsite_cr,prefix=subsite_cr,min=2, max=4,reflev=1,missing=.,usemiss=0,
          label1='proximal colon',
          label2='distal colon',
          label3='rectum',
          label4='Unspecified');

      %indic3(vbl=stage_cr,prefix=stage_cr,min=2, max=4,reflev=1,missing=.,usemiss=1,
          label1='stage I',
          label2='stage II',
          label3='stage III',
          label4='stage IV');

      %indic3(vbl=grade_cr,prefix=grade_cr,min=2, max=4,reflev=1,missing=.,usemiss=0,
          label1='Well Differentiated',
          label2='Grade 2 Moderately Diff',
          label3='Grade 3 Poorly Diff',
          label4='Unspecified');
      
/* Family history of diabetes, MI, stroke, and cancer */
    if dbfh89=1 or dbfh97=1  or dbfh01=1  or dbfh05=1 or dbfh09=1 or dbfh13=1  then dbfh=1;else dbfh=0; /* last dbfh 13*/
    if mifh89=1 or mifh93=1 or mifh97=1 or mifh01=1 or mifh05=1  or mifh09=1 or mifh13=1 then mifh=1; else mifh=0;
    if strfh97=1 or strfh01=1 or strfh05=1 or strfh09=1 or strfh13=1 then strfh=1; else strfh=0;
    if cafh89=1 or cafh97=1 or cafh01=1 or cafh05=1 or cafh09=1 or cafh13=1 then cafh=1; else cafh=0;  


/* chronic disease status (flag for stopping diet updates) */
    array canc_prev {*} can89  can91  can93  can95  can97  can99  can01  can03  can05  can07  can09  can11 can13  can15 can17 can19;
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
        if mi(i) =1 or strk(i)=1 or cabg(i)=1  or diab(i)=1  /*or canc_prev(i)=1*/   then chrodis(i)=1; else chrodis(i)=0; 
    end;

    do i = 2 to dim(chrodis); if chrodis{i-1} = 1 then chrodis{i} = 1; 
    end; /* Make chrodis cumulative */

*** nSES;
if nSES_89 = . then nSES_89m = -0.426; else nSES_89m = nSES_89;
if nSES_91 = . then nSES_91m = nSES_89m; else nSES_91m = nSES_91;

/***************** Diet scores: carry forward + STOP updating after cancer/CVD *******************/
/* original, simple updated scores */ 
array irt   {*} retmo91  retmo95  retmo99  retmo03  retmo07 retmo11 retmo15 retmo19;
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
array alcon {*} alco91n alco95n alco99n alco03n alco07n  alco11n  alco15n alco19n;

array    carbp  {*}   carbo91p         carbo95p       carbo99p         carbo03p       carbo07p         carbo11p       carbo15p         carbo19p   ;
array    hcarbp  {*}   hcarbo91p        hcarbo95p      hcarbo99p        hcarbo03p      hcarbo07p        hcarbo11p      hcarbo15p        hcarbo19p  ;
array sumref_carbp {*}  sumref_carb91p   sumref_carb95p sumref_carb99p   sumref_carb03p sumref_carb07p   sumref_carb11p  sumref_carb15p   sumref_carb19p      ;
array     protp  {*}   prot91p          prot95p        prot99p          prot03p        prot07p          prot11p        prot15p          prot19p;
array    aprotp  {*}   aprot91p         aprot95p       aprot99p         aprot03p       aprot07p         aprot11p       aprot15p         aprot19p ;
array    vprotp  {*}   vprot91p         vprot95p       vprot99p         vprot03p       vprot07p         vprot11p       vprot15p         vprot19p ;
array    tfatp   {*}   tfat91p          tfat95p        tfat99p          tfat03p        tfat07p          tfat11p        tfat15p          tfat19p;
array    afatp   {*}   afat91p          afat95p        afat99p          afat03p        afat07p          afat11p        afat15p          afat19p;
array    vfatp   {*}   vfat91p          vfat95p        vfat99p          vfat03p        vfat07p          vfat11p        vfat15p          vfat19p;

/* stop updated scores */
array calors {*} calor91ns calor95ns calor99ns calor03ns calor07ns calor11ns calor15ns calor19ns; 
array aheims {*} ahei91s ahei95s ahei99s ahei03s ahei07s ahei11s ahei15s ahei19s;
array lcss  {*} lcs91s lcs95s lcs99s lcs03s lcs07s lcs11s lcs15s lcs19s;
array alcss {*} alcs91s alcs95s alcs99s alcs03s alcs07s alcs11s alcs15s alcs19s;
array vlcss {*} vlcs91s vlcs95s vlcs99s vlcs03s vlcs07s vlcs11s vlcs15s vlcs19s;
array ulcss {*} ulcs91s ulcs95s ulcs99s ulcs03s ulcs07s ulcs11s ulcs15s ulcs19s;
array hlcss {*} hlcs91s hlcs95s hlcs99s hlcs03s hlcs07s hlcs11s hlcs15s hlcs19s;
array lfss  {*} lfs91s lfs95s lfs99s lfs03s lfs07s lfs11s lfs15s lfs19s;
array ulfss {*} ulfs91s ulfs95s ulfs99s ulfs03s ulfs07s ulfs11s ulfs15s ulfs19s;
array hlfss {*} hlfs91s hlfs95s hlfs99s hlfs03s hlfs07s hlfs11s hlfs15s hlfs19s;
array alfss {*} alfs91s alfs95s alfs99s alfs03s alfs07s alfs11s alfs15s alfs19s;
array vlfss {*} vlfs91s vlfs95s vlfs99s vlfs03s vlfs07s vlfs11s vlfs15s vlfs19s;
array alcons {*} alco91ns alco95ns alco99ns alco03ns alco07ns  alco11ns  alco15ns alco19ns;

array    carbs  {*}   carbo91s         carbo95s       carbo99s         carbo03s       carbo07s         carbo11s       carbo15s         carbo19s   ;
array    hcarbs  {*}   hcarbo91s        hcarbo95s      hcarbo99s        hcarbo03s      hcarbo07s        hcarbo11s      hcarbo15s        hcarbo19s  ;
array sumref_carbs {*}  sumref_carb91s   sumref_carb95s sumref_carb99s   sumref_carb03s sumref_carb07s   sumref_carb11s  sumref_carb15s   sumref_carb19s      ;
array     prots  {*}   prot91s          prot95s        prot99s          prot03s        prot07s          prot11s        prot15s          prot19s; 
array    aprots  {*}   aprot91s         aprot95s       aprot99s         aprot03s       aprot07s         aprot11s       aprot15s         aprot19s ;
array    vprots  {*}   vprot91s         vprot95s       vprot99s         vprot03s       vprot07s         vprot11s       vprot15s         vprot19s ;
array    tfats   {*}   tfat91s          tfat95s        tfat99s          tfat03s        tfat07s          tfat11s        tfat15s          tfat19s; 
array    afats   {*}   afat91s          afat95s        afat99s          afat03s        afat07s          afat11s        afat15s          afat19s; 
array    vfats   {*}   vfat91s          vfat95s        vfat99s          vfat03s        vfat07s          vfat11s        vfat15s          vfat19s; 

/* nominal FFQ dates (June of each FFQ year, months since 1900), used when irt is missing */
array ffqdt {8} _temporary_ (1098 1146 1194 1242 1290 1338 1386 1434);   /* 1991 1995 1999 2003 2007 2011 2015 2019 */

/******** STEP 1. Post-diagnosis inputs: RAW values, post-diagnosis FFQs only ********/
/* Must run BEFORE the carry-forward in Step 2, so a missing post-dx FFQ is not      */
/* filled with a pre-diagnosis value.                                                */
do i = 1 to dim(calor);
    calors{i} = calor{i};   aheims{i} = aheim{i};
    lcss{i}   = lcs{i};     alcss{i}  = alcs{i};  vlcss{i} = vlcs{i};  ulcss{i} = ulcs{i};  hlcss{i} = hlcs{i};
    lfss{i}   = lfs{i};     alfss{i}  = alfs{i};  vlfss{i} = vlfs{i};  ulfss{i} = ulfs{i};  hlfss{i} = hlfs{i};
    alcons{i} = alcon{i};

    carbs{i}  = carbp{i};   hcarbs{i} = hcarbp{i}; sumref_carbs{i} = sumref_carbp{i};
    prots{i}  = protp{i};   aprots{i} = aprotp{i}; vprots{i} = vprotp{i};
    tfats{i}  = tfatp{i};   afats{i}  = afatp{i};  vfats{i}  = vfatp{i};

    /* keep cycle i only if its FFQ was completed after diagnosis */
    ffq_date = ifn(irt{i} > 0, irt{i}, ffqdt{i});
    if not (0 < dtdxca < ffq_date) then do;          /* FFQ on/before diagnosis, or no diagnosis */
        calors{i} = .;  aheims{i} = .;
        lcss{i} = .;  alcss{i} = .;  vlcss{i} = .;  ulcss{i} = .;  hlcss{i} = .;
        lfss{i} = .;  alfss{i} = .;  vlfss{i} = .;  ulfss{i} = .;  hlfss{i} = .;
        alcons{i} = .;
        carbs{i} = .; hcarbs{i} = .; sumref_carbs{i} = .;
        prots{i} = .; aprots{i} = .; vprots{i} = .;
        tfats{i} = .; afats{i}  = .; vfats{i}  = .;
    end;
end;
drop ffq_date;


/******** STEP 2. Carry forward the simple updated (pre-dx) scores -- unchanged ********/
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
if alcon{i}=.  then alcon{i}=alcon{i-1};

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


     /* pre-dx cumulative average up to pre FFQ (1 to i-1) */
array diet{8,13}
    calor91n ahei91 lcs91  alcs91 vlcs91  ulcs91  hlcs91 lfs91    ulfs91 hlfs91 alfs91 vlfs91  alco91n
    calor95n ahei95 lcs95  alcs95 vlcs95  ulcs95  hlcs95 lfs95    ulfs95 hlfs95 alfs95 vlfs95  alco95n
    calor99n ahei99 lcs99  alcs99 vlcs99  ulcs99  hlcs99 lfs99    ulfs99 hlfs99 alfs99 vlfs99  alco99n
    calor03n ahei03 lcs03  alcs03 vlcs03  ulcs03  hlcs03 lfs03    ulfs03 hlfs03 alfs03 vlfs03  alco03n
    calor07n ahei07 lcs07  alcs07 vlcs07  ulcs07  hlcs07 lfs07    ulfs07 hlfs07 alfs07 vlfs07  alco07n
    calor11n ahei11 lcs11  alcs11 vlcs11  ulcs11  hlcs11 lfs11    ulfs11 hlfs11 alfs11 vlfs11  alco11n
    calor15n ahei15 lcs15  alcs15 vlcs15  ulcs15  hlcs15 lfs15    ulfs15 hlfs15 alfs15 vlfs15  alco15n
    calor19n ahei19 lcs19  alcs19 vlcs19  ulcs19  hlcs19 lfs19    ulfs19 hlfs19 alfs19 vlfs19  alco19n
    ;

array avg_pre{13} avg_calor_pre avg_aheim_pre avg_lcs_pre avg_alcs_pre avg_vlcs_pre avg_ulcs_pre avg_hlcs_pre avg_lfs_pre avg_ulfs_pre avg_hlfs_pre avg_alfs_pre avg_vlfs_pre avg_alcon_pre ;
array sum_pre{13} _temporary_;
array n_pre  {13} _temporary_;

do k = 1 to 13;
    avg_pre{k} = .; sum_pre{k} = 0; n_pre{k}   = 0;
end;

if dtdxca > 0 then do;
    pre_idx = .;
    do j = 1 to dim(irt); if 0 < irt{j} < dtdxca then pre_idx = j; end;
    if pre_idx > 0 then do;
        do j = 1 to pre_idx;          /* loop over FFQ cycles      */
            do k = 1 to 13;           /* loop over diet indices     */
                if not missing(diet{j,k}) then do;
                    sum_pre{k} = sum_pre{k} + diet{j,k};
                    n_pre{k}   = n_pre{k} + 1;
                end;end;end;
        do k = 1 to 13;
            if n_pre{k} > 0 then avg_pre{k} = sum_pre{k} / n_pre{k};
            else avg_pre{k} = .;
        end;end;end;



/******** Post-diagnosis cumulative averages restarting at post_cycle ********/

array dietin{8,13}
calor91ns ahei91s lcs91s alcs91s vlcs91s ulcs91s hlcs91s lfs91s ulfs91s hlfs91s alfs91s vlfs91s alco91ns
calor95ns ahei95s lcs95s alcs95s vlcs95s ulcs95s hlcs95s lfs95s ulfs95s hlfs95s alfs95s vlfs95s alco95ns
calor99ns ahei99s lcs99s alcs99s vlcs99s ulcs99s hlcs99s lfs99s ulfs99s hlfs99s alfs99s vlfs99s alco99ns
calor03ns ahei03s lcs03s alcs03s vlcs03s ulcs03s hlcs03s lfs03s ulfs03s hlfs03s alfs03s vlfs03s alco03ns
calor07ns ahei07s lcs07s alcs07s vlcs07s ulcs07s hlcs07s lfs07s ulfs07s hlfs07s alfs07s vlfs07s alco07ns
calor11ns ahei11s lcs11s alcs11s vlcs11s ulcs11s hlcs11s lfs11s ulfs11s hlfs11s alfs11s vlfs11s alco11ns
calor15ns ahei15s lcs15s alcs15s vlcs15s ulcs15s hlcs15s lfs15s ulfs15s hlfs15s alfs15s vlfs15s alco15ns
calor19ns ahei19s lcs19s alcs19s vlcs19s ulcs19s hlcs19s lfs19s ulfs19s hlfs19s alfs19s vlfs19s alco19ns
;
array dietout{8,13}
calor91nsv ahei91sv lcs91sv alcs91sv vlcs91sv ulcs91sv hlcs91sv lfs91sv ulfs91sv hlfs91sv alfs91sv vlfs91sv alco91sv
calor95nsv ahei95sv lcs95sv alcs95sv vlcs95sv ulcs95sv hlcs95sv lfs95sv ulfs95sv hlfs95sv alfs95sv vlfs95sv alco95sv
calor99nsv ahei99sv lcs99sv alcs99sv vlcs99sv ulcs99sv hlcs99sv lfs99sv ulfs99sv hlfs99sv alfs99sv vlfs99sv alco99sv
calor03nsv ahei03sv lcs03sv alcs03sv vlcs03sv ulcs03sv hlcs03sv lfs03sv ulfs03sv hlfs03sv alfs03sv vlfs03sv alco03sv
calor07nsv ahei07sv lcs07sv alcs07sv vlcs07sv ulcs07sv hlcs07sv lfs07sv ulfs07sv hlfs07sv alfs07sv vlfs07sv alco07sv
calor11nsv ahei11sv lcs11sv alcs11sv vlcs11sv ulcs11sv hlcs11sv lfs11sv ulfs11sv hlfs11sv alfs11sv vlfs11sv alco11sv
calor15nsv ahei15sv lcs15sv alcs15sv vlcs15sv ulcs15sv hlcs15sv lfs15sv ulfs15sv hlfs15sv alfs15sv vlfs15sv alco15sv
calor19nsv ahei19sv lcs19sv alcs19sv vlcs19sv ulcs19sv hlcs19sv lfs19sv ulfs19sv hlfs19sv alfs19sv vlfs19sv alco19sv
;

/* calculate cumulative average intake */
do j=1 to 13;
  dietout{1,j}=dietin{1,j};
 do i=2 to 8;
  sumvar=0;
   n=0;
     do k=1 to i;
      if (dietin{k,j} ne .) then do;
      n=n+1; sumvar=sumvar+dietin{k,j};
      end;
     end;
    if n = 0 then dietout{i,j} = .;          /* no post-dx FFQ yet -> missing */
   else dietout{i,j}=sumvar/n;
  end;
end;


/*calculate cumulative average of updated score**/
array macron{8,9}
carbo91ns ahei91s  sumref_carb91s prot91s aprot91s vprot91s tfat91s  afat91s vfat91s
carbo95ns ahei95s  sumref_carb95s prot95s aprot95s vprot95s tfat95s  afat95s vfat95s
carbo99ns ahei99s  sumref_carb99s prot99s aprot99s vprot99s tfat99s  afat99s vfat99s
carbo03ns ahei03s  sumref_carb03s prot03s aprot03s vprot03s tfat03s  afat03s vfat03s
carbo07ns ahei07s  sumref_carb07s prot07s aprot07s vprot07s tfat07s  afat07s vfat07s
carbo11ns ahei11s  sumref_carb11s prot11s aprot11s vprot11s tfat11s  afat11s vfat11s
carbo15ns ahei15s  sumref_carb15s prot15s aprot15s vprot15s tfat15s  afat15s vfat15s
carbo19ns ahei19s  sumref_carb19s prot19s aprot19s vprot19s tfat19s  afat19s vfat19s
;
array macronv{8,9}
carbo91nsv ahei91sv  sumref_carb91sv prot91sv aprot91sv vprot91sv tfat91sv  afat91sv vfat91sv
carbo95nsv ahei95sv  sumref_carb95sv prot95sv aprot95sv vprot95sv tfat95sv  afat95sv vfat95sv
carbo99nsv ahei99sv  sumref_carb99sv prot99sv aprot99sv vprot99sv tfat99sv  afat99sv vfat99sv
carbo03nsv ahei03sv  sumref_carb03sv prot03sv aprot03sv vprot03sv tfat03sv  afat03sv vfat03sv
carbo07nsv ahei07sv  sumref_carb07sv prot07sv aprot07sv vprot07sv tfat07sv  afat07sv vfat07sv
carbo11nsv ahei11sv  sumref_carb11sv prot11sv aprot11sv vprot11sv tfat11sv  afat11sv vfat11sv
carbo15nsv ahei15sv  sumref_carb15sv prot15sv aprot15sv vprot15sv tfat15sv  afat15sv vfat15sv
carbo19nsv ahei19sv  sumref_carb19sv prot19sv aprot19sv vprot19sv tfat19sv  afat19sv vfat19sv
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
   if n=0 then macronv{i,j}=.;
   else macronv{i,j}=sumvar/n;
  end;
end;
run;


PROC DATASETS;
delete    disease crc_nhs2 ahei2010_9119  der8919  act8917 supp8915  meds8919  oc8911
          nur89 nur91 nur93 nur95 nur97 nur99  nur01 nur03 nur05 nur07 nur09 nur11 nur13 nur15 nur17 nur19
          n91_nts n95_nts n99_nts n03_nts n07_nts n11_nts n15_nts n19_nts 
          ;
RUN;

/*************************************************************************************/
/*                          Large section of arrays                                  */
/*************************************************************************************/
data nhs2;
     set base end=_end_;

cutoff=1434;    /*** Set cutoff at 2019.6 ***/
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

array   conf1119a {15};
array   alcq    {15}    alcc91n   alcc91n   alcc95n     alcc95n   alcc99n   alcc99n    alcc03n   alcc03n   alcc07n   alcc07n   alcc11n   alcc11n alcc15n alcc15n alcc19n;
array   act     {15}    act91m    act91m    act91m      act97m    act97m     act01m    act01m    act05m    act05m    act09m    act09m    act13m   act13m  act17m act17m;
array   mvit    {15}    mvitu91  mvitu93   mvitu95    mvitu97   mvitu99  mvitu01    mvitu03   mvitu05   mvitu07   mvitu09  mvitu11   mvitu13  mvitu15  mvitu15 mvitu15;
array   aspi    {15}    aspu91   aspu93    aspu95     aspu97    aspu99   aspu01     aspu03    aspu05    aspu07    aspu09   aspu11    aspu13   aspu15 aspu17 aspu19;
array   ocuse   {15}    nocu91   nocu93    nocu95     nocu97    nocu99   nocu01     nocu03    nocu05    nocu07    nocu09   nocu11    nocu11   nocu11 nocu11 nocu11;
array   hbp     {15}    hbp91  hbp93  hbp95  hbp97  hbp99  hbp01  hbp03  hbp05  hbp07  hbp09  hbp11 hbp13  hbp15 hbp17  hbp19;
array   chol    {15}    chol91 chol93 chol95 chol97 chol99 chol01 chol03 chol05 chol07 chol09 chol11 chol13 chol15 chol17 chol19;
array   diab    {15}    db91   db93   db95   db97   db99   db01   db03   db05   db07   db09   db11 db13  db15 db17  db19;
array   canc    {15}    can91  can93  can95  can97  can99  can01  can03  can05  can07  can09  can11 can13  can15 can17 can19;


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

/* There are a few cases with date of diagnosis after date of death.  Since were more sure about date of death, use that date before we assign period of diagnosis */
    if 0<dtdxca<9999 and 0<dtdth<9999 then do;
        if dtdxca>dtdth then dtdxca=dtdth;
    end;

do j=2 to DIM(bmi);
      if age{j}=.      then age{j}=age{j-1}+2;
      *if bmi{j}=.      then bmi{j}=bmi{j-1};
      if phms{j}=.      then phms{j}=phms{j-1}; 
      if hor{j}=.      then hor{j}=hor{j-1}; 
      *if smok{j}=.     then smok{j}=smok{j-1}; /* smoking status range 0-15 0.missing 1.never 2-8.past 9-15.curr*/
      *if alcon{j}=.    then alcon{j}=alcon{j-1};
      if act{j}=.      then act{j}=act{j-1};
      if mvit{j}=.     then mvit{j}=mvit{j-1};
      if aspi{j}=.     then aspi{j}=aspi{j-1};
      if ocuse{j}=.     then ocuse{j}=ocuse{j-1};
      *if marriage{j}=. then marriage{j}=marriage{j-1};
      *if alone{j}=.    then alone{j}=alone{j-1};
      if hbp{j}=.     then hbp{j}=hbp{j-1};
      if chol{j}=.    then chol{j}=chol{j-1};
      if nSESArr{j}=.   then  nSESArr{j}= nSESArr{j-1};
end;

do i=2 to DIM(bmi);
if irt{i-1}<=dtdxca and dtdxca <irt{i} then do;
  *age_dx=age{i};
  age_dx=age{i-1} + (dtdxca - irt{i-1})/12; /* age at diagnosis*/
  bmi_dx=bmi{i-1}; /* BMI prior to diagnosis*/
  *bmi_dx=bmi{i};
  highbp_dx=hbp{i};
  hightc_dx=chol{i};
  db_dx=diab{i};
end;end;
if age_dx = . and dtdxca >= irt{15} then do;
  age_dx=age{15};
  bmi_dx=bmi{15};
  highbp_dx=hbp{15};
  hightc_dx=chol{15};
  db_dx=diab{15};
end;



/* -original  s-stop updated  sv-postdx averaged   pre-predx   post-postdx single */
/* -original*/
array calor {*} calor91n calor91n calor95n calor95n calor99n calor99n calor03n calor03n calor07n calor07n calor11n calor11n calor15n calor15n calor19n; 
array aheim {*} ahei91 ahei91 ahei95 ahei95 ahei99 ahei99 ahei03 ahei03 ahei07 ahei07 ahei11 ahei11 ahei15 ahei15 ahei19;
array lcs   {*} lcs91 lcs91 lcs95 lcs95 lcs99 lcs99 lcs03 lcs03 lcs07 lcs07 lcs11 lcs11 lcs15 lcs15 lcs19;
array alcs  {*} alcs91 alcs91 alcs95 alcs95 alcs99 alcs99 alcs03 alcs03 alcs07 alcs07 alcs11 alcs11 alcs15 alcs15 alcs19;
array vlcs  {*} vlcs91 vlcs91 vlcs95 vlcs95 vlcs99 vlcs99 vlcs03 vlcs03 vlcs07 vlcs07 vlcs11 vlcs11 vlcs15 vlcs15 vlcs19;
array ulcs  {*} ulcs91 ulcs91 ulcs95 ulcs95 ulcs99 ulcs99 ulcs03 ulcs03 ulcs07 ulcs07 ulcs11 ulcs11 ulcs15 ulcs15 ulcs19;
array hlcs  {*} hlcs91 hlcs91 hlcs95 hlcs95 hlcs99 hlcs99 hlcs03 hlcs03 hlcs07 hlcs07 hlcs11 hlcs11 hlcs15 hlcs15 hlcs19;
array lfs   {*} lfs91  lfs91  lfs95  lfs95  lfs99  lfs99  lfs03  lfs03  lfs07  lfs07  lfs11  lfs11  lfs15  lfs15  lfs19;
array ulfs  {*} ulfs91 ulfs91 ulfs95 ulfs95 ulfs99 ulfs99 ulfs03 ulfs03 ulfs07 ulfs07 ulfs11 ulfs11 ulfs15 ulfs15 ulfs19;
array hlfs  {*} hlfs91 hlfs91 hlfs95 hlfs95 hlfs99 hlfs99 hlfs03 hlfs03 hlfs07 hlfs07 hlfs11 hlfs11 hlfs15 hlfs15 hlfs19;
array alfs  {*} alfs91 alfs91 alfs95 alfs95 alfs99 alfs99 alfs03 alfs03 alfs07 alfs07 alfs11 alfs11 alfs15 alfs15 alfs19;
array vlfs  {*} vlfs91 vlfs91 vlfs95 vlfs95 vlfs99 vlfs99 vlfs03 vlfs03 vlfs07 vlfs07 vlfs11 vlfs11 vlfs15 vlfs15 vlfs19;

/* -s stop updated */
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

array    carbs  {*}   carbo91s  carbo91s carbo95s  carbo95s carbo99s  carbo99s carbo03s  carbo03s carbo07s  carbo07s   carbo11s  carbo11s carbo15s  carbo15s  carbo19s   ;
array    hcarbs  {*}   hcarbo91s  hcarbo91s hcarbo95s  hcarbo95s hcarbo99s  hcarbo99s hcarbo03s  hcarbo03s hcarbo07s  hcarbo07s   hcarbo11s  hcarbo11s hcarbo15s  hcarbo15s  hcarbo19s   ;
array sumref_carbs {*}  sumref_carb91s  sumref_carb91s sumref_carb95s  sumref_carb95s sumref_carb99s  sumref_carb99s sumref_carb03s  sumref_carb03s sumref_carb07s  sumref_carb07s   sumref_carb11s  sumref_carb11s sumref_carb15s  sumref_carb15s  sumref_carb19s   ;
array     prots  {*}   prot91s  prot91s prot95s  prot95s prot99s  prot99s prot03s  prot03s prot07s  prot07s   prot11s  prot11s prot15s  prot15s  prot19s   ;
array    aprots  {*}   aprot91s  aprot91s aprot95s  aprot95s aprot99s  aprot99s aprot03s  aprot03s aprot07s  aprot07s   aprot11s  aprot11s aprot15s  aprot15s  aprot19s   ;
array    vprots  {*}   vprot91s  vprot91s vprot95s  vprot95s vprot99s  vprot99s vprot03s  vprot03s vprot07s  vprot07s   vprot11s  vprot11s vprot15s  vprot15s  vprot19s   ;
array    tfats   {*}   tfat91s  tfat91s tfat95s  tfat95s tfat99s  tfat99s tfat03s  tfat03s tfat07s  tfat07s   tfat11s  tfat11s tfat15s  tfat15s  tfat19s   ;
array    afats   {*}   afat91s  afat91s afat95s  afat95s afat99s  afat99s afat03s  afat03s afat07s  afat07s   afat11s  afat11s afat15s  afat15s  afat19s   ;
array    vfats   {*}   vfat91s  vfat91s vfat95s  vfat95s vfat99s  vfat99s vfat03s  vfat03s vfat07s  vfat07s   vfat11s  vfat11s vfat15s  vfat15s  vfat19s   ;
array alcons {15} alco91ns alco91ns alco95ns alco95ns alco99ns alco99ns alco03ns alco03ns alco07ns alco07ns alco11ns alco11ns alco15ns alco15ns alco19ns;

/* -sv  post-dx cumulative average */
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
array alcosv {15} alco91sv alco91sv alco95sv alco95sv alco99sv alco99sv alco03sv alco03sv alco07sv alco07sv alco11sv alco11sv alco15sv alco15sv alco19sv;

array carbsv {15} carbo91nsv carbo91nsv carbo95nsv carbo95nsv carbo99nsv carbo99nsv carbo03nsv carbo03nsv carbo07nsv carbo07nsv carbo11nsv carbo11nsv carbo15nsv carbo15nsv carbo19nsv;
array hcarbsv {15} hcarbo91sv hcarbo91sv hcarbo95sv hcarbo95sv hcarbo99sv hcarbo99sv hcarbo03sv hcarbo03sv hcarbo07sv hcarbo07sv hcarbo11sv hcarbo11sv hcarbo15sv hcarbo15sv hcarbo19sv; 
array sumref_carbsv {15} sumref_carb91sv sumref_carb91sv sumref_carb95sv sumref_carb95sv sumref_carb99sv sumref_carb99sv sumref_carb03sv sumref_carb03sv sumref_carb07sv sumref_carb07sv sumref_carb11sv sumref_carb11sv sumref_carb15sv sumref_carb15sv sumref_carb19sv;
array protsv {15} prot91sv prot91sv prot95sv prot95sv prot99sv prot99sv prot03sv prot03sv prot07sv prot07sv prot11sv prot11sv prot15sv prot15sv prot19sv;
array aprotsv {15} aprot91sv aprot91sv aprot95sv aprot95sv aprot99sv aprot99sv aprot03sv aprot03sv aprot07sv aprot07sv aprot11sv aprot11sv aprot15sv aprot15sv aprot19sv;
array vprotsv {15} vprot91sv vprot91sv vprot95sv vprot95sv vprot99sv vprot99sv vprot03sv vprot03sv vprot07sv vprot07sv vprot11sv vprot11sv vprot15sv vprot15sv vprot19sv;
array tfatsv {15} tfat91sv tfat91sv tfat95sv tfat95sv tfat99sv tfat99sv tfat03sv tfat03sv tfat07sv tfat07sv tfat11sv tfat11sv tfat15sv tfat15sv tfat19sv;
array afatsv {15} afat91sv afat91sv afat95sv afat95sv afat99sv afat99sv afat03sv afat03sv afat07sv afat07sv afat11sv afat11sv afat15sv afat15sv afat19sv;
array vfatsv {15} vfat91sv vfat91sv vfat95sv vfat95sv vfat99sv vfat99sv vfat03sv vfat03sv vfat07sv vfat07sv vfat11sv vfat11sv vfat15sv vfat15sv vfat19sv;

/* carry forward for updated  scores */
do i=2 to dim(calorsv);
if calorsv{i}=.   then calorsv{i}=calorsv{i-1};
if aheimsv{i}=.   then aheimsv{i}=aheimsv{i-1};
if lcssv{i}=.     then lcssv{i}=lcssv{i-1};
if alcssv{i}=.    then alcssv{i}=alcssv{i-1};
if vlcssv{i}=.    then vlcssv{i}=vlcssv{i-1};
if ulcssv{i}=.    then ulcssv{i}=ulcssv{i-1};
if hlcssv{i}=.    then hlcssv{i}=hlcssv{i-1};
if lfssv{i}=.     then lfssv{i}=lfssv{i-1};
if ulfssv{i}=.    then ulfssv{i}=ulfssv{i-1};
if hlfssv{i}=.    then hlfssv{i}=hlfssv{i-1};
if alfssv{i}=.   then alfssv{i}=alfssv{i-1};
if vlfssv{i}=.   then vlfssv{i}=vlfssv{i-1};
if alcosv{i}=.   then alcosv{i}=alcosv{i-1};

if carbsv{i}=.     then carbsv{i}=carbsv{i-1};
if hcarbsv{i}=.    then hcarbsv{i}=hcarbsv{i-1};
if sumref_carbsv{i}=.    then sumref_carbsv{i}=sumref_carbsv{i-1};
if protsv{i}=.   then protsv{i}=protsv{i-1};
if aprotsv{i}=.   then aprotsv{i}=aprotsv{i-1};
if vprotsv{i}=.     then vprotsv{i}=vprotsv{i-1};
if tfatsv{i}=.    then tfatsv{i}=tfatsv{i-1};
if afatsv{i}=.    then afatsv{i}=afatsv{i-1};
if vfatsv{i}=.   then vfatsv{i}=vfatsv{i-1};
end;


/*-pre pre-dx */
array calor_pre {15}  calor91p_pre calor91p_pre calor95p_pre calor95p_pre calor99p_pre calor99p_pre calor03p_pre calor03p_pre calor07p_pre calor07p_pre calor11p_pre calor11p_pre calor15p_pre calor15p_pre calor19p_pre ;
array aheim_pre {15}  ahei91_pre ahei91_pre ahei95_pre ahei95_pre ahei99_pre ahei99_pre ahei03_pre ahei03_pre ahei07_pre ahei07_pre ahei11_pre ahei11_pre ahei15_pre ahei15_pre ahei19_pre;   
array lcs_pre   {15}  lcs91_pre    lcs91_pre    lcs95_pre    lcs95_pre    lcs99_pre    lcs99_pre    lcs03_pre    lcs03_pre    lcs07_pre    lcs07_pre    lcs11_pre    lcs11_pre    lcs15_pre   lcs15_pre   lcs19_pre    ;
array vlcs_pre  {15}  vlcs91_pre   vlcs91_pre   vlcs95_pre   vlcs95_pre   vlcs99_pre   vlcs99_pre   vlcs03_pre   vlcs03_pre   vlcs07_pre   vlcs07_pre   vlcs11_pre   vlcs11_pre   vlcs15_pre  vlcs15_pre  vlcs19_pre   ;
array alcs_pre  {15}  alcs91_pre   alcs91_pre   alcs95_pre   alcs95_pre   alcs99_pre   alcs99_pre   alcs03_pre   alcs03_pre   alcs07_pre   alcs07_pre   alcs11_pre   alcs11_pre   alcs15_pre  alcs15_pre  alcs19_pre   ;
array ulcs_pre  {15}  ulcs91_pre   ulcs91_pre   ulcs95_pre   ulcs95_pre   ulcs99_pre   ulcs99_pre   ulcs03_pre   ulcs03_pre   ulcs07_pre   ulcs07_pre   ulcs11_pre   ulcs11_pre   ulcs15_pre  ulcs15_pre  ulcs19_pre   ; 
array hlcs_pre  {15}  hlcs91_pre   hlcs91_pre   hlcs95_pre   hlcs95_pre   hlcs99_pre   hlcs99_pre   hlcs03_pre   hlcs03_pre   hlcs07_pre   hlcs07_pre   hlcs11_pre   hlcs11_pre   hlcs15_pre  hlcs15_pre  hlcs19_pre   ;
array lfs_pre   {15}  lfs91_pre    lfs91_pre    lfs95_pre    lfs95_pre    lfs99_pre    lfs99_pre    lfs03_pre    lfs03_pre    lfs07_pre    lfs07_pre    lfs11_pre    lfs11_pre    lfs15_pre   lfs15_pre   lfs19_pre    ;
array vlfs_pre  {15}  vlfs91_pre   vlfs91_pre   vlfs95_pre   vlfs95_pre   vlfs99_pre   vlfs99_pre   vlfs03_pre   vlfs03_pre   vlfs07_pre   vlfs07_pre   vlfs11_pre   vlfs11_pre   vlfs15_pre  vlfs15_pre  vlfs19_pre   ;
array alfs_pre  {15}  alfs91_pre   alfs91_pre   alfs95_pre   alfs95_pre   alfs99_pre   alfs99_pre   alfs03_pre   alfs03_pre   alfs07_pre   alfs07_pre   alfs11_pre   alfs11_pre   alfs15_pre  alfs15_pre  alfs19_pre   ;
array ulfs_pre  {15}  ulfs91_pre   ulfs91_pre   ulfs95_pre   ulfs95_pre   ulfs99_pre   ulfs99_pre   ulfs03_pre   ulfs03_pre   ulfs07_pre   ulfs07_pre   ulfs11_pre   ulfs11_pre   ulfs15_pre  ulfs15_pre  ulfs19_pre   ; 
array hlfs_pre  {15}  hlfs91_pre   hlfs91_pre   hlfs95_pre   hlfs95_pre   hlfs99_pre   hlfs99_pre   hlfs03_pre   hlfs03_pre   hlfs07_pre   hlfs07_pre   hlfs11_pre   hlfs11_pre   hlfs15_pre  hlfs15_pre  hlfs19_pre   ;
array alco_pre  {15}  alco91_pre   alco91_pre   alco95_pre   alco95_pre   alco99_pre   alco99_pre   alco03_pre   alco03_pre   alco07_pre   alco07_pre   alco11_pre   alco11_pre   alco15_pre  alco15_pre  alco19_pre   ;

/*Post-dx, donnot use the carry forward value to ensure the cumulaitve averaged validity*/
array calor_post {15}  calor91p_post calor91p_post calor95p_post calor95p_post calor99p_post calor99p_post calor03p_post calor03p_post calor07p_post calor07p_post calor11p_post calor11p_post calor15p_post calor15p_post calor19p_post;
array aheim_post {15}  ahei91_post ahei91_post ahei95_post ahei95_post ahei99_post ahei99_post ahei03_post ahei03_post ahei07_post ahei07_post ahei11_post ahei11_post ahei15_post ahei15_post ahei19_post;   
array lcs_post   {15}  lcs91_post    lcs91_post    lcs95_post    lcs95_post    lcs99_post    lcs99_post    lcs03_post    lcs03_post    lcs07_post    lcs07_post    lcs11_post    lcs11_post    lcs15_post  lcs15_post  lcs19_post    ;
array vlcs_post  {15}  vlcs91_post   vlcs91_post   vlcs95_post   vlcs95_post   vlcs99_post   vlcs99_post   vlcs03_post   vlcs03_post   vlcs07_post   vlcs07_post   vlcs11_post   vlcs11_post   vlcs15_post vlcs15_post   vlcs19_post   ;
array alcs_post  {15}  alcs91_post   alcs91_post   alcs95_post   alcs95_post   alcs99_post   alcs99_post   alcs03_post   alcs03_post   alcs07_post   alcs07_post   alcs11_post   alcs11_post   alcs15_post alcs15_post   alcs19_post   ;
array ulcs_post  {15}  ulcs91_post   ulcs91_post   ulcs95_post   ulcs95_post   ulcs99_post   ulcs99_post   ulcs03_post   ulcs03_post   ulcs07_post   ulcs07_post   ulcs11_post   ulcs11_post   ulcs15_post ulcs15_post   ulcs19_post   ; 
array hlcs_post  {15}  hlcs91_post   hlcs91_post   hlcs95_post   hlcs95_post   hlcs99_post   hlcs99_post   hlcs03_post   hlcs03_post   hlcs07_post   hlcs07_post   hlcs11_post   hlcs11_post   hlcs15_post hlcs15_post   hlcs19_post   ;
array lfs_post   {15}  lfs91_post    lfs91_post    lfs95_post    lfs95_post    lfs99_post    lfs99_post    lfs03_post    lfs03_post    lfs07_post    lfs07_post    lfs11_post    lfs11_post    lfs15_post  lfs15_post    lfs19_post   ;
array vlfs_post  {15}  vlfs91_post   vlfs91_post   vlfs95_post   vlfs95_post   vlfs99_post   vlfs99_post   vlfs03_post   vlfs03_post   vlfs07_post   vlfs07_post   vlfs11_post   vlfs11_post   vlfs15_post vlfs15_post   vlfs19_post   ;
array alfs_post  {15}  alfs91_post   alfs91_post   alfs95_post   alfs95_post   alfs99_post   alfs99_post   alfs03_post   alfs03_post   alfs07_post   alfs07_post   alfs11_post   alfs11_post   alfs15_post alfs15_post   alfs19_post   ;
array ulfs_post  {15}  ulfs91_post   ulfs91_post   ulfs95_post   ulfs95_post   ulfs99_post   ulfs99_post   ulfs03_post   ulfs03_post   ulfs07_post   ulfs07_post   ulfs11_post   ulfs11_post   ulfs15_post ulfs15_post   ulfs19_post   ; 
array hlfs_post  {15}  hlfs91_post   hlfs91_post   hlfs95_post   hlfs95_post   hlfs99_post   hlfs99_post   hlfs03_post   hlfs03_post   hlfs07_post   hlfs07_post   hlfs11_post   hlfs11_post   hlfs15_post hlfs15_post   hlfs19_post   ;
array alco_post  {15}  alco91_post   alco91_post   alco95_post   alco95_post   alco99_post   alco99_post   alco03_post   alco03_post   alco07_post   alco07_post   alco11_post   alco11_post   alco15_post  alco15_post  alco19_post   ;


/*-crc Post-dx, within 6 month to 4 years */
array calor_crc {15}  calor91p_crc calor91p_crc calor95p_crc calor95p_crc calor99p_crc calor99p_crc calor03p_crc calor03p_crc calor07p_crc calor07p_crc calor11p_crc calor11p_crc calor15p_crc calor15p_crc calor19p_crc ;
array aheim_crc {15}  ahei91_crc ahei91_crc ahei95_crc ahei95_crc ahei99_crc ahei99_crc ahei03_crc ahei03_crc ahei07_crc ahei07_crc ahei11_crc ahei11_crc ahei15_crc ahei15_crc ahei19_crc;   
array lcs_crc   {15}  lcs91_crc    lcs91_crc    lcs95_crc    lcs95_crc    lcs99_crc    lcs99_crc    lcs03_crc    lcs03_crc    lcs07_crc    lcs07_crc    lcs11_crc    lcs11_crc    lcs15_crc    lcs15_crc  lcs19_crc    ;
array vlcs_crc  {15}  vlcs91_crc   vlcs91_crc   vlcs95_crc   vlcs95_crc   vlcs99_crc   vlcs99_crc   vlcs03_crc   vlcs03_crc   vlcs07_crc   vlcs07_crc   vlcs11_crc   vlcs11_crc   vlcs15_crc   vlcs15_crc vlcs19_crc   ;
array alcs_crc  {15}  alcs91_crc   alcs91_crc   alcs95_crc   alcs95_crc   alcs99_crc   alcs99_crc   alcs03_crc   alcs03_crc   alcs07_crc   alcs07_crc   alcs11_crc   alcs11_crc   alcs15_crc   alcs15_crc alcs19_crc   ;
array ulcs_crc  {15}  ulcs91_crc   ulcs91_crc   ulcs95_crc   ulcs95_crc   ulcs99_crc   ulcs99_crc   ulcs03_crc   ulcs03_crc   ulcs07_crc   ulcs07_crc   ulcs11_crc   ulcs11_crc   ulcs15_crc   ulcs15_crc ulcs19_crc   ; 
array hlcs_crc  {15}  hlcs91_crc   hlcs91_crc   hlcs95_crc   hlcs95_crc   hlcs99_crc   hlcs99_crc   hlcs03_crc   hlcs03_crc   hlcs07_crc   hlcs07_crc   hlcs11_crc   hlcs11_crc   hlcs15_crc   hlcs15_crc hlcs19_crc   ;
array lfs_crc   {15}  lfs91_crc    lfs91_crc    lfs95_crc    lfs95_crc    lfs99_crc    lfs99_crc    lfs03_crc    lfs03_crc    lfs07_crc    lfs07_crc    lfs11_crc    lfs11_crc    lfs15_crc    lfs15_crc  lfs19_crc    ;
array vlfs_crc  {15}  vlfs91_crc   vlfs91_crc   vlfs95_crc   vlfs95_crc   vlfs99_crc   vlfs99_crc   vlfs03_crc   vlfs03_crc   vlfs07_crc   vlfs07_crc   vlfs11_crc   vlfs11_crc   vlfs15_crc   vlfs15_crc vlfs19_crc   ;
array alfs_crc  {15}  alfs91_crc   alfs91_crc   alfs95_crc   alfs95_crc   alfs99_crc   alfs99_crc   alfs03_crc   alfs03_crc   alfs07_crc   alfs07_crc   alfs11_crc   alfs11_crc   alfs15_crc   alfs15_crc alfs19_crc   ;
array ulfs_crc  {15}  ulfs91_crc   ulfs91_crc   ulfs95_crc   ulfs95_crc   ulfs99_crc   ulfs99_crc   ulfs03_crc   ulfs03_crc   ulfs07_crc   ulfs07_crc   ulfs11_crc   ulfs11_crc   ulfs15_crc   ulfs15_crc ulfs19_crc   ; 
array hlfs_crc  {15}  hlfs91_crc   hlfs91_crc   hlfs95_crc   hlfs95_crc   hlfs99_crc   hlfs99_crc   hlfs03_crc   hlfs03_crc   hlfs07_crc   hlfs07_crc   hlfs11_crc   hlfs11_crc   hlfs15_crc   hlfs15_crc hlfs19_crc   ;
array alco_crc  {15}  alco91_crc   alco91_crc   alco95_crc   alco95_crc   alco99_crc   alco99_crc   alco03_crc   alco03_crc   alco07_crc   alco07_crc   alco11_crc   alco11_crc   alco15_crc  alco15_crc  alco19_crc   ;


/* Compute pre, post, and chg */
do i=2 to DIM(calor)+1;
    if irt{i-1}<dtdxca and dtdxca<=irt{i} then do;
        /* pre-dx (one time): FFQ just before dx */
        pre_calor_fix = calor{i-1};
        pre_aheim_fix = aheim{i-1};
        pre_lcs_fix   = lcs{i-1};
        pre_alcs_fix  = alcs{i-1};
        pre_vlcs_fix  = vlcs{i-1};
        pre_ulcs_fix  = ulcs{i-1};
        pre_hlcs_fix  = hlcs{i-1};
        pre_lfs_fix   = lfs{i-1};
        pre_ulfs_fix  = ulfs{i-1};
        pre_hlfs_fix  = hlfs{i-1};
        pre_alfs_fix  = alfs{i-1};
        pre_vlfs_fix  = vlfs{i-1};
        pre_alco_fix  = alcon{i-1};  
        /*Set all pre value the same after dx*/
        do j=i-1 to DIM(calor);
          calor_pre{j}=pre_calor_fix;
          aheim_pre{j}=pre_aheim_fix;
          lcs_pre{j} =pre_lcs_fix;
          alcs_pre{j}=pre_alcs_fix;
          vlcs_pre{j}=pre_vlcs_fix;
          ulcs_pre{j}=pre_ulcs_fix;
          hlcs_pre{j}=pre_hlcs_fix;
          lfs_pre{j} =pre_lfs_fix;
          alfs_pre{j}=pre_alfs_fix;
          vlfs_pre{j}=pre_vlfs_fix;
          ulfs_pre{j}=pre_ulfs_fix;
          hlfs_pre{j}=pre_hlfs_fix;
          alco_pre{j}=pre_alco_fix;
        end;
    end;
end;

        /* If missing, impute with pre-dx cumulative average */

do i=1 to DIM(calor);
    if calor_pre{i} = . then calor_pre{i} = avg_calor_pre;
    if aheim_pre{i} = . then aheim_pre{i} = avg_aheim_pre;
    if lcs_pre{i}   = . then lcs_pre{i}   = avg_lcs_pre;
    if alcs_pre{i}  = . then alcs_pre{i}  = avg_alcs_pre;
    if vlcs_pre{i}  = . then vlcs_pre{i}  = avg_vlcs_pre;
    if ulcs_pre{i}  = . then ulcs_pre{i}  = avg_ulcs_pre;
    if hlcs_pre{i}  = . then hlcs_pre{i}  = avg_hlcs_pre;
    if lfs_pre{i}   = . then lfs_pre{i}   = avg_lfs_pre;
    if ulfs_pre{i}  = . then ulfs_pre{i}  = avg_ulfs_pre;
    if hlfs_pre{i}  = . then hlfs_pre{i}  = avg_hlfs_pre;
    if alfs_pre{i}  = . then alfs_pre{i}  = avg_alfs_pre;
    if vlfs_pre{i}  = . then vlfs_pre{i}  = avg_vlfs_pre;
    if alco_pre{i}  = . then alco_pre{i}  = avg_alcon_pre;
end;

/* FFQ positions in the 17-element period arrays: 91 95 99 03 07 11 15 19  */
array ffqk {8} _temporary_ (1 3 5 7 9 11 13 15 );

/*---- Post-diagnosis, single time: first valid FFQ completed after diagnosis ----*/
post_calor_fix=.; post_aheim_fix=.; post_lcs_fix=.; post_alcs_fix=.; post_vlcs_fix=.; post_ulcs_fix=.;
post_hlcs_fix=.;  post_lfs_fix=.;   post_alfs_fix=.; post_vlfs_fix=.; post_ulfs_fix=.;  post_hlfs_fix=.;
if dtdxca > 0 then do f = 1 to dim(ffqk);
    k = ffqk{f};
    if post_lcs_fix = . and lcss{k} ne . then do;       /* *s = raw, post-dx only, first valid post-diagnosis FFQ */
        post_calor_fix = calors{k}; post_aheim_fix = aheims{k};
        post_lcs_fix  = lcss{k};  post_alcs_fix = alcss{k}; post_vlcs_fix = vlcss{k};
        post_ulcs_fix = ulcss{k}; post_hlcs_fix = hlcss{k};
        post_lfs_fix  = lfss{k};  post_alfs_fix = alfss{k}; post_vlfs_fix = vlfss{k};
        post_ulfs_fix = ulfss{k}; post_hlfs_fix = hlfss{k};
        post_alco_fix = alcons{k};  post_ffq_dt = irt{k};
    end;
end;
if post_lcs_fix ne . then do j = 1 to dim(calor);
    if irt{j+1} > dtdxca then do;                       /* periods from diagnosis onward */
        calor_post{j}=post_calor_fix; aheim_post{j}=post_aheim_fix;
        lcs_post{j}=post_lcs_fix;   alcs_post{j}=post_alcs_fix; vlcs_post{j}=post_vlcs_fix;
        ulcs_post{j}=post_ulcs_fix; hlcs_post{j}=post_hlcs_fix;
        lfs_post{j}=post_lfs_fix;   alfs_post{j}=post_alfs_fix; vlfs_post{j}=post_vlfs_fix;
        ulfs_post{j}=post_ulfs_fix; hlfs_post{j}=post_hlfs_fix;
    end;
end;

/*---- Post-diagnosis, 6 months to 4 years: first valid FFQ in that window ----*/
crc_calor_fix=.; crc_aheim_fix=.; crc_lcs_fix=.; crc_alcs_fix=.; crc_vlcs_fix=.; crc_ulcs_fix=.;
crc_hlcs_fix=.;  crc_lfs_fix=.;   crc_alfs_fix=.; crc_vlfs_fix=.; crc_ulfs_fix=.;  crc_hlfs_fix=.;
if dtdxca > 0 then do f = 1 to dim(ffqk);
    k = ffqk{f};
    if crc_lcs_fix = . and lcss{k} ne . and 6 <= irt{k} - dtdxca <= 48 then do;
        crc_calor_fix = calors{k}; crc_aheim_fix = aheims{k};
        crc_lcs_fix  = lcss{k};  crc_alcs_fix = alcss{k}; crc_vlcs_fix = vlcss{k};
        crc_ulcs_fix = ulcss{k}; crc_hlcs_fix = hlcss{k};
        crc_lfs_fix  = lfss{k};  crc_alfs_fix = alfss{k}; crc_vlfs_fix = vlfss{k};
        crc_ulfs_fix = ulfss{k}; crc_hlfs_fix = hlfss{k};
        crc_alco_fix  = alcons{k};  crc_ffq_dt  = irt{k};
    end;
end;

if crc_calor_fix ne . then do i=1 to DIM(calor);
    if irt{i} > dtdxca then do;
        calor_crc{i} = crc_calor_fix;
        aheim_crc{i} = crc_aheim_fix;
        lcs_crc{i}   = crc_lcs_fix;
        alcs_crc{i}  = crc_alcs_fix;
        vlcs_crc{i}  = crc_vlcs_fix;
        ulcs_crc{i}  = crc_ulcs_fix;
        hlcs_crc{i}  = crc_hlcs_fix;
        lfs_crc{i}   = crc_lfs_fix;
        alfs_crc{i}  = crc_alfs_fix;
        vlfs_crc{i}  = crc_vlfs_fix;
        ulfs_crc{i}  = crc_ulfs_fix;
        hlfs_crc{i}  = crc_hlfs_fix;
    end;
end;

        elig_any = 0;
        if dtdxca > 0 then do f = 1 to dim(ffqk);  k = ffqk{f};   /* -- FFQ positions */
           nomd = 1074+24*k;       /* -- NHS ; NHSII 1074+24*k ; HPFS 1009+24*k */
           if nomd > dtdxca and nomd < cutoff and (dtdth in (.,0) or dtdth > nomd)
              then elig_any = 1;      /* -- alive when an FFQ was mailed after dx */
        end;



%beginex();
/*************************************************************/
/*                DO-LOOP OVER TIME PERIODS                  */

   
    do i=1 to DIM(irt)-1;     
      
    period=i;
      do j=1 to DIM(tvar);
      tvar{j}=0;
      end;
    tvar{i}=1; 

    conf1119a{i}=0;
    /* Confirmed invasive breast cancer case if conf1119=1 */
    if irt{i}<dtdxca<=irt{i+1} and conf1119=1 then do;
        conf1119a{i}=1;
    end;
    /* Carry forward conf1119=1 indefinitely so that cases stay in analysis */
    if i>1 then do;
        if conf1119a{i-1}=1 then conf1119a{i}=1;
    end;


/****** OUTCOME ******/
    /*** death ***/
    dead = 0; tdead = irt{i+1} - irt{i};
    if irt{i} < dtdth <= irt{i+1} then do; dead = 1; tdead = dtdth - irt{i}; end;
    
    /*death within 1 year of cancer diagnose*/
    d1y=0; if 0 <= (dtdth - dtdxca)/12 and (dtdth - dtdxca)/12 <= 1 then d1y = 1;

    /* Cause-specific death */
    dead_can   = 0; tdead_can=irt{i+1}-irt{i}; 
    if irt{i}< dtdth <=irt{i+1} and (140<=newicda3<=208) then do; dead_can=1; tdead_can=dtdth-irt{i}; end;

    dead_nocan=0; tdead_nocan=irt{i+1}-irt{i}; dtdx_dead_nocan=.;
    if irt{i}<dtdth<=irt{i+1} and dead_can ne 1 then do; dead_nocan=1; tdead_nocan=dtdth-irt{i}; dtdx_dead_nocan=dtdth; end;

    dead_brc   = 0; tdead_brc=irt{i+1}-irt{i}; dtdx_dead_brc=.; /*breast cancer death*/
    if irt{i}< dtdth <=irt{i+1} and (newicda3=174) then do; dead_brc=1; tdead_brc=dtdth-irt{i}; dtdx_dead_brc=dtdth; end;

    dead_crc   = 0; tdead_crc=irt{i+1}-irt{i}; dtdx_dead_crc=.; /*colorectal cancer death*/
    if irt{i}< dtdth <=irt{i+1} and newicda3 in (153,154) then do; dead_crc=1; tdead_crc=dtdth-irt{i}; dtdx_dead_crc=dtdth; end;

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
    
    /** cancer dignosis **/
   cancase=0;tcan=irt{i+1}-irt{i};
     if irt{i}<dtdxca<=irt{i+1} then do;
        cancase=1; tcan=dtdxca-irt{i};
        end;

   canpre=0;
   if 0<dtdxca<=irt{i+1} then do;
   canpre=1;
   end; 

  /****** crc duration ******/
      crcduration = (irt{i} - dtdxca)/12;
      if crcduration < 1 then crcduration = 1;
      *Define crc duration category; 
        if 0<crcduration<=5 then crcduration_c=1;
        else if  5<crcduration<=10 then crcduration_c=2;
        else if 10<crcduration<=15 then crcduration_c=3;
        else if 15<crcduration then crcduration_c=4;
      %indic3(vbl=crcduration_c, prefix=crcduration_c, min=2, max=4, reflev=1, usemiss=0, missing=.,
           label1='0-5 years', label2='5-10 years',label3='10-15 yaers',label4='>15 years');

   /* Year of diagnosis */
      if dtdxca=9999 then yeardx_crc=.;
      else yeardx_crc=round(((dtdxca/12)+1900-0.5)); /* Calendar year of diagnosis */


/***************time since diagnosis (post-dx survival time scale) ******************/
            
      if dtdxca > .  and conf1119 = 1 then do;
            end_cal = min(irt{i+1}, cutoff);         /*calendar end of the interval*/
            if end_cal <= dtdxca then do;
                tsd         = .;
                time_begin  = .;
                time_end    = .;
                person_time = 0;
            end;
            else do; 

                if irt{i} le dtdxca lt irt{i+1} then do; /* Interval contains the diagnosis date */
                    time_begin=0;
                end;          
                else do;                                 /* Intervals after diagnosis */
                time_begin = irt{i} - dtdxca; end;
                tsd = time_begin;

                if 0 < dtdth < end_cal then end_cal = dtdth;
                time_end = end_cal - dtdxca;
                if time_end < time_begin then time_end = time_begin;
                person_time = time_end - time_begin; 
            end;  
        end;
        else do;
        tsd        = .;
        time_begin = .;
        time_end   = .;
        person_time = 0;
    end;


/******EXPOSURE******/
    /* -single updated */
    calorcon=calor{i};
    aheimcon=aheim{i};
    lcscon=lcs{i};
    alcscon=alcs{i};
    vlcscon=vlcs{i};
    ulcscon=ulcs{i};
    hlcscon=hlcs{i};
    lfscon=lfs{i};
    ulfscon=ulfs{i};
    hlfscon=hlfs{i};
    alfscon=alfs{i};
    vlfscon=vlfs{i};

    carb=carbs {i};
    hcarb=hcarbs {i};
    sumref_carb=sumref_carbs {i};
    prot=prots {i};
    aprot=aprots {i};
    vprot=vprots {i};
    tfat=tfats {i};
    afat=afats {i};
    vfat=vfats {i};
    /* 2. Cumulatively averaged post-dx score (restarting at post_cycle) */
    calorsvcon=calorsv{i};   * post dx TEI, cumulative from i+1 ;
    aheimsvcon=aheimsv{i};   * post dx AHEI, cumulative from i+1;
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

    carb_svcon=carbsv {i};
    hcarb_svcon=hcarbsv {i};
    sumref_carb_svcon=sumref_carbsv {i};
    prot_svcon=protsv {i};
    aprot_svcon=aprotsv {i};
    vprot_svcon=vprotsv {i};
    tfat_svcon=tfatsv {i};
    afat_svcon=afatsv {i};
    vfat_svcon=vfatsv {i};
    
    /* 3. pre-dx score */

    /* 4. post-dx single updated score */

    /* 5. post-dx 6 months-4 years score */

    /* 6. -svchg: change from pre-dx: cumulative post - single pre*/
    lcs_svchg = .;  if lcs_pre{i} ne . and lcssvcon ne . then lcs_svchg = lcssvcon - lcs_pre{i};
    alcs_svchg = .; if alcs_pre{i} ne . and alcssvcon ne . then alcs_svchg = alcssvcon - alcs_pre{i};
    vlcs_svchg = .; if vlcs_pre{i} ne . and vlcssvcon ne . then vlcs_svchg = vlcssvcon - vlcs_pre{i};
    ulcs_svchg= .; if ulcs_pre{i} ne . and ulcssvcon ne . then ulcs_svchg = ulcssvcon - ulcs_pre{i};
    hlcs_svchg = .; if hlcs_pre{i} ne . and hlcssvcon ne . then hlcs_svchg = hlcssvcon - hlcs_pre{i};
    lfs_svchg = .;  if lfs_pre{i} ne . and lfssvcon ne . then lfs_svchg = lfssvcon - lfs_pre{i};
    alfs_svchg = .; if alfs_pre{i} ne . and alfssvcon ne . then alfs_svchg = alfssvcon - alfs_pre{i};
    vlfs_svchg = .; if vlfs_pre{i} ne . and vlfssvcon ne . then vlfs_svchg = vlfssvcon - vlfs_pre{i};
    ulfs_svchg = .; if ulfs_pre{i} ne . and ulfssvcon ne . then ulfs_svchg = ulfssvcon - ulfs_pre{i};
    hlfs_svchg = .; if hlfs_pre{i} ne . and hlfssvcon ne . then hlfs_svchg = hlfssvcon - hlfs_pre{i};

    /* 7. change from pre-dx: single post - single pre*/
    lcs_change = .;  if lcs_pre{i} ne . and lcs_post{i} ne . then lcs_change = lcs_post{i} - lcs_pre{i};
    alcs_change = .; if alcs_pre{i} ne . and alcs_post{i} ne . then alcs_change = alcs_post{i} - alcs_pre{i};
    vlcs_change = .; if vlcs_pre{i} ne . and vlcs_post{i} ne . then vlcs_change = vlcs_post{i} - vlcs_pre{i};
    ulcs_change= .; if ulcs_pre{i} ne . and ulcs_post{i} ne . then ulcs_change = ulcs_post{i} - ulcs_pre{i};
    hlcs_change = .; if hlcs_pre{i} ne . and hlcs_post{i} ne . then hlcs_change = hlcs_post{i} - hlcs_pre{i};
    lfs_change = .;  if lfs_pre{i} ne . and lfs_post{i} ne . then lfs_change = lfs_post{i} - lfs_pre{i};
    alfs_change = .; if alfs_pre{i} ne . and alfs_post{i} ne . then alfs_change = alfs_post{i} - alfs_pre{i};
    vlfs_change = .; if vlfs_pre{i} ne . and vlfs_post{i} ne . then vlfs_change = vlfs_post{i} - vlfs_pre{i};
    ulfs_change = .; if ulfs_pre{i} ne . and ulfs_post{i} ne . then ulfs_change = ulfs_post{i} - ulfs_pre{i};
    hlfs_change = .; if hlfs_pre{i} ne . and hlfs_post{i} ne . then hlfs_change = hlfs_post{i} - hlfs_pre{i};


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
array   smkka    {15}    smkk91   smkk93   smkk95   smkk97   smkk99    smkk01   smkk03   smkk05   smkk07   smkk09   smkk11   smkk13  smkk15 smkk17  smkk19;
array   smksta   {15}    smksta91   smksta93   smksta95   smksta97   smksta99    smksta01  smksta03   smksta05   smksta07   smksta09   smksta11   smksta13  smksta15 smksta17  smksta19;

length smkk 3 smkstatus 3;
        if smok{i}=1 then smkk=1;  *never;
        else if smok{i} in(2,3,4,5,6,7,8) then smkk=2; *past smoker;
        else if smok{i} in (9,10) then smkk=3; *current, 1-14 cigs ;
        else if smok{i}=11 then smkk=4; *current, 15-24 cigs ;
        else if smok{i} in(12,13,14) then smkk=5; *current, 25+ cigs ;
        else smkk=.; 
    /* store time-varying 5-level category */
    smkka{i} = smkk;
     %indic3(vbl=smkk, prefix=smkk, min=2,max=5, reflev=1,missing=.,  usemiss=1,
              label1='never smker',
              label2='past smker',
              label3='curr 1-14 smker',
              label4='curr 15+ smker',
              label5='curr 25+ smker');

      /* 3-level: 1 never 2 past 3 current */
      if smkk=1 then smkstatus=1;
      else if smkk=2 then smkstatus=2;
      else if smkk>=3 then smkstatus=3;
      else smkstatus=.;
      /* store time-varying 3-level category */
      smksta{i} = smkstatus;
      %indic3(vbl=smkstatus, prefix=smkstatus, min=2, max=3, reflev=1, missing=., usemiss=1);

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

            if alcon{i}=0.0                          then alcq{i}=1;
            else if alcon{i}>0.0  and alcon{i}<5.0   then alcq{i}=2;
            else if alcon{i}>=5.0 and alcon{i}<10.0  then alcq{i}=3;
            else if alcon{i}>=10.0 and alcon{i}<15.0 then alcq{i}=4;
            else if alcon{i}>=15.0 and alcon{i}<30.0 then alcq{i}=5;
            else if alcon{i}>=30.0                   then alcq{i}=6;
            else                                          alcq{i}=9;

       alc=alcq{i};
       %indic3(vbl=alc, prefix=alc, min=2, max=6, reflev=1, missing=9, usemiss=1,
                label1='0 g/d',
                label2='0.1-4.9 g/d',
                label3='5.0-9.9 g/d',
                label4='10.0-14.9 g/d',
                label5='15.0-29.9 g/d',
                label6='30+ g/d');

        if alc =1 then neverdrinker=1; else neverdrinker=0;

        /* post dx alcohol intake, cumulative from i+1*/
        alcosvcon =alcosv{i};

              if      alcosv{i}=0       then alcovc_post=1;
            else if 0 <  alcosv{i}<5    then alcovc_post=2;
            else if 5<=  alcosv{i}<15   then alcovc_post=3;
            else if 15<=  alcosv{i}<30  then alcovc_post=4;
            else if 30<= alcosv{i}      then alcovc_post=5;
            else                             alcovc_post=.;
    
          label alcovc_post='Post-dx Cumulative average daily alcohol intake, g/d, categorical';
          %indic3(vbl=alcovc_post, prefix=alcovc_post, min=2, max=5, reflev=1, missing=., usemiss=0,
              label1='none',
              label2='0.1-4.9 g/d',
              label3='5.0-14.9 g/d',
              label4='15.0-29.9 g/d',
              label5='30+ g/d');

/*** Indicator for BMI ***/
         bmicon=bmi{i};

         bmicon=bmi{i};
         if      0< bmicon <25     then bmic3cat=1;
         else if 25<= bmicon < 30  then bmic3cat=2;
         else if 30<= bmicon       then bmic3cat=3;
         else                           bmic3cat=.;
         %indic3(vbl= bmic3cat, prefix= bmic3cat, min=2, max=3, reflev=1,  usemiss=1);

         if bmicon>=25 then bmic25=1;  else  bmic25=0;
         if bmicon>=30 then bmic30=1;  else  bmic30=0;

         bmibcon=bmi{1};
         if      0<   bmibcon <25 then bmib3cat=1;
         else if 25<= bmibcon <30 then bmib3cat=2;
         else if 30<= bmibcon     then bmib3cat=3;
         else                          bmib3cat=.;
         %indic3(vbl= bmib3cat, prefix= bmib3cat, min=2, max=3, reflev=1,  usemiss=1);

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

        if       0< bmi_dx<21    then bmid9=1;
        else if 21<=bmi_dx<23    then bmid9=2;
        else if 23<=bmi_dx<25    then bmid9=3;
        else if 25<=bmi_dx<27    then bmid9=4;
        else if 27<=bmi_dx<29    then bmid9=5;
        else if 29<=bmi_dx<31    then bmid9=6;
        else if 31<=bmi_dx<33    then bmid9=7;
        else if 33<=bmi_dx<35    then bmid9=8;
        else if     bmi_dx>=35   then bmid9=9;
        else bmid9=.;
        %indic3(vbl=bmid9, prefix=bmid9, min=2, max=9, reflev=1, usemiss=1);

        if      0<   bmi_dx <25 then bmid3cat=1;
        else if 25<= bmi_dx <30 then bmid3cat=2;
        else if 30<= bmi_dx     then bmid3cat=3;
        else                         bmid3cat=.;
        %indic3(vbl=bmid3cat, prefix=bmid3cat, min=2, max=3, reflev=1, usemiss=1);


/*** Indicator for OC use ******/
        if ocuse{i}=1 then ocu=1;
        else if ocuse{i}=2 then ocu=2;
        else if ocuse{i}=3 then ocu=3;
        else                 ocu=1;

%indic3(vbl=ocu, prefix=ocu, min=2, max=3, reflev=1,
                missing=., usemiss=0,
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
            if highbp_dx=1 then diaghbp=1; else diaghbp=0;
         
/****** Indicator for TC lowering ******
            select(antitc{i});
                when (1)     antihtc=1;
                otherwise    antihtc=0;
            end; */
            if chol{i}=1 then hightc=1; else hightc=0;            
            if chol{1}=1 then basehightc=1; else basehightc=0;
            if hightc_dx=1 then diaghtc=1; else diaghtc=0;

/****** Indicator for History of diabetes *******/
          if diab{i}=1 then diabetes=1; else diabetes=0;

/** nSES ***/
nSES=nSESArr{i};

/* Covariates at diagnosis, pre and post-diagnosis */
/*
    (D) Diagnosis: = from questionnaire in which cancer is reported
    (P) Pre: = from questionnaire preceding the cycle in which cancer is reported
    (A) After: = from questionnaire after the cycle in which cancer is reported
    (CH) Change:    */
    

/* Covariates at diagnosis */
 do b=2 to dim(irt)-1; if irt{b-1}<=dtdxca<irt{b} then do; 
            bmicond=bmi{b}; agecond=age{b};  actcond=act{b}; smkkd=smkka{b}; smkstd=smksta{b};
            mvitd=mvit{b}; aspd=aspi{b}; phmsd=phms{b}; hord=hor{b}; 
            hbpd=hbp{b}; chold=chol{b}; diabd=diab{b};
            alcd=alcq{b};
            nSESd=nSESArr{b}; 
            aheid=aheim{b}; lcsd=lcs{b}; alcsd=alcs{b}; vlcsd=vlcs{b}; ulcsd=ulcs{b}; hlcsd=hlcs{b};
            calord=calor{b}; lfsd=lfs{b}; alfsd=alfs{b}; vlfsd=vlfs{b}; ulfsd=ulfs{b}; hlfsd=hlfs{b};
            end; 
          end;
 
 do b=2 to dim(irt)-2; if irt{b-1}<=dtdxca<irt{b} then do; 
 /* Covariates pre diagnosis */
        bmiconp=bmi{b-1}; ageconp=age{b-1};  actconp=act{b-1}; smkkp=smkka{b-1}; smkstp=smksta{b-1};
        mvitp=mvit{b-1}; aspp=aspi{b-1}; phmsp=phms{b-1}; horp=hor{b-1}; 
        hbpp=hbp{b-1}; cholp=chol{b-1}; diabp=diab{b-1};
        alcp=alcq{b-1};
        nSESp=nSESArr{b-1}; 
        aheip=aheim{b-1}; /*lcsp=lcs{b-1}; alcsp=alcs{b-1}; vlcsp=vlcs{b-1}; ulcsp=ulcs{b-1}; hlcsp=hlcs{b-1};*/
        calorp=calor{b-1}; /*lfsp=lfs{b-1}; alfsp=alfs{b-1}; vlfsp=vlfs{b-1}; ulfsp=ulfs{b-1}; hlfsp=hlfs{b-1};*/
 /* Covariates post diagnosis */
        bmicona=bmi{b+1}; agecona=age{b+1}; actcona=act{b+1}; smkkaa=smkka{b+1}; smkstaa=smksta{b+1};
        mvita=mvit{b+1}; aspa=aspi{b+1}; phmsa=phms{b+1}; hora=hor{b+1}; 
        hbpa=hbp{b+1}; chola=chol{b+1}; diaba=diab{b+1};
        alca=alcq{b+1};
        nSESa=nSESArr{b+1};
        aheia=aheim{b+1}; /*lcsa=lcs{b+1}; alcsa=alcs{b+1}; vlcsa=vlcs{b+1}; ulcsa=ulcs{b+1}; hlcsa=hlcs{b+1};*/
        calora=calor{b+1}; /*lfsa=lfs{b+1}; alfsa=alfs{b+1}; vlfsa=vlfs{b+1}; ulfsa=ulfs{b+1}; hlfsa=hlfs{b+1};*/
 /* change  */
        bmiconch=bmi{b+1}-bmi{b-1}; ageconch=age{b+1}-age{b-1}; 
        actconch=act{b+1}-act{b-1}; 
        smkkch=smkka{b+1}-smkka{b-1}; smkstch=smksta{b+1}-smksta{b-1};
        mvitch=mvit{b+1}-mvit{b-1};   aspch=aspi{b+1}-aspi{b-1}; phmsch=phms{b+1}-phms{b-1}; horch=hor{b+1}-hor{b-1}; 
        hbpch=hbp{b+1}-hbp{b-1}; cholch=chol{b+1}-chol{b-1}; diabch=diab{b+1}-diab{b-1};
        nSESch=nSESArr{b+1}-nSESArr{b-1};
        aheich=aheim{b+1}-aheim{b-1}; calorch=calor{b+1}-calor{b-1};

          end;
        end;
        
if 0 <= actcona and actcona < 7  then postact=1;
      else if actcona>=7 and actcona<15 then postact=2;
      else if actcona>=15 and actcona<25 then postact=3;
      else if actcona>=25  then postact=4;
      else  postact=2;

%indic3(vbl=postact, prefix=postact, min=2, max=4, reflev=1, missing=., usemiss=0,
        label1='<7 total mets/wk',
        label2='7 to < 15 mets/wk',
        label3='15 to < 25 mets/wk',
        label4='25+ mets/wk');


%indic3(vbl=smkstaa, prefix=smkstaa, min=2, max=3, reflev=1, missing=., usemiss=1);
if aspa = . then aspa = 0;
if mvita = . then mvita = 0;


/****************  BASELINE EXCLUSIONS ********************************************/

if i=1 then do;
   %exclude(exrec eq 1);      *before 1991;
   %exclude(crc_cases eq 1);  * not in CRC file ;
   %exclude(conf1119 ne 1);   *no comfired cases;

   %exclude(0 lt dtdth   le irt{i});   
   %exclude(0 lt dtdxca  le irt{i});    
   %exclude(dtdxca ge cutoff); *cases diagnosed after the end of follow-up;

   %exclude(ltf eq irt{1});
   %exclude(birthday eq .);
   %exclude(age91 eq .);
   %exclude(stage_cr eq 4);    *stage IV CRC;
   %exclude(stage_cr eq .);    *stage missing CRC;
   %exclude(canpre eq 0, skip=t);  * skip individuals diagnosed after current time period; 
   %output();
  end;

/* EXCLUSIONS DURING FOLLOW-UP */
   if i>1 then do;
   %exclude(irt{i-1} lt dtdth le irt{i} );
   %exclude(stage_cr eq 4);    *stage IV CRC;   
   %exclude(stage_cr eq .);    *stage missing CRC;
   %exclude(canpre eq 0, skip=t);  * skip individuals diagnosed after current time period;   
   %output();
  end;
end;        

/* END OF DO-LOOP OVER TIME PERIODs */
%endex();








