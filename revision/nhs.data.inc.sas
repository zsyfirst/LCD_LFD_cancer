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
/* death */
data dead;
infile 'deadff.current.nhs1' lrecl=91 recfm=d;
    input
    @1              id         6.
    @8              icda  $    5.
    @15             mod        2.
    @17             yod        2.
    @21             icda10 $   5.    
    @28             qyr        2.
    @32             state  $   2.
    @34             source     1.
    @46             conf       2.
    @49             modx       2.
    @51             yodx       2.
    @54             nhsicda $  4.
    @88             dxmonthdd  4.
    ;
    label
    id              = 'id number'
    icda            = 'icda code for death'
    mod             = 'reported month of death'
    yod             = 'reported year of death'
    conf            = 'confirmation code'
    modx            = 'month of DIAGNOSIS'
    yodx            =  'year of DIAGNOSIS'
    nhsicda         = 'nhs icda code'
    dxmonthdd       = 'Date of death'
    ;
   dtdth=9999;
   if dxmonthdd > 0 then dtdth=dxmonthdd;
    newicda=compress(nhsicda, 'E');    
    newicda3=substr(newicda,1,3);
    /*if dtdth=9999 then delete*/
    proc sort; by id;
    keep id dtdth newicda newicda3;
run;

 
/* cancer 2022 */
data cancer;
infile 'canc7622.082823.cases' lrecl=53 recfm=d;

input
@1      id              6.
@7      cd              $ 1.
@10     cohort          1.
@13     icda            5.1
@20     qyr             2.
@26     mdx1            2.
@28     ydx1            2.
@32     conf            2.
@36     nhsicda         5.1
@43     mdx             2.
@45     ydx             2.
@49     dxmonth         4.
;

label
        id      = 'ID'
        cd      = 'Check Digit'
        cohort      = 'Cohort'
        icda        = 'Reported ICDA'
        qyr     = 'Questionnaire cycle year when reported'
        mdx1        = 'Reported month of diagnosis'
        ydx1        = 'Reported year of diagnosis'
        conf        = 'Confirmation Code'
        nhsicda     = 'NHS reviewed and confirmed ICDA code'
        mdx     = 'Confirmed Month of Diagnosis'
        ydx     = 'Confirmed Year of Diagnosis'
        dxmonth     = 'Date of dx months since 1900'
;

icda2=nhsicda;

if conf in (11,12,13,14,15,16,17,18,19) then ca=1; /* ca=1: confirmed cancer ; exclude in situ ca by taking out conf code 25*/
if ca=1 and icda2 ne . then caflag=1;
if dxmonth>0 and ca=1 and icda2 ne 173.2 and icda2 ne 173.1 then dtdxca=dxmonth;

        dtdx_allca=.;
        if dtdxca>0 and ca=1 and (140<=icda2<=172 or 174<=icda2<=207 or icda2 =225) then dtdx_allca=dtdxca;

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
     if obeca_nobcop=1 and dtdxca>0 and ca=1 and icda2 in (150,151, /*153,154,*/ 155,156, /*157,*/  /*174,*/  /*182, 183,*/ 189,193,203,225) then dtdx_obeca_nobcop=dtdxca;

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


/*  diabetes 2024 */
data diabetes;
infile 'db7624.110123.cases' lrecl=183 recfm=d;

input
@1		id		    6.
@20     qyr             2.
@32     conf            2.
@49     dxmonth         4.
@66     type            1.
@67     probabil        1.
@69     gestdm          1.
;

type2db=0;
if type=2 and probabil=1 then type2db=1; *type 2 diabetes;
type1db=0;
if type eq 1 or type eq 5 then type1db=1; /* definite or probable type 1 diabetes */
dtdxdb2=.;
if 0< dxmonth <9999 and type2db=1 then dtdxdb2=dxmonth;

keep id dtdxdb2; run;


/*  CABG data  */
*We do not follow CABG, All cases in this file (1976-2012) are reported;
data cabg;
    %include 'cabg7612.060816.nodups.input'; 

if dxmonth>0 then dtdxcabg=dxmonth;

keep id dtdxcabg;
run;

*Read reproted cabg from %nurXX;
%nur14(keep=cabg14 cabgd14 );
run;

%nur16(keep=cabg16 cabgd16 );
run;

*We did not sent out questionnaires in 2018. Instead, the questionnaires were sent out in 2019. The variables were labled as 20; 
%nur20(keep=cabg20  cabgd20);
run;

data update_cabg;
  merge cabg nur14 nur16 nur20;
  by id;

if dtdxcabg=. then do;
                              if cabg20=1 then do;
                                                 if cabgd20=1 then dtdxcabg=1404;
                                            else if cabgd20=2 then dtdxcabg=1422;
                                            else if cabgd20=3 then dtdxcabg=1434;
                                            else                   dtdxcabg=1422;
                                                end;

                              if cabg16=1 then do;
                                                 if cabgd16=1 then dtdxcabg=1374;
                                            else if cabgd16=2 then dtdxcabg=1386;
                                            else if cabgd16=3 then dtdxcabg=1398;
                                            else                   dtdxcabg=1386;
						end;

                              if cabg14=1 then do;
                                                 if cabgd14=1 then dtdxcabg=1350;
                                            else if cabgd14=2 then dtdxcabg=1362;
                                            else if cabgd14=3 then dtdxcabg=1374;
                                            else                   dtdxcabg=1362;
                                              end;
end;

if dtdxcabg>0 ;
cabg=1;

keep id dtdxcabg;
run;

/* stroke 2020 */
data tot_str;
    infile 'str7620.031621.cases' lrecl=81 recfm=d;

    input
     @1      id          6.
     @7      cd      $   1.
     @10     cohort          1.
     @13     icda    $        5.
     @20     qyr         2.
     @26     mdx1            2.
     @28     ydx1            2.
     @32     strconf            2.
     @36     nhsicda         5.
     @43     mdx         2.
     @45     ydx         2.
     @49     dxmonth         4.
     @55     rc_record_id        6.
     @63     recur           1.
     @66     death_record        1.
     @69     str_type        1.
     @70     lvl_conf        1.
     @71     primorsec       1.
     @72     thromb          1.
     @73     embolic         1.
     @74     vasc_abn        1.
     @75     ct_perf         1.
     @76     mri_perf        1.
     @77     autopsy         1.
     @78     lart_infr       1.
     @79     intra_hem       1.
     @80     iph_loc         1.
    ;

    label
        id      = 'ID'
        cd      = 'Check Digit'
        cohort      = 'Cohort'
        icda        = 'Reported ICDA'
        qyr     = 'Questionnaire cycle year when reported'
        mdx1        = 'Reported month of diagnosis'
        ydx1        = 'Reported year of diagnosis'
        strconf        = 'Confirmation Code'
        nhsicda     = 'NHS reviewed and confirmed ICDA code'
        mdx     = 'Confirmed Month of Diagnosis'
        ydx     = 'Confirmed Year of Diagnosis'
        dxmonth     = 'Date of dx months since 1900'
        rc_record_id    = 'REDCap Record ID'
        recur       = 'Possible recurrence'
        death_record    = 'From the death file?'
        str_type    = 'Stroke subtype'
        lvl_conf    = 'Level of confirmation'
        primorsec   = 'Primary or secondary'
        thromb      = 'Thrombotic distribution (thrombotic only)'
        embolic     = 'Embolic sources (embolic only)'
        vasc_abn    = 'Vascular abnormality (hemorrhagic only)'
        ct_perf     = 'CT performed (all stroke subtypes)'
        mri_perf    = 'MRI performed (all stroke subtypes)'
        autopsy     = 'Autopsy performed (all stroke subtypes)'
        lart_infr   = 'Large artery infarction territory (large artery infarction only)'
        intra_hem   = 'Intraventricular hemorrhage (IPH only)'
        iph_loc     = 'IPH location (IPH only)'
    ;

    fst=0;
    if 11<=strconf<=19 & death_record=1 then fst=1; 

    dtdxnfst=.;
    if dxmonth> 0 and fst=0 then dtdxnfst=dxmonth;

    dtdxfst=.;
    if dxmonth> 0 and fst=1 then dtdxfst=dxmonth;

    dtdxstr=.;
    if 0<dtdxnfst<9999 or 0<dtdxfst<9999 then dtdxstr=min(dtdxnfst,dtdxfst, .);

    *type of stroke cases;
    if str_type in (1,2) then hem_str=1; else hem_str=0;   /*hemorrhagic stroke cases*/
    if str_type in (3,4,5) then isc_str=1; else isc_str=0; /*ischemic stroke cases*/
    if str_type=6 then unk_str=1; else unk_str=0;          /*unknown type of stroke cases*/
    
    keep id dtdxstr;    
run;

/** MI 2020 **/
data mi;
    infile 'mi7620.073120.cases' lrecl=117 recfm=d;

    input
     @1      id          6.
     @7      cd      $   1.
     @10     cohort          1.
     @13     icda       $     5.
     @20     qyr         2.
     @26     mdx1            2.
     @28     ydx1            2.
     @32     conf        $    2.
     @36     nhsicda         5.1
     @43     mdx         2.
     @45     ydx         2.
     @49     dxmonth         4.
     @55     rc_record_id        6.
     @63     recur           1.
     @66     dth_subcode     1.
     @69     misite          1.
     @70     acang           1.
     @71     acthromb        1.
     @72     no_therapy      1.
     @73     ekgavail        1.
     @74     st          1.
     @75     extent          1.
     @76     enzavail        1.
     @77     peakck          5.
     @82     hck         4.
     @86     peakmb          6.1
     @92     hmb         2.
     @94     peaktrpn        6.1
     @100        htrpn           3.1
     @103        creat           4.1
     @107        creatna         1.
     @108        lvef            2.
     @110        lvefna          1.
     @111        hrtfail         1.
     @112        dismed          1.
     @113        ace         1.
     @114        betablk         1.
     @115        antiplt         1.
     @116        statin          1.
    ;

    label
        id      = 'ID'
        cd      = 'Check Digit'
        cohort      = 'Cohort'
        icda        = 'Reported ICDA'
        qyr     = 'Questionnaire cycle year when reported'
        mdx1        = 'Reported month of diagnosis'
        ydx1        = 'Reported year of diagnosis'
        conf        = 'Confirmation Code'
        nhsicda     = 'NHS reviewed and confirmed ICDA code'
        mdx     = 'Confirmed Month of Diagnosis'
        ydx     = 'Confirmed Year of Diagnosis'
        dxmonth     = 'Date of dx months since 1900'
        rc_record_id    = 'REDCap Record ID'
        recur       = 'Possible recurrence'
        dth_subcode = 'Confirmation of mi/chd'
        misite      = 'MI Site'
        acang       = 'Acute Therapy: Angioplasty'
        acthromb    = 'Acute Therapy: Thrombolytics'
        no_therapy  = 'Acute Therapy: Neither'
        ekgavail    = 'EKG Available'
        st      = 'ST'
        extent      = 'Extent'
        enzavail    = 'Enzymes Available'
        peakck      = 'Peak CK'
        hck     = 'Hospital ULN CK'
        peakmb      = 'Peak MB'
        hmb     = 'Hospital ULN MB'
        peaktrpn    = 'Peak Troponin'
        htrpn       = 'Hospital ULN Troponin'
        creat       = 'Initial creatinine'
        creatna     = 'initial creatinine: not available'
        lvef        = 'Left Ventricular Ejection Fraction % (last)'
        lvefna      = 'Left Ventricular Ejection Fraction'
        hrtfail     = 'Heart Failure during hospitalization'
        dismed      = 'Discharge medications available'
        ace     = 'ACE Inhibitors'
        betablk     = 'Beta-Blockers'
        antiplt     = 'Antiplatelet agents'
        statin      = 'Statins'
    ;

    if dxmonth> 0 then dtdxmi=dxmonth;
    
    if mod(qyr,2)=1 then fmi=1;else fmi=0;
    if mod(qyr,2)=0 then nfmi=1;else nfmi=0;
    if fmi=1 then dtdxfmi=dxmonth; 
    if nfmi=1 then dtdxnfmi=dxmonth; 
    keep id dtdxmi;    
run;   
   
/* disease outcomes */    
    data disease;
    merge dead cancer diabetes update_cabg mi tot_str;
    by id;

    /* CVD */
    dtdxchd=.;
    if 0<dtdxmi<9999 or  0<dtdxcabg<9999 then dtdxchd=min(dtdxmi,dtdxcabg, .);
    
    dtdxcvd=.;
    if 0<dtdxchd<9999 or  0<dtdxstr<9999 then dtdxcvd=min(dtdxchd,dtdxstr,.);

    keep id dtdth dtdxca dtdx_: newicda newicda3 dtdxdb2 dtdxcvd dtdxchd dtdxstr dtdxmi obeca_site ;
    run; proc sort; by id; run;

    PROC DATASETS;
    delete diabetes dead mi cabg update_cabg nur14 nur16 nur20 tot_str cancer; run;


*----------------2.Covariables-----------------------------;
** Endoscopy **;
*%include 'endoscopy_nhs_hpfs.sas';

/*  AHEI  */
%ahei2010_8410(keep= 
    ahei2010_noETOH84 ahei2010_noETOH86 ahei2010_noETOH90 ahei2010_noETOH94 ahei2010_noETOH98 ahei2010_noETOH02 ahei2010_noETOH06 ahei2010_noETOH10 
    ahei2010_frtI84 ahei2010_ssbI84 ahei2010_whgrnI84 ahei2010_nutI84 ahei2010_rmtI84 ahei2010_ptranI84 ahei2010_omegaI84 ahei2010_polyI84 ahei2010_naI84
    ahei2010_frtI86 ahei2010_ssbI86 ahei2010_whgrnI86 ahei2010_nutI86 ahei2010_rmtI86 ahei2010_ptranI86 ahei2010_omegaI86 ahei2010_polyI86 ahei2010_naI86
    ahei2010_frtI90 ahei2010_ssbI90 ahei2010_whgrnI90 ahei2010_nutI90 ahei2010_rmtI90 ahei2010_ptranI90 ahei2010_omegaI90 ahei2010_polyI90 ahei2010_naI90
    ahei2010_frtI94 ahei2010_ssbI94 ahei2010_whgrnI94 ahei2010_nutI94 ahei2010_rmtI94 ahei2010_ptranI94 ahei2010_omegaI94 ahei2010_polyI94 ahei2010_naI94
    ahei2010_frtI98 ahei2010_ssbI98 ahei2010_whgrnI98 ahei2010_nutI98 ahei2010_rmtI98 ahei2010_ptranI98 ahei2010_omegaI98 ahei2010_polyI98 ahei2010_naI98
    ahei2010_frtI02 ahei2010_ssbI02 ahei2010_whgrnI02 ahei2010_nutI02 ahei2010_rmtI02 ahei2010_ptranI02 ahei2010_omegaI02 ahei2010_polyI02 ahei2010_naI02
    ahei2010_frtI06 ahei2010_ssbI06 ahei2010_whgrnI06 ahei2010_nutI06 ahei2010_rmtI06 ahei2010_ptranI06 ahei2010_omegaI06 ahei2010_polyI06 ahei2010_naI06
    ahei2010_frtI10 ahei2010_ssbI10 ahei2010_whgrnI10 ahei2010_nutI10 ahei2010_rmtI10 ahei2010_ptranI10 ahei2010_omegaI10 ahei2010_polyI10 ahei2010_naI10
    ahei84 ahei86 ahei90 ahei94 ahei98 ahei02 ahei06 ahei10);

    ahei84=ahei2010_noETOH84; /*0-100*/
    ahei86=ahei2010_noETOH86;
    ahei90=ahei2010_noETOH90;
    ahei94=ahei2010_noETOH94;
    ahei98=ahei2010_noETOH98;
    ahei02=ahei2010_noETOH02;
    ahei06=ahei2010_noETOH06;
    ahei10=ahei2010_noETOH10;
run;

/* Physical activity - METs/week - 90 02 06 10 not calculated */
%act8614(keep= act86m  act88m  act92m  act94m  act96m act98m  act00m  act04m  act08m act12m act14m);

array actm{*} act86m act88m act92m act94m act96m act98m act00m act04m act08m act12m act14m;
    do i=1 to dim(actm); if actm{i} >900 then actm{i}=.; end;
    do i=2 to dim(actm);
        if actm{i}=.  then actm{i}=actm{i-1};
    end;
run;

/* family history of MI */
    %n767880(keep=  mmi76 fmi76 ammi76 afmi76
                    mothmi76 fathmi76 mifh76);
      if mmi76=1 and ammi76<=60 then mothmi76=1;
      if fmi76=1 and afmi76<=60 then fathmi76=1;
      if mothmi76=1 or fathmi76=1 then mifh76=1; else mifh76=0;
    %nur84(keep= id mmi84 ammi84 fmi84 afmi84  mothmi84 fathmi84 mifh84);
      if mmi84=2 and  1<=ammi84<=2 then mothmi84=1;
      if fmi84=2 and  1<=afmi84<=2 then fathmi84=1;
      if mothmi84=1 or fathmi84=1 then mifh84=1; else mifh84=0;

/*  Family hostory of cancer, diabetes */
    %nur82(keep= fdb82 sdb82 mdb82 bdb82 dbfh82 
                 mclc82 fclc82 sclc82 bclc82 mbrcn82 sbrcn82 mmel82 fmel82 smel82 bmel82 cafh82);
                 /* family history of cancer   */
                 cafh82=0; if mclc82=1 or fclc82=1 or sclc82=1 or bclc82=1 or mbrcn82=1 or sbrcn82=1 then cafh82=1;
                 /* family history of diabetes */
                 dbfh82=0; if fdb82=1 or mdb82=1 or sdb82=1 or bdb82=1 then dbfh82=1;
                 run;
    %nur88(keep= fdb88 sdb88 mdb88 bdb88 dbfh88
                 mclc88 fclc88 bclc88 sclc88 mbrcn88 sbrcn88 mcanc88 fcanc88 bcanc88 scanc88 cafh88
                 );
                 if fdb88=1 or mdb88=1 or sdb88=1 or bdb88=1 then dbfh88=1; else dbfh88=0;
                 if mclc88=1 or fclc88=1 or bclc88=1 or sclc88=1 or mbrcn88=1 or sbrcn88=1 or mcanc88=1 or fcanc88=1 or bcanc88=1 or scanc88=1 then cafh88=1; else cafh88=0;
                 run;

/*  Family history of cancer, diabetes & hypertension */
    %nur92(keep= fdb92 sdb92 mdb92 dbfh92 fhbp92 shbp92 mhbp92 hbpfh92
                 mbrcn92 mov92 mcolc92 mmel92 fbrcn92 fcolc92 fmel92 sbrcn92 sov92 scolc92 smel92 cafh92);
                 if fdb92=1 or mdb92=1 or sdb92=1 then dbfh92=1;
                 else dbfh92=0;
                 if fhbp92=1 or mhbp92=1 or shbp92=1 then hbpfh92=1;
                 else hbpfh92=0;
                 if mbrcn92=1 or mov92=1 or mcolc92=1 or mmel92=1 or fbrcn92=1 or fcolc92=1 or fmel92=1 or sbrcn92=1 or sov92=1 or scolc92=1 or smel92=1 then cafh92=1;
                 else cafh92=0;
                 run;

/* family history of MI, stroke, cancer */
    %nur96(keep= brmi96 brmid96 smi96 smid96 mifh96 
                 bstr96 bstrd96 sstr96 sstrd96 strfh96
                 mov96 sov96 mbrcn96 sbrc196 sbrc296 mut96 sut96 pclc96 sclc196 sclc296 fpro96 bpro196
                 bpro296 ppan96 span96 pmel96 smel96 cafh96);
      if (brmi96=1 and brmid96 in (1,2)) or (smi96=1 and smid96 in (1,2)) then mifh96=1; else mifh96=0;
      if (bstr96=1 and bstrd96 in (1,2)) or (sstr96=1 and sstrd96 in (1,2)) then strfh96=1; else strfh96=0;
      if pclc96=1  or sclc196=1 or sclc296=1 or
      fpro96=1  or bpro196=1 or bpro296=1 or mov96=1  or sov96=1  or 
      mbrcn96=1 or sbrc196=1 or sbrc296=1 or 
      mut96=1   or sut96=1 or ppan96=1  or span96=1  or
      pmel96=1  or smel96=1 
      then cafh96=1; else cafh96=0;
      run;

/* family history of MI, stroke, cancer */
   %nur08(keep= brmi08 brmid08 smi08 smid08 omi08 omid08 mifh08 
                mstr08 mstrd08 fstr08 fstrd08 sstr08 sstrd08 strfh08
                mov08 sov08 sbrc108 sbrc208 sclc108 sclc208 ppan08 span08 pmel08 smel08 mut08 sut08 pkdc08 skidc08 cafh08);
      if (brmi08=1 and brmid08 in (1,2)) or (smi08=1 and smid08 in (1,2)) or (omi08=1 and omid08 in (1,2)) then mifh08=1; else mifh08=0;
      if (mstr08=1 and mstrd08 in (1,2)) or (fstr08=1 and fstrd08 in (1,2)) or (sstr08=1 and sstrd08 in (1,2)) then strfh08=1; else strfh08=0;
      if mov08=1 or sov08=1 or sbrc108=1 or sbrc208=1 or sclc108=1 or sclc208=1 or ppan08=1 or span08=1 or pmel08=1 or smel08=1 or mut08=1 or sut08=1 or pkdc08=1 or skidc08=1 then cafh08=1; else cafh08=0;run;

data familyhis;
      merge  n767880 nur82 nur84 nur88  nur92 nur96 nur08;
      by id;
      dbfh=0;  if dbfh82=1 or dbfh88=1 or dbfh92=1                 then dbfh=1;
      hbpfh=0; if hbpfh92=1                                        then hbpfh=1;
      mifh=0;  if mifh76=1 or mifh84=1 or mifh96=1 or mifh08=1     then mifh=1; /* note: family history of MI before 60 years, ref: Fung T AJCN 2009paper */
      strfh=0; if strfh96=1 or strfh08=1                           then strfh=1;
      cafh=0;  if cafh82=1 or cafh88=1 or cafh92=1 or cafh96=1 or cafh08=1 then cafh=1;     
      keep id dbfh mifh strfh cafh;
      proc sort; by id; run;

    PROC DATASETS;
    delete n767880 nur82 nur84 nur88  nur92 nur96 nur08;
    RUN;

/*   other covariables   */
%der7620(keep= mobf                                      /*month of birth*/
               yobf                                      /*year  of birth*/
               
               irt76 irt78 irt80 irt82 irt84 irt86 irt88 irt90 irt92 irt94 irt96 
               irt98 irt00 irt02 irt04 irt06 irt08 irt10 irt12 irt14 irt16 /* returned date */
               
               age76 age78 age80 age82 age84 age86 age88 age90 age92 
               age94 age96 age98 age00 age02 age04 age06 age08 age10 age12 age14 age16

               bmiage18 bmi76 bmi78 bmi80 bmi82 bmi84 bmi86 bmi88 bmi90 bmi92 
               bmi94 bmi96 bmi98 bmi00 bmi02 bmi04 bmi06 bmi08 bmi10 bmi12 bmi14 bmi16
               
               race9204  /* 1.White; 2.Black; 3.American Indian; 4.Asian; 5.Hawaiian */
               mrace9204 /* 1.White;\ 2.Black;\ 3.American Indian;\4.Asian;\5.Hawaiian;\6.Other/unknown;\7.Multi-racial */
               eth9204   /*1.Hispanic;\ 2.Not Hispanic */    
               race mrace white eth3g

               nhor76 nhor78 nhor80 nhor82 nhor84 nhor86 nhor88 nhor90 nhor92 nhor94 
               nhor96 nhor98 nhor00 nhor02 nhor04 nhor06 nhor08 nhor10 nhor12 nhor14 nhor16
               /* New PMH use, $label 1.Pre/Missing Mnp; 2.Never Used; 3.Current User; 4.Past User; 5.Status Unknown; */

               men76 men78 men80 men82 men84 men86 men88 men90 men92 men94 men96 men98 men00 men02 men04 
               /* Simple menopause status,  $label 1.Premenopause; 2.Post Menopause; 3.Missing */

               dmnp76 dmnp78 dmnp80 dmnp82 dmnp84 dmnp86 dmnp88 dmnp90 dmnp92 dmnp94 dmnp96 dmnp98 dmnp00 dmnp02 dmnp04
               /* Dubious menopause $label  1.Premenopause; 2.Post Menopause; 3.Dubious Menopause; 4.Missing */

               ocu76 ocu78 ocu80 ocu82 ocu84 /* Oral contraceptive use $lable 1. current; 2.past; 3.never; 4.missing */
               
               can76 can78 can80 can82 can84 can86 can88 can90 can92 can94 can96 can98 can00 
               can02 can04 can06 can08 can10 can12 can14 can16 /* 1=yes */

               hrt76 hrt78 hrt80 hrt82 hrt84 hrt86 hrt88 hrt90 hrt92 hrt94 hrt96 hrt98 hrt00 
               hrt02 hrt04 hrt06 hrt08 hrt10 hrt12 hrt14 hrt16 /* 1.MI; 2.Angina; 3.MI+Angina */

               smkdr76 smkdr78 smkdr80 smkdr82 smkdr84 smkdr86 smkdr88 smkdr90 smkdr92 smkdr94 smkdr96 
               smkdr98 smkdr00 smkdr02 smkdr04 smkdr06 smkdr08 smkdr10 smkdr12 smkdr14 smkdr16
               /* smoking status range 0-15 0.missing 1.never 2-8.past 9-15.curr*/
               );

             race=.; /* 1. white; 2.asian; 3.black; 4.hispanic*/
             if  eth9204 =1 then race=4;
      else   if  race9204=1 then race=1;
      else   if  race9204=4 then race=2;
      else   if  race9204=2 then race=3;

             mrace=.; /* 1. white; 2.black; 3.hispanic; 4.others*/
             if  eth9204 =1  then mrace=3;
      else   if  mrace9204=1 then mrace=1;
      else   if  mrace9204=2 then mrace=2;
      else   mrace=4;

             eth3g=.; /* 1. white; 2.black; 3.others */
             if  race9204=1 then eth3g=1;
      else   if  race9204=2 then eth3g=2;
      else  eth3g=3;

      if     race9204=1 then white=1; else white=0;
      proc sort; by id; run;

%n767880 (keep= db76 hbp76 chol76
                db78 hbp78 chol78
                db80 hbp80 chol80);

%nur82(keep=db82 hbp82 chol82 str82 ); /* cabg82 not asked */

%nur84(keep=  db84 hbp84 chol84 cabg84 str84); /* CODING CHANGE: 0=blank, 1=no, 2=yes, range=0-2 */
       array fix db84 hbp84 chol84 cabg84 str84;
       do over fix; fix=fix-1; end;

%nur86(keep=db86 hbp86 chol86 cabg86 str86 );

%n86_dt(keep=  waist86  hip86 nblnk86d  cerbr86d);
%serv86 (keep=
/* whole grains */   cer_s86  ckcer_s86 oat_s86 dkbr_s86 brice_s86 otgrn_s86 bran_s86 wgerm_s86 popc_s86  
/* refined grains */ whbr_s86 engl_s86 muff_s86 wrice_s86 pasta_s86 pcake_s86  
/* fruits */         rais_s86 prune_s86 ban_s86 cant_s86 h2om_s86 avo_s86 appl_s86 apsau_s86 pear_s86 cpear_s86 peach_s86 cpeac_s86 fcock_s86 oran_s86 grfr_s86 straw_s86 blueb_s86  
/* vegetables */     tom_s86 toj_s86 tosau_s86 chsau_s86 brocc_s86 kraut_s86 ccabb_s86 slaw_s86 cauli_s86 bruss_s86 rcar_s86 ccar_s86 mixv_s86 yam_s86 ysqua_s86 eggpl_s86 kale_s86 cspin_s86 rspin_s86 ilett_s86 rlett_s86 cel_s86 mush_s86 alfsp_s86 grpep_s86 cuke_s86 corn_s86
/* legumes */        sbean_s86 tofu_s86 bean_s86 peas_s86
/* tea-coffee */     tea_s86 dcaf_s86 coff_s86 

/* potatoes */  fries_s86 pot_s86 pchip_s86
/* dairy */ skim_s86 whole_s86 sherb_s86 icecr_s86 yog_s86 cream_s86 sour_s86 cotch_s86 crmch_s86 otch_s86
/* eggs */ egg_s86
/* fish */ ctuna_s86 dkfsh_s86 ofish_s86 shrim_s86
/* fruit juices */ aj_s86 oj_s86 grj_s86 othj_s86  

/* ssb */  soda_s86 punch_s86   
/* sweets */ choco_s86 cdyw_s86 cdywo_s86 cokh_s86 cokr_s86 brwn_s86 donut_s86 cakh_s86 cakr_s86 pieh_s86 pier_s86 sroll_s86 jam_s86
/* diet bev */ local_s86
/* alcoholic bev */ beer_s86 wwine_s86 rwine_s86 liq_s86
/* meat */ chwi_s86 chwo_s86 bacon_s86 dog_s86 procm_s86 hamb_s86 bmix_s86 bmain_s86 livb_s86 livc_s86
/* rpmeats_s86 ssb_s86 swtsvg_s86 swtsnochoc_s86 coff_s86 whgrn_s86 frui_s86 vege_s86 fv_s86 */
);
data serv86;
  merge serv86 n86_dt;
  by id;

if cerbr86d in ( 0, 4, 11, 12, 13, 15, 16, 17, 19, 20, 21, 22, 23, 25, 29, 30, 31, 33, 35,37,38, 42, 43, 45, 46, 47,
59, 60,61,62,63,64,65,66,68,73,74,75,76,79,80,81,82,83, 84,85,86,88,95,96,97,98,
101,109,111, 112,113,114, 115,119, 123,124, 127, 128, 129,130,131, 132, 133, 134, 146, 147, 148,149,150,
151, 152,153, 154, 155, 157,158, 164,165, 166, 167,170,173, 176, 177,182, 184, 185, 186, 188/*,
190,197,198,997*/)
   then rcer_s86=cer_s86;
   else wcer_s86=cer_s86;
   if rcer_s86=. then rcer_s86=0;
   if wcer_s86=. then wcer_s86=0;

/* take sum of caffeinated and decaffeinated tea & coffee */
tea_s86  = SUM(tea_s86);
coff_s86  = sum(coff_s86, dcaf_s86);
/** whole grains **/
whgrn_s86 =sum(wcer_s86, oat_s86, ckcer_s86, dkbr_s86, brice_s86, otgrn_s86, bran_s86, wgerm_s86, popc_s86);
/** refined grains **/
rgrain_s86 =sum(rcer_s86, whbr_s86, engl_s86, muff_s86, wrice_s86, pasta_s86, pcake_s86);
/* fruit and vegetable*/
frui_s86 =sum(rais_s86,  avo_s86,  prune_s86,  ban_s86,  cant_s86,  h2om_s86,  appl_s86,  apsau_s86,  pear_s86,  
         cpear_s86,  peach_s86,  cpeac_s86,  fcock_s86,  oran_s86,  grfr_s86,  straw_s86,  blueb_s86);
vege_s86 =sum(sbean_s86,brocc_s86,  kraut_s86,    slaw_s86,  ccabb_s86,  cauli_s86,  bruss_s86,  rcar_s86,     ccar_s86,    /*corn_s86,*/   peas_s86,   mixv_s86,     
                bean_s86,  alfsp_s86,    cel_s86,    mush_s86,  ysqua_s86,  eggpl_s86,    yam_s86, cspin_s86,    rspin_s86,    kale_s86,  ilett_s86,  rlett_s86,   
                grpep_s86,   cuke_s86,    tom_s86,     toj_s86,  tosau_s86,  chsau_s86,   tofu_s86 );
fv_s86 =sum(frui_s86, vege_s86);

/* potatoes */
pota_s86=sum(fries_s86, pot_s86, pchip_s86);
/*legumes */
leg_s86=sum(sbean_s86, tofu_s86, bean_s86, peas_s86);
/* fruit juices */
fj_s86=sum(aj_s86, oj_s86, grj_s86, othj_s86  );
/* vegetables minus potatoes and legumes */
veg_noleg_s86=vege_s86-leg_s86;
/* dairy and fish */
dairy_s86=sum(skim_s86, whole_s86, sherb_s86, icecr_s86, yog_s86, cream_s86, sour_s86 ,cotch_s86 ,crmch_s86, otch_s86);
fish_s86=sum(ctuna_s86, dkfsh_s86, ofish_s86, shrim_s86);
/* red and processed meat */
rpmeats_s86= sum(bacon_s86, dog_s86, procm_s86, hamb_s86, bmix_s86, bmain_s86);
/* poultry */
poul_s86=sum(chwi_s86, chwo_s86);
/* ssb*/
ssb_s86  = sum( soda_s86, punch_s86);
/*** SWEETS & DESERTS  ***/
swtsvg_s86  = sum(choco_s86, cdyw_s86, cdywo_s86, cokh_s86, cokr_s86, brwn_s86, donut_s86, cakh_s86, cakr_s86, pieh_s86, pier_s86, sroll_s86, jam_s86);
swtsnochoc_s86  = sum( cdyw_s86, cdywo_s86, cokh_s86, cokr_s86, brwn_s86, donut_s86, cakh_s86, cakr_s86, pieh_s86, pier_s86, sroll_s86, jam_s86 );
run; 

%nur88(keep=db88 hbp88 chol88 str88  cabg88
            thiaz88 ccblo88 betab88 ace88 bprx88 /* medications asked */
            clrx88
            insul88 odbrx88 /* Hypoglycemic */
            batch88 sigm88 
            mamm88 ymamm88 
            );

%nur90(keep=db90 hbp90 chol90 cabg90 str90
            batch90 sigm90 sigmr90 sigm190 vbld90 diarr90 fecal90 
            hxclc90 abdpn90 rout90 q8b88pt 
            mamm90 );

%serv90(keep= 
/* whole grains */  cer_s90 ckcer_s90 oat_s90 dkbr_s90 brice_s90 otgrn_s90 bran_s90 oatbr_s90 wgerm_s90 popc_s90  
/* fruits */ rais_s90 prune_s90 ban_s90 cant_s90 h2om_s90 appl_s90 apsau_s90 peach_s90 oran_s90 grfr_s90 straw_s90 blueb_s90  
/* vegetables */ tom_s90 toj_s90 tosau_s90 chsau_s90 brocc_s90 cabb_s90 cauli_s90 bruss_s90 rcar_s90 ccar_s90 mixv_s90 yam_s90 osqua_s90 eggpl_s90 kale_s90 cspin_s90 rspin_s90 ilett_s90 rlett_s90 cel_s90 beet_s90 oniog_s90 oniov_s90 corn_s90
/* legumes */ sbean_s90 tofu_s90 bean_s90 peas_s90  
/* tea-coffee */ tea_s90 dcaf_s90 coff_s90 

/* potatoes */  fries_s90 pot_s90 pchip_s90
/* dairy */ skim_s90 whole_s90 sherb_s90 icecr_s90 yog_s90 cream_s90 sour_s90 cotch_s90 crmch_s90 otch_s90
/* eggs */ egg_s90
/* fish */ ctuna_s90 dkfsh_s90 ofish_s90 shrim_s90
/* fruit juices */ aj_s90 oj_s90 grj_s90 othj_s90  

/* ssb */  cola_s90 cnoc_s90 otsug_s90 punch_s90   
/* sweets */ choco_s90 cdyw_s90 cdywo_s90 cokh_s90 cokr_s90 brwn_s90 donut_s90 cakh_s90 cakr_s90 pieh_s90 pier_s90 srolr_s90 srolh_s90 jam_s90
/* diet bev */ lccaf_s90 lcnoc_s90 lcoth_s90
/* alcoholic bev */ beer_s90 wwine_s90 rwine_s90 liq_s90
/* fish */ ctuna_s90 dkfsh_s90 ofish_s90 shrim_s90
/* meat */ chwi_s90 chwo_s90 bacon_s90 dog_s90 procm_s90 hamb_s90 bmix_s90 pmain_s90 bmain_s90 livb_s90 livc_s90
/* rpmeats_s90 ssb_s90 swtsvg_s90 swtsnochoc_s90 whgrn_s90 frui_s90 vege_s90 fv_s90 */
/* vege_s90 frui_s90 fv_s90 rpmeats_s90 whgrn_s90 cerbr_s90 ssb_s90 wcer_s90 swtsvg_s90 */
);

%n90_dt(keep= cerbr90d );

data serv90;
  merge serv90 n90_dt;
  by id;

if cerbr90d in ( 0, 4, 11, 12, 13, 15, 16, 17, 19, 20, 21, 22, 23, 25, 29, 30, 31, 33, 35,37,38, 42, 43, 45, 46, 47,
59, 60,61,62,63,64,65,66,68,73,74,75,76,79,80,81,82,83, 84,85,86,88,95,96,97,98,
101,109,111, 112,113,114, 115,119, 123,124, 127, 128, 129,130,131, 132, 133, 134, 146, 147, 148,149,150,
151, 152,153, 154, 155, 157,158, 164,165, 166, 167,170,173, 176, 177,182, 184, 185, 186, 188/*,
190,197,198,997*/)
   then rcer_s90=cer_s90;
   else wcer_s90=cer_s90;
   if rcer_s90=. then rcer_s90=0;
   if wcer_s90=. then wcer_s90=0;

/* take sum of caffeinated and decaffeinated tea & coffee */
tea_s90  = SUM(tea_s90);
coff_s90  = sum(coff_s90, dcaf_s90);
/** whole grains **/
whgrn_s90 =sum(wcer_s90, oat_s90, ckcer_s90, dkbr_s90, brice_s90, otgrn_s90, oatbr_s90, bran_s90, wgerm_s90);
/* fruit and vegetable*/
frui_s90 =sum(rais_s90, prune_s90, ban_s90, cant_s90, h2om_s90, appl_s90, apsau_s90, peach_s90, oran_s90, grfr_s90, straw_s90, blueb_s90);
vege_s90 =sum(tom_s90, toj_s90, tosau_s90, chsau_s90, tofu_s90,sbean_s90,brocc_s90, cabb_s90,cauli_s90,bruss_s90, rcar_s90,ccar_s90, /*corn_s90,*/
                 peas_s90, mixv_s90, bean_s90,osqua_s90,eggpl_s90,  yam_s90,cspin_s90,rspin_s90,kale_s90,ilett_s90,  
                rlett_s90,  cel_s90, beet_s90,oniog_s90,oniov_s90);
fv_s90 =sum(frui_s90, vege_s90);
/* potatoes */
pota_s90=sum(fries_s90, pot_s90, pchip_s90);
/*legumes */
leg_s90=sum(sbean_s90, tofu_s90, bean_s90, peas_s90);
/* fruit juices */
fj_s90=sum(aj_s90, oj_s90, grj_s90, othj_s90  );
/* vegetables minus potatoes and legumes */
veg_noleg_s90=vege_s90-leg_s90;
/* dairy and fish */
dairy_s90=sum(skim_s90, whole_s90, sherb_s90, icecr_s90, yog_s90, cream_s90, sour_s90 ,cotch_s90 ,crmch_s90, otch_s90);
fish_s90=sum(ctuna_s90, dkfsh_s90, ofish_s90, shrim_s90);
/* red and processed meat */
rpmeats_s90= sum(bacon_s90, dog_s90, procm_s90, hamb_s90, bmix_s90, pmain_s90, bmain_s90);
/* poultry */
poul_s90=sum(chwi_s90, chwo_s90);
/* ssb*/
ssb_s90  = sum( cola_s90, cnoc_s90, otsug_s90, punch_s90);
/*** SWEETS & DESERTS  ***/
swtsvg_s90  = sum(choco_s90, cdyw_s90, cdywo_s90, cokh_s90, cokr_s90, brwn_s90, donut_s90, cakh_s90, cakr_s90, pieh_s90, pier_s90, srolr_s90, srolh_s90, jam_s90);
swtsnochoc_s90  = sum( cdyw_s90, cdywo_s90, cokh_s90, cokr_s90, brwn_s90, donut_s90, cakh_s90, cakr_s90, pieh_s90, pier_s90, srolr_s90, srolh_s90, jam_s90);
run; 

%nur92(keep=db92 hbp92 chol92 cabg92 str92 alone92
            batch92 sigm92  mamm92 ymamm92 );

%nur94(keep=db94 hbp94 chol94 cabg94 str94 
            thiaz94 lasix94 ccblo94 betab94 bprx94 
            couma94 digox94 
            clrx94 
            tag94 zantc94 /* PPIs */
            insul94 ohypo94 /* Hypoglycemic */
            batch94 sigm94 vbld94 hxclc94 fecal94 abdpn94 diarr94 rout94 q11ap94 
            mam8894 mam8994 mam9094 mam9194 mam9294 mam9394 mam9494 tblg94 db94 wt94 fbbd94 dbd94 
            ndmam94 mab8894 
            );

%serv94(keep= 
/* whole grains */  cer_s94 ckcer_s94 oat_s94 dkbr_s94 brice_s94 otgrn_s94 bran_s94 oatbr_s94 wgerm_s94 popc_s94 tort_s94  
/* fruits */ rais_s94 prune_s94 ban_s94 cant_s94 avo_s94 appl_s94 peach_s94 oran_s94 grfr_s94 straw_s94 blueb_s94  
/* vegetables */ tom_s94 toj_s94 tosau_s94 salsa_s94 brocc_s94 cabb_s94 cauli_s94 bruss_s94 rcar_s94 ccar_s94 mixv_s94 yam_s94 osqua_s94 eggpl_s94 kale_s94 cspin_s94 rspin_s94 ilett_s94 rlett_s94 cel_s94 grpep_s94 oniog_s94 oniov_s94 garl_s94 corn_s94
/* legumes */ sbean_s94 tofu_s94 bean_s94 peas_s94  

/* potatoes */  fries_s94 pot_s94 pchip_s94
/* dairy */ skim_s94 m1or2_s94 whole_s94 sherb_s94 icecr_s94 flyog_s94 plyog_s94 cream_s94 cotch_s94 crmch_s94 otch_s94
/* eggs */ egg_s94
/* fish */ ctuna_s94 dkfsh_s94 ofish_s94 shrim_s94
/* fruit juices */ aj_s94 oj_s94 grj_s94 othj_s94  

/* tea-coffee */ tea_s94 dcaf_s94 coff_s94 
/* ssb */  cola_s94 cnoc_s94 otsug_s94 punch_s94   
/* sweets */ choco_s94 cdyw_s94 cdywo_s94 cokh_s94 cokr_s94 brwn_s94 donut_s94 cakh_s94 cakr_s94 pieh_s94 pier_s94 srolr_s94 srolh_s94 jam_s94
/* diet bev */ lccaf_s94 lcnoc_s94 lcoth_s94
/* alcoholic bev */ beer_s94 lbeer_s94 wwine_s94 rwine_s94 liq_s94
/* meat */ chwi_s94 chwo_s94 bacon_s94 bfdog_s94 ctdog_s94 procm_s94 hamb_s94 hambl_s94 bmix_s94 pmain_s94 bmain_s94 livb_s94 livc_s94
/* rpmeats_s94 ssb_s94 swtsvg_s94 swtsnochoc_s94 whgrn_s94 frui_s94 vege_s94 fv_s94  */
/*vege_s94 frui_s94 fv_s94 rpmeats_s94 whgrn_s94 cerbr_s94 ssb_s94 wcer_s94 swtsvg_s94 */
);

%n94_dt(keep=cerbr94d);

data serv94;
  merge serv94 n94_dt;
  by id;

if cerbr94d in ( 0, 4, 11, 12, 13, 15, 16, 17, 19, 20, 21, 22, 23, 25, 29, 30, 31, 33, 35,37,38, 42, 43, 45, 46, 47,
59, 60,61,62,63,64,65,66,68,73,74,75,76,79,80,81,82,83, 84,85,86,88,95,96,97,98,
101,109,111, 112,113,114, 115,119, 123,124, 127, 128, 129,130,131, 132, 133, 134, 146, 147, 148,149,150,
151, 152,153, 154, 155, 157,158, 164,165, 166, 167,170,173, 176, 177,182, 184, 185, 186, 188/*,
190,197,198,997*/)
   then rcer_s94=cer_s94;
   else wcer_s94=cer_s94;
   if rcer_s94=. then rcer_s94=0;
   if wcer_s94=. then wcer_s94=0;

/* take sum of caffeinated and decaffeinated tea & coffee */
tea_s94  = SUM(tea_s94);
coff_s94  = sum(coff_s94, dcaf_s94);
/** whole grains **/
whgrn_s94 =sum(wcer_s94, oat_s94,ckcer_s94, dkbr_s94,brice_s94,otgrn_s94,oatbr_s94,bran_s94,wgerm_s94);

/* fruit and vegetable*/
frui_s94 =sum(rais_s94, prune_s94, ban_s94, cant_s94, avo_s94, appl_s94, peach_s94, oran_s94, grfr_s94, straw_s94, blueb_s94);
vege_s94 =sum(tom_s94  ,   toj_s94  ,   tosau_s94,   salsa_s94,   tofu_s94 ,   sbean_s94,   brocc_s94,   
               cabb_s94 ,   cauli_s94,   bruss_s94,   rcar_s94 ,   ccar_s94 ,   /*corn_s94 ,*/   peas_s94 ,   
               mixv_s94 ,   bean_s94 ,   osqua_s94,   eggpl_s94,   yam_s94  ,   cspin_s94,   rspin_s94,   
               kale_s94 ,   ilett_s94,   rlett_s94,   cel_s94  ,   grpep_s94,   oniog_s94,   oniov_s94 );
fv_s94 =sum(frui_s94, vege_s94);

/* potatoes */
pota_s94=sum(fries_s94, pot_s94, pchip_s94);
/*legumes */
leg_s94=sum(sbean_s94, tofu_s94, bean_s94, peas_s94);
/* fruit juices */
fj_s94=sum(aj_s94, oj_s94, grj_s94, othj_s94  );
/* vegetables minus potatoes and legumes */
veg_noleg_s94=vege_s94-leg_s94;
/* dairy and fish */
dairy_s94=sum(skim_s94, m1or2_s94, whole_s94, sherb_s94, icecr_s94, flyog_s94, plyog_s94, cream_s94, cotch_s94, crmch_s94, otch_s94);
fish_s94=sum(ctuna_s94, dkfsh_s94, ofish_s94, shrim_s94);
/* poultry */
poul_s94=sum(chwi_s94, chwo_s94, ctdog_s94 );
/* red and processed meat */
rpmeats_s94= sum(bacon_s94, bfdog_s94, procm_s94, hamb_s94, hambl_s94, bmix_s94, pmain_s94, bmain_s94);
/* ssb*/
ssb_s94  = sum( cola_s94, cnoc_s94, otsug_s94, punch_s94);
/*** SWEETS & DESERTS  ***/
swtsvg_s94  = sum(choco_s94, cdyw_s94, cdywo_s94, cokh_s94, cokr_s94, brwn_s94, donut_s94, cakh_s94, cakr_s94, pieh_s94, pier_s94, srolr_s94, srolh_s94, jam_s94);
swtsnochoc_s94  = sum( cdyw_s94, cdywo_s94, cokh_s94, cokr_s94, brwn_s94, donut_s94, cakh_s94, cakr_s94, pieh_s94, pier_s94, srolr_s94, srolh_s94, jam_s94);
run; 


%nur96(keep=db96 hbp96 chol96 cabg96 str96  alone96  marry96 
            thiaz96  lasix96  ccblo96  betab96  ace96  bprx96  
            couma96  digox96  antia96 /* anticoagulant  heart-failure antiarrhythmic */
            chrx96
            antid96 /* antidepressant */
            tag96 h2blk96 /* PPIs */
            insul96 ohypo96 /* Hypoglycemic */
            batch96 sigm96 vbld96 hxclc96 fecal96 abdpn96 diarr96 rout96 q11ap96 
            mamno96 mamsc96 mamsy96 mampt96 
            );

%nur98(keep=db98 hbp98 chol98 cabg98 str98  
            thiaz98  lasix98  ccblo98  betab98  ace98 bprx98 
            couma98   digox98   antia98 
            chrx98
            antid98 /* antidepressant */
            tag98 h2blk98 /* PPIs */
            insul98 ohypo98 /* Hypoglycemic */
            q98 sigm98 vbld98 hxclc98 fecal98 abdpn98 diarr98 rout98 bariu98 q13ap98 
            mamno98 mamsc98 mamsy98 mampt98 
            );

%serv98(keep= 
/* whole grains */  cer_s98 ckcer_s98 oat_s98 dkbr_s98 brice_s98 otgrn_s98 bran_s98 oatbr_s98 wgerm_s98 popc_s98 tort_s98 
/* fruits */ rais_s98 prune_s98 ban_s98 cant_s98 avo_s98 appl_s98 peach_s98 oran_s98 grfr_s98 straw_s98 blueb_s98  
/* vegetables */ tom_s98 toj_s98 tosau_s98 salsa_s98 brocc_s98 cabb_s98 cauli_s98 bruss_s98 rcar_s98 ccar_s98 mixv_s98 yam_s98 osqua_s98 eggpl_s98 kale_s98 cspin_s98 rspin_s98 ilett_s98 rlett_s98 cel_s98 grpep_s98 oniog_s98 oniov_s98 garl_s98 corn_s98
/* legumes */ sbean_s98 tofu_s98 bean_s98 peas_s98 soy_s98  

/* potatoes */  fries_s98 pot_s98 pchip_s98
/* dairy */ skim_s98 m2_s98 whole_s98 sherb_s98 icecr_s98 flyog_s98 plyog_s98 cream_s98 cotch_s98 crmch_s98 otch_s98
/* eggs */ egg_s98
/* fish */ ctuna_s98 dkfsh_s98 ofish_s98 shrim_s98
/* fruit juices */ aj_s98 oj_s98 ojca_s98 grj_s98 othj_s98  

/* tea-coffee */ tea_s98 dtea_s98 dcaf_s98 coff_s98 
/* ssb */  cola_s98 otsnc_s98 otsug_s98 punch_s98   
/* sweets */ choco_s98 cdyw_s98 cdywo_s98 coknf_s98 cokh_s98 cokr_s98 brwn_s98 donut_s98 cakr_s98 cakh_s98 pieh_s98 srolf_s98 srolr_s98 srolh_s98 jam_s98
/* alcoholic bev */ beer_s98 lbeer_s98 wwine_s98 rwine_s98 liq_s98
/* meat */ chksa_s98 chwi_s98 chwo_s98 bacon_s98 bfdog_s98 ctdog_s98 pmsan_s98 procm_s98 hamb_s98 hambl_s98 bmix_s98 pmain_s98 bmain_s98 livb_s98 livc_s98
/* rpmeats_s98 ssb_s98 swtsvg_s98 swtsnochoc_s98 whgrn_s98 frui_s98 vege_s98 fv_s98 */
/*vege_s98 frui_s98 fv_s98 rpmeats_s98 whgrn_s98 cerbr_s98 ssb_s98 wcer_s98 swtsvg_s98 */

);

%n98_dt(keep=cerbr98d);

data serv98;
  merge serv98 n98_dt;
  by id;

if cerbr98d in ( 0, 4, 11, 12, 13, 15, 16, 17, 19, 20, 21, 22, 23, 25, 29, 30, 31, 33, 35,37,38, 42, 43, 45, 46, 47,
59, 60,61,62,63,64,65,66,68,73,74,75,76,79,80,81,82,83, 84,85,86,88,95,96,97,98,
101,109,111, 112,113,114, 115,119, 123,124, 127, 128, 129,130,131, 132, 133, 134, 146, 147, 148,149,150,
151, 152,153, 154, 155, 157,158, 164,165, 166, 167,170,173, 176, 177,182, 184, 185, 186, 188/*,
190,197,198,997*/)
   then rcer_s98=cer_s98;
   else wcer_s98=cer_s98;
   if rcer_s98=. then rcer_s98=0;
   if wcer_s98=. then wcer_s98=0;

/* take sum of caffeinated and decaffeinated tea & coffee */
tea_s98  = SUM(tea_s98, dtea_s98);
coff_s98  = sum(dcaf_s98, coff_s98);
/** whole grains **/
whgrn_s98 =sum(wcer_s98, oat_s98, ckcer_s98, dkbr_s98, brice_s98, otgrn_s98, oatbr_s98, bran_s98, wgerm_s98);
/* fruit and vegetable*/
frui_s98 =sum(rais_s98, prune_s98, ban_s98, cant_s98, avo_s98, appl_s98, peach_s98, oran_s98, grfr_s98, straw_s98, blueb_s98);
vege_s98 =sum(tom_s98   ,  toj_s98   ,  tosau_s98 ,  salsa_s98 ,  tofu_s98  ,  sbean_s98 , soy_s98, brocc_s98 ,  cabb_s98  ,  
               cauli_s98 ,  bruss_s98 ,  rcar_s98  ,  ccar_s98  ,  /*corn_s98,*/  peas_s98  ,  mixv_s98  ,  bean_s98  ,  
                yam_s98  ,  osqua_s98 ,  eggpl_s98 ,  kale_s98  ,  cspin_s98 ,  rspin_s98 ,  ilett_s98 ,  rlett_s98 ,  
                 cel_s98   ,  grpep_s98 ,  oniog_s98 ,  oniov_s98);
fv_s98 =sum(frui_s98, vege_s98);
/* potatoes */
pota_s98=sum(fries_s98, pot_s98, pchip_s98);
/*legumes */
leg_s98=sum(sbean_s98, tofu_s98, bean_s98, peas_s98, soy_s98);
/* fruit juices */
fj_s98=sum(aj_s98, oj_s98, ojca_s98, grj_s98, othj_s98  );
/* vegetables minus potatoes and legumes */
veg_noleg_s98=vege_s98-leg_s98;
/* dairy and fish */
dairy_s98=sum(skim_s98, m2_s98, whole_s98, sherb_s98, icecr_s98, flyog_s98, plyog_s98, cream_s98, cotch_s98, crmch_s98, otch_s98);
fish_s98=sum(ctuna_s98, dkfsh_s98, ofish_s98, shrim_s98);
/* poultry */
poul_s98=sum(chksa_s98, chwi_s98, chwo_s98, ctdog_s98 );
/* red and processed meat */
rpmeats_s98= sum(bacon_s98, bfdog_s98, pmsan_s98, procm_s98, hamb_s98, hambl_s98, bmix_s98, pmain_s98, bmain_s98);
/* ssb*/
ssb_s98  = sum( cola_s98, otsnc_s98, otsug_s98, punch_s98);
/*** SWEETS & DESERTS  ***/
swtsvg_s98  = sum(choco_s98, cdyw_s98, cdywo_s98, coknf_s98, cokh_s98, cokr_s98, brwn_s98, donut_s98, cakr_s98, cakh_s98, pieh_s98, srolf_s98, srolr_s98, srolh_s98, jam_s98);
swtsnochoc_s98  = sum(cdyw_s98, cdywo_s98, coknf_s98, cokh_s98, cokr_s98, brwn_s98, donut_s98, cakr_s98, cakh_s98, pieh_s98, srolf_s98, srolr_s98, srolh_s98, jam_s98);
run; 

%nur00(keep=db00 hbp00 chol00 cabg00 str00 alone00 marry00 
            thiaz00  lasix00  ccblo00  betab00  ace00 bprx00   
            couma00   digox00   antia00
            chrx00 oclrx00 /* Statin  other cholesterol-lowering */
            przc00 zol00 paxil00 celex00 antid00 /* antidepressant */
            tag00 h2blk00 prilo00 /* PPIs */
            insul00 ohypo00 /* Hypoglycemic */
            q00 sigm00 vbld00 hxclc00 fecal00 abdpn00 diarr00 rout00 bariu00 q15ap00 
            mamno00 mamsc00 mamsy00 mampt00
            );

%nur02(keep=db02 hbp02 chol02 cabg02 str02 
            thiaz02 lasix02 ace02 ccblo02 betab02 bprx02 
            przc02 zol02 paxil02 celex02 antid02 
            couma02   digox02   antia02 
            chrx02 oclrx02
            h2blk02 prilo02 /* PPIs */
            insul02 ohypo02
            cerbr02d
            q02 sigm02 vbld02 hxclc02 fecal02 abdpn02 diarr02 rout02 bariu02 q16apt02 
            mamno02 mamsc02 mamsy02 mampt02 );
%serv02(keep= 
/* whole grains */  cer_s02 ckcer_s02 oat_s02 dkbr_s02 ryebr_s02 brice_s02 oatbr_s02 bran_s02 wgerm_s02 popc_s02 ffpop_s02 tort_s02 
/* fruits */ rais_s02 prune_s02 ban_s02 cant_s02 appl_s02 apsau_s02 peach_s02 oran_s02 grfr_s02 straw_s02 blueb_s02  
/* vegetables */ tom_s02 toj_s02 tosau_s02 salsa_s02 brocc_s02 cauli_s02 cabb_s02 bruss_s02 rcar_s02 ccar_s02 mixv_s02 yam_s02 osqua_s02 eggpl_s02 kale_s02 cspin_s02 rspin_s02 ilett_s02 rlett_s02 cel_s02 grpep_s02 oniog_s02 oniov_s02 garl_s02 corn_s02
/* legumes */ sbean_s02 tofu_s02 bean_s02 peas_s02 soy_s02 

/* potatoes */  fries_s02 pot_s02 pchip_s02
/* dairy */ skim_s02 m1or2_s02 whole_s02 sherb_s02 icecr_s02 plyog_s02 flyog_s02 cream_s02 cotch_s02 crmch_s02 otch_s02
/* eggs */ eggwh_s02 egg_s02
/* fish */ ctuna_s02 bfsh_s02 dkfsh_s02 ofish_s02 shrim_s02
/* fruit juices */ aj_s02 oj_s02 ojca_s02 prunj_s02 othj_s02  


/* tea-coffee */ tea_s02 dtea_s02 dcaf_s02 coff_s02 
/* ssb */  cola_s02 otsug_s02 punch_s02   
/* sweets */ choco_s02 cdyw_s02 cdywo_s02 coknf_s02 cokr_s02 cokh_s02 brwn_s02 donut_s02 cakr_s02 cakh_s02 pieh_s02 srolf_s02 srolr_s02 srolh_s02 jam_s02
/* alcoholic bev */ beer_s02 lbeer_s02 rwine_s02 wwine_s02 liq_s02
/* meat */ chksa_s02 chwi_s02 chwo_s02 bacon_s02 bfdog_s02 ctdog_s02 pmsan_s02 procm_s02 hamb_s02 hambl_s02 bmix_s02 pmain_s02 bmain_s02 livb_s02 livc_s02
/* rpmeats_s02 ssb_s02 swtsvg_s02 swtsnochoc_s02 whgrn_s02 frui_s02 vege_s02 fv_s02 */
/*vege_s02 frui_s02 fv_s02 rpmeats_s02 whgrn_s02 cerbr_s02 ssb_s02 wcer_s02 swtsvg_s02 */

);

data serv02; 
  merge serv02 nur02 (keep=id cerbr02d);
  by id;

if cerbr02d in ( 0, 4, 11, 12, 13, 15, 16, 17, 19, 20, 21, 22, 23, 25, 29, 30, 31, 33, 35,37,38, 42, 43, 45, 46, 47,
59, 60,61,62,63,64,65,66,68,73,74,75,76,79,80,81,82,83, 84,85,86,88,95,96,97,98,
101,109,111, 112,113,114, 115,119, 123,124, 127, 128, 129,130,131, 132, 133, 134, 146, 147, 148,149,150,
151, 152,153, 154, 155, 157,158, 164,165, 166, 167,170,173, 176, 177,182, 184, 185, 186, 188/*,
190,197,198,997*/)
   then rcer_s02=cer_s02;
   else wcer_s02=cer_s02;
   if rcer_s02=. then rcer_s02=0;
   if wcer_s02=. then wcer_s02=0;


/* take sum of caffeinated and decaffeinated tea & coffee */
tea_s02  = SUM(tea_s02, dtea_s02);
coff_s02  = sum(dcaf_s02, coff_s02);
/** whole grains **/
whgrn_s02 =sum(cer_s02,  oat_s02,  ckcer_s02,  dkbr_s02,  brice_s02,  oatbr_s02,  bran_s02,  wgerm_s02);

/* fruit and vegetable*/
frui_s02 =sum(rais_s02, prune_s02, ban_s02, cant_s02, appl_s02, apsau_s02, peach_s02, oran_s02, grfr_s02, straw_s02, blueb_s02);
vege_s02 =sum(tom_s02,  toj_s02,   tosau_s02,  salsa_s02,   tofu_s02, sbean_s02,  bean_s02, soy_s02,   peas_s02,  
                brocc_s02,cauli_s02,    cabb_s02,  bruss_s02,   rcar_s02,  ccar_s02,  /*corn_s02,*/    mixv_s02,  
                 yam_s02,osqua_s02,   eggpl_s02,   kale_s02,  cspin_s02, rspin_s02, ilett_s02,   rlett_s02,  
                  cel_s02,grpep_s02,   oniog_s02,  oniov_s02);
fv_s02 =sum(frui_s02, vege_s02);
/* potatoes */
pota_s02=sum(fries_s02, pot_s02, pchip_s02);
/*legumes */
leg_s02=sum(sbean_s02, tofu_s02, bean_s02, peas_s02, soy_s02);
/* fruit juices */
fj_s02=sum(aj_s02, oj_s02, ojca_s02, prunj_s02, othj_s02  );
/* vegetables minus potatoes and legumes */
veg_noleg_s02=vege_s02-leg_s02;
/* dairy and fish */
dairy_s02=sum(skim_s02, m1or2_s02, whole_s02, sherb_s02, icecr_s02, plyog_s02, flyog_s02, cream_s02, cotch_s02, crmch_s02, otch_s02);
fish_s02=sum(ctuna_s02, bfsh_s02, dkfsh_s02, ofish_s02, shrim_s02);
egg_s02=sum(eggwh_s02, egg_s02);
/* red and processed meat */
rpmeats_s02= sum(bacon_s02, bfdog_s02, pmsan_s02, procm_s02, hamb_s02, hambl_s02, bmix_s02, pmain_s02, bmain_s02);
/* poultry */
poul_s02=sum(chksa_s02, chwi_s02, chwo_s02, ctdog_s02 );
/* ssb*/
ssb_s02  = sum( cola_s02, otsug_s02, punch_s02);
/*** SWEETS & DESERTS  ***/
swtsvg_s02  = sum(choco_s02, cdyw_s02, cdywo_s02, coknf_s02, cokh_s02, cokr_s02, brwn_s02, donut_s02, cakr_s02, cakh_s02, pieh_s02, srolf_s02, srolr_s02, srolh_s02, jam_s02);
swtsnochoc_s02  = sum( cdyw_s02, cdywo_s02, coknf_s02, cokh_s02, cokr_s02, brwn_s02, donut_s02, cakr_s02, cakh_s02, pieh_s02, srolf_s02, srolr_s02, srolh_s02, jam_s02);
run; 

%nur04(keep=db04 hbp04 chol04 cabg04 str04 alone04 marry04
            lasix04 thiaz04 ace04 ccblo04 betab04 k04 bprex04 
            mev04 zoc04 crest04 prav04 lip04 lesc04 oclrx04
            insul04 ohypo04 /* Hypoglycemic */
            ssri04 antid04  /* antidepressant */
            h2blk04 prilo04
            q04 virtc04 colsc04 sigsc04 sigm04 
            vbld04 diarr04 bariu04 cpolp04 fecal04 hxclc04 virt04 rout04 abdpn04
            q18apt04 q19pt04
            signo04 sig9304 sig9404 sig9604 sig9804 sig0004 sig0204 sigpt04
            colno04 col9304 col9404 col9604 col9804 col0004 col0204 colpt04
            mamno04 mamsc04 mamsy04 mampt04 
            );

%nur06(keep=db06 hbp06 chol06 str06 cabg06 
            antid06 ssri06 
            bprx06 lasix06 thiaz06 ace06 ccblo06 betab06 k06
            mev06 zoc06 crest06 prav06 lip06 ostat06 oclrx06
            insul06 ohypo06 /* Hypoglycemic */
            ssri06 antid06  /* antidepressant */
            prilo06 h2blk06 cerbr06d
            q06 virtc06 colsc06 sigsc06 q15apt06 
            vbld06 diarr06 bariu06 cpolp06 fecal06 hxclc06 virt06 rout06 abdpn06 q15apt06
            mamno06 mamsc06 mamsy06 mampt06 
            );
%serv06(keep= 
mchoc_s06 /*milk chocolate*/
dchoc_s06 /*dark chocolate*/
cdyw_s06 /*Candy bar*/
cdywo_s06/*Candy wo Chocolate*/
coknf_s06/*Cookies, Fat Free*/
cokr_s06/*Cookies, Commercial*/
cokh_s06/*Cookies, Home Baked*/
brwn_s06 /*Brownies*/
donut_s06/*Doughnuts*/
cake_s06/*Cake, Homemade*/
piehr_s06/*Pie, Ready-Made*/
jam_s06  /*Jams*/
pbut_s06 /*peanut butter*/
ffpop_s06/*Popcorn: Fat Free*/
popc_s06/*Popcorn: Regular*/
srolf_s06/*Sweetroll Fat Free*/
srolr_s06/*Sweetroll, Coffee Cake-Ready Made*/
srolh_s06/*Sweetroll, Coffee Cake-Home Made*/
pretz_s06/*Pretzels*/
pnut_s06/*Peanuts*/
wnut_s06/*Walnuts*/
onut_s06/*Other Nuts*/

pchip_s06 /*chips*/
icecr_s06 /*Ice cream*/
plyog_s06 /*DAIRY - Yogurt-plain, lo-carb*/
flyog_s06 /*Yogurt*/

cola_s06   otsug_s06  punch_s06
rais_s06 prune_s06 ban_s06 cant_s06   avo_s06   appl_s06   oran_s06  grfr_s06  straw_s06   blueb_s06   peach_s06 apric_s06 
tom_s06   toj_s06    tosau_s06   salsa_s06    tofu_s06  sbean_s06   bean_s06     peas_s06    soy_s06
brocc_s06 cauli_s06     cabb_s06   bruss_s06    rcar_s06   ccar_s06   corn_s06     mixv_s06   
yam_s06 osqua_s06    eggpl_s06    kale_s06   cspin_s06  rspin_s06  ilett_s06    rlett_s06   
cel_s06 grpep_s06    oniog_s06   oniov_s06  
hamb_s06 hambl_s06 bfdog_s06 procm_s06 bacon_s06 bmix_s06 bmain_s06 pmain_s06 pmsan_s06
whbr_s06 dkbr_s06  ryebr_s06 cer_s06 bran_s06 oatbr_s06 oat_s06 brice_s06 ckcer_s06  popc_s06 ffpop_s06
tea_s06 dtea_s06 coff_s06 dcaf_s06 rwine_s06 wwine_s06
/* potatoes */  fries_s06 pot_s06 pchip_s06
/* dairy */ skim_s06 m1or2_s06 whole_s06 sherb_s06 icecr_s06 plyog_s06 flyog_s06 cream_s06 cotch_s06 crmch_s06 otch_s06
/* eggs */ eggom_s06 egg_s06
/* fish */ ctuna_s06 bfsh_s06 dkfsh_s06 ofish_s06 shrim_s06
/* fruit juices */ aj_s06 oj_s06 ojca_s06 prunj_s06 othj_s06  
chksa_s06 chwi_s06 chwo_s06 ctdog_s06 
);

data serv06;
  merge serv06 nur06 (keep=id cerbr06d);
  by id;

if cerbr06d in ( 0, 4, 11, 12, 13, 15, 16, 17, 19, 20, 21, 22, 23, 25, 29, 30, 31, 33, 35,37,38, 42, 43, 45, 46, 47,
59, 60,61,62,63,64,65,66,68,73,74,75,76,79,80,81,82,83, 84,85,86,88,95,96,97,98,
101,109,111, 112,113,114, 115,119, 123,124, 127, 128, 129,130,131, 132, 133, 134, 146, 147, 148,149,150,
151, 152,153, 154, 155, 157,158, 164,165, 166, 167,170,173, 176, 177,182, 184, 185, 186, 188/*,
190,197,198,997*/)
then rcer_s06=cer_s06;
else wcer_s06=cer_s06;
if rcer_s06=. then rcer_s06=0;
if wcer_s06=. then wcer_s06=0;

choco_s06=sum(mchoc_s06 , dchoc_s06);
/* take sum of caffeinated and decaffeinated tea & coffee */
tea06d  = SUM(tea_s06, dtea_s06);
coff06d  = sum(coff_s06, dcaf_s06);
/** whole grains **/
whgrn_s06  = sum(wcer_s06  ,  oat_s06  ,  ckcer_s06,  dkbr_s06 ,  brice_s06,  oatbr_s06,  bran_s06, popc_s06  );
/* fruit and vegetables */ 
vege_s06  = sum(  tom_s06,  toj_s06,   tosau_s06,  salsa_s06,   tofu_s06, sbean_s06,  bean_s06,    peas_s06,  soy_s06,
                brocc_s06,cauli_s06,    cabb_s06,  bruss_s06,   rcar_s06,  ccar_s06,  corn_s06,    mixv_s06,  
                 yam_s06,osqua_s06,   eggpl_s06,   kale_s06,  cspin_s06, rspin_s06, ilett_s06,   rlett_s06,  
                  cel_s06,grpep_s06,   oniog_s06,  oniov_s06       );   
frui_s06  = sum(rais_s06 ,  prune_s06, ban_s06, cant_s06 , avo_s06,  appl_s06 ,  
                oran_s06 ,grfr_s06 ,straw_s06,  blueb_s06,  peach_s06,apric_s06  );
fv_s06=sum (frui_s06, vege_s06);

/* potatoes */
pota_s06=sum(fries_s06, pot_s06, pchip_s06);
/*legumes */
leg_s06=sum(sbean_s06, tofu_s06, bean_s06, peas_s06, soy_s06);
/* fruit juices */
fj_s06=sum(aj_s06, oj_s06, ojca_s06, prunj_s06, othj_s06  );
/* vegetables minus potatoes and legumes */
veg_noleg_s06=vege_s06-leg_s06;
/* dairy and fish */
dairy_s06=sum(skim_s06, m1or2_s06, whole_s06, sherb_s06, icecr_s06, plyog_s06, flyog_s06, cream_s06, cotch_s06, crmch_s06, otch_s06);
fish_s06=sum(ctuna_s06, bfsh_s06, dkfsh_s06, ofish_s06, shrim_s06);
egg_s06=sum(eggom_s06, egg_s06);
/* poultry */
poul_s06=sum(chksa_s06, chwi_s06, chwo_s06, ctdog_s06 );
/* red and processed meat */
rpmeats_s06= sum(hamb_s06, hambl_s06, bfdog_s06, procm_s06, bacon_s06,  bmix_s06, bmain_s06, pmain_s06, pmsan_s06);
/* ssb*/
ssb_s06  = sum( cola_s06, otsug_s06, punch_s06);
/*** SWEETS & DESERTS  ***/
swtsvg_s06  = sum(mchoc_s06, dchoc_s06 ,  cdyw_s06 ,  cdywo_s06,  coknf_s06,  cokr_s06 ,  cokh_s06 ,  brwn_s06 ,  donut_s06,  
                 cake_s06 ,  piehr_s06 ,  srolf_s06,  srolr_s06,  srolh_s06 ,icecr_s06 );
swtsnochoc_s06  = sum(cdyw_s06 ,  cdywo_s06,  coknf_s06,  cokr_s06 ,  cokh_s06 ,  brwn_s06 ,  donut_s06,  
                 cake_s06 ,  piehr_s06 ,  srolf_s06,  srolr_s06,  srolh_s06 ,icecr_s06 );
run; 

%nur08(keep=db08 hbp08 chol08 cabg08 str08 alone08 marry08 
            antid08 ssri08 
            bprx08 lasix08 thiaz08 ace08 ccblo08 betab08 k08 
            mev08 zoc08 crest08 prav08 lip08 ostat08 oclrx08
            insul08 metf08 ohypo08
            ssri08 antid08  /* antidepressant */
            prilo08 h2blk08
            q08 virtc08 colsc08 sigsc08  q16apt08
            vbld08 diarr08 bariu08 cpolp08 fecal08 hxclc08 virt08 rout08 abdpn08 q16apt08
            mamno08 mamsc08 mamsy08 mampt08 
            );

%nur10(keep=db10 hbp10 chol10 str10 cabg10 
            antid10 ssri10 
            bprx10 lasix10 thiaz10 ace10 ccblo10 betab10 angio10 k10
            stat10 mev10 prav10 zoc10 lip10 crest10 ostat10 oclrx10
            insul10 metf10 avand10 ohypo10 
            prilo10 h2blk10
            cerbr10d
            q10 virtc10 colsc10 sigsc10 q16apt10 
            vbld10 diarr10 bariu10 cpolp10 fecal10 hxclc10 virt10 rout10 abdpn10 q16apt10
            mamsc10 mamsy10 mamno10 mampt10 
            );

%serv10(keep= 
mchoc_s10 /*milk chocolate*/
dchoc_s10 /*dark chocolate*/
cdyw_s10 /*Candy bar*/
cdywo_s10/*Candy wo Chocolate*/
coknf_s10/*Cookies, Fat Free*/
cokr_s10/*Cookies, Commercial*/
cokh_s10/*Cookies, Home Baked*/
donut_s10/*Doughnuts*/
cake_s10/*Cake, Homemade*/
piehr_s10/*Pie, Ready-Made*/
jam_s10  /*Jams*/
pbut_s10 /*peanut butter*/
ffpop_s10/*Popcorn: Fat Free*/
popc_s10/*Popcorn: Regular*/
sroll_s10 /*Sweetroll*/
pretz_s10/*Pretzels*/
pnut_s10/*Peanuts*/
wnut_s10/*Walnuts*/
onut_s10/*Other Nuts*/

pchip_s10 /*chips*/
icecr_s10 /*Ice cream*/
plyog_s10 /*DAIRY - Yogurt-plain, lo-carb*/
flyog_s10 /*Yogurt*/

cola_s10   otsug_s10  punch_s10

rais_s10  prune_s10  prunj_s10  ban_s10  cant_s10  appl_s10  aj_s10 oran_s10  ojca_s10  oj_s10  grfr_s10  othj_s10  straw_s10  blueb_s10  peach_s10
tom_s10 toj_s10  tosau_s10  salsa_s10  tofu_s10  sbean_s10  bean_s10  peas_s10 soy_s10 brocc_s10  cauli_s10 cabb_s10  bruss_s10  rcar_s10  ccar_s10  corn_s10  mixv_s10  yam_s10  osqua_s10  eggpl_s10 
kale_s10  cspin_s10  rspin_s10  ilett_s10  rlett_s10  cel_s10  grpep_s10  oniov_s10  
hamb_s10 hambl_s10 bfdog_s10 procm_s10 bacon_s10 bmix_s10 bmain_s10 pmain_s10 pmsan_s10
whbr_s10 dkbr_s10 wgerm_s10 cer_s10 ryebr_s10  oatbr_s10 oat_s10 brice_s10 ckcer_s10  popc_s10 ffpop_s10
tea_s10 dtea_s10 coff_s10 dcaf_s10 rwine_s10 wwine_s10
/* potatoes */  fries_s10 pot_s10 pchip_s10
/* dairy */ skim_s10 m1or2_s10 whole_s10 sherb_s10 icecr_s10 plyog_s10 flyog_s10 cream_s10 cotch_s10 crmch_s10 otch_s10
/* eggs */ eggom_s10 egg_s10
/* fish */ ctuna_s10 bfsh_s10 dkfsh_s10 ofish_s10 shrim_s10
/* fruit juices */ aj_s10 oj_s10 ojca_s10 prunj_s10 othj_s10  
chksa_s10 chwi_s10 chwo_s10 ctdog_s10 
 );


data serv10;
 merge serv10 nur10 (keep=id cerbr10d);
 by id;

if cerbr10d in ( 0, 4, 11, 12, 13, 15, 16, 17, 19, 20, 21, 22, 23, 25, 29, 30, 31, 33, 35,37,38, 42, 43, 45, 46, 47,
59, 60,61,62,63,64,65,66,68,73,74,75,76,79,80,81,82,83, 84,85,86,88,95,96,97,98,
101,109,111, 112,113,114, 115,119, 123,124, 127, 128, 129,130,131, 132, 133, 134, 146, 147, 148,149,150,
151, 152,153, 154, 155, 157,158, 164,165, 166, 167,170,173, 176, 177,182, 184, 185, 186, 188/*,
190,197,198,997*/)
    then rcer_s10=cer_s10;
    else wcer_s10=cer_s10;
    if rcer_s10=. then rcer_s10=0;
    if wcer_s10=. then wcer_s10=0;

choco_s10=sum(mchoc_s10, dchoc_s10);

/* take sum of caffeinated and decaffeinated tea & coffee */
tea10d  = SUM(tea_s10, dtea_s10);
coff10d  = sum(coff_s10, dcaf_s10);
/** whole grains **/
whgrn_s10  = sum(wcer_s10  ,  oat_s10  ,  ckcer_s10,  dkbr_s10 ,  brice_s10,  oatbr_s10, popc_s10, wgerm_s10 );
/* fruit and vegetables */ 
frui_s10=sum (rais_s10, prune_s10, prunj_s10, ban_s10, cant_s10 , appl_s10, /* aj_s10, oran_s10, ojca_s10, oj_s10, othj_s10, */ grfr_s10,  straw_s10, blueb_s10, peach_s10 );
vege_s10=sum (tom_s10, toj_s10, tosau_s10, salsa_s10, tofu_s10, sbean_s10, bean_s10, peas_s10, brocc_s10, cauli_s10, cabb_s10, bruss_s10, rcar_s10, ccar_s10, corn_s10, mixv_s10, yam_s10, osqua_s10, eggpl_s10, 
             kale_s10, cspin_s10, rspin_s10, ilett_s10, rlett_s10, cel_s10, grpep_s10, oniov_s10);
fv_s10=sum (frui_s10, vege_s10);

/* potatoes */
pota_s10=sum(fries_s10, pot_s10, pchip_s10);
/*legumes */
leg_s10=sum(sbean_s10, tofu_s10, bean_s10, peas_s10, soy_s10);
/* fruit juices */
fj_s10=sum(aj_s10, oj_s10, ojca_s10, prunj_s10, othj_s10  );
/* vegetables minus potatoes and legumes */
veg_noleg_s10=vege_s10-leg_s10;
/* dairy and fish */
dairy_s10=sum(skim_s10, m1or2_s10, whole_s10, sherb_s10, icecr_s10, plyog_s10, flyog_s10, cream_s10, cotch_s10, crmch_s10, otch_s10);
fish_s10=sum(ctuna_s10, bfsh_s10, dkfsh_s10, ofish_s10, shrim_s10);
egg_s10=sum(eggom_s10, egg_s10);
/* poultry */
poul_s10=sum(chksa_s10, chwi_s10, chwo_s10, ctdog_s10 );
/* red and processed meat */
rpmeats_s10= sum(hamb_s10, hambl_s10, bfdog_s10, procm_s10, bacon_s10,  bmix_s10, bmain_s10, pmain_s10, pmsan_s10);
/* ssb*/
ssb_s10  = sum( cola_s10, otsug_s10, punch_s10);
/*** SWEETS & DESERTS  ***/
swtsvg_s10  = sum(mchoc_s10, dchoc_s10 ,  cdyw_s10 ,  cdywo_s10,  coknf_s10,  cokr_s10 ,  cokh_s10 ,   donut_s10,  
                 cake_s10 ,  piehr_s10 ,  sroll_s10 ,icecr_s10 );
swtsnochoc_s10  = sum(cdyw_s10 ,  cdywo_s10,  coknf_s10,  cokr_s10 ,  cokh_s10 ,   donut_s10,  
                 cake_s10 ,  piehr_s10 ,  sroll_s10 ,icecr_s10 );
run; 

/*** NSAIDs ***/
%meds8020(keep=aspu84 aspu86 aspu88 aspu90 aspu92 aspu94 aspu96 aspu98 aspu00 aspu02 aspu04 aspu06 aspu08 aspu10
               asp84 asp86 asp88 asp90 asp92 asp94 asp96 asp98 asp00 asp02 asp04 asp06 asp08 asp10);
            asp84=0; if aspu84=1 then asp84=1; /* $label 0.nonuser;\1.current user;\9.unknown status */
            asp86=0; if aspu86=1 then asp86=1;
            asp88=0; if aspu88=1 then asp88=1;
            asp90=0; if aspu90=1 then asp90=1;
            asp92=0; if aspu92=1 then asp92=1;
            asp94=0; if aspu94=1 then asp94=1;
            asp96=0; if aspu96=1 then asp96=1;
            asp98=0; if aspu98=1 then asp98=1;
            asp00=0; if aspu00=1 then asp00=1;
            asp02=0; if aspu02=1 then asp02=1;
            asp04=0; if aspu04=1 then asp04=1;
            asp06=0; if aspu06=1 then asp06=1;
            asp08=0; if aspu08=1 then asp08=1;
            asp10=0; if aspu10=1 then asp10=1;

%supp8016(keep= mvitu84 mvitu86 mvitu88 mvitu90 mvitu92 mvitu94 mvitu96 mvitu98 mvitu00 mvitu02 mvitu04 mvitu06 mvitu08 mvitu10
                mvt84  mvt86  mvt88  mvt90  mvt92  mvt94  mvt96  mvt98  mvt00  mvt02  mvt04  mvt06  mvt08  mvt10);
            mvt84=0; if mvitu84=1 then mvt84=1; /* $label   0.nonuser;\ 1.current user;\ 9.unknown status */
            mvt86=0; if mvitu86=1 then mvt86=1;
            mvt88=0; if mvitu88=1 then mvt88=1;
            mvt90=0; if mvitu90=1 then mvt90=1;
            mvt92=0; if mvitu92=1 then mvt92=1;
            mvt94=0; if mvitu94=1 then mvt94=1;
            mvt96=0; if mvitu96=1 then mvt96=1;
            mvt98=0; if mvitu98=1 then mvt98=1;
            mvt00=0; if mvitu00=1 then mvt00=1;
            mvt02=0; if mvitu02=1 then mvt02=1;
            mvt04=0; if mvitu04=1 then mvt04=1;
            mvt06=0; if mvitu06=1 then mvt06=1;
            mvt08=0; if mvitu08=1 then mvt08=1;
            mvt10=0; if mvitu10=1 then mvt10=1; 


/* census tract-level neighborhood socioeconomic (SES) status */
%nses7616(keep=
id 
nSES_86 nSES_88 nSES_90 nSES_92 nSES_94 nSES_96 nSES_98 nSES_00 nSES_02 nSES_04 nSES_06 nSES_08 nSES_10 nSES_12 nSES_14 nSES_16
);



/*** Nutrient ***/
%n86_nts(keep= alco86n calor86n na86n
  carbo86n 
  pot_carb86n/*Carbohydrate from potato*/
  leg_carb86n/*Carbohdyrate from legumes*/
  frtcb_njs86n/*Fruit carbohydrate minus juice sugar*/
  vegcb_npl86n/*Vegetable Carbohydrate minus potato and legume carbohydrate*/
  sumref_carb86n/*Sum of carbohydrate from potato, added sugar in foods, and refined grain, foods only*/
  sum_whgcarb86n/*Sum Carbohydrate From Wholegrain*/
  sumoth_crb_njs86n/*Sum of carbohydrate from fruit minus juice sugar, plus carb from vegetable, legume, and lact*/
  mill_grn86n/*Carbohydrate from milled whole grain*/
  int_grn86n/*Carbohydrate from intact whole grain*/
  ref_grn86n/*Carbohydrate from refined grain*/
  adds_cer86n/********Added sugar minus cereal fruit for Carb Quality Variables only********/
  adds_fds86n/********Added Sugar From Foods for Carb Quality Variables only********/

  tfat86n
  afat86n/*Animal Fat*/
  vfat86n/*Vegetable Fat*/
  dfat86n/*Dairy Fat*/
  sat86n/*Total Saturated Fat*/
  mon86n/*Total Monounsaturated Fat*/
  mufa_pl86n/*Plant Monounsaturated Fat*/
  mufa_a86n/*Animal Monounsaturated Fat*/
  poly86n/*Total Polyunsaturated Fat*/
  omega86n/*Omega 3 (f205+f226) gm, no alpha 18:3*/
  trnss86n/*Total Trans Fat*/
 
  prot86n
  aprot86n/*Animal Protein*/
  vprot86n/*Vegetable Protein*/
  dprot86n/*Dairy Protein*/
  ); 

data n86_nts;
set n86_nts;

array gram_cp carbo86n pot_carb86n leg_carb86n frtcb_njs86n vegcb_npl86n sumref_carb86n sum_whgcarb86n sumoth_crb_njs86n mill_grn86n int_grn86n ref_grn86n prot86n aprot86n vprot86n dprot86n adds_cer86n adds_fds86n;
array perc_cp carbo86p pot_carb86p leg_carb86p frtcb_njs86p vegcb_npl86p sumref_carb86p sum_whgcarb86p sumoth_crb_njs86p mill_grn86p int_grn86p ref_grn86p prot86p aprot86p vprot86p dprot86p adds_cer86p adds_fds86p;

array gram_f  tfat86n afat86n vfat86n dfat86n sat86n mon86n mufa_pl86n mufa_a86n poly86n omega86n trnss86n;
array perc_f  tfat86p afat86p vfat86p dfat86p sat86p mon86p mufa_pl86p mufa_a86p poly86p omega86p trnss86p;

do over gram_cp;
perc_cp=(4*gram_cp)/calor86n*100;
end;

do over gram_f;
perc_f=(9*gram_f)/calor86n*100;
end;

hfat86p=sum(mufa_pl86p,poly86p);     /*healthy fat as % of energy*/
ufat86p=sum(sat86p,mufa_a86p,trnss86p);     /*unhealthy fat as % of energy*/

hcarbo86p=sum(leg_carb86p,frtcb_njs86p,vegcb_npl86p,sum_whgcarb86p); /*healthy carb as % of energy*/
run;

proc rank data=n86_nts groups=11 out=n86_nts;  * assign 11 groups based on these values (high prot=10, low prot = 0);
var    prot86p  aprot86p  vprot86p  tfat86p  afat86p  vfat86p  sum_whgcarb86p  sumoth_crb_njs86p  hfat86p  ufat86p;
ranks  prot86pd aprot86pd vprot86pd tfat86pd afat86pd vfat86pd sum_whgcarb86pd sumoth_crb_njs86pd hfat86pd ufat86pd;
run;

proc rank data=n86_nts descending groups=11 out=n86_nts;  * high carb = 0, low carb = 10;
var    carbo86p  sumref_carb86p  adds_cer86p  adds_fds86p  hcarbo86p;
ranks  carbo86pd sumref_carb86pd adds_cer86pd adds_fds86pd hcarbo86pd;
run;

data n86_nts;
set n86_nts;
/* low carb score: more prot/fat (higher group ranking), less carb (higher group ranking)*/
lcs86=sum(carbo86pd,prot86pd,tfat86pd);     * low carb score = sum of total carb, prot, fat;
alcs86=sum(carbo86pd,aprot86pd,afat86pd);   * animal-based low carb score = sum or carb, animal protein, animal fat;
vlcs86=sum(carbo86pd,vprot86pd,vfat86pd);   * vege(plant-based) low carb score = sum of carb, vege protein, vege fat;

hlcs86=sum(sumref_carb86pd,vprot86pd,vfat86pd); * healthy(high-quality) low carb score = sum of a reversed rank of unhealthful carb, vege protein, vege fat;
ulcs86=sum(hcarbo86pd     ,aprot86pd,afat86pd); * unhealthy(low-quality) low carb score = sum of a reversed rank of healthy carb, animal protein, animal fat;

cqs86=sum(sum_whgcarb86pd,sumoth_crb_njs86pd,sumref_carb86pd);
sug_cqs86=sum(adds_cer86pd,adds_fds86pd);
run;

proc rank data=n86_nts groups=11 out=n86_nts;
var    carbo86p  sumref_carb86p  hcarbo86p   prot86p  aprot86p  vprot86p   ;
ranks  carbo86pdd sumref_carb86pdd hcarbo86pdd  prot86pdd aprot86pdd vprot86pdd  ;
run;

proc rank data=n86_nts descending groups=11 out=n86_nts ;
var    tfat86p  afat86p  vfat86p  ;
ranks  tfat86pdd afat86pdd vfat86pdd ;
run;

data n86_nts;
set n86_nts;

lfs86=sum(carbo86pdd,prot86pdd,tfat86pdd);
alfs86=sum(carbo86pdd,aprot86pdd,vfat86pdd);
vlfs86=sum(carbo86pdd,vprot86pdd,afat86pdd);

hlfs86=sum(hcarbo86pdd, vprot86pdd, afat86pdd);
ulfs86=sum(sumref_carb86pdd,aprot86pdd,vfat86pdd);
run;


%n90_nts(keep=alco90n calor90n na90n

  carbo90n 
  pot_carb90n/*Carbohydrate from potato*/
  leg_carb90n/*Carbohdyrate from legumes*/
  frtcb_njs90n/*Fruit carbohydrate minus juice sugar*/
  vegcb_npl90n/*Vegetable Carbohydrate minus potato and legume carbohydrate*/
  sumref_carb90n/*Sum of carbohydrate from potato, added sugar in foods, and refined grain, foods only*/
  sum_whgcarb90n/*Sum Carbohydrate From Wholegrain*/
  sumoth_crb_njs90n/*Sum of carbohydrate from fruit minus juice sugar, plus carb from vegetable, legume, and lact*/
  mill_grn90n/*Carbohydrate from milled whole grain*/
  int_grn90n/*Carbohydrate from intact whole grain*/
  ref_grn90n/*Carbohydrate from refined grain*/
  adds_cer90n/********Added sugar minus cereal fruit for Carb Quality Variables only********/
  adds_fds90n/********Added Sugar From Foods for Carb Quality Variables only********/

  tfat90n
  afat90n/*Animal Fat*/
  vfat90n/*Vegetable Fat*/
  dfat90n/*Dairy Fat*/
  sat90n/*Total Saturated Fat*/
  mon90n/*Total Monounsaturated Fat*/
  mufa_pl90n/*Plant Monounsaturated Fat*/
  mufa_a90n/*Animal Monounsaturated Fat*/
  poly90n/*Total Polyunsaturated Fat*/
  omega90n/*Omega 3 (f205+f226) gm, no alpha 18:3*/
  trnss90n/*Total Trans Fat*/
 
  prot90n
  aprot90n/*Animal Protein*/
  vprot90n/*Vegetable Protein*/
  dprot90n/*Dairy Protein*/
  ); 

data n90_nts;
set n90_nts;

array gram_cp carbo90n pot_carb90n leg_carb90n frtcb_njs90n vegcb_npl90n sumref_carb90n sum_whgcarb90n sumoth_crb_njs90n mill_grn90n int_grn90n ref_grn90n prot90n aprot90n vprot90n dprot90n adds_cer90n adds_fds90n;
array perc_cp carbo90p pot_carb90p leg_carb90p frtcb_njs90p vegcb_npl90p sumref_carb90p sum_whgcarb90p sumoth_crb_njs90p mill_grn90p int_grn90p ref_grn90p prot90p aprot90p vprot90p dprot90p adds_cer90p adds_fds90p;

array gram_f  tfat90n afat90n vfat90n dfat90n sat90n mon90n mufa_pl90n mufa_a90n poly90n omega90n trnss90n;
array perc_f  tfat90p afat90p vfat90p dfat90p sat90p mon90p mufa_pl90p mufa_a90p poly90p omega90p trnss90p;

do over gram_cp;
perc_cp=(4*gram_cp)/calor90n*100;
end;

do over gram_f;
perc_f=(9*gram_f)/calor90n*100;
end;

hfat90p=sum(mufa_pl90p,poly90p);
ufat90p=sum(sat90p,mufa_a90p,trnss90p);

hcarbo90p=sum(leg_carb90p,frtcb_njs90p,vegcb_npl90p,sum_whgcarb90p);

run;

proc rank data=n90_nts groups=11 out=n90_nts;
var    prot90p  aprot90p  vprot90p  tfat90p  afat90p  vfat90p  sum_whgcarb90p  sumoth_crb_njs90p  hfat90p  ufat90p;
ranks  prot90pd aprot90pd vprot90pd tfat90pd afat90pd vfat90pd sum_whgcarb90pd sumoth_crb_njs90pd hfat90pd ufat90pd;
run;

proc rank data=n90_nts descending groups=11 out=n90_nts;
var    carbo90p  sumref_carb90p  adds_cer90p  adds_fds90p  hcarbo90p;
ranks  carbo90pd sumref_carb90pd adds_cer90pd adds_fds90pd hcarbo90pd;
run;

data n90_nts;
set n90_nts;

lcs90=sum(carbo90pd,prot90pd,tfat90pd);
alcs90=sum(carbo90pd,aprot90pd,afat90pd);
vlcs90=sum(carbo90pd,vprot90pd,vfat90pd);

hlcs90=sum(sumref_carb90pd,vprot90pd,vfat90pd);
ulcs90=sum(hcarbo90pd     ,aprot90pd,afat90pd);

cqs90=sum(sum_whgcarb90pd,sumoth_crb_njs90pd,sumref_carb90pd);
sug_cqs90=sum(adds_cer90pd,adds_fds90pd);
run;

proc rank data=n90_nts groups=11 out=n90_nts;
var    carbo90p  sumref_carb90p  hcarbo90p   prot90p  aprot90p  vprot90p   ;
ranks  carbo90pdd sumref_carb90pdd hcarbo90pdd  prot90pdd aprot90pdd vprot90pdd  ;
run;

proc rank data=n90_nts descending groups=11 out=n90_nts ;
var    tfat90p  afat90p  vfat90p  ;
ranks  tfat90pdd afat90pdd vfat90pdd ;
run;

data n90_nts;
set n90_nts;

lfs90=sum(carbo90pdd,prot90pdd,tfat90pdd);
alfs90=sum(carbo90pdd,aprot90pdd,vfat90pdd);
vlfs90=sum(carbo90pdd,vprot90pdd,afat90pdd);

hlfs90=sum(hcarbo90pdd, vprot90pdd, afat90pdd);
ulfs90=sum(sumref_carb90pdd,aprot90pdd,vfat90pdd);
run;

%n94_nts(keep=alco94n calor94n na94n

  carbo94n 
  pot_carb94n/*Carbohydrate from potato*/
  leg_carb94n/*Carbohdyrate from legumes*/
  frtcb_njs94n/*Fruit carbohydrate minus juice sugar*/
  vegcb_npl94n/*Vegetable Carbohydrate minus potato and legume carbohydrate*/
  sumref_carb94n/*Sum of carbohydrate from potato, added sugar in foods, and refined grain, foods only*/
  sum_whgcarb94n/*Sum Carbohydrate From Wholegrain*/
  sumoth_crb_njs94n/*Sum of carbohydrate from fruit minus juice sugar, plus carb from vegetable, legume, and lact*/
  mill_grn94n/*Carbohydrate from milled whole grain*/
  int_grn94n/*Carbohydrate from intact whole grain*/
  ref_grn94n/*Carbohydrate from refined grain*/
  adds_cer94n/********Added sugar minus cereal fruit for Carb Quality Variables only********/
  adds_fds94n/********Added Sugar From Foods for Carb Quality Variables only********/

  tfat94n
  afat94n/*Animal Fat*/
  vfat94n/*Vegetable Fat*/
  dfat94n/*Dairy Fat*/
  sat94n/*Total Saturated Fat*/
  mon94n/*Total Monounsaturated Fat*/
  mufa_pl94n/*Plant Monounsaturated Fat*/
  mufa_a94n/*Animal Monounsaturated Fat*/
  poly94n/*Total Polyunsaturated Fat*/
  omega94n/*Omega 3 (f205+f226) gm, no alpha 18:3*/
  trnss94n/*Total Trans Fat*/
 
  prot94n
  aprot94n/*Animal Protein*/
  vprot94n/*Vegetable Protein*/
  dprot94n/*Dairy Protein*/
  ); 

data n94_nts;
set n94_nts;

array gram_cp carbo94n pot_carb94n leg_carb94n frtcb_njs94n vegcb_npl94n sumref_carb94n sum_whgcarb94n sumoth_crb_njs94n mill_grn94n int_grn94n ref_grn94n prot94n aprot94n vprot94n dprot94n adds_cer94n adds_fds94n;
array perc_cp carbo94p pot_carb94p leg_carb94p frtcb_njs94p vegcb_npl94p sumref_carb94p sum_whgcarb94p sumoth_crb_njs94p mill_grn94p int_grn94p ref_grn94p prot94p aprot94p vprot94p dprot94p adds_cer94p adds_fds94p;

array gram_f  tfat94n afat94n vfat94n dfat94n sat94n mon94n mufa_pl94n mufa_a94n poly94n omega94n trnss94n;
array perc_f  tfat94p afat94p vfat94p dfat94p sat94p mon94p mufa_pl94p mufa_a94p poly94p omega94p trnss94p;

do over gram_cp;
perc_cp=(4*gram_cp)/calor94n*100;
end;

do over gram_f;
perc_f=(9*gram_f)/calor94n*100;
end;

hfat94p=sum(mufa_pl94p,poly94p);
ufat94p=sum(sat94p,mufa_a94p,trnss94p);

hcarbo94p=sum(leg_carb94p,frtcb_njs94p,vegcb_npl94p,sum_whgcarb94p);

run;

proc rank data=n94_nts groups=11 out=n94_nts;
var    prot94p  aprot94p  vprot94p  tfat94p  afat94p  vfat94p  sum_whgcarb94p  sumoth_crb_njs94p  hfat94p  ufat94p;
ranks  prot94pd aprot94pd vprot94pd tfat94pd afat94pd vfat94pd sum_whgcarb94pd sumoth_crb_njs94pd hfat94pd ufat94pd;
run;

proc rank data=n94_nts descending groups=11 out=n94_nts;
var    carbo94p  sumref_carb94p  adds_cer94p  adds_fds94p  hcarbo94p;
ranks  carbo94pd sumref_carb94pd adds_cer94pd adds_fds94pd hcarbo94pd;
run;

data n94_nts;
set n94_nts;

lcs94=sum(carbo94pd,prot94pd,tfat94pd);
alcs94=sum(carbo94pd,aprot94pd,afat94pd);
vlcs94=sum(carbo94pd,vprot94pd,vfat94pd);

hlcs94=sum(sumref_carb94pd,vprot94pd,vfat94pd);
ulcs94=sum(hcarbo94pd     ,aprot94pd,afat94pd);

cqs94=sum(sum_whgcarb94pd,sumoth_crb_njs94pd,sumref_carb94pd);
sug_cqs94=sum(adds_cer94pd,adds_fds94pd);
run;

proc rank data=n94_nts groups=11 out=n94_nts;
var    carbo94p  sumref_carb94p  hcarbo94p   prot94p  aprot94p  vprot94p   ;
ranks  carbo94pdd sumref_carb94pdd hcarbo94pdd  prot94pdd aprot94pdd vprot94pdd  ;
run;

proc rank data=n94_nts descending groups=11 out=n94_nts ;
var    tfat94p  afat94p  vfat94p  ;
ranks  tfat94pdd afat94pdd vfat94pdd ;
run;

data n94_nts;
set n94_nts;

lfs94=sum(carbo94pdd,prot94pdd,tfat94pdd);
alfs94=sum(carbo94pdd,aprot94pdd,vfat94pdd);
vlfs94=sum(carbo94pdd,vprot94pdd,afat94pdd);
hlfs94=sum(hcarbo94pdd, vprot94pdd, afat94pdd);
ulfs94=sum(sumref_carb94pdd,aprot94pdd,vfat94pdd);
run;

%n98_nts(keep=alco98n calor98n na98n

  carbo98n 
  pot_carb98n/*Carbohydrate from potato*/
  leg_carb98n/*Carbohdyrate from legumes*/
  frtcb_njs98n/*Fruit carbohydrate minus juice sugar*/
  vegcb_npl98n/*Vegetable Carbohydrate minus potato and legume carbohydrate*/
  sumref_carb98n/*Sum of carbohydrate from potato, added sugar in foods, and refined grain, foods only*/
  sum_whgcarb98n/*Sum Carbohydrate From Wholegrain*/
  sumoth_crb_njs98n/*Sum of carbohydrate from fruit minus juice sugar, plus carb from vegetable, legume, and lact*/
  mill_grn98n/*Carbohydrate from milled whole grain*/
  int_grn98n/*Carbohydrate from intact whole grain*/
  ref_grn98n/*Carbohydrate from refined grain*/
  adds_cer98n/********Added sugar minus cereal fruit for Carb Quality Variables only********/
  adds_fds98n/********Added Sugar From Foods for Carb Quality Variables only********/

  tfat98n
  afat98n/*Animal Fat*/
  vfat98n/*Vegetable Fat*/
  dfat98n/*Dairy Fat*/
  sat98n/*Total Saturated Fat*/
  mon98n/*Total Monounsaturated Fat*/
  mufa_pl98n/*Plant Monounsaturated Fat*/
  mufa_a98n/*Animal Monounsaturated Fat*/
  poly98n/*Total Polyunsaturated Fat*/
  omega98n/*Omega 3 (f205+f226) gm, no alpha 18:3*/
  trn0098n/*Total Trans Fat*/
 
  prot98n
  aprot98n/*Animal Protein*/
  vprot98n/*Vegetable Protein*/
  dprot98n/*Dairy Protein*/
  ); 

data n98_nts;
set n98_nts;

array gram_cp carbo98n pot_carb98n leg_carb98n frtcb_njs98n vegcb_npl98n sumref_carb98n sum_whgcarb98n sumoth_crb_njs98n mill_grn98n int_grn98n ref_grn98n prot98n aprot98n vprot98n dprot98n adds_cer98n adds_fds98n;
array perc_cp carbo98p pot_carb98p leg_carb98p frtcb_njs98p vegcb_npl98p sumref_carb98p sum_whgcarb98p sumoth_crb_njs98p mill_grn98p int_grn98p ref_grn98p prot98p aprot98p vprot98p dprot98p adds_cer98p adds_fds98p;

array gram_f  tfat98n afat98n vfat98n dfat98n sat98n mon98n mufa_pl98n mufa_a98n poly98n omega98n trn0098n;
array perc_f  tfat98p afat98p vfat98p dfat98p sat98p mon98p mufa_pl98p mufa_a98p poly98p omega98p trnss98p;

do over gram_cp;
perc_cp=(4*gram_cp)/calor98n*100;
end;

do over gram_f;
perc_f=(9*gram_f)/calor98n*100;
end;

hfat98p=sum(mufa_pl98p,poly98p);
ufat98p=sum(sat98p,mufa_a98p,trnss98p);

hcarbo98p=sum(leg_carb98p,frtcb_njs98p,vegcb_npl98p,sum_whgcarb98p);

run;

proc rank data=n98_nts groups=11 out=n98_nts;
var    prot98p  aprot98p  vprot98p  tfat98p  afat98p  vfat98p  sum_whgcarb98p  sumoth_crb_njs98p  hfat98p  ufat98p;
ranks  prot98pd aprot98pd vprot98pd tfat98pd afat98pd vfat98pd sum_whgcarb98pd sumoth_crb_njs98pd hfat98pd ufat98pd;
run;

proc rank data=n98_nts descending groups=11 out=n98_nts;
var    carbo98p  sumref_carb98p  adds_cer98p  adds_fds98p  hcarbo98p;
ranks  carbo98pd sumref_carb98pd adds_cer98pd adds_fds98pd hcarbo98pd;
run;

data n98_nts;
set n98_nts;

lcs98=sum(carbo98pd,prot98pd,tfat98pd);
alcs98=sum(carbo98pd,aprot98pd,afat98pd);
vlcs98=sum(carbo98pd,vprot98pd,vfat98pd);

hlcs98=sum(sumref_carb98pd,vprot98pd,vfat98pd);
ulcs98=sum(hcarbo98pd     ,aprot98pd,afat98pd);

cqs98=sum(sum_whgcarb98pd,sumoth_crb_njs98pd,sumref_carb98pd);
sug_cqs98=sum(adds_cer98pd,adds_fds98pd);
run;

proc rank data=n98_nts groups=11 out=n98_nts;
var    carbo98p  sumref_carb98p  hcarbo98p   prot98p  aprot98p  vprot98p   ;
ranks  carbo98pdd sumref_carb98pdd hcarbo98pdd  prot98pdd aprot98pdd vprot98pdd  ;
run;

proc rank data=n98_nts descending groups=11 out=n98_nts ;
var    tfat98p  afat98p  vfat98p  ;
ranks  tfat98pdd afat98pdd vfat98pdd ;
run;

data n98_nts;
set n98_nts;

lfs98=sum(carbo98pdd,prot98pdd,tfat98pdd);
alfs98=sum(carbo98pdd,aprot98pdd,vfat98pdd);
vlfs98=sum(carbo98pdd,vprot98pdd,afat98pdd);
hlfs98=sum(hcarbo98pdd, vprot98pdd, afat98pdd);
ulfs98=sum(sumref_carb98pdd,aprot98pdd,vfat98pdd);
run;


%n02_nts(keep=alco02n calor02n na02n

  carbo02n 
  pot_carb02n/*Carbohydrate from potato*/
  leg_carb02n/*Carbohdyrate from legumes*/
  frtcb_njs02n/*Fruit carbohydrate minus juice sugar*/
  vegcb_npl02n/*Vegetable Carbohydrate minus potato and legume carbohydrate*/
  sumref_carb02n/*Sum of carbohydrate from potato, added sugar in foods, and refined grain, foods only*/
  sum_whgcarb02n/*Sum Carbohydrate From Wholegrain*/
  sumoth_crb_njs02n/*Sum of carbohydrate from fruit minus juice sugar, plus carb from vegetable, legume, and lact*/
  mill_grn02n/*Carbohydrate from milled whole grain*/
  int_grn02n/*Carbohydrate from intact whole grain*/
  ref_grn02n/*Carbohydrate from refined grain*/
  adds_cer02n/********Added sugar minus cereal fruit for Carb Quality Variables only********/
  adds_fds02n/********Added Sugar From Foods for Carb Quality Variables only********/

  tfat02n
  afat02n/*Animal Fat*/
  vfat02n/*Vegetable Fat*/
  dfat02n/*Dairy Fat*/
  sat02n/*Total Saturated Fat*/
  mon02n/*Total Monounsaturated Fat*/
  mufa_pl02n/*Plant Monounsaturated Fat*/
  mufa_a02n/*Animal Monounsaturated Fat*/
  poly02n/*Total Polyunsaturated Fat*/
  omega02n/*Omega 3 (f205+f226) gm, no alpha 18:3*/
  trn0202n/*Total Trans Fat*/
 
  prot02n
  aprot02n/*Animal Protein*/
  vprot02n/*Vegetable Protein*/
  dprot02n/*Dairy Protein*/
  ); 

data n02_nts;
set n02_nts;

array gram_cp carbo02n pot_carb02n leg_carb02n frtcb_njs02n vegcb_npl02n sumref_carb02n sum_whgcarb02n sumoth_crb_njs02n mill_grn02n int_grn02n ref_grn02n prot02n aprot02n vprot02n dprot02n adds_cer02n adds_fds02n;
array perc_cp carbo02p pot_carb02p leg_carb02p frtcb_njs02p vegcb_npl02p sumref_carb02p sum_whgcarb02p sumoth_crb_njs02p mill_grn02p int_grn02p ref_grn02p prot02p aprot02p vprot02p dprot02p adds_cer02p adds_fds02p;

array gram_f  tfat02n afat02n vfat02n dfat02n sat02n mon02n mufa_pl02n mufa_a02n poly02n omega02n trn0202n;
array perc_f  tfat02p afat02p vfat02p dfat02p sat02p mon02p mufa_pl02p mufa_a02p poly02p omega02p trnss02p;

do over gram_cp;
perc_cp=(4*gram_cp)/calor02n*100;
end;

do over gram_f;
perc_f=(9*gram_f)/calor02n*100;
end;

hfat02p=sum(mufa_pl02p,poly02p);
ufat02p=sum(sat02p,mufa_a02p,trnss02p);

hcarbo02p=sum(leg_carb02p,frtcb_njs02p,vegcb_npl02p,sum_whgcarb02p);

run;

proc rank data=n02_nts groups=11 out=n02_nts;
var    prot02p  aprot02p  vprot02p  tfat02p  afat02p  vfat02p  sum_whgcarb02p  sumoth_crb_njs02p  hfat02p  ufat02p;
ranks  prot02pd aprot02pd vprot02pd tfat02pd afat02pd vfat02pd sum_whgcarb02pd sumoth_crb_njs02pd hfat02pd ufat02pd;
run;

proc rank data=n02_nts descending groups=11 out=n02_nts;
var    carbo02p  sumref_carb02p  adds_cer02p  adds_fds02p  hcarbo02p;
ranks  carbo02pd sumref_carb02pd adds_cer02pd adds_fds02pd hcarbo02pd;
run;

data n02_nts;
set n02_nts;

lcs02=sum(carbo02pd,prot02pd,tfat02pd);
alcs02=sum(carbo02pd,aprot02pd,afat02pd);
vlcs02=sum(carbo02pd,vprot02pd,vfat02pd);

hlcs02=sum(sumref_carb02pd,vprot02pd,vfat02pd);
ulcs02=sum(hcarbo02pd     ,aprot02pd,afat02pd);

cqs02=sum(sum_whgcarb02pd,sumoth_crb_njs02pd,sumref_carb02pd);
sug_cqs02=sum(adds_cer02pd,adds_fds02pd);
run;

proc rank data=n02_nts groups=11 out=n02_nts;
var    carbo02p  sumref_carb02p  hcarbo02p   prot02p  aprot02p  vprot02p   ;
ranks  carbo02pdd sumref_carb02pdd hcarbo02pdd  prot02pdd aprot02pdd vprot02pdd  ;
run;

proc rank data=n02_nts descending groups=11 out=n02_nts ;
var    tfat02p  afat02p  vfat02p  ;
ranks  tfat02pdd afat02pdd vfat02pdd ;
run;

data n02_nts;
set n02_nts;

lfs02=sum(carbo02pdd,prot02pdd,tfat02pdd);
alfs02=sum(carbo02pdd,aprot02pdd,vfat02pdd);
vlfs02=sum(carbo02pdd,vprot02pdd,afat02pdd);
hlfs02=sum(hcarbo02pdd, vprot02pdd, afat02pdd);
ulfs02=sum(sumref_carb02pdd,aprot02pdd,vfat02pdd);
run;


%n06_nts(keep=alco06n calor06n na06n
  carbo06n 
  pot_carb06n/*Carbohydrate from potato*/
  leg_carb06n/*Carbohdyrate from legumes*/
  frtcb_njs06n/*Fruit carbohydrate minus juice sugar*/
  vegcb_npl06n/*Vegetable Carbohydrate minus potato and legume carbohydrate*/
  sumref_carb06n/*Sum of carbohydrate from potato, added sugar in foods, and refined grain, foods only*/
  sum_whgcarb06n/*Sum Carbohydrate From Wholegrain*/
  sumoth_crb_njs06n/*Sum of carbohydrate from fruit minus juice sugar, plus carb from vegetable, legume, and lact*/
  mill_grn06n/*Carbohydrate from milled whole grain*/
  int_grn06n/*Carbohydrate from intact whole grain*/
  ref_grn06n/*Carbohydrate from refined grain*/
  adds_cer06n/********Added sugar minus cereal fruit for Carb Quality Variables only********/
  adds_fds06n/********Added Sugar From Foods for Carb Quality Variables only********/

  tfat06n
  afat06n/*Animal Fat*/
  vfat06n/*Vegetable Fat*/
  dfat06n/*Dairy Fat*/
  sat06n/*Total Saturated Fat*/
  mon06n/*Total Monounsaturated Fat*/
  mufa_pl06n/*Plant Monounsaturated Fat*/
  mufa_a06n/*Animal Monounsaturated Fat*/
  poly06n/*Total Polyunsaturated Fat*/
  omega06n/*Omega 3 (f205+f226) gm, no alpha 18:3*/
  trn0706n/*Total Trans Fat*/
 
  prot06n
  aprot06n/*Animal Protein*/
  vprot06n/*Vegetable Protein*/
  dprot06n/*Dairy Protein*/
  ); 

data n06_nts;
set n06_nts;

array gram_cp carbo06n pot_carb06n leg_carb06n frtcb_njs06n vegcb_npl06n sumref_carb06n sum_whgcarb06n sumoth_crb_njs06n mill_grn06n int_grn06n ref_grn06n prot06n aprot06n vprot06n dprot06n adds_cer06n adds_fds06n;
array perc_cp carbo06p pot_carb06p leg_carb06p frtcb_njs06p vegcb_npl06p sumref_carb06p sum_whgcarb06p sumoth_crb_njs06p mill_grn06p int_grn06p ref_grn06p prot06p aprot06p vprot06p dprot06p adds_cer06p adds_fds06p;

array gram_f  tfat06n afat06n vfat06n dfat06n sat06n mon06n mufa_pl06n mufa_a06n poly06n omega06n trn0706n;
array perc_f  tfat06p afat06p vfat06p dfat06p sat06p mon06p mufa_pl06p mufa_a06p poly06p omega06p trnss06p;

do over gram_cp;
perc_cp=(4*gram_cp)/calor06n*100;
end;

do over gram_f;
perc_f=(9*gram_f)/calor06n*100;
end;

hfat06p=sum(mufa_pl06p,poly06p);
ufat06p=sum(sat06p,mufa_a06p,trnss06p);

hcarbo06p=sum(leg_carb06p,frtcb_njs06p,vegcb_npl06p,sum_whgcarb06p);

run;

proc rank data=n06_nts groups=11 out=n06_nts;
var    prot06p  aprot06p  vprot06p  tfat06p  afat06p  vfat06p  sum_whgcarb06p  sumoth_crb_njs06p  hfat06p  ufat06p;
ranks  prot06pd aprot06pd vprot06pd tfat06pd afat06pd vfat06pd sum_whgcarb06pd sumoth_crb_njs06pd hfat06pd ufat06pd;
run;

proc rank data=n06_nts descending groups=11 out=n06_nts;
var    carbo06p  sumref_carb06p  adds_cer06p  adds_fds06p  hcarbo06p;
ranks  carbo06pd sumref_carb06pd adds_cer06pd adds_fds06pd hcarbo06pd;
run;

data n06_nts;
set n06_nts;

lcs06=sum(carbo06pd,prot06pd,tfat06pd);
alcs06=sum(carbo06pd,aprot06pd,afat06pd);
vlcs06=sum(carbo06pd,vprot06pd,vfat06pd);

hlcs06=sum(sumref_carb06pd,vprot06pd,vfat06pd);
ulcs06=sum(hcarbo06pd     ,aprot06pd,afat06pd);

cqs06=sum(sum_whgcarb06pd,sumoth_crb_njs06pd,sumref_carb06pd);
sug_cqs06=sum(adds_cer06pd,adds_fds06pd);
run;

proc rank data=n06_nts groups=11 out=n06_nts;
var    carbo06p  sumref_carb06p  hcarbo06p   prot06p  aprot06p  vprot06p   ;
ranks  carbo06pdd sumref_carb06pdd hcarbo06pdd  prot06pdd aprot06pdd vprot06pdd  ;
run;

proc rank data=n06_nts descending groups=11 out=n06_nts ;
var    tfat06p  afat06p  vfat06p  ;
ranks  tfat06pdd afat06pdd vfat06pdd ;
run;

data n06_nts;
set n06_nts;

lfs06=sum(carbo06pdd,prot06pdd,tfat06pdd);
alfs06=sum(carbo06pdd,aprot06pdd,vfat06pdd);
vlfs06=sum(carbo06pdd,vprot06pdd,afat06pdd);
hlfs06=sum(hcarbo06pdd, vprot06pdd, afat06pdd);
ulfs06=sum(sumref_carb06pdd,aprot06pdd,vfat06pdd);
run;


%n10_nts(keep=id alco10n calor10n na10n

  carbo10n 
  pot_carb10n/*Carbohydrate from potato*/
  leg_carb10n/*Carbohdyrate from legumes*/
  frtcb_njs10n/*Fruit carbohydrate minus juice sugar*/
  vegcb_npl10n/*Vegetable Carbohydrate minus potato and legume carbohydrate*/
  sumref_carb10n/*Sum of carbohydrate from potato, added sugar in foods, and refined grain, foods only*/
  sum_whgcarb10n/*Sum Carbohydrate From Wholegrain*/
  sumoth_crb_njs10n/*Sum of carbohydrate from fruit minus juice sugar, plus carb from vegetable, legume, and lact*/
  mill_grn10n/*Carbohydrate from milled whole grain*/
  int_grn10n/*Carbohydrate from intact whole grain*/
  ref_grn10n/*Carbohydrate from refined grain*/
  adds_cer10n/********Added sugar minus cereal fruit for Carb Quality Variables only********/
  adds_fds10n/********Added Sugar From Foods for Carb Quality Variables only********/

  tfat10n
  afat10n/*Animal Fat*/
  vfat10n/*Vegetable Fat*/
  dfat10n/*Dairy Fat*/
  sat10n/*Total Saturated Fat*/
  mon10n/*Total Monounsaturated Fat*/
  mufa_pl10n/*Plant Monounsaturated Fat*/
  mufa_a10n/*Animal Monounsaturated Fat*/
  poly10n/*Total Polyunsaturated Fat*/
  omega10n/*Omega 3 (f205+f226) gm, no alpha 18:3*/
  trn1110n/*Total Trans Fat*/
 
  prot10n
  aprot10n/*Animal Protein*/
  vprot10n/*Vegetable Protein*/
  dprot10n/*Dairy Protein*/
  ); 

data n10_nts;
set n10_nts;

array gram_cp carbo10n pot_carb10n leg_carb10n frtcb_njs10n vegcb_npl10n sumref_carb10n sum_whgcarb10n sumoth_crb_njs10n mill_grn10n int_grn10n ref_grn10n prot10n aprot10n vprot10n dprot10n adds_cer10n adds_fds10n;
array perc_cp carbo10p pot_carb10p leg_carb10p frtcb_njs10p vegcb_npl10p sumref_carb10p sum_whgcarb10p sumoth_crb_njs10p mill_grn10p int_grn10p ref_grn10p prot10p aprot10p vprot10p dprot10p adds_cer10p adds_fds10p;

array gram_f  tfat10n afat10n vfat10n dfat10n sat10n mon10n mufa_pl10n mufa_a10n poly10n omega10n trn1110n;
array perc_f  tfat10p afat10p vfat10p dfat10p sat10p mon10p mufa_pl10p mufa_a10p poly10p omega10p trnss10p;

do over gram_cp;
perc_cp=(4*gram_cp)/calor10n*100;
end;

do over gram_f;
perc_f=(9*gram_f)/calor10n*100;
end;

hfat10p=sum(mufa_pl10p,poly10p);
ufat10p=sum(sat10p,mufa_a10p,trnss10p);

hcarbo10p=sum(leg_carb10p,frtcb_njs10p,vegcb_npl10p,sum_whgcarb10p);

run;

proc rank data=n10_nts groups=11 out=n10_nts;
var    prot10p  aprot10p  vprot10p  tfat10p  afat10p  vfat10p  sum_whgcarb10p  sumoth_crb_njs10p  hfat10p  ufat10p;
ranks  prot10pd aprot10pd vprot10pd tfat10pd afat10pd vfat10pd sum_whgcarb10pd sumoth_crb_njs10pd hfat10pd ufat10pd;
run;

proc rank data=n10_nts descending groups=11 out=n10_nts;
var    carbo10p  sumref_carb10p  adds_cer10p  adds_fds10p  hcarbo10p;
ranks  carbo10pd sumref_carb10pd adds_cer10pd adds_fds10pd hcarbo10pd;
run;

data n10_nts;
set n10_nts;

lcs10=sum(carbo10pd,prot10pd,tfat10pd);
alcs10=sum(carbo10pd,aprot10pd,afat10pd);
vlcs10=sum(carbo10pd,vprot10pd,vfat10pd);

hlcs10=sum(sumref_carb10pd,vprot10pd,vfat10pd);
ulcs10=sum(hcarbo10pd     ,aprot10pd,afat10pd);

cqs10=sum(sum_whgcarb10pd,sumoth_crb_njs10pd,sumref_carb10pd);
sug_cqs10=sum(adds_cer10pd,adds_fds10pd);
run;

proc rank data=n10_nts groups=11 out=n10_nts;
var    carbo10p  sumref_carb10p  hcarbo10p   prot10p  aprot10p  vprot10p   ;
ranks  carbo10pdd sumref_carb10pdd hcarbo10pdd  prot10pdd aprot10pdd vprot10pdd  ;
run;

proc rank data=n10_nts descending groups=11 out=n10_nts ;
var    tfat10p  afat10p  vfat10p  ;
ranks  tfat10pdd afat10pdd vfat10pdd ;
run;

data n10_nts;
set n10_nts;

lfs10=sum(carbo10pdd,prot10pdd,tfat10pdd);
alfs10=sum(carbo10pdd,aprot10pdd,vfat10pdd);
vlfs10=sum(carbo10pdd,vprot10pdd,afat10pdd);
hlfs10=sum(hcarbo10pdd, vprot10pdd, afat10pdd);
ulfs10=sum(sumref_carb10pdd,aprot10pdd,vfat10pdd);
run;

*----------------3.Merge-----------------------------;
    data base;
    merge der7620(in=mst) disease familyhis  
          n86_nts(in=ffq) n90_nts n94_nts  n98_nts  n02_nts  n06_nts n10_nts
          n767880 nur82 nur84 nur86 nur88 nur90 nur92 nur94 nur96 nur98 nur00 nur02 nur04 nur06 nur08 nur10
          ahei2010_8410 meds8020 supp8016 act8614
          serv86 serv90 serv94 serv98 serv02 serv06 serv10 
          nses7616 end=_end_
          ;
     by id;
     exrec=1; if first.id and mst then exrec=0; 
     baseffq=1; if  ffq then baseffq=0; 
    
          /* ocever */
          if (ocu76 in (1,2) or ocu78 in (1,2) or ocu80 in (1,2) or ocu82 in (1,2) or ocu84 in (1,2)) then ocever=1;   else ocever=0;
   
          /* antihypertensive */
           if (thiaz88=1 or ccblo88=1 or betab88=1 or ace88=1 or bprx88=1) then htnrx88=1;   else htnrx88=0;
           if (thiaz94=1 or lasix94 =1 or ccblo94 =1 or betab94 =1 or bprx94=1 ) then htnrx94=1;   else htnrx94=0;
           if (thiaz96=1 or lasix96 =1 or ccblo96 =1 or betab96 =1 or ace96=1 or bprx96=1 ) then htnrx96=1;   else htnrx96=0;
           if (thiaz98=1 or lasix98 =1 or ccblo98 =1 or betab98 =1 or ace98=1 or bprx98=1 ) then htnrx98=1;   else htnrx98=0;
           if (thiaz00=1 or lasix00 =1 or ccblo00 =1 or betab00 =1 or ace00=1 or bprx00=1 ) then htnrx00=1;   else htnrx00=0;
           if (thiaz02=1 or lasix02 =1 or ccblo02 =1 or betab02 =1 or ace02=1 or bprx02=1 ) then htnrx02=1;   else htnrx02=0;
           if (bprex04=1 or lasix04=1 or  thiaz04=1 or ace04=1 or ccblo04=1 or  betab04=1 or  k04=1) then htnrx04=1;   else htnrx04=0;
           if (bprx06=1 or lasix06=1 or thiaz06=1 or ace06=1 or ccblo06=1 or betab06=1 or k06=1) then htnrx06=1;   else htnrx06=0;
           if (bprx08=1 or lasix08=1 or thiaz08=1 or ace08=1 or ccblo08=1 or betab08=1 or k08=1) then htnrx08=1;   else htnrx08=0;
           if (bprx10=1 or lasix10=1 or thiaz10=1 or ace10=1 or ccblo10=1 or betab10=1 or angio10=1 or k10=1) then htnrx10=1;  else htnrx10=0;

          /* heart diseases medications */
           if (couma94=1 or digox94=1 or betab94 =1 or ccblo94 =1) then chdrx94=1; else chdrx94=0;/***digox94 inclde digoxin or antiarrhythmic***/
           if (couma96=1 or digox96=1 or antia96=1 or betab96 =1 or ccblo96 =1) then chdrx96=1; else chdrx96=0;
           if (couma98=1 or digox98=1 or antia98=1 or betab98 =1 or ccblo98 =1) then chdrx98=1; else chdrx98=0;
           if (couma00=1 or digox00=1 or antia00=1 or betab00 =1 or ccblo00 =1) then chdrx00=1; else chdrx00=0;
           if (couma02=1 or digox02=1 or antia02=1 or betab02 =1 or ccblo02 =1) then chdrx02=1; else chdrx02=0;

          /* antihtc */
           if clrx88=1 then hchtx88=1; else hchtx88=0;
           if clrx94=1 then hchtx94=1; else hchtx94=0;
           if chrx96=1 then hchtx96=1; else hchtx96=0;
           if chrx98=1 then hchtx98=1; else hchtx98=0;
           if chrx00=1 or oclrx00=1 then hchtx00=1; else hchtx00=0;
           if chrx02=1 or oclrx02=1 then hchtx02=1; else hchtx02=0;
           if mev04=1 or zoc04=1 or crest04=1 or prav04=1 or lip04=1 or lesc04=1 or oclrx04=1 then hchtx04=1; else hchtx04=0;
           if mev06=1 or zoc06=1 or crest06=1 or prav06=1 or lip06=1 or ostat06=1 or oclrx06=1 then hchtx06=1; else hchtx06=0;
           if mev08=1 or zoc08=1 or crest08=1 or prav08=1 or lip08=1 or ostat08=1 or oclrx08=1 then hchtx08=1; else hchtx08=0;
           if stat10=1 or mev10=1 or prav10=1 or zoc10=1 or lip10=1 or crest10=1 or ostat10=1 or oclrx10=1 then hchtx10=1; else hchtx10=0;


*** nSES - impute to missing;
if nSES_86 = . then nSES_86m = -0.464;   else nSES_86m = nSES_86;
if nSES_88 = . then nSES_88m = nSES_86m; else nSES_88m = nSES_88;
if nSES_90 = . then nSES_90m = nSES_88m; else nSES_90m = nSES_90;
if nSES_92 = . then nSES_92m = nSES_90m; else nSES_92m = nSES_92;
if nSES_94 = . then nSES_94m = nSES_92m; else nSES_94m = nSES_94;
if nSES_96 = . then nSES_96m = nSES_94m; else nSES_96m = nSES_96;
if nSES_98 = . then nSES_98m = nSES_96m; else nSES_98m = nSES_98;
if nSES_00 = . then nSES_00m = nSES_98m; else nSES_00m = nSES_00;
if nSES_02 = . then nSES_02m = nSES_00m; else nSES_02m = nSES_02;
if nSES_04 = . then nSES_04m = nSES_02m; else nSES_04m = nSES_04;
if nSES_06 = . then nSES_06m = nSES_04m; else nSES_06m = nSES_06;
if nSES_08 = . then nSES_08m = nSES_06m; else nSES_08m = nSES_08;
if nSES_10 = . then nSES_10m = nSES_08m; else nSES_10m = nSES_10;


/************** cancer: any confirmed cancer (from registry and self-reported) *****************/
  array canc_prev {*} can76  can78  can80  can82  can84 can86  can88  can90  can92  can94 can96  can98  can00  can02  can04 can06  can08  can10;
  array canc_new  {*} cann76  cann78  cann80  cann82  cann84 cann86 cann88 cann90 cann92 cann94 cann96 cann98 cann00 cann02 cann04 cann06 cann08 cann10;

   /* set first cancer */
   if dtdx_allca > 0 then first_dtdx_ca = dtdx_allca; /* any confirmed cancer */
   else if dtdxca > 0 then first_dtdx_ca = dtdxca;
   else first_dtdx_ca = .;

    array qqyear {*} irt76 irt78   irt80   irt82   irt84 irt86   irt88   irt90   irt92   irt94   irt96   irt98   irt00   irt02   irt04   irt06   irt08   irt10; 

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
    *make intermediate outcomes cumulative varibles: i.e. once cancer, always cancer;
    /*array canc_prev{18} can76  can78  can80  can82  can84 can86  can88  can90  can92  can94 
                        can96  can98  can00  can02  can04 can06  can08  can10;*/
    array hrt {18} hrt76  hrt78  hrt80  hrt82  hrt84 hrt86  hrt88  hrt90  hrt92  hrt94
                   hrt96 hrt98 hrt00 hrt02 hrt04 hrt06 hrt08 hrt10;
    array ang {18} ang76  ang78  ang80  ang82  ang84 ang86  ang88  ang90  ang92  ang94
                   ang96  ang98  ang00  ang02  ang04 ang06  ang08  ang10;
    array mi  {18} mi76  mi78  mi80  mi82  mi84 mi86  mi88  mi90  mi92  mi94
                   mi96  mi98  mi00  mi02  mi04 mi06  mi08  mi10;
    array hbp {18} hbp76 hbp78 hbp80 hbp82 hbp84 hbp86 hbp88 hbp90 hbp92 hbp94
                   hbp96 hbp98 hbp00 hbp02 hbp04 hbp06 hbp08 hbp10;
    array chol{18} chol76 chol78 chol80 chol82 chol84 chol86 chol88 chol90 chol92 chol94
                   chol96 chol98 chol00 chol02 chol04 chol06 chol08 chol10;
    array diab{18} db76 db78 db80 db82 db84 db86 db88 db90 db92 db94 
                   db96 db98 db00 db02 db04 db06 db08 db10;
    array strk{18} str82 str82 str82 str82 str84 str86 str88 str90 str92 str94
                   str96 str98 str00 str02 str04 str06 str08 str10;
    array cabg{18} cabg84 cabg84 cabg84 cabg84 cabg84 cabg86 cabg88 cabg90 cabg92 cabg94
                   cabg96 cabg98 cabg00 cabg02 cabg04 cabg06 cabg08 cabg10;
    array chrodis {18} chrodis76 chrodis78 chrodis80 chrodis82 chrodis84 chrodis86 chrodis88 chrodis90 chrodis92 
                       chrodis94 chrodis96 chrodis98 chrodis00 chrodis02 chrodis04 chrodis06 chrodis08 chrodis10;

    do i=1 to 18;
        if hrt{i} in (1,3) then mi  {i} = 1;
        if hrt{i} in (2,3) then ang {i} = 1; end;
    
    do i=2 to 18;
        if canc_prev {i-1} = 1   then canc_prev {i} = 1;
        *if ang {i-1} = 1   then ang  {i} = 1;
        if mi   {i-1} = 1   then mi   {i} = 1;
        if hbp  {i-1} = 1   then hbp  {i} = 1;
        if chol {i-1} = 1   then chol {i} = 1;
        if diab {i-1} = 1   then diab {i} = 1;
        if strk {i-1} = 1   then strk {i} = 1;
        if cabg {i-1} = 1   then cabg {i} = 1; 
    end;

    do i=1 to 18;
        if canc_prev {i} ne 1    then canc_prev {i} = 0;
        *if ang {i} ne 1    then ang  {i} = 0;
        if mi   {i} ne 1    then mi   {i} = 0;
        if hbp  {i} ne 1    then hbp  {i} = 0;
        if chol {i} ne 1    then chol {i} = 0;
        if diab {i} ne 1    then diab {i} = 0;
        if strk {i} ne 1    then strk {i} = 0;
        if cabg {i} ne 1    then cabg {i} = 0;
        if mi(i) =1 or strk(i)=1 or cabg(i)=1  or diab(i)=1 or canc_prev(i)=1  then chrodis(i)=1; else chrodis(i)=0; 
    end;

    do i = 2 to dim(chrodis); if chrodis{i-1} = 1 then chrodis{i} = 1; 
    end; /* Make chrodis cumulative */

/***************** Diet scores: carry forward + STOP updating after cancer/CVD *******************/
/* updated  scores */  
array dis    {*}  chrodis86 chrodis90 chrodis94 chrodis98 chrodis02 chrodis06 chrodis10;
array calor  {*}  calor86n  calor90n  calor94n  calor98n  calor02n  calor06n  calor10n;   
array aheim  {*}  ahei86  ahei90  ahei94  ahei98  ahei02   ahei06  ahei10;
array lcs    {*}  lcs86    lcs90    lcs94    lcs98    lcs02    lcs06    lcs10;
array alcs   {*}  alcs86   alcs90   alcs94   alcs98   alcs02   alcs06   alcs10;
array vlcs   {*}  vlcs86   vlcs90   vlcs94   vlcs98   vlcs02   vlcs06   vlcs10;
array ulcs   {*}  ulcs86   ulcs90   ulcs94   ulcs98   ulcs02   ulcs06   ulcs10;
array hlcs   {*}  hlcs86   hlcs90   hlcs94   hlcs98   hlcs02   hlcs06   hlcs10;
array lfs    {*}  lfs86    lfs90    lfs94    lfs98    lfs02    lfs06    lfs10;
array ulfs   {*}  ulfs86   ulfs90   ulfs94   ulfs98   ulfs02   ulfs06   ulfs10;
array hlfs   {*}  hlfs86   hlfs90   hlfs94   hlfs98   hlfs02   hlfs06   hlfs10;
array alfs  {*}  alfs86   alfs90   alfs94   alfs98   alfs02   alfs06   alfs10 ;
array vlfs  {*}  vlfs86   vlfs90   vlfs94   vlfs98   vlfs02   vlfs06   vlfs10 ;


array    carbp  {*}   carbo86p         carbo90p        carbo94p         carbo98p         carbo02p        carbo06p      carbo10p  ;
array    hcarbp  {*}   hcarbo86p        hcarbo90p       hcarbo94p        hcarbo98p        hcarbo02p       hcarbo06p     hcarbo10p   ;
array sumref_carbp {*}  sumref_carb86p   sumref_carb90p  sumref_carb94p   sumref_carb98p   sumref_carb02p  sumref_carb06p  sumref_carb10p ;
array     protp  {*}   prot86p          prot90p         prot94p          prot98p          prot02p         prot06p       prot10p ;
array    aprotp  {*}   aprot86p         aprot90p        aprot94p         aprot98p         aprot02p        aprot06p      aprot10p;
array    vprotp  {*}   vprot86p         vprot90p        vprot94p         vprot98p         vprot02p        vprot06p      vprot10p;
array    tfatp   {*}   tfat86p          tfat90p         tfat94p          tfat98p          tfat02p         tfat06p       tfat10p;
array    afatp   {*}   afat86p          afat90p         afat94p          afat98p          afat02p         afat06p       afat10p;
array    vfatp   {*}   vfat86p          vfat90p         vfat94p          vfat98p          vfat02p         vfat06p       vfat10p;

array    carbs  {*}   carbo86s         carbo90s        carbo94s         carbo98s         carbo02s        carbo06s         carbo10s         ;
array    hcarbs  {*}   hcarbo86s        hcarbo90s       hcarbo94s        hcarbo98s        hcarbo02s       hcarbo06s        hcarbo10s        ;
array sumref_carbs {*}  sumref_carb86s   sumref_carb90s  sumref_carb94s   sumref_carb98s   sumref_carb02s  sumref_carb06s   sumref_carb10s   ;
array     prots  {*}   prot86s          prot90s         prot94s          prot98s          prot02s         prot06s          prot10s           ;
array    aprots  {*}   aprot86s         aprot90s        aprot94s         aprot98s         aprot02s        aprot06s         aprot10s          ;
array    vprots  {*}   vprot86s         vprot90s        vprot94s         vprot98s         vprot02s        vprot06s         vprot10s          ;
array    tfats   {*}   tfat86s          tfat90s         tfat94s          tfat98s          tfat02s         tfat06s          tfat10s           ;
array    afats   {*}   afat86s          afat90s         afat94s          afat98s          afat02s         afat06s          afat10s          ;
array    vfats   {*}   vfat86s          vfat90s         vfat94s          vfat98s          vfat02s         vfat06s          vfat10s          ;


array    whgrn_o  {*}   whgrn_s86       whgrn_s90      whgrn_s94       whgrn_s98       whgrn_s02      whgrn_s06    whgrn_s10;
array    rgrain_o {*}   rgrain_s86      rgrain_s90     rgrain_s94      rgrain_s98      rgrain_s02     rgrain_s06   rgrain_s10 ;
array    frui_o   {*}   frui_s86    frui_s90    frui_s94 frui_s98 frui_s02 frui_s06 frui_s10 ;
array     leg_o   {*}   leg_s86        leg_s90         leg_s94        leg_s98        leg_s02       leg_s06     leg_s10  ;
array    pota_o   {*}   pota_s86       pota_s90        pota_s94       pota_s98       pota_s02      pota_s06    pota_s10 ;
array  veg_noleg_o{*}   veg_noleg_s86   veg_noleg_s90  veg_noleg_s94       veg_noleg_s98       veg_noleg_s02      veg_noleg_s06    veg_noleg_s10  ;
array    fish_o   {*}   fish_s86        fish_s90       fish_s94        fish_s98        fish_s02       fish_s06     fish_s10 ;
array   rpmeats_o {*}   rpmeats_s86     rpmeats_s90    rpmeats_s94        rpmeats_s98        rpmeats_s02       rpmeats_s06     rpmeats_s10  ;
array    ssb_o    {*}   ssb_s86         ssb_s90        ssb_s94        ssb_s98        ssb_s02       ssb_s06     ssb_s10  ;
array    dairy_o   {*}   dairy_s86        dairy_s90       dairy_s94        dairy_s98        dairy_s02       dairy_s06     dairy_s10 ;
array swtsnochoc_o {*}   swtsnochoc_s86     swtsnochoc_s90    swtsnochoc_s94        swtsnochoc_s98        swtsnochoc_s02       swtsnochoc_s06     swtsnochoc_s10  ;
array    poul_o    {*}   poul_s86         poul_s90        poul_s94        poul_s98        poul_s02       poul_s06     poul_s10  ;


/* stop updated scores */
array calors {*}  calor86ns    calor90ns    calor94ns    calor98ns    calor02ns    calor06ns    calor10ns;   
array aheims {*}  ahei86s  ahei90s  ahei94s  ahei98s  ahei02s   ahei06s  ahei10s;
array lcss   {*}  lcs86s    lcs90s    lcs94s    lcs98s    lcs02s    lcs06s    lcs10s;
array alcss  {*}  alcs86s   alcs90s   alcs94s   alcs98s   alcs02s   alcs06s   alcs10s;
array vlcss  {*}  vlcs86s   vlcs90s   vlcs94s   vlcs98s   vlcs02s   vlcs06s   vlcs10s;
array ulcss  {*}  ulcs86s   ulcs90s   ulcs94s   ulcs98s   ulcs02s   ulcs06s   ulcs10s;
array hlcss  {*}  hlcs86s   hlcs90s   hlcs94s   hlcs98s   hlcs02s   hlcs06s   hlcs10s;
array lfss   {*}  lfs86s    lfs90s    lfs94s    lfs98s    lfs02s    lfs06s    lfs10s;
array ulfss  {*}  ulfs86s   ulfs90s   ulfs94s   ulfs98s   ulfs02s   ulfs06s   ulfs10s;
array hlfss  {*}  hlfs86s   hlfs90s   hlfs94s   hlfs98s   hlfs02s   hlfs06s   hlfs10s;
array alfss {*} alfs86s alfs90s alfs94s alfs98s alfs02s alfs06s alfs10s ;
array vlfss {*} vlfs86s vlfs90s vlfs94s vlfs98s vlfs02s vlfs06s vlfs10s ;

array    whgrn_os  {*}   whgrn_s86s       whgrn_s90s      whgrn_s94s       whgrn_s98s       whgrn_s02s      whgrn_s06s    whgrn_s10s;
array    rgrain_os {*}   rgrain_s86s      rgrain_s90s     rgrain_s94s      rgrain_s98s      rgrain_s02s     rgrain_s06s   rgrain_s10s ;
array    frui_os   {*}   frui_s86s    frui_s90s    frui_s94s frui_s98s frui_s02s frui_s06s frui_s10s ;
array     leg_os   {*}   leg_s86s        leg_s90s         leg_s94s        leg_s98s        leg_s02s       leg_s06s     leg_s10s  ;
array    pota_os   {*}   pota_s86s       pota_s90s        pota_s94s       pota_s98s       pota_s02s      pota_s06s    pota_s10s ;
array  veg_noleg_os{*}   veg_noleg_s86s   veg_noleg_s90s  veg_noleg_s94s       veg_noleg_s98s       veg_noleg_s02s      veg_noleg_s06s    veg_noleg_s10s  ;
array    fish_os   {*}   fish_s86s        fish_s90s       fish_s94s        fish_s98s        fish_s02s       fish_s06s     fish_s10s ;
array   rpmeats_os {*}   rpmeats_s86s     rpmeats_s90s    rpmeats_s94s        rpmeats_s98s        rpmeats_s02s       rpmeats_s06s     rpmeats_s10s  ;
array    ssb_os    {*}   ssb_s86s         ssb_s90s        ssb_s94s        ssb_s98s        ssb_s02s       ssb_s06s     ssb_s10s  ;
array    dairy_os   {*}   dairy_s86s        dairy_s90s       dairy_s94s        dairy_s98s        dairy_s02s       dairy_s06s     dairy_s10s ;
array swtsnochoc_os {*}   swtsnochoc_s86s     swtsnochoc_s90s    swtsnochoc_s94s        swtsnochoc_s98s        swtsnochoc_s02s       swtsnochoc_s06s     swtsnochoc_s10s  ;
array    poul_os    {*}   poul_s86s         poul_s90s        poul_s94s        poul_s98s        poul_s02s       poul_s06s     poul_s10s  ;


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

if whgrn_o{i}=.   then whgrn_o{i}=whgrn_o{i-1};
if rgrain_o{i}=.   then rgrain_o{i}=rgrain_o{i-1};
if frui_o{i}=.     then frui_o{i}=frui_o{i-1};
if leg_o{i}=.    then leg_o{i}=leg_o{i-1};
if pota_o{i}=.    then pota_o{i}=pota_o{i-1};
if veg_noleg_o{i}=.    then veg_noleg_o{i}=veg_noleg_o{i-1};
if fish_o{i}=.    then fish_o{i}=fish_o{i-1};
if rpmeats_o{i}=.     then rpmeats_o{i}=rpmeats_o{i-1};
if ssb_o{i}=.    then ssb_o{i}=ssb_o{i-1};
if dairy_o{i}=.    then dairy_o{i}=dairy_o{i-1};
if swtsnochoc_o{i}=.   then swtsnochoc_o{i}=swtsnochoc_o{i-1};
if poul_o{i}=.   then poul_o{i}=poul_o{i-1};
end;

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

  carbs{i}=carbp{i};  
  hcarbs{i}=hcarbp{i};
  sumref_carbs{i}=sumref_carbp{i};  
  prots{i}=protp{i};
  aprots{i}=aprotp{i};
  vprots{i}=vprotp{i};
  tfats{i}=tfatp{i};
  afats{i}=afatp{i};
  vfats{i}=vfatp{i};

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


do i=2 to DIM(carbp); /* STOP updating scores after cancer OR CVD */
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

/*calculate cumulative average of updated score**/
array nut{7,12}
calor86ns  ahei86s  lcs86s  alcs86s vlcs86s  ulcs86s  hlcs86s  lfs86s  ulfs86s  hlfs86s alfs86s  vlfs86s
calor90ns  ahei90s  lcs90s  alcs90s vlcs90s  ulcs90s  hlcs90s  lfs90s  ulfs90s  hlfs90s alfs90s  vlfs90s
calor94ns  ahei94s  lcs94s  alcs94s vlcs94s  ulcs94s  hlcs94s  lfs94s  ulfs94s  hlfs94s alfs94s  vlfs94s
calor98ns  ahei98s  lcs98s  alcs98s vlcs98s  ulcs98s  hlcs98s  lfs98s  ulfs98s  hlfs98s alfs98s  vlfs98s
calor02ns  ahei02s  lcs02s  alcs02s vlcs02s  ulcs02s  hlcs02s  lfs02s  ulfs02s  hlfs02s alfs02s  vlfs02s
calor06ns  ahei06s  lcs06s  alcs06s vlcs06s  ulcs06s  hlcs06s  lfs06s  ulfs06s  hlfs06s alfs06s  vlfs06s
calor10ns  ahei10s  lcs10s  alcs10s vlcs10s  ulcs10s  hlcs10s  lfs10s  ulfs10s  hlfs10s alfs10s  vlfs10s
;
array nutv{7,12}
calor86nsv  ahei86sv  lcs86sv  alcs86sv vlcs86sv  ulcs86sv  hlcs86sv  lfs86sv  ulfs86sv  hlfs86sv alfs86sv  vlfs86sv
calor90nsv  ahei90sv  lcs90sv  alcs90sv vlcs90sv  ulcs90sv  hlcs90sv  lfs90sv  ulfs90sv  hlfs90sv alfs90sv  vlfs90sv
calor94nsv  ahei94sv  lcs94sv  alcs94sv vlcs94sv  ulcs94sv  hlcs94sv  lfs94sv  ulfs94sv  hlfs94sv alfs94sv  vlfs94sv
calor98nsv  ahei98sv  lcs98sv  alcs98sv vlcs98sv  ulcs98sv  hlcs98sv  lfs98sv  ulfs98sv  hlfs98sv alfs98sv  vlfs98sv
calor02nsv  ahei02sv  lcs02sv  alcs02sv vlcs02sv  ulcs02sv  hlcs02sv  lfs02sv  ulfs02sv  hlfs02sv alfs02sv  vlfs02sv
calor06nsv  ahei06sv  lcs06sv  alcs06sv vlcs06sv  ulcs06sv  hlcs06sv  lfs06sv  ulfs06sv  hlfs06sv alfs06sv  vlfs06sv
calor10nsv  ahei10sv  lcs10sv  alcs10sv vlcs10sv  ulcs10sv  hlcs10sv  lfs10sv  ulfs10sv  hlfs10sv alfs10sv  vlfs10sv
;

/* calculate cumulative average intake */
do j=1 to 12;
  nutv{1,j}=nut{1,j};
 do i=2 to 7;
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

/* calculate cumulative average intake */
array foodlist{7,12}
whgrn_s86s  rgrain_s86s  frui_s86s  leg_s86s pota_s86s  veg_noleg_s86s  fish_s86s  rpmeats_s86s  ssb_s86s  dairy_s86s swtsnochoc_s86s  poul_s86s
whgrn_s90s  rgrain_s90s  frui_s90s  leg_s90s pota_s90s  veg_noleg_s90s  fish_s90s  rpmeats_s90s  ssb_s90s  dairy_s90s swtsnochoc_s90s  poul_s90s
whgrn_s94s  rgrain_s94s  frui_s94s  leg_s94s pota_s94s  veg_noleg_s94s  fish_s94s  rpmeats_s94s  ssb_s94s  dairy_s94s swtsnochoc_s94s  poul_s94s
whgrn_s98s  rgrain_s98s  frui_s98s  leg_s98s pota_s98s  veg_noleg_s98s  fish_s98s  rpmeats_s98s  ssb_s98s  dairy_s98s swtsnochoc_s98s  poul_s98s
whgrn_s02s  rgrain_s02s  frui_s02s  leg_s02s pota_s02s  veg_noleg_s02s  fish_s02s  rpmeats_s02s  ssb_s02s  dairy_s02s swtsnochoc_s02s  poul_s02s
whgrn_s06s  rgrain_s06s  frui_s06s  leg_s06s pota_s06s  veg_noleg_s06s  fish_s06s  rpmeats_s06s  ssb_s06s  dairy_s06s swtsnochoc_s06s  poul_s06s
whgrn_s10s  rgrain_s10s  frui_s10s  leg_s10s pota_s10s  veg_noleg_s10s  fish_s10s  rpmeats_s10s  ssb_s10s  dairy_s10s swtsnochoc_s10s  poul_s10s
;
array foodlistv{7,12}
whgrn_s86sv  rgrain_s86sv  frui_s86sv  leg_s86sv pota_s86sv  veg_noleg_s86sv  fish_s86sv  rpmeats_s86sv  ssb_s86sv  dairy_s86sv swtsnochoc_s86sv  poul_s86sv
whgrn_s90sv  rgrain_s90sv  frui_s90sv  leg_s90sv pota_s90sv  veg_noleg_s90sv  fish_s90sv  rpmeats_s90sv  ssb_s90sv  dairy_s90sv swtsnochoc_s90sv  poul_s90sv
whgrn_s94sv  rgrain_s94sv  frui_s94sv  leg_s94sv pota_s94sv  veg_noleg_s94sv  fish_s94sv  rpmeats_s94sv  ssb_s94sv  dairy_s94sv swtsnochoc_s94sv  poul_s94sv
whgrn_s98sv  rgrain_s98sv  frui_s98sv  leg_s98sv pota_s98sv  veg_noleg_s98sv  fish_s98sv  rpmeats_s98sv  ssb_s98sv  dairy_s98sv swtsnochoc_s98sv  poul_s98sv
whgrn_s02sv  rgrain_s02sv  frui_s02sv  leg_s02sv pota_s02sv  veg_noleg_s02sv  fish_s02sv  rpmeats_s02sv  ssb_s02sv  dairy_s02sv swtsnochoc_s02sv  poul_s02sv
whgrn_s06sv  rgrain_s06sv  frui_s06sv  leg_s06sv pota_s06sv  veg_noleg_s06sv  fish_s06sv  rpmeats_s06sv  ssb_s06sv  dairy_s06sv swtsnochoc_s06sv  poul_s06sv
whgrn_s10sv  rgrain_s10sv  frui_s10sv  leg_s10sv pota_s10sv  veg_noleg_s10sv  fish_s10sv  rpmeats_s10sv  ssb_s10sv  dairy_s10sv swtsnochoc_s10sv  poul_s10sv
;

do j=1 to 12;
  foodlistv{1,j}=foodlist{1,j};
 do i=2 to 7;
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


/* calculate cumulative average */
array macron{7,9}
carbo86s hcarbo86s sumref_carb86s  prot86s aprot86s  vprot86s  tfat86s  afat86s  vfat86s
carbo90s hcarbo90s sumref_carb90s  prot90s aprot90s  vprot90s  tfat90s  afat90s  vfat90s
carbo94s hcarbo94s sumref_carb94s  prot94s aprot94s  vprot94s  tfat94s  afat94s  vfat94s
carbo98s hcarbo98s sumref_carb98s  prot98s aprot98s  vprot98s  tfat98s  afat98s  vfat98s
carbo02s hcarbo02s sumref_carb02s  prot02s aprot02s  vprot02s  tfat02s  afat02s  vfat02s
carbo06s hcarbo06s sumref_carb06s  prot06s aprot06s  vprot06s  tfat06s  afat06s  vfat06s
carbo10s hcarbo10s sumref_carb10s  prot10s aprot10s  vprot10s  tfat10s  afat10s  vfat10s 
;
array macronv{7,9}
carbo86sv  hcarbo86sv  sumref_carb86sv prot86sv aprot86sv  vprot86sv  tfat86sv  afat86sv  vfat86sv
carbo90sv  hcarbo90sv  sumref_carb90sv prot90sv aprot90sv  vprot90sv  tfat90sv  afat90sv  vfat90sv
carbo94sv  hcarbo94sv  sumref_carb94sv prot94sv aprot94sv  vprot94sv  tfat94sv  afat94sv  vfat94sv
carbo98sv  hcarbo98sv  sumref_carb98sv prot98sv aprot98sv  vprot98sv  tfat98sv  afat98sv  vfat98sv
carbo02sv  hcarbo02sv  sumref_carb02sv prot02sv aprot02sv  vprot02sv  tfat02sv  afat02sv  vfat02sv
carbo06sv  hcarbo06sv  sumref_carb06sv prot06sv aprot06sv  vprot06sv  tfat06sv  afat06sv  vfat06sv
carbo10sv  hcarbo10sv  sumref_carb10sv prot10sv aprot10sv  vprot10sv  tfat10sv  afat10sv  vfat10sv 
;

/* calculate cumulative average intake */
do j=1 to 9;
  macronv{1,j}=macron{1,j};
 do i=2 to 7;
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


/****** MAMMOGRAPHY USE *******/
* use one i+1 mem memography status *;
/*Mammogram in last 2 yrs - yes/no/missing*/

if mamm90 in (2,3) or mam8994=1 or mam9094=1 then m88=1;
if mamm90=1 then m88=0; 
if mamm92 in (2,3) or mam9194=1 or mam9294=1 then m90=1;
if mamm92=1 then m90=0; 
if mam9394=1 or mam9494=1 then m92=1;
if ndmam94=1 then m92=0;
if mamsc96=1 or mamsy96=1 then m94=1;
if mamno96=1 then m94=0;
if mamsc98=1 or mamsy98=1 then m96=1;
if mamno98=1 then m96=0;
if mamsc00=1 or mamsy00=1 then m98=1;
if mamno00=1 then m98=0;
if mamsc02=1 or mamsy02=1 then m00=1;
if mamno02=1 then m00=0;
if mamsc04=1 or mamsy04=1 then m02=1; 
if mamno04=1 then m02=0;
if mamsc06=1 or mamsy06=1 then m04=1; 
if mamno06=1 then m04=0;
if mamsc08=1 or mamsy08=1 then m06=1; 
if mamno08=1 then m06=0;
if mamsc10=1 or mamsy10=1 then m08=1; 
if mamno10=1 then m08=0;



/***********/
/*endoscopy*/
/***********/
if sigm88=1 then endo88=1; /*endo no*/
else if sigm88 in (2,3) then endo88=2; /*endo yes*/

if sigm90=1 then endo90=2; /*endo yes*/
else if sigm90=2 then endo90=1; /*endo no*/

if sigm92=1 then endo92=1; /*endo no*/
else if sigm92 in (2,3) then endo92=2; /*endo yes*/

if sigm94=1 then endo94=1; /*endo no*/
else if sigm94=2 then endo94=2; /*endo yes*/

if sigm96=1 then endo96=1; /*endo no*/
else if sigm96=2 then endo96=2; /*endo yes*/

if sigm98=1 then endo98=1; /*endo no*/
else if sigm98=2 then endo98=2; /*endo yes*/

if sigm00=1 then endo00=2; /*endo yes*/
else if sigm00=2 then endo00=1; /*endo no*/

if sigm02=1 then endo02=1; /*endo no*/
else if sigm02=2 then endo02=2; /*endo yes*/

if colsc04=2 or sigsc04=2 or virtc04=2 or sigm04=2 or sig0204=1 or col0204=1 then endo04=2; /*endo yes*/
else if colsc04=1 or sigsc04=1 or virtc04=1 then endo04=1; /*endo no*/

if colsc06=2 or sigsc06=2 or virtc06=2 then endo06=2; /*endo yes*/
else if colsc06=1 or sigsc06=1 or virtc06=1 then endo06=1; /*endo no*/

if colsc08=2 or sigsc08=2 or virtc08=2 then endo08=2; /*endo yes*/
else if colsc08=1 or sigsc08=1 or virtc08=1 then endo08=1; /*endo no*/

if colsc10=2 or sigsc10=2 or virtc10=2 then endo10=2; /*endo yes*/
else if colsc10=1 or sigsc10=1 or virtc10=1 then endo10=1; /*endo no*/

*-----------*
 Use 1990 Qx
*-----------*;
/*endo76-endo86 derived from 1990 info.*/
/*endo76-endo86 indicate experience of endoscopy before 1980 (1=BEFORE 1980)*/
if sigm190 in (1) or sigmr90 in (1) then do; endo80=2; endo78=2; endo76=2; end;
else if sigm190 in (2,3,4,5,6,7,8) then do;  endo80=1; endo78=1; endo76=1; end;

if sigm190 in (2) or sigmr90 in (2) then do; endo82=2; endo84=2; end; /*we will clean successive endo later*/
else if sigm190 in (3,4,5,6,7,8) or sigmr90 in (3) then do; endo82=1; endo84=1; end; /*if recent endo is 1984-85*/

if sigm190 in (3) or sigmr90 in (3) then endo86=2;
else if (sigm190 in (4,5,6,7,8) or sigmr90 in (4,5,6)) then endo86=1; /*if recent endo is within 5 years*/

/*when endo88-endo90 were missing, 1990 info was used to fill the answer.
 We don't impute endo88 if endo90=yes becase person shouldn't take
 endoscopy twice within 4 years. Same for endo90*/
if ((endo90 ne 2) and (sigm190 in (4,5) or sigmr90 in (4,5)) and endo88=. )
   then endo88=2;
else if (endo88=. and sigm190 in (6,7,8)) then endo88=1;
else if (endo88=. and sigmr90 in (3,6,7,8)) then endo88=1;

if ((endo88 ne 2) and (sigm190 in (6,7,8) or sigmr90 in (6,7,8)) and endo90=. )
   then endo90=2;
else if (endo90=. and sigmr90 in (2,3,4,5)) then endo90=1;

*-----------*
 Use 2004 Qx
*-----------*;
/*For only missing variables of the past Qx, we use 2004 info*/
/*We do not impute endo[i] if person already endo[i-1]=yes or endo[i+1]=yes*/
/*never*/
if (endo76=2 or endo78=2 or endo80=2 or endo82=2 or endo84=2 or
    endo86=2 or endo88=2 or endo90=2 or endo92=2 or endo94=2 or
    endo96=2 or endo98=2 or endo00=2 or endo02=2 or endo04=2 )
   then atleast=1; /*yes edo at least once in some year. They cannot be "Never endo"*/
array arendo[*] endo76 endo78 endo80 endo82 endo84 endo86 endo88 endo90
                endo92 endo94 endo96 endo98 endo00 endo02 endo04;
if (atleast ne 1 and signo04=1 and colno04=1) then do; /*never in 2004 qx*/
   do i=1 to dim(arendo);
      if arendo[i]=. then arendo[i]=1; /*no*/
   end;
end;

/*we fill from the latest information*/
/*2002/Later endo can fill endo04*/
if (endo02 ne 2 and endo04=. and (sig0204=1 or col0204=1))
   then endo04=2; /*yes*/
else if (endo04=. and sig0204 ne 1 and col0204 ne 1 and q04 in (1,2,3)/*all04miss ne 1 and q19pt04 ne 1*/)
   then endo04=1; /*no*/                               /*q04=1,2,3 includes endoscopy question*/
/*2000-2001 endo can fill endo02*/
if (endo00 ne 2 and endo04 ne 2 and endo02=. and (sig0004=1 or col0004=1))
   then endo02=2; /*yes*/
else if (endo02=. and sig0004 ne 1 and col0004 ne 1 and q04 in (1,2,3)/*all04miss ne 1 and q19pt04 ne 1*/)
   then endo02=1; /*no*/
/*1998-1999 endo can fill endo00*/
if (endo98 ne 2 and endo02 ne 2 and endo00=. and (sig9804=1 or col9804=1))
   then endo00=2; /*yes*/
else if (endo00=. and sig9804 ne 1 and col9804 ne 1 and q04 in (1,2,3)/*all04miss ne 1 and q19pt04 ne 1*/)
   then endo00=1; /*no*/
/*1996-1997 endo can fill endo98*/
if (endo96 ne 2 and endo00 ne 2 and endo98=. and (sig9604=1 or col9604=1))
   then endo98=2; /*yes*/
else if (endo98=. and sig9604 ne 1 and col9604 ne 1 and q04 in (1,2,3)/*all04miss ne 1 and q19pt04 ne 1*/)
   then endo98=1; /*no*/
/*1994-1995 endo can fill endo96*/
if (endo94 ne 2 and endo98 ne 2 and endo96=. and (sig9404=1 or col9404=1))
   then endo96=2; /*yes*/
else if (endo96=. and sig9404 ne 1 and col9404 ne 1 and q04 in (1,2,3)/*all04miss ne 1 and q19pt04 ne 1*/)
   then endo96=1; /*no*/
/*1993/Earlier, we just do not know exact year. but can fill no-endo*/
*if (endo92 ne 2 and endo96 ne 2 and endo94=. and (sig9304=1 or col9304=1))
   then endo94=2; /*yes*/                                     /*q19pt04=1 is all path throu*/
if (endo94=. and sig9304 ne 1 and col9304 ne 1 and q04 in (1,2,3)/*all04miss ne 1 and q19pt04 ne 1*/)
   then endo94=1; /*no*/
else if (endo92=. and sig9304 ne 1 and col9304 ne 1 and q04 in (1,2,3)) then endo92=1; /*no*/
else if (endo90=. and sig9304 ne 1 and col9304 ne 1 and q04 in (1,2,3)) then endo90=1; /*no*/
else if (endo88=. and sig9304 ne 1 and col9304 ne 1 and q04 in (1,2,3)) then endo88=1; /*no*/
else if (endo86=. and sig9304 ne 1 and col9304 ne 1 and q04 in (1,2,3)) then endo86=1; /*no*/
else if (endo84=. and sig9304 ne 1 and col9304 ne 1 and q04 in (1,2,3)) then endo84=1; /*no*/
else if (endo82=. and sig9304 ne 1 and col9304 ne 1 and q04 in (1,2,3)) then endo82=1; /*no*/
else if (endo80=. and sig9304 ne 1 and col9304 ne 1 and q04 in (1,2,3)) then endo80=1; /*no*/
else if (endo78=. and sig9304 ne 1 and col9304 ne 1 and q04 in (1,2,3)) then endo78=1; /*no*/
else if (endo76=. and sig9304 ne 1 and col9304 ne 1 and q04 in (1,2,3)) then endo76=1; /*no*/

*****fill in if missing*******;
if endo88=. and .<batch88<800 then endo88=1; /*one answered long qx but missing, we assume no endo*/
if endo90=. and .<batch90<800 then endo90=1; /*one answered long qx but missing, we assume no endo*/
if endo00=. and q00 in (1,2,3,4) then endo00=1; /*one answered long qx but missing, we assume no endo*/
if endo92=. and .<batch92<800 then endo92=1; /*one answered long qx but missing, we assume no endo*/
if endo94=. and .<batch94<800 then endo94=1; /*one answered long qx but missing, we assume no endo*/
if endo96=. and .<batch96<800 then endo96=1; /*one answered long qx but missing, we assume no endo*/
if endo98=. and q98 in (1,2) then endo98=1;
if endo02=. and q02 in (1,2,3) then endo02=1;
if endo04=. and q04 in (1,2,3) then endo04=1;
if endo06=. and q06 in (1,2,3,4) then endo06=1;
if endo08=. and q08 in (1,2,3) then endo08=1;
if endo10=. and q10 in (1,2,3) then endo10=1;

*-----------*
 Use 1990 Qx
*-----------*;
if endo80=. and endo78=. and endo76=. and .<batch90<800 and sigm190 in (9,.) and sigmr90 in (9,.) then do; endo80=1 ; endo78=1 ; endo76=1 ; end;
       /*when she answered long Qx, impute pass thr to no endo*/
if endo82=. and endo84=. and .<batch90<800 and sigm190 in (9,.) and sigmr90 in (9,.) then do; endo82=1; endo84=1; end;
       /*when she answered long Qx, impute pass thr to no endo*/
if endo86=. and .<batch90<800 and sigm190 in (9,.) and sigmr90 in (9,.) then endo86=1;
/*when endo88-endo90 were missing, 1990 info was used to fill the answer.
 We don't impute endo88 if endo90=yes becase person shouldn't take
 endoscopy twice within 4 years. Same for endo90*/
if endo88=. and .<batch90<800 then endo88=1; /*long qx*/
if endo90=. and .<batch90<800 then endo90=1; /*long qx*/

*--------------------------------------*
 reason for endoscopy, make dichotomous
*--------------------------------------*;

   /*1980-1986, derived from 1990 Qx which asked why did you have the first colono or sigmo?*/
   /*reason for endo in 1980-1986 can be used for the reason of first endoscopy*/
   if (vbld90=1 or fecal90=1 or abdpn90=1 or diarr90=1) then reasoned80=1; /*symptom*/
   else if hxclc90=1 or rout90=1 then reasoned80=2;          /*screeting or family hx*/

   if (vbld90=1 or fecal90=1 or abdpn90=1 or diarr90=1) then reasoned82=1; /*symptom*/
   else if hxclc90=1 or rout90=1 then reasoned82=2;          /*screeting or family hx*/

   if (vbld90=1 or fecal90=1 or abdpn90=1 or diarr90=1) then reasoned84=1; /*symptom*/
   else if hxclc90=1 or rout90=1 then reasoned84=2;          /*screeting or family hx*/

   if (vbld90=1 or fecal90=1 or abdpn90=1 or diarr90=1) then reasoned86=1; /*symptom*/
   else if hxclc90=1 or rout90=1 then reasoned86=2;          /*screeting or family hx*/

   /*1988-*/
   if sigm88=2 then reasoned88=1; /*symptom. NOTE, Precise reasons were not asked in 1988 Qx*/
   if sigm88=3 then reasoned88=2; /*screening*/

   if (vbld90=1 or fecal90=1 or abdpn90=1 or diarr90=1) then reasoned90=1; /*symptom*/
   else if hxclc90=1 or rout90=1 then reasoned90=2;                 /*screeting or family hx*/

   if sigm92=2 then reasoned92=2; /*screening*/
   if sigm92=3 then reasoned92=1; /*symptom. NOTE, Precise reasons were not asked in 1992 Qx*/

   if (vbld94=1 or fecal94=1 or abdpn94=1 or diarr94=1) then reasoned94=1; /*symptom*/
   else if hxclc94=1 or rout94=1 then reasoned94=2;                 /*screeting or family hx*/

   if (vbld96=1 or fecal96=1 or abdpn96=1 or diarr96=1) then reasoned96=1; /*symptom*/
   else if hxclc96=1 or rout96=1 then reasoned96=2;                 /*screeting or family hx*/

   if (vbld98=1 or fecal98=1 or abdpn98=1 or diarr98=1 or bariu98=1) then reasoned98=1; /*symptom*/
   else if hxclc98=1 or rout98=1 then reasoned98=2;                 /*screeting or family hx*/

   if (vbld00=1 or fecal00=1 or abdpn00=1 or diarr00=1 or bariu00=1) then reasoned00=1; /*symptom*/
   else if hxclc00=1 or rout00=1 then reasoned00=2;                 /*screeting or family hx*/

   if (vbld02=1 or fecal02=1 or abdpn02=1 or diarr02=1 or bariu02=1) then reasoned02=1; /*symptom*/
   else if hxclc02=1 or rout02=1 then reasoned02=2;                 /*screeting or family hx*/

   if (vbld04=1 or fecal04=1 or abdpn04=1 or diarr04=1 or bariu04=1) then reasoned04=1; /*symptom*/
   else if hxclc04=1 or rout04=1 then reasoned04=2;                 /*screeting or family hx*/

   if (vbld06=1 or fecal06=1 or abdpn06=1 or diarr06=1 or bariu06=1) then reasoned06=1; /*symptom*/
   else if (hxclc06=1 or rout06=1 or cpolp06=1 or virt06=1) then reasoned06=2; /*screeting or family hx*/

   if (vbld08=1 or fecal08=1 or abdpn08=1 or diarr08=1 or bariu08=1) then reasoned08=1; /*symptom*/
   else if (hxclc08=1 or rout08=1 or cpolp08=1 or virt08=1) then reasoned08=2; /*screeting or family hx*/

   if (vbld10=1 or fecal10=1 or abdpn10=1 or diarr10=1 or bariu10=1) then reasoned10=1; /*symptom*/
   else if (hxclc10=1 or rout10=1 or cpolp10=1 or virt10=1) then reasoned10=2; /*screeting or family hx*/

array endoxx [*] endo86 endo88 endo90 endo92 endo94 endo96 endo98 endo00 endo02 endo04 endo06 endo08 endo10   ; /*1=no endo, 2=yes endo*/
array reasonxx[*] reasoned86 reasoned88 reasoned90 reasoned92 reasoned94 reasoned96
                  reasoned98 reasoned00 reasoned02 reasoned04 reasoned06 reasoned08 reasoned10 ;/*At the time of qx 1=symptom,2=screening*/
array endoreason[*] endoreason86 endoreason88 endoreason90 endoreason92 endoreason94 endoreason96
                    endoreason98 endoreason00 endoreason02 endoreason04 endoreason06 endoreason08 endoreason10 ;/*Carry forwarded status 1=symptom,2=screening*/

do i=1 to dim(endoreason);
   if endoxx[i]=2 and reasonxx[i]=. then do; /*missing*/
      do j=i to dim(endoreason);
         endoreason[j]=.;
      end;
   end;
   /*overwrite if status changes*/
   if endoxx[i]=2 and reasonxx[i]=1 then do; /*symptom*/
      do j=i to dim(endoreason);
         endoreason[j]=1;
      end;
   end;
   /*overwrite if status changes*/
   if endoxx[i]=2 and reasonxx[i]=2 then do; /*screening*/
      do j=i to dim(endoreason);
         endoreason[j]=2;
      end;
   end;
end;


/***1. time varying, updated every 2 years, no endo, endo for symptoms, screening endo***/
array endo3c [*]    endo3c86    endo3c88    endo3c90    endo3c92    endo3c94    endo3c96    endo3c98    endo3c00    endo3c02    endo3c04    endo3c06    endo3c08    endo3c10 ;

do i=1 to dim(endo3c);
   if endoxx[i]=1 then endo3c[i]=1; *no endo;
   else if endoxx[i]=2 and reasonxx[i]=1 then endo3c[i]=2; *endo for symptoms;
   else if endoxx[i]=2 and reasonxx[i]=2 then endo3c[i]=3; *endo for screening;
   else endo3c[i]=.; *missing;
end;

/***2. ever endo3c, carry forward*****************************/
array eendo3c {*}   eendo3c86    eendo3c88    eendo3c90    eendo3c92    eendo3c94    eendo3c96    eendo3c98    eendo3c00    eendo3c02    eendo3c04    eendo3c06    eendo3c08    eendo3c10 ;

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


/*** 3. ever screening endo, carry forward
      create screening endo vs no screening endo (including endo for symptoms).
      1=no screening endo, 2=screening endo ****/
array eendo2c {*}  eendo2c86    eendo2c88    eendo2c90    eendo2c92    eendo2c94    eendo2c96    eendo2c98    eendo2c00    eendo2c02    eendo2c04    eendo2c06    eendo2c08    eendo2c10 ;

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
delete    disease familyhis der7620
          n86_nts n90_nts n94_nts  n98_nts  n02_nts  n06_nts n10_nts
          n767880 nur82 nur84 nur86 nur88 nur90 nur92 nur94 nur96 nur98 nur00 nur02 nur04 nur06 nur08 nur10
          ahei2010_8410 meds8020 supp8016 act8614 
          serv86 serv90 serv94 serv98 serv02 serv06 serv10
          ;
RUN;

data nhs1;
     set base end=_end_;

/***PPP, for all cancer file, the suggested cutoff is:
NHS1 June of 2018 (June of 2020 acceptable)
NHS2 June of 2019 (June of 2021 acceptable)
HPFS January of 2018 (January 2020 acceptable)
Please check and modify accordingly for nhs1, nhs2 and hpfs.
*/

/*Siyu, revise the cutoff: NHS1 2020.6, NHS2 2021.6, HPFS 2020.1 */
cutoff=1446;    /*** Set cutoff at 2020.6 ***/
array   irt     {14}    irt86   irt88   irt90   irt92   irt94   irt96   irt98   irt00   irt02   irt04   irt06   irt08   irt10   cutoff;
array   tvar    {13}    t86     t88     t90     t92     t94     t96     t98     t00     t02     t04     t06     t08     t10 ;
array   bmi     {13}    bmi86   bmi88   bmi90   bmi92   bmi94   bmi96   bmi98   bmi00   bmi02   bmi04   bmi06   bmi08   bmi10   ;
array   age     {13}    age86   age88   age90   age92   age94   age96   age98   age00   age02   age04   age06   age08   age10   ;
array   hor     {13}    nhor86  nhor88  nhor90  nhor92  nhor94  nhor96  nhor98  nhor00  nhor02  nhor04  nhor06  nhor08  nhor10  ;
array   phms    {13}    phms86  phms88  phms90  phms92  phms94  phms96  phms98  phms00  phms02  phms04  phms06  phms08  phms10  ;
array   smok    {13}    smkdr86 smkdr88 smkdr90 smkdr92 smkdr94 smkdr96 smkdr98 smkdr00 smkdr02 smkdr04 smkdr06 smkdr08 smkdr10 ;
array   alcon   {13}    alco86n alco86n alco90n alco90n alco94n alco94n alco98n alco98n alco02n alco02n alco06n alco06n alco10n ;
array   alcq    {13}    alcc86  alcc88  alcc90  alcc92  alcc94  alcc96  alcc98  alcc00  alcc02  alcc04  alcc06  alcc08  alcc10  ;
array nSESArr   {13}    nSES_86m nSES_88m nSES_90m nSES_92m nSES_94m nSES_96m nSES_98m nSES_00m nSES_02m nSES_04m nSES_06m nSES_08m nSES_10m ;

array   act     {13}    act86m  act88m  act88m  act92m act94m   act96m   act98m act00m act00m   act04m act04m   act08m act08m  ;
*array  alone   {13}    alone92 alone92 alone92 alone92 alone96 alone96 alone00 alone00 alone04 alone04 alone08 alone08 alone08 ;
*array  marriage{13}    marry96 marry96 marry96 marry96 marry96 marry96 marry00 marry00 marry04 marry04 marry08 marry08 marry08 ;
array   mvit    {13}    mvt86  mvt88  mvt90  mvt92  mvt94  mvt96  mvt98  mvt00  mvt02  mvt04  mvt06  mvt08  mvt10  ;
array   aspi    {13}    asp86 asp88 asp90 asp92 asp94 asp96 asp98 asp00 asp02 asp04 asp06 asp08 asp10 ;
array   hbp     {13}    hbp86  hbp88  hbp90  hbp92  hbp94  hbp96  hbp98  hbp00  hbp02  hbp04  hbp06  hbp08  hbp10;
array   chol    {13}    chol86 chol88 chol90 chol92 chol94 chol96 chol98 chol00 chol02 chol04 chol06 chol08 chol10;
array   diab    {13}    db86    db88    db90    db92     db94    db96    db98    db00    db02    db04    db06    db08    db10;
array   antihp  {13}    htnrx88 htnrx88 htnrx88 htnrx88 htnrx94 htnrx96 htnrx98 htnrx00 htnrx02 htnrx04 htnrx06 htnrx08 htnrx10 ;
array   antitc  {13}    hchtx88 hchtx88 hchtx88 hchtx88 hchtx94 hchtx96 hchtx98 hchtx00 hchtx02 hchtx04 hchtx06 hchtx08 hchtx10 ;
*array  strk    {13}    str86   str88   str90   str92   str94   str96   str98   str00   str02   str04   str06   str08   str10   ;
array   canc    {13}    can86   can88   can90   can92   can94   can96   can98   can00   can02   can04   can06   can08   can10;

array    carbp  {*}   carbo86p    carbo86p     carbo90p    carbo90p    carbo94p    carbo94p     carbo98p   carbo98p      carbo02p   carbo02p     carbo06p   carbo06p   carbo10p  ;
array    hcarbp  {*}   hcarbo86p   hcarbo86p    hcarbo90p   hcarbo90p   hcarbo94p   hcarbo94p    hcarbo98p  hcarbo98p     hcarbo02p  hcarbo02p    hcarbo06p  hcarbo06p  hcarbo10p   ;
array sumref_carbp {*}  sumref_carb86p   sumref_carb86p sumref_carb90p  sumref_carb90p sumref_carb94p  sumref_carb94p sumref_carb98p  sumref_carb98p sumref_carb02p  sumref_carb02p sumref_carb06p  sumref_carb06p sumref_carb10p ;
array     protp  {*}   prot86p     prot86p      prot90p     prot90p     prot94p     prot94p      prot98p     prot98p     prot02p     prot02p    prot06p    prot06p   prot10p ;
array    aprotp  {*}   aprot86p     aprot86p      aprot90p     aprot90p     aprot94p     aprot94p      aprot98p     aprot98p     aprot02p     aprot02p    aprot06p    aprot06p   aprot10p;
array    vprotp  {*}   vprot86p     vprot86p      vprot90p     vprot90p     vprot94p     vprot94p      vprot98p     vprot98p     vprot02p     vprot02p    vprot06p    vprot06p   vprot10p ;
array    tfatp   {*}   tfat86p     tfat86p      tfat90p     tfat90p     tfat94p     tfat94p      tfat98p     tfat98p     tfat02p     tfat02p    tfat06p    tfat06p   tfat10p ;
array    afatp   {*}   afat86p     afat86p      afat90p     afat90p     afat94p     afat94p      afat98p     afat98p     afat02p     afat02p    afat06p    afat06p   afat10p ;
array    vfatp   {*}   vfat86p     vfat86p      vfat90p     vfat90p     vfat94p     vfat94p      vfat98p     vfat98p     vfat02p     vfat02p    vfat06p    vfat06p   vfat10p ;

array eendo2c    {13} eendo2c86 eendo2c88 eendo2c90 eendo2c92 eendo2c94 eendo2c96 eendo2c98 eendo2c00 eendo2c02 eendo2c04 eendo2c06 eendo2c08 eendo2c10 ;
array eendo3c    {13} eendo3c86 eendo3c88    eendo3c90    eendo3c92    eendo3c94    eendo3c96    eendo3c98    eendo3c00 eendo3c02
                      eendo3c04    eendo3c06    eendo3c08    eendo3c10 ;
       
array mama {13} m88 m88 m90 m92 m94 m96 m98 m00 m02 m04 m06 m08 m10  ;


do j=2 to DIM(bmi);
      if age{j}=.      then age{j}=age{j-1}+2;
      if bmi{j}=.      then bmi{j}=bmi{j-1};
      if hor{j}=.      then hor{j}=hor{j-1}; /* New PMH use, $label 1.Pre/Missing Mnp; 2.Never Used; 3.Current User; 4.Past User; 5.Status Unknown; 6.Missing */
      if smok{j}=.     then smok{j}=smok{j-1}; /* smoking status range 0-15 0.missing 1.never 2-8.past 9-15.curr*/
      if alcon{j}=.    then alcon{j}=alcon{j-1};
      if act{j}=.      then act{j}=act{j-1};
      if mvit{j}=.     then mvit{j}=mvit{j-1};
      if aspi{j}=.     then aspi{j}=aspi{j-1};
      *if marriage{j}=. then marriage{j}=marriage{j-1};
      *if alone{j}=.    then alone{j}=alone{j-1};
      if antihp{j}=.   then antihp{j}=antihp{j-1};
      if antitc{j}=.   then antitc{j}=antitc{j-1};
      if canc{j}=.     then canc{j}=canc{j-1};
      if nSESArr{j}=.  then nSESArr{j}=nSESArr{j-1};
end;

/* -single  s-stop updated  v-average  sv-stop updated average*/
/* -s */
array calors {13}  calor86ns calor86ns calor90ns calor90ns calor94ns calor94ns calor98ns calor98ns calor02ns calor02ns calor06ns calor06ns calor10ns ;
array aheims {13}  ahei86s ahei86s ahei90s ahei90s ahei94s ahei94s ahei98s ahei98s ahei02s ahei02s ahei06s ahei06s ahei10s ;   
array lcss   {13}  lcs86s    lcs86s    lcs90s    lcs90s    lcs94s    lcs94s    lcs98s    lcs98s    lcs02s    lcs02s    lcs06s    lcs06s    lcs10s;
array alcss  {13}  alcs86s   alcs86s   alcs90s   alcs90s   alcs94s   alcs94s   alcs98s   alcs98s   alcs02s   alcs02s   alcs06s   alcs06s   alcs10s;
array vlcss  {13}  vlcs86s   vlcs86s   vlcs90s   vlcs90s   vlcs94s   vlcs94s   vlcs98s   vlcs98s   vlcs02s   vlcs02s   vlcs06s   vlcs06s   vlcs10s;
array ulcss  {13}  ulcs86s   ulcs86s   ulcs90s   ulcs90s   ulcs94s   ulcs94s   ulcs98s   ulcs98s   ulcs02s   ulcs02s   ulcs06s   ulcs06s   ulcs10s;
array hlcss  {13}  hlcs86s   hlcs86s   hlcs90s   hlcs90s   hlcs94s   hlcs94s   hlcs98s   hlcs98s   hlcs02s   hlcs02s   hlcs06s   hlcs06s   hlcs10s;
array lfss   {13}  lfs86s    lfs86s    lfs90s    lfs90s    lfs94s    lfs94s    lfs98s    lfs98s    lfs02s    lfs02s    lfs06s    lfs06s    lfs10s;
array ulfss  {13}  ulfs86s   ulfs86s   ulfs90s   ulfs90s   ulfs94s   ulfs94s   ulfs98s   ulfs98s   ulfs02s   ulfs02s   ulfs06s   ulfs06s   ulfs10s;
array hlfss  {13}  hlfs86s   hlfs86s   hlfs90s   hlfs90s   hlfs94s   hlfs94s   hlfs98s   hlfs98s   hlfs02s   hlfs02s   hlfs06s   hlfs06s   hlfs10s;

array alfss {13} alfs86s alfs86s alfs90s alfs90s alfs94s alfs94s alfs98s alfs98s alfs02s alfs02s alfs06s alfs06s alfs10s ;
array vlfss {13} vlfs86s vlfs86s vlfs90s vlfs90s vlfs94s vlfs94s vlfs98s vlfs98s vlfs02s vlfs02s vlfs06s vlfs06s vlfs10s ;

/* -sv */
array calorsv {13}  calor86nsv calor86nsv calor90nsv calor90nsv calor94nsv calor94nsv calor98nsv calor98nsv calor02nsv calor02nsv calor06nsv calor06nsv calor10nsv ;
array aheimsv {13}  ahei86sv ahei86sv ahei90sv ahei90sv ahei94sv ahei94sv ahei98sv ahei98sv ahei02sv ahei02sv ahei06sv ahei06sv ahei10sv ;   
array lcssv   {13}  lcs86sv    lcs86sv    lcs90sv    lcs90sv    lcs94sv    lcs94sv    lcs98sv    lcs98sv    lcs02sv    lcs02sv    lcs06sv    lcs06sv    lcs10sv;
array alcssv  {13}  alcs86sv   alcs86sv   alcs90sv   alcs90sv   alcs94sv   alcs94sv   alcs98sv   alcs98sv   alcs02sv   alcs02sv   alcs06sv   alcs06sv   alcs10sv;
array vlcssv  {13}  vlcs86sv   vlcs86sv   vlcs90sv   vlcs90sv   vlcs94sv   vlcs94sv   vlcs98sv   vlcs98sv   vlcs02sv   vlcs02sv   vlcs06sv   vlcs06sv   vlcs10sv;
array ulcssv  {13}  ulcs86sv   ulcs86sv   ulcs90sv   ulcs90sv   ulcs94sv   ulcs94sv   ulcs98sv   ulcs98sv   ulcs02sv   ulcs02sv   ulcs06sv   ulcs06sv   ulcs10sv;
array hlcssv  {13}  hlcs86sv   hlcs86sv   hlcs90sv   hlcs90sv   hlcs94sv   hlcs94sv   hlcs98sv   hlcs98sv   hlcs02sv   hlcs02sv   hlcs06sv   hlcs06sv   hlcs10sv;
array lfssv   {13}  lfs86sv    lfs86sv    lfs90sv    lfs90sv    lfs94sv    lfs94sv    lfs98sv    lfs98sv    lfs02sv    lfs02sv    lfs06sv    lfs06sv    lfs10sv;
array ulfssv  {13}  ulfs86sv   ulfs86sv   ulfs90sv   ulfs90sv   ulfs94sv   ulfs94sv   ulfs98sv   ulfs98sv   ulfs02sv   ulfs02sv   ulfs06sv   ulfs06sv   ulfs10sv;
array hlfssv  {13}  hlfs86sv   hlfs86sv   hlfs90sv   hlfs90sv   hlfs94sv   hlfs94sv   hlfs98sv   hlfs98sv   hlfs02sv   hlfs02sv   hlfs06sv   hlfs06sv   hlfs10sv;
array alfssv {13} alfs86sv alfs86sv alfs90sv alfs90sv alfs94sv alfs94sv alfs98sv alfs98sv alfs02sv alfs02sv alfs06sv alfs06sv alfs10sv ;
array vlfssv {13} vlfs86sv vlfs86sv vlfs90sv vlfs90sv vlfs94sv vlfs94sv vlfs98sv vlfs98sv vlfs02sv vlfs02sv vlfs06sv vlfs06sv vlfs10sv ;


array whgrn_sv {13}  whgrn_s86sv whgrn_s86sv whgrn_s90sv whgrn_s90sv whgrn_s94sv whgrn_s94sv whgrn_s98sv whgrn_s98sv whgrn_s02sv whgrn_s02sv whgrn_s06sv whgrn_s06sv whgrn_s10sv ;
array rgrain_sv {13}  rgrain_s86sv rgrain_s86sv rgrain_s90sv rgrain_s90sv rgrain_s94sv rgrain_s94sv rgrain_s98sv rgrain_s98sv rgrain_s02sv rgrain_s02sv rgrain_s06sv rgrain_s06sv rgrain_s10sv ;   
array frui_sv   {13}  frui_s86sv    frui_s86sv    frui_s90sv    frui_s90sv    frui_s94sv    frui_s94sv    frui_s98sv    frui_s98sv    frui_s02sv    frui_s02sv    frui_s06sv    frui_s06sv    frui_s10sv;
array leg_sv  {13}  leg_s86sv   leg_s86sv   leg_s90sv   leg_s90sv   leg_s94sv   leg_s94sv   leg_s98sv   leg_s98sv   leg_s02sv   leg_s02sv   leg_s06sv   leg_s06sv   leg_s10sv;
array pota_sv  {13}  pota_s86sv   pota_s86sv   pota_s90sv   pota_s90sv   pota_s94sv   pota_s94sv   pota_s98sv   pota_s98sv   pota_s02sv   pota_s02sv   pota_s06sv   pota_s06sv   pota_s10sv;
array veg_noleg_sv  {13}  veg_noleg_s86sv   veg_noleg_s86sv   veg_noleg_s90sv   veg_noleg_s90sv   veg_noleg_s94sv   veg_noleg_s94sv   veg_noleg_s98sv   veg_noleg_s98sv   veg_noleg_s02sv   veg_noleg_s02sv   veg_noleg_s06sv   veg_noleg_s06sv   veg_noleg_s10sv;
array fish_sv  {13}  fish_s86sv   fish_s86sv   fish_s90sv   fish_s90sv   fish_s94sv   fish_s94sv   fish_s98sv   fish_s98sv   fish_s02sv   fish_s02sv   fish_s06sv   fish_s06sv   fish_s10sv;
array rpmeats_sv   {13}  rpmeats_s86sv    rpmeats_s86sv    rpmeats_s90sv    rpmeats_s90sv    rpmeats_s94sv    rpmeats_s94sv    rpmeats_s98sv    rpmeats_s98sv    rpmeats_s02sv    rpmeats_s02sv    rpmeats_s06sv    rpmeats_s06sv    rpmeats_s10sv;
array ssb_sv  {13}  ssb_s86sv   ssb_s86sv   ssb_s90sv   ssb_s90sv   ssb_s94sv   ssb_s94sv   ssb_s98sv   ssb_s98sv   ssb_s02sv   ssb_s02sv   ssb_s06sv   ssb_s06sv   ssb_s10sv;
array dairy_sv  {13}  dairy_s86sv   dairy_s86sv   dairy_s90sv   dairy_s90sv   dairy_s94sv   dairy_s94sv   dairy_s98sv   dairy_s98sv   dairy_s02sv   dairy_s02sv   dairy_s06sv   dairy_s06sv   dairy_s10sv;
array swtsnochoc_sv {13} swtsnochoc_s86sv swtsnochoc_s86sv swtsnochoc_s90sv swtsnochoc_s90sv swtsnochoc_s94sv swtsnochoc_s94sv swtsnochoc_s98sv swtsnochoc_s98sv swtsnochoc_s02sv swtsnochoc_s02sv swtsnochoc_s06sv swtsnochoc_s06sv swtsnochoc_s10sv ;
array poul_sv {13} poul_s86sv poul_s86sv poul_s90sv poul_s90sv poul_s94sv poul_s94sv poul_s98sv poul_s98sv poul_s02sv poul_s02sv poul_s06sv poul_s06sv poul_s10sv ;

array carbsv {13}  carbo86sv carbo86sv carbo90sv carbo90sv carbo94sv carbo94sv carbo98sv carbo98sv carbo02sv carbo02sv carbo06sv carbo06sv carbo10sv  ;
array hcarbsv {13}  hcarbo86sv hcarbo86sv hcarbo90sv hcarbo90sv hcarbo94sv hcarbo94sv hcarbo98sv hcarbo98sv hcarbo02sv hcarbo02sv hcarbo06sv hcarbo06sv hcarbo10sv  ;   
array sumref_carbsv   {13}  sumref_carb86sv    sumref_carb86sv    sumref_carb90sv    sumref_carb90sv    sumref_carb94sv    sumref_carb94sv    sumref_carb98sv    sumref_carb98sv    sumref_carb02sv    sumref_carb02sv    sumref_carb06sv    sumref_carb06sv    sumref_carb10sv  ;
array protsv  {13}  prot86sv   prot86sv   prot90sv   prot90sv   prot94sv   prot94sv   prot98sv   prot98sv   prot02sv   prot02sv   prot06sv   prot06sv   prot10sv  ;
array aprotsv  {13}  aprot86sv   aprot86sv   aprot90sv   aprot90sv   aprot94sv   aprot94sv   aprot98sv   aprot98sv   aprot02sv   aprot02sv   aprot06sv   aprot06sv   aprot10sv  ;
array vprotsv  {13}  vprot86sv   vprot86sv   vprot90sv   vprot90sv   vprot94sv   vprot94sv   vprot98sv   vprot98sv   vprot02sv   vprot02sv   vprot06sv   vprot06sv   vprot10sv  ;
array tfatsv  {13}  tfat86sv   tfat86sv   tfat90sv   tfat90sv   tfat94sv   tfat94sv   tfat98sv   tfat98sv   tfat02sv   tfat02sv   tfat06sv   tfat06sv   tfat10sv  ;
array afatsv {13} afat86sv afat86sv afat90sv afat90sv afat94sv afat94sv afat98sv afat98sv afat02sv afat02sv afat06sv afat06sv afat10sv  ;
array vfatsv {13} vfat86sv vfat86sv vfat90sv vfat90sv vfat94sv vfat94sv vfat98sv vfat98sv vfat02sv vfat02sv vfat06sv vfat06sv vfat10sv  ;

**************************************************************************
***** If an irt date is before June of that qq year or after or equal ****
***** to the next qq year it is incorrect and should be defaulted to  ****
***** June of that qq year.    Make time period indicator tvar=0.       ****
**************************************************************************;

/* from 1984 */
/* do i=1 to DIM(irt)-1;
  if (irt{i}<(990+24*i) | irt{i}>=(1014+24*i)) then irt{i}=990+24*i;
  age{i}=int((irt{i}-bdt)/12);
end; */

/* From 1986: create ltf using the original irt values */
ltf=.;
do i=1 to dim(irt)-1; if irt{i}>0 then ltf=irt{i}; end;
if ltf=irt{dim(irt)-1} then ltf=.;

/* Correct irt values */
do i=1 to dim(irt)-1;
    if (irt{i} <  (1014+24*i)) or (irt{i} >= (1038+24*i)) then irt{i}=1014+24*i;
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
    dead_can   = 0; tdead_can=irt{i+1}-irt{i}; dtdx_dead_can=.;
    if irt{i}< dtdth <=irt{i+1} and (140<=newicda3<=208) then do; dead_can=1; tdead_can=dtdth-irt{i}; dtdx_dead_can=dtdth; end;

    dead_cvd   = 0; tdead_cvd=irt{i+1}-irt{i}; dtdx_dead_cvd=.;
    if irt{i}< dtdth <=irt{i+1} and (390<=newicda3<=458) then do; dead_cvd=1; tdead_cvd=dtdth-irt{i}; dtdx_dead_cvd=dtdth; end;

    dead_resp=0;  tdead_resp=irt{i+1}-irt{i}; dtdx_dead_resp=.;
    if irt{i}<dtdth<=irt{i+1} and (460<=newicda3<520) then do; dead_resp=1; tdead_resp=dtdth-irt{i}; dtdx_dead_resp=dtdth; end;

    dead_oth=0;tdead_oth=irt{i+1}-irt{i}; dtdx_dead_oth=.;
    if irt{i}<dtdth<=irt{i+1} and (dead_can=0 and dead_cvd=0 and dead_resp=0 ) then do; dead_oth=1; tdead_oth=dtdth-irt{i}; dtdx_dead_oth=dtdth; end;


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

/* race */
%indic3(vbl=eth3g, prefix=eth3g, reflev=1, min=2, max=3, missing=., usemiss=0,
           label1='causcasian',
           label2='african',
           label3='asian & others');

/* AGE GROUP */
     agecon=age{i};
     agegp=int((age{i}-30)/5);     **** Define the agegp in the i-th period;
     if agegp<=0 then agegp=1;
     if agegp>8 then agegp=8;

     %indic3(vbl=agegp, prefix=agegp, reflev=1, min=2, max=8, missing=., usemiss=0,
                     label1='age<40',
                     label2='age40-44',
                     label3='age45-49',
                     label4='age50-54',
                     label5='age55-59',
                     label6='age60-64',
                     label7='age65-69',
                     label8='age70+');
     if agecon>65 then agesub=1; else agesub=0;

/*** phms ***/
         if hor{i} in (2,5) then phms{i}=2;
         else if hor{i} =3  then phms{i}=3;
         else if hor{i} =4  then phms{i}=4;
         else phms{i}=1;
        
         phmsstatus=phms{i};
         %indic3(vbl=phmsstatus, prefix=phmsstatus, reflev=1, missing=., usemiss=0, min=2, max=4,
                     label1='pre mnp',
                     label2='never/unknown pmhuser, post mnp',
                     label3='curr pmh user and post mnp',
                     label4='past pmhuser and post mnp');

if phmsstatus=1 then postpau=0; else postpau=1;
if phmsstatus in (1,2) then pmhever=0; else pmhever=1; 

/** smoking status **/
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
       
       if 0 <= actcon and actcon < 3  then actc=1;
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

                                                     alcq{i}=1;
            if alcon{i}=0.0                     then alcq{i}=1;
            if alcon{i}>0.0  and alcon{i}<5.0   then alcq{i}=2;
            if alcon{i}>=5.0 and alcon{i}<10.0  then alcq{i}=3;
            if alcon{i}>=10.0 and alcon{i}<15.0 then alcq{i}=4;
            if alcon{i}>=15.0 and alcon{i}<30.0 then alcq{i}=5;
            if alcon{i}>=30.0                   then alcq{i}=6;

       alc=alcq{i};
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

         if      25<= bmicon and bmicon < 30  then bmic3cat=2;
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

         if      25=< bmibcon <30 then bmib3cat=2;
         else if 30=< bmibcon     then bmib3cat=3;
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



/*** multiple vitamin ***/
             select(mvit{i});
                when (1) mv=1;
                otherwise  mv=0;
            end;

/*** aspirin ***/
            select(aspi{i});
                when (1)     asp=1;
                otherwise    asp=0;
            end;

/***family history ***/
       if mifh ne 1 then mifh=0;
       if dbfh ne 1 then dbfh=0;
       if cafh ne 1 then cafh=0;
       
/****** Indicator for History of antihypertensive user *******/
            select(antihp{i});
                when (1)     antihbp=1;
                otherwise    antihbp=0;
            end;

            if hbp{i}=1  then highbp=1; else highbp=0;
            if hbp{1}=1 then basehighbp=1; else basehighbp=0;
     

/****** Indicator for TC lowering *******/
            select(antitc{i});
                when (1)     antihtc=1;
                otherwise    antihtc=0;
            end;

           if chol{i}=1 then hightc=1; else hightc=0;
           if chol{1}=1 then basehightc=1; else basehightc=0;


/****** Indicator for History of diabetes *******/
           if diab{i}=1 then diabetes=1; else diabetes=0;

/****** endoscopy for screening purpose ******/
    endo2c=eendo2c{i}; /* ever screening endo, 1=no screening endo, 2=screening endo */
    endo3c=eendo3c{i}; /* ever endo3c, 1=no endo, 2=endo for symptoms, 3=screening endo */
    if endo2c = . then endo2c = 1;
/****** Mammogram in last 2 yrs *******/
   mam=mama{i};                     
   if mam = . then mam = 0;

/** Indicator for nSES **/
nSES = nSESArr{i};

/* Baseline food */
basewhgrn=whgrn_s86;
basergrain=rgrain_s86;
basefrui=frui_s86;
baseleg=leg_s86;
basepota=pota_s86;
baseveg_noleg=veg_noleg_s86;
basefish=fish_s86;
baserpmeats=rpmeats_s86;
basessb=ssb_s86;
basedairy=dairy_s86;
basesweets=swtsnochoc_s86;
basepoultry=poul_s86;

/****************  BASELINE EXCLUSIONS ********************************************/

if i=1 then do;
   
   %exclude(exrec eq 1);                 *multiple records and not in master file;
   %exclude(baseffq eq 1);  
   %exclude(lcs86 eq .);
   %exclude(alcs86 eq .);
   %exclude(vlcs86 eq .);
   %exclude(ulcs86 eq .);
   %exclude(hlcs86 eq .);
   %exclude(lfs86 eq .);
   %exclude(ulfs86 eq .);
   %exclude(hlfs86 eq .);

   %exclude(ltf eq irt{1});     * Only completed baseline questionnaire;

   %exclude(0 lt dtdth   le irt{i});
   %exclude(0 lt dtdx_allca  le irt{i}); 
   %exclude(can86   eq 1);
   %exclude(0 lt dtdxcvd le irt{i}); 
   %exclude(mi86    eq 1);
   %exclude(str86   eq 1);   
   %exclude(cabg86  eq 1); 

   %exclude(0 lt dtdxdb2 le irt{i});
   %exclude(db86    eq 1); 

   %exclude(yobf le 20);
   %exclude(age86 eq .);

   %exclude(calor86n eq .);
   %exclude(calor86n < 500);
   %exclude(calor86n > 3500);
   %exclude(ahei86 eq .);
   
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

proc sort data=nhs1;by period;run;

proc rank data=nhs1 out=nhs1 group=5;
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

data nhs1;
set  nhs1 end=_end_;

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
     irt86 irt88 irt90 irt92 irt94 irt96 irt98 irt00 irt02 irt04 irt06 irt08 irt10 cutoff
     t86 t88 t90 t92 t94 t96 t98 t00 t02 t04 t06 t08 t10
     dtdxchd dtdxmi dtdth dtdxca  dtdx_allca dtdx_obeca  dtdx_nonobeca  fup_allca
     /** Outcome **/
     dead tdead 
     dead_can dead_cvd dead_resp dead_oth tdead_can tdead_cvd tdead_resp tdead_oth 
     can_obese tobeca can_nonobese tnonobeca 
     cvd tcvd chd tchd str tstr  
     /** Site-specific cancer **/
     can_: t: dtdx_:      obeca_site can_obese_site         
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
     age86   age88   age90   age92   age94   age96   age98   age00   age02   age04   age06   age08   age10
     white eth3g dbfh mifh cafh basehighbp basehightc highbp hightc diabetes mv asp 
     agecon agegp  agesub  bmicon alcocon actcon 
     
     bmic3cat &bmic3cat_ bmib3cat &bmib3cat_ bmic25 bmic30  bmic9 &bmic9_ bmib9 &bmib9_ 
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
proc means data=nhs1 mean std median n nmiss min max;
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