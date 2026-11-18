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
/** Death **/
%hp_dead (keep= mmdth yydth icda newicda dtdth);

   if 0<=yydth<50  then yydth=yydth+100;

   newicda=compress(icda, 'E');
   if newicda='V812' then newicda=812;

   if mmdth<=0 and yydth>0 then mmdth=6;
   dtdth=(yydth*12)+mmdth;

   if dtdth=. then delete;
run;

/** Cancer **/
data cancer;
infile 'Cancer2020_March2024_withduplicate.dat' lrecl=200;         
        input
        hpfsid 1-6
        dateca 8-11  /* date of cancer */
        dod 12-15  /* date of death for all deaths */
        cancer   16 /* all cancers 1,0 */
        mel    18 /* melanoma */
        lal    19 /* lymphoma + leukemia */
        other   20 /* other cancers */
        colorect 21  /* colorectal */
        lung 22  /* */
        carcin  24  /* carcinoma = cancers without lal,brain,sarcoma */
        newag 25  /* prostate cancer aggressive only */
        cancerp 26 /* total cancers including ONLY aggressive prostate ca */
        carcinp 27 /* total carconoma including ONLY aggressive prostate ca */
        prosnoa1 29 /* prostate no A1 stage*/
        acan 30  /* icd code is: 141,143,144,145,146,148,149,150,153,154 */
        bladder 31
        pancreas 32
        kidney 33
        stomach 34 /*icd code 151*/
        liver 35
        brain 36
        sarcoma 37
        oralc 38 /*mouth only*/
        status 39  /* prostate, no A1 stage, from prostate cancer file */
        esoph 40
        leuk 41
        lymp_nh 42
        lymp_h 43
        myeloma 44
        pharyn 45 /* pharynx separate */
        oroph  46
        conf 47-48   /* confirmation code */
        cancerdth 49
        icdx      52-54   /* first report cancer ICD*/
        ;

        id=hpfsid; 
        if conf in (11,12,13,14,15,16,17,18,19,25) then ca=1;  /* confirmed cancer */

        if ca=1 and icdx ne . then caflag=1;

        dtdx_allca=.;
        if cancerp=1 and dateca>0 and ca=1 and (140<=icdx<=172 or 174<=icdx<=207 or icdx =225) then dtdx_allca=dateca;

        if dateca>0 and ca=1 and icdx ne 173.2 and icdx ne 173.1 then dtdxca=dateca;


    /* ---- Site groups ---- */
      /* Digestive */
      if 140 <= icdx <= 149 and ca=1 then do; oral_ca=1;      dtdx_oral      = dtdxca; end;
      if icdx=150 and ca=1 then do; esoph_ca=1;               dtdx_esoph     = dtdxca; end;
      if icdx=151 and ca=1 then do; stomach_ca=1;             dtdx_stomach   = dtdxca; end;
      if icdx=152 and ca=1 then do; smallintest_ca=1;         dtdx_smallintest=dtdxca; end;
      if icdx in (153,154) and ca=1 then do; colorect_ca=1;   dtdx_colorect  = dtdxca; end;
      if icdx=155 and ca=1 then do; liver_ca=1;               dtdx_liver     = dtdxca; end;
      if icdx=156 and ca=1 then do; gallbladder_ca=1;         dtdx_gallbladder=dtdxca; end;
      if icdx=157 and ca=1 then do; panc_ca=1;                dtdx_panc      = dtdxca; end;
      if 140 <= icdx <= 159 and ca=1 then do; diges_ca=1;     dtdx_diges     = dtdxca; end;

      /* Respiratory */
      if icdx=160 and ca=1 then do; nose_ca=1;                dtdx_nose      = dtdxca; end;
      if icdx=161 and ca=1 then do; larynx_ca=1;              dtdx_larynx    = dtdxca; end;
      if icdx=162 and ca=1 then do; lung_ca=1;                dtdx_lung      = dtdxca; end;
      if 160 <= icdx <= 165 and ca=1 then do; resp_ca=1;      dtdx_resp      = dtdxca; end;

      /* Skin / other selected */
      if icdx=170 and ca=1 then do; bone_ca=1;                dtdx_bone      = dtdxca; end;
      if icdx=172 and ca=1 then do; melanoma_ca=1;            dtdx_melanoma  = dtdxca; end;
      if icdx=174 and ca=1 then do; brca_ca=1;                dtdx_brca      = dtdxca; end;

      /* Genito-urinary Cancers */
      if icdx=180 and ca=1 then do; cervical=1;               dtdx_cervical  = dtdxca; end;
      if icdx=182 and ca=1 then do; endo_ca=1;                dtdx_endo      = dtdxca; end;
      if icdx=183 and ca=1 then do; ovca_ca=1;                dtdx_ovca      = dtdxca; end;
      if icdx=188 and ca=1 then do; bladder_ca=1;             dtdx_bladder   = dtdxca; end;
      if icdx=189 and ca=1 then do; kidney_ca=1;              dtdx_kidney    = dtdxca; end;
      if icdx in (185,186,187,188,189) and ca=1 then do; uri_ca=1; dtdx_uri  = dtdxca; end;

      /* Other sites */
      if icdx=191 and ca=1 then do; brain_ca=1;               dtdx_brain     = dtdxca; end;
      if icdx=193 and ca=1 then do; thyroid_ca=1;             dtdx_thyroid   = dtdxca; end;

      /* Lymphatic / blood */
      if icdx=203 and ca=1 then do; myeloma_ca=1;             dtdx_myeloma   = dtdxca; end;
      if icdx=205 and ca=1 then do; myeleuk_ca=1;             dtdx_myeleuk   = dtdxca; end;

      /* Meningioma */
      if icdx=225 and ca=1 then do; menin_ca=1;               dtdx_menin     = dtdxca; end;

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
     if obese_ca=1 and dtdxca>0 and ca=1 and icdx in (150,151,153,154,155,156,157,174,182,183,189,193,203,225) then dtdx_obeca=dtdxca;
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
     if obeca_nobc=1 and dtdxca>0 and ca=1 and icdx in (150,151, /*153,154,*/ 155,156,157, /*174,*/ 182,183,189,193,203,225) then dtdx_obeca_nobc=dtdxca;
     
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
     if obeca_nobcop=1 and dtdxca>0 and ca=1 and icdx in (150,151, /*153,154,*/ 155,156, /*157,*/ /*174,*/ /*182, 183,*/ 189,193,203,225) then dtdx_obeca_nobcop=dtdxca;

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
     if obese_ca=1 and dtdxca>0 and ca=1 and icdx in (150,  151, 155,156,  203 ) then dtdx_othsite=dtdxca;
    
/*==============================================================*/
/* Obesity-related cancer subtype: 13 cancer sites              */
/*==============================================================*/

obeca_site = .;

if ca=1 and dtdxca>0 then do;
    select (icdx);
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



    /** CVD **/
filename cardio16 'cardio2016.02062020';   /* most updated */               

data cvd16;
infile cardio16 lrecl=700;
        input
        id 1-6
        mi_nf 7  /*0=noncase, 1=definite, 2=probable, 3=selfreport, 4=rejected*/
        mnyr_nmi 8-11
        fat_chd 27 /* fatal mi; 0=noncase, 1=def, 2=presumed, 3=sudden*/
        mnyr_fat 28-31
        cabg 12 /*same as nf mi*/
        mnyr_ca 13-16
        angina 17  /*same*/
        mnyr_an 18-21
        tot_chd 32 /*nf, fatal mi + cabg*/
        mnyr_chd 33-36
        tot_cvd 37 /* 0: no, 1: yes */
        mnyr_cvd 38-41
        nf_str 22 /*0=noncase, 1=definite, 2=probable, 3=selfreport, 4=rejected**/
        fat_str 42 /* fatal mi; 0=noncase, 1=def, 2=presumed, 3=sudden*/
        tot_str 43 /*0: no 1: yes */
        mnyr_nfst 23-26
        mnyr_tst 47-50 /** p-t total stroke **/
        s1 $ 44
         
         /*
         Stroke Subcode 1, type of stroke
         $label 1.430: subarachnoid hemorrhage;\
         2.icda 431:intracerebral hemorrhage;\
         3.432-4,6:thrombotic stroke;\
         4.432-4,6:embolic stroke;\
         5.432-4,6:thrombotic/embolic stroke;\
         6.436-8:-unknown type;\
         */
        ;

    /** MI **/
        nfmi=0;
        if mi_nf in (1,2) then nfmi=1;
        dtdxnfmi=.;
        if mnyr_nmi > 0 and nfmi=1 then dtdxnfmi=mnyr_nmi;

        fmi=0;
        if fat_chd in (1,2,3) then fmi=1;
        dtdxfmi=.;
        if mnyr_fat > 0 and nfmi=1 then dtdxfmi=mnyr_fat;

        tot_mi = sum (nfmi,fmi);
        totmi=0;
        if tot_mi in (1,2) then totmi=1;

        dtdxmi=.;
        if 0<dtdxnfmi<9999 or 0<dtdxfmi<9999 then dtdxmi=min(dtdxnfmi,dtdxfmi,.);

        dtdxcvd=.;
        if 0<mnyr_cvd<9999 and tot_cvd=1 then dtdxcvd=mnyr_cvd;

        dtdxchd=.;
        if 0<mnyr_chd<9999 and tot_chd=1 then dtdxchd=mnyr_chd;

        dtdxstr=.;
        if 0<mnyr_tst<9999 and tot_str=1 then dtdxstr=mnyr_tst;

    /** stroke **/
        if s1 in (1,2) then hem_str=1; else hem_str=0;
        if s1 in (3,4,5) then isc_str=1; else isc_str=0;
        if s1=6 then unk_str=1; else unk_str=0;


        nfstr=0;
        if nf_str in (1,2) then nfstr=1;

        fatstr=0;
        if fat_str in (1,2) then fatstr=1;

        keep id dtdxchd dtdxcvd dtdxstr dtdxmi;
        run;

    /** diabetes **/
data diabetes;
infile 'diab8622.05MAY23';
        input
        @1 id        6.
        @7 cd        2.
        @10 yrdx     2.
        @13 dbconf   2.
        @16 diag_mn  2.
        @19 diag_yr  2.
        @30 type     1.
        @32 prob     1.
        @34 dbcase   1.
        @36 dtdx     4.
        @43 fu     $ 4.
        ;

        dtdxdb=.;
        if 0< dtdx <9999 and dbcase=1 then dtdxdb=dtdx;          

        type2db=0;          
        if type=2 and prob=1 then type2db=1; *type 2 diabees;          

        type1db=0;          
        if type eq 1 or type eq 5 then type1db=1;          

        dtdxdb2=.;          
        if 0< dtdx <9999 and type2db=1 then dtdxdb2=dtdx;  
          
        keep  id  dtdxdb2; /* dtdxdb2 is the time variable */          
        run;                                        



*----------------2.Exposure: Low-Carb Score-----------------------------;

/*** nutrient 86, 90, 94, 98, 02, 06, 10, 14, 18 ***/
/*  nts86 */
  %h86_nts (keep=alco86n calor86n na86n
  carbo86n/*Carbohydrate gm*/
  pot_carb86n/*Carbohydrate from potato*/
  sumref_carb86n/********Sum of carbohydrate from potato, added sugar in foods, and refined grain, foods only********/
  sumrf_nfr86n/*Sum of added sugar minus fruit carb, plus carb from refined grain, cereals only*/ 
  adds_cer86n/********Added sugar minus cereal fruit for Carb Quality Variables only********/
  adds_fds86n/********Added Sugar From Foods for Carb Quality Variables only********/
  ref_grn86n/*Carbohydrate from refined grain*/

  leg_carb86n/*Carbohdyrate from legumes*/
  frtcb_njs86n/*Fruit carbohydrate minus juice sugar*/
  vegcb_npl86n/*Vegetable Carbohydrate minus potato and legume carbohydrate*/
  sum_whgcarb86n/********Sum Carbohydrate From Wholegrain********/
  sumoth_crb_njs86n/********Sum of carbohydrate from fruit minus juice sugar, plus carb from vegetable, legume, and lact********/
  mill_grn86n/*Carbohydrate from milled whole grain*/
  int_grn86n/*Carbohydrate from intact whole grain*/

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

data h86_nts;
set h86_nts;

array gram_cp carbo86n 
              pot_carb86n sumref_carb86n sumrf_nfr86n adds_cer86n adds_fds86n ref_grn86n 
              leg_carb86n frtcb_njs86n vegcb_npl86n sum_whgcarb86n sumoth_crb_njs86n mill_grn86n int_grn86n 
              prot86n aprot86n vprot86n dprot86n ;
array perc_cp carbo86p 
              pot_carb86p sumref_carb86p sumrf_nfr86p adds_cer86p adds_fds86p ref_grn86p 
              leg_carb86p frtcb_njs86p vegcb_npl86p sum_whgcarb86p sumoth_crb_njs86p mill_grn86p int_grn86p  
              prot86p aprot86p vprot86p dprot86p ;

array gram_f  tfat86n afat86n vfat86n dfat86n sat86n mon86n mufa_pl86n mufa_a86n poly86n omega86n trnss86n;
array perc_f  tfat86p afat86p vfat86p dfat86p sat86p mon86p mufa_pl86p mufa_a86p poly86p omega86p trnss86p;

do over gram_cp;
perc_cp=(4*gram_cp)/calor86n*100;
end;

do over gram_f;
perc_f=(9*gram_f)/calor86n*100;
end;

hfat86p=sum(mufa_pl86p,poly86p);
ufat86p=sum(sat86p,mufa_a86p,trnss86p);

hcarbo86p=sum(leg_carb86p,frtcb_njs86p,vegcb_npl86p,sum_whgcarb86p);
run;

proc rank data=h86_nts groups=11 out=h86_nts;
var    prot86p  aprot86p  vprot86p  tfat86p  afat86p  vfat86p  hfat86p  ufat86p;
ranks  prot86pd aprot86pd vprot86pd tfat86pd afat86pd vfat86pd hfat86pd ufat86pd;
run;
 
proc rank data=h86_nts descending groups=11 out=h86_nts;
var    carbo86p  sumref_carb86p  hcarbo86p;
ranks  carbo86pd sumref_carb86pd hcarbo86pd;
run;

data h86_nts;
set h86_nts;

lcs86=sum(carbo86pd,prot86pd,tfat86pd);
alcs86=sum(carbo86pd,aprot86pd,afat86pd);
vlcs86=sum(carbo86pd,vprot86pd,vfat86pd);

hlcs86=sum(sumref_carb86pd,vprot86pd,vfat86pd);
ulcs86=sum(hcarbo86pd     ,aprot86pd,afat86pd);
run;

proc   rank data=h86_nts groups=11 out=h86_nts;
var    carbo86p  sumref_carb86p  hcarbo86p   prot86p  aprot86p  vprot86p   ;
ranks  carbo86pdd sumref_carb86pdd hcarbo86pdd  prot86pdd aprot86pdd vprot86pdd  ;
run;

proc   rank data=h86_nts descending groups=11 out=h86_nts ;
var    tfat86p  afat86p  vfat86p  ;
ranks  tfat86pdd afat86pdd vfat86pdd ;
run;

data h86_nts;
set  h86_nts;

lfs86=sum(carbo86pdd,prot86pdd,tfat86pdd);
alfs86=sum(carbo86pdd,aprot86pdd,vfat86pdd); /*updated: from tfat86pdd to vfat86pdd */
vlfs86=sum(carbo86pdd,vprot86pdd,afat86pdd); /*updated: from tfat86pdd to afat86pdd */

hlfs86=sum(hcarbo86pdd, vprot86pdd, afat86pdd);
ulfs86=sum(sumref_carb86pdd,aprot86pdd,vfat86pdd);
run;

/* nts90 */
%h90_nts(keep=id alco90n calor90n na90n

  carbo90n 
  pot_carb90n/*Carbohydrate from potato*/
  sumref_carb90n/********Sum of carbohydrate from potato, added sugar in foods, and refined grain, foods only********/
  sumrf_nfr90n/*Sum of added sugar minus fruit carb, plus carb from refined grain, cereals only*/ 
  adds_cer90n/********Added sugar minus cereal fruit for Carb Quality Variables only********/
  adds_fds90n/********Added Sugar From Foods for Carb Quality Variables only********/
  ref_grn90n/*Carbohydrate from refined grain*/

  leg_carb90n/*Carbohdyrate from legumes*/
  frtcb_njs90n/*Fruit carbohydrate minus juice sugar*/
  vegcb_npl90n/*Vegetable Carbohydrate minus potato and legume carbohydrate*/
  sum_whgcarb90n/********Sum Carbohydrate From Wholegrain********/
  sumoth_crb_njs90n/********Sum of carbohydrate from fruit minus juice sugar, plus carb from vegetable, legume, and lact********/
  mill_grn90n/*Carbohydrate from milled whole grain*/
  int_grn90n/*Carbohydrate from intact whole grain*/

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

data h90_nts;
set  h90_nts;

array gram_cp carbo90n 
              pot_carb90n sumref_carb90n sumrf_nfr90n adds_cer90n adds_fds90n ref_grn90n 
              leg_carb90n frtcb_njs90n vegcb_npl90n sum_whgcarb90n sumoth_crb_njs90n mill_grn90n int_grn90n 
              prot90n aprot90n vprot90n dprot90n ;
array perc_cp carbo90p 
              pot_carb90p sumref_carb90p sumrf_nfr90p adds_cer90p adds_fds90p ref_grn90p 
              leg_carb90p frtcb_njs90p vegcb_npl90p sum_whgcarb90p sumoth_crb_njs90p mill_grn90p int_grn90p  
              prot90p aprot90p vprot90p dprot90p ;

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

proc   rank data=h90_nts groups=11 out=h90_nts;
var    prot90p  aprot90p  vprot90p  tfat90p  afat90p  vfat90p  hfat90p  ufat90p;
ranks  prot90pd aprot90pd vprot90pd tfat90pd afat90pd vfat90pd hfat90pd ufat90pd;
run;

proc   rank data=h90_nts descending groups=11 out=h90_nts;
var    carbo90p  sumref_carb90p  hcarbo90p;
ranks  carbo90pd sumref_carb90pd hcarbo90pd;
run;

data h90_nts;
set h90_nts;

lcs90=sum(carbo90pd,prot90pd,tfat90pd);
alcs90=sum(carbo90pd,aprot90pd,afat90pd);
vlcs90=sum(carbo90pd,vprot90pd,vfat90pd);

hlcs90=sum(sumref_carb90pd,vprot90pd,vfat90pd);
ulcs90=sum(hcarbo90pd     ,aprot90pd,afat90pd);

run;

proc rank data=h90_nts groups=11 out=h90_nts;
var    carbo90p  sumref_carb90p  hcarbo90p   prot90p  aprot90p  vprot90p   ;
ranks  carbo90pdd sumref_carb90pdd hcarbo90pdd  prot90pdd aprot90pdd vprot90pdd  ;
run;

proc rank data=h90_nts descending groups=11 out=h90_nts ;
var    tfat90p  afat90p  vfat90p  ;
ranks  tfat90pdd afat90pdd vfat90pdd ;
run;

data h90_nts;
set h90_nts;

lfs90=sum(carbo90pdd,prot90pdd,tfat90pdd);
alfs90=sum(carbo90pdd,aprot90pdd,vfat90pdd);
vlfs90=sum(carbo90pdd,vprot90pdd,afat90pdd);
hlfs90=sum(hcarbo90pdd, vprot90pdd, afat90pdd);
ulfs90=sum(sumref_carb90pdd,aprot90pdd,vfat90pdd);
run;

/* nts94 */
%h94_nts(keep=id alco94n calor94n na94n

  carbo94n 
  pot_carb94n/*Carbohydrate from potato*/
  sumref_carb94n/********Sum of carbohydrate from potato, added sugar in foods, and refined grain, foods only********/
  sumrf_nfr94n/*Sum of added sugar minus fruit carb, plus carb from refined grain, cereals only*/ 
  adds_cer94n/********Added sugar minus cereal fruit for Carb Quality Variables only********/
  adds_fds94n/********Added Sugar From Foods for Carb Quality Variables only********/
  ref_grn94n/*Carbohydrate from refined grain*/

  leg_carb94n/*Carbohdyrate from legumes*/
  frtcb_njs94n/*Fruit carbohydrate minus juice sugar*/
  vegcb_npl94n/*Vegetable Carbohydrate minus potato and legume carbohydrate*/
  sum_whgcarb94n/********Sum Carbohydrate From Wholegrain********/
  sumoth_crb_njs94n/********Sum of carbohydrate from fruit minus juice sugar, plus carb from vegetable, legume, and lact********/
  mill_grn94n/*Carbohydrate from milled whole grain*/
  int_grn94n/*Carbohydrate from intact whole grain*/

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

data h94_nts;
set h94_nts;

array gram_cp carbo94n 
              pot_carb94n sumref_carb94n sumrf_nfr94n adds_cer94n adds_fds94n ref_grn94n 
              leg_carb94n frtcb_njs94n vegcb_npl94n sum_whgcarb94n sumoth_crb_njs94n mill_grn94n int_grn94n 
              prot94n aprot94n vprot94n dprot94n ;
array perc_cp carbo94p 
              pot_carb94p sumref_carb94p sumrf_nfr94p adds_cer94p adds_fds94p ref_grn94p 
              leg_carb94p frtcb_njs94p vegcb_npl94p sum_whgcarb94p sumoth_crb_njs94p mill_grn94p int_grn94p  
              prot94p aprot94p vprot94p dprot94p ;

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

proc rank data=h94_nts groups=11 out=h94_nts;
var    prot94p  aprot94p  vprot94p  tfat94p  afat94p  vfat94p  hfat94p  ufat94p;
ranks  prot94pd aprot94pd vprot94pd tfat94pd afat94pd vfat94pd hfat94pd ufat94pd;
run;

proc rank data=h94_nts descending groups=11 out=h94_nts;
var    carbo94p  sumref_carb94p  hcarbo94p;
ranks  carbo94pd sumref_carb94pd hcarbo94pd;
run;

data h94_nts;
set h94_nts;

lcs94=sum(carbo94pd,prot94pd,tfat94pd);
alcs94=sum(carbo94pd,aprot94pd,afat94pd);
vlcs94=sum(carbo94pd,vprot94pd,vfat94pd);

hlcs94=sum(sumref_carb94pd,vprot94pd,vfat94pd);
ulcs94=sum(hcarbo94pd     ,aprot94pd,afat94pd);

run;

proc rank data=h94_nts groups=11 out=h94_nts;
var    carbo94p  sumref_carb94p  hcarbo94p   prot94p  aprot94p  vprot94p   ;
ranks  carbo94pdd sumref_carb94pdd hcarbo94pdd  prot94pdd aprot94pdd vprot94pdd  ;
run;

proc rank data=h94_nts descending groups=11 out=h94_nts ;
var    tfat94p  afat94p  vfat94p  ;
ranks  tfat94pdd afat94pdd vfat94pdd ;
run;

data h94_nts;
set h94_nts;

lfs94=sum(carbo94pdd,prot94pdd,tfat94pdd);
alfs94=sum(carbo94pdd,aprot94pdd,vfat94pdd);
vlfs94=sum(carbo94pdd,vprot94pdd,afat94pdd);
hlfs94=sum(hcarbo94pdd, vprot94pdd, afat94pdd);
ulfs94=sum(sumref_carb94pdd,aprot94pdd,vfat94pdd);
run;

/* nts98 */
%h98_nts(keep=id alco98n calor98n na98n

  carbo98n 
  pot_carb98n/*Carbohydrate from potato*/
  sumref_carb98n/********Sum of carbohydrate from potato, added sugar in foods, and refined grain, foods only********/
  sumrf_nfr98n/*Sum of added sugar minus fruit carb, plus carb from refined grain, cereals only*/ 
  adds_cer98n/********Added sugar minus cereal fruit for Carb Quality Variables only********/
  adds_fds98n/********Added Sugar From Foods for Carb Quality Variables only********/
  ref_grn98n/*Carbohydrate from refined grain*/

  leg_carb98n/*Carbohdyrate from legumes*/
  frtcb_njs98n/*Fruit carbohydrate minus juice sugar*/
  vegcb_npl98n/*Vegetable Carbohydrate minus potato and legume carbohydrate*/
  sum_whgcarb98n/********Sum Carbohydrate From Wholegrain********/
  sumoth_crb_njs98n/********Sum of carbohydrate from fruit minus juice sugar, plus carb from vegetable, legume, and lact********/
  mill_grn98n/*Carbohydrate from milled whole grain*/
  int_grn98n/*Carbohydrate from intact whole grain*/

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
  ttrn098n/*Total Trans Fat*/
 
  prot98n
  aprot98n/*Animal Protein*/
  vprot98n/*Vegetable Protein*/
  dprot98n/*Dairy Protein*/
  ); 

data h98_nts;
set h98_nts;

array gram_cp carbo98n 
              pot_carb98n sumref_carb98n sumrf_nfr98n adds_cer98n adds_fds98n ref_grn98n 
              leg_carb98n frtcb_njs98n vegcb_npl98n sum_whgcarb98n sumoth_crb_njs98n mill_grn98n int_grn98n 
              prot98n aprot98n vprot98n dprot98n ;
array perc_cp carbo98p 
              pot_carb98p sumref_carb98p sumrf_nfr98p adds_cer98p adds_fds98p ref_grn98p 
              leg_carb98p frtcb_njs98p vegcb_npl98p sum_whgcarb98p sumoth_crb_njs98p mill_grn98p int_grn98p  
              prot98p aprot98p vprot98p dprot98p ;

array gram_f  tfat98n afat98n vfat98n dfat98n sat98n mon98n mufa_pl98n mufa_a98n poly98n omega98n ttrn098n;
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

proc rank data=h98_nts groups=11 out=h98_nts;
var    prot98p  aprot98p  vprot98p  tfat98p  afat98p  vfat98p  hfat98p  ufat98p;
ranks  prot98pd aprot98pd vprot98pd tfat98pd afat98pd vfat98pd hfat98pd ufat98pd;
run;

proc rank data=h98_nts descending groups=11 out=h98_nts;
var    carbo98p  sumref_carb98p  hcarbo98p;
ranks  carbo98pd sumref_carb98pd hcarbo98pd;
run;

data h98_nts;
set h98_nts;

lcs98=sum(carbo98pd,prot98pd,tfat98pd);
alcs98=sum(carbo98pd,aprot98pd,afat98pd);
vlcs98=sum(carbo98pd,vprot98pd,vfat98pd);

hlcs98=sum(sumref_carb98pd,vprot98pd,vfat98pd);
ulcs98=sum(hcarbo98pd     ,aprot98pd,afat98pd);
run;

proc rank data=h98_nts groups=11 out=h98_nts;
var    carbo98p  sumref_carb98p  hcarbo98p   prot98p  aprot98p  vprot98p   ;
ranks  carbo98pdd sumref_carb98pdd hcarbo98pdd  prot98pdd aprot98pdd vprot98pdd  ;
run;

proc rank data=h98_nts descending groups=11 out=h98_nts ;
var    tfat98p  afat98p  vfat98p  ;
ranks  tfat98pdd afat98pdd vfat98pdd ;
run;

data h98_nts;
set h98_nts;

lfs98=sum(carbo98pdd,prot98pdd,tfat98pdd);
alfs98=sum(carbo98pdd,aprot98pdd,vfat98pdd);
vlfs98=sum(carbo98pdd,vprot98pdd,afat98pdd);
hlfs98=sum(hcarbo98pdd, vprot98pdd, afat98pdd);
ulfs98=sum(sumref_carb98pdd,aprot98pdd,vfat98pdd);
run;

/* nts02 */
%h02_nts(keep=id alco02n calor02n na02n

  carbo02n 
  pot_carb02n/*Carbohydrate from potato*/
  sumref_carb02n/********Sum of carbohydrate from potato, added sugar in foods, and refined grain, foods only********/
  sumrf_nfr02n/*Sum of added sugar minus fruit carb, plus carb from refined grain, cereals only*/ 
  adds_cer02n/********Added sugar minus cereal fruit for Carb Quality Variables only********/
  adds_fds02n/********Added Sugar From Foods for Carb Quality Variables only********/
  ref_grn02n/*Carbohydrate from refined grain*/

  leg_carb02n/*Carbohdyrate from legumes*/
  frtcb_njs02n/*Fruit carbohydrate minus juice sugar*/
  vegcb_npl02n/*Vegetable Carbohydrate minus potato and legume carbohydrate*/
  sum_whgcarb02n/********Sum Carbohydrate From Wholegrain********/
  sumoth_crb_njs02n/********Sum of carbohydrate from fruit minus juice sugar, plus carb from vegetable, legume, and lact********/
  mill_grn02n/*Carbohydrate from milled whole grain*/
  int_grn02n/*Carbohydrate from intact whole grain*/

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
  ttrn0202n/*Total Trans Fat*/
 
  prot02n
  aprot02n/*Animal Protein*/
  vprot02n/*Vegetable Protein*/
  dprot02n/*Dairy Protein*/
  ); 

data h02_nts;
set h02_nts;

array gram_cp carbo02n 
              pot_carb02n sumref_carb02n sumrf_nfr02n adds_cer02n adds_fds02n ref_grn02n 
              leg_carb02n frtcb_njs02n vegcb_npl02n sum_whgcarb02n sumoth_crb_njs02n mill_grn02n int_grn02n 
              prot02n aprot02n vprot02n dprot02n ;
array perc_cp carbo02p 
              pot_carb02p sumref_carb02p sumrf_nfr02p adds_cer02p adds_fds02p ref_grn02p 
              leg_carb02p frtcb_njs02p vegcb_npl02p sum_whgcarb02p sumoth_crb_njs02p mill_grn02p int_grn02p  
              prot02p aprot02p vprot02p dprot02p ;

array gram_f  tfat02n afat02n vfat02n dfat02n sat02n mon02n mufa_pl02n mufa_a02n poly02n omega02n ttrn0202n;
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

proc rank data=h02_nts groups=11 out=h02_nts;
var    prot02p  aprot02p  vprot02p  tfat02p  afat02p  vfat02p  hfat02p  ufat02p;
ranks  prot02pd aprot02pd vprot02pd tfat02pd afat02pd vfat02pd hfat02pd ufat02pd;
run;

proc rank data=h02_nts descending groups=11 out=h02_nts;
var    carbo02p  sumref_carb02p  hcarbo02p;
ranks  carbo02pd sumref_carb02pd hcarbo02pd;
run;

data h02_nts;
set h02_nts;

lcs02=sum(carbo02pd,prot02pd,tfat02pd);
alcs02=sum(carbo02pd,aprot02pd,afat02pd);
vlcs02=sum(carbo02pd,vprot02pd,vfat02pd);

hlcs02=sum(sumref_carb02pd,vprot02pd,vfat02pd);
ulcs02=sum(hcarbo02pd     ,aprot02pd,afat02pd);
run;

proc rank data=h02_nts groups=11 out=h02_nts;
var    carbo02p  sumref_carb02p  hcarbo02p   prot02p  aprot02p  vprot02p   ;
ranks  carbo02pdd sumref_carb02pdd hcarbo02pdd  prot02pdd aprot02pdd vprot02pdd  ;
run;

proc rank data=h02_nts descending groups=11 out=h02_nts ;
var    tfat02p  afat02p  vfat02p  ;
ranks  tfat02pdd afat02pdd vfat02pdd ;
run;

data h02_nts;
set h02_nts;

lfs02=sum(carbo02pdd,prot02pdd,tfat02pdd);
alfs02=sum(carbo02pdd,aprot02pdd,vfat02pdd);
vlfs02=sum(carbo02pdd,vprot02pdd,afat02pdd);
hlfs02=sum(hcarbo02pdd, vprot02pdd, afat02pdd);
ulfs02=sum(sumref_carb02pdd,aprot02pdd,vfat02pdd);
run;

/* nts06 */
%h06_nts(keep=id alco06n calor06n na06n

  carbo06n 
  pot_carb06n/*Carbohydrate from potato*/
  sumref_carb06n/********Sum of carbohydrate from potato, added sugar in foods, and refined grain, foods only********/
  sumrf_nfr06n/*Sum of added sugar minus fruit carb, plus carb from refined grain, cereals only*/ 
  adds_cer06n/********Added sugar minus cereal fruit for Carb Quality Variables only********/
  adds_fds06n/********Added Sugar From Foods for Carb Quality Variables only********/
  ref_grn06n/*Carbohydrate from refined grain*/

  leg_carb06n/*Carbohdyrate from legumes*/
  frtcb_njs06n/*Fruit carbohydrate minus juice sugar*/
  vegcb_npl06n/*Vegetable Carbohydrate minus potato and legume carbohydrate*/
  sum_whgcarb06n/********Sum Carbohydrate From Wholegrain********/
  sumoth_crb_njs06n/********Sum of carbohydrate from fruit minus juice sugar, plus carb from vegetable, legume, and lact********/
  mill_grn06n/*Carbohydrate from milled whole grain*/
  int_grn06n/*Carbohydrate from intact whole grain*/

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
  ttrn0706n/*Total Trans Fat*/
 
  prot06n
  aprot06n/*Animal Protein*/
  vprot06n/*Vegetable Protein*/
  dprot06n/*Dairy Protein*/
  ); 

data h06_nts;
set h06_nts;

array gram_cp carbo06n 
              pot_carb06n sumref_carb06n sumrf_nfr06n adds_cer06n adds_fds06n ref_grn06n 
              leg_carb06n frtcb_njs06n vegcb_npl06n sum_whgcarb06n sumoth_crb_njs06n mill_grn06n int_grn06n 
              prot06n aprot06n vprot06n dprot06n ;
array perc_cp carbo06p 
              pot_carb06p sumref_carb06p sumrf_nfr06p adds_cer06p adds_fds06p ref_grn06p 
              leg_carb06p frtcb_njs06p vegcb_npl06p sum_whgcarb06p sumoth_crb_njs06p mill_grn06p int_grn06p  
              prot06p aprot06p vprot06p dprot06p ;

array gram_f  tfat06n afat06n vfat06n dfat06n sat06n mon06n mufa_pl06n mufa_a06n poly06n omega06n ttrn0706n;
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

proc rank data=h06_nts groups=11 out=h06_nts;
var    prot06p  aprot06p  vprot06p  tfat06p  afat06p  vfat06p  hfat06p  ufat06p;
ranks  prot06pd aprot06pd vprot06pd tfat06pd afat06pd vfat06pd hfat06pd ufat06pd;
run;

proc rank data=h06_nts descending groups=11 out=h06_nts;
var    carbo06p  sumref_carb06p  hcarbo06p;
ranks  carbo06pd sumref_carb06pd hcarbo06pd;
run;

data h06_nts;
set h06_nts;

lcs06=sum(carbo06pd,prot06pd,tfat06pd);
alcs06=sum(carbo06pd,aprot06pd,afat06pd);
vlcs06=sum(carbo06pd,vprot06pd,vfat06pd);

hlcs06=sum(sumref_carb06pd,vprot06pd,vfat06pd);
ulcs06=sum(hcarbo06pd     ,aprot06pd,afat06pd);

run;

proc rank data=h06_nts groups=11 out=h06_nts;
var    carbo06p  sumref_carb06p  hcarbo06p   prot06p  aprot06p  vprot06p   ;
ranks  carbo06pdd sumref_carb06pdd hcarbo06pdd  prot06pdd aprot06pdd vprot06pdd  ;
run;

proc rank data=h06_nts descending groups=11 out=h06_nts ;
var    tfat06p  afat06p  vfat06p  ;
ranks  tfat06pdd afat06pdd vfat06pdd ;
run;

data h06_nts;
set h06_nts;

lfs06=sum(carbo06pdd,prot06pdd,tfat06pdd);
alfs06=sum(carbo06pdd,aprot06pdd,vfat06pdd);
vlfs06=sum(carbo06pdd,vprot06pdd,afat06pdd);
hlfs06=sum(hcarbo06pdd, vprot06pdd, afat06pdd);
ulfs06=sum(sumref_carb06pdd,aprot06pdd,vfat06pdd);
run;

/* nts10 */
%h10_nts(keep=id alco10n calor10n na10n

  carbo10n 
  pot_carb10n/*Carbohydrate from potato*/
  sumref_carb10n/********Sum of carbohydrate from potato, added sugar in foods, and refined grain, foods only********/
  sumrf_nfr10n/*Sum of added sugar minus fruit carb, plus carb from refined grain, cereals only*/ 
  adds_cer10n/********Added sugar minus cereal fruit for Carb Quality Variables only********/
  adds_fds10n/********Added Sugar From Foods for Carb Quality Variables only********/
  ref_grn10n/*Carbohydrate from refined grain*/

  leg_carb10n/*Carbohdyrate from legumes*/
  frtcb_njs10n/*Fruit carbohydrate minus juice sugar*/
  vegcb_npl10n/*Vegetable Carbohydrate minus potato and legume carbohydrate*/
  sum_whgcarb10n/********Sum Carbohydrate From Wholegrain********/
  sumoth_crb_njs10n/********Sum of carbohydrate from fruit minus juice sugar, plus carb from vegetable, legume, and lact********/
  mill_grn10n/*Carbohydrate from milled whole grain*/
  int_grn10n/*Carbohydrate from intact whole grain*/

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
  ttrn1110n/*Total Trans Fat*/
 
  prot10n
  aprot10n/*Animal Protein*/
  vprot10n/*Vegetable Protein*/
  dprot10n/*Dairy Protein*/
  ); 

data h10_nts;
set h10_nts;

array gram_cp carbo10n 
              pot_carb10n sumref_carb10n sumrf_nfr10n adds_cer10n adds_fds10n ref_grn10n 
              leg_carb10n frtcb_njs10n vegcb_npl10n sum_whgcarb10n sumoth_crb_njs10n mill_grn10n int_grn10n 
              prot10n aprot10n vprot10n dprot10n ;
array perc_cp carbo10p 
              pot_carb10p sumref_carb10p sumrf_nfr10p adds_cer10p adds_fds10p ref_grn10p 
              leg_carb10p frtcb_njs10p vegcb_npl10p sum_whgcarb10p sumoth_crb_njs10p mill_grn10p int_grn10p  
              prot10p aprot10p vprot10p dprot10p ;

array gram_f  tfat10n afat10n vfat10n dfat10n sat10n mon10n mufa_pl10n mufa_a10n poly10n omega10n ttrn1110n;
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

proc rank data=h10_nts groups=11 out=h10_nts;
var    prot10p  aprot10p  vprot10p  tfat10p  afat10p  vfat10p  hfat10p  ufat10p;
ranks  prot10pd aprot10pd vprot10pd tfat10pd afat10pd vfat10pd hfat10pd ufat10pd;
run;

proc rank data=h10_nts descending groups=11 out=h10_nts;
var    carbo10p  sumref_carb10p  hcarbo10p;
ranks  carbo10pd sumref_carb10pd hcarbo10pd;
run;

data h10_nts;
set h10_nts;

lcs10=sum(carbo10pd,prot10pd,tfat10pd);
alcs10=sum(carbo10pd,aprot10pd,afat10pd);
vlcs10=sum(carbo10pd,vprot10pd,vfat10pd);

hlcs10=sum(sumref_carb10pd,vprot10pd,vfat10pd);
ulcs10=sum(hcarbo10pd     ,aprot10pd,afat10pd);

run;

proc rank data=h10_nts groups=11 out=h10_nts;
var    carbo10p  sumref_carb10p  hcarbo10p   prot10p  aprot10p  vprot10p   ;
ranks  carbo10pdd sumref_carb10pdd hcarbo10pdd  prot10pdd aprot10pdd vprot10pdd  ;
run;

proc rank data=h10_nts descending groups=11 out=h10_nts ;
var    tfat10p  afat10p  vfat10p  ;
ranks  tfat10pdd afat10pdd vfat10pdd ;
run;

data h10_nts;
set h10_nts;

lfs10=sum(carbo10pdd,prot10pdd,tfat10pdd);
alfs10=sum(carbo10pdd,aprot10pdd,vfat10pdd);
vlfs10=sum(carbo10pdd,vprot10pdd,afat10pdd);
hlfs10=sum(hcarbo10pdd, vprot10pdd, afat10pdd);
ulfs10=sum(sumref_carb10pdd,aprot10pdd,vfat10pdd);
run;

/* nts14 */
%h14_nts(keep=id alco14n calor14n na14n

  carbo14n 
  pot_carb14n/*Carbohydrate from potato*/
  sumref_carb14n/********Sum of carbohydrate from potato, added sugar in foods, and refined grain, foods only********/
  sumrf_nfr14n/*Sum of added sugar minus fruit carb, plus carb from refined grain, cereals only*/ 
  adds_cer14n/********Added sugar minus cereal fruit for Carb Quality Variables only********/
  adds_fds14n/********Added Sugar From Foods for Carb Quality Variables only********/
  ref_grn14n/*Carbohydrate from refined grain*/

  leg_carb14n/*Carbohdyrate from legumes*/
  frtcb_njs14n/*Fruit carbohydrate minus juice sugar*/
  vegcb_npl14n/*Vegetable Carbohydrate minus potato and legume carbohydrate*/
  sum_whgcarb14n/********Sum Carbohydrate From Wholegrain********/
  sumoth_crb_njs14n/********Sum of carbohydrate from fruit minus juice sugar, plus carb from vegetable, legume, and lact********/
  mill_grn14n/*Carbohydrate from milled whole grain*/
  int_grn14n/*Carbohydrate from intact whole grain*/

  tfat14n
  afat14n/*Animal Fat*/
  vfat14n/*Vegetable Fat*/
  dfat14n/*Dairy Fat*/
  sat14n/*Total Saturated Fat*/
  mon14n/*Total Monounsaturated Fat*/
  mufa_pl14n/*Plant Monounsaturated Fat*/
  mufa_a14n/*Animal Monounsaturated Fat*/
  poly14n/*Total Polyunsaturated Fat*/
  omega14n/*Omega 3 (f205+f226) gm, no alpha 18:3*/
  ttrn1114n/*Total Trans Fat*/
 
  prot14n
  aprot14n/*Animal Protein*/
  vprot14n/*Vegetable Protein*/
  dprot14n/*Dairy Protein*/
  ); 

data h14_nts;
set h14_nts;

array gram_cp carbo14n 
              pot_carb14n sumref_carb14n sumrf_nfr14n adds_cer14n adds_fds14n ref_grn14n 
              leg_carb14n frtcb_njs14n vegcb_npl14n sum_whgcarb14n sumoth_crb_njs14n mill_grn14n int_grn14n 
              prot14n aprot14n vprot14n dprot14n ;
array perc_cp carbo14p 
              pot_carb14p sumref_carb14p sumrf_nfr14p adds_cer14p adds_fds14p ref_grn14p 
              leg_carb14p frtcb_njs14p vegcb_npl14p sum_whgcarb14p sumoth_crb_njs14p mill_grn14p int_grn14p  
              prot14p aprot14p vprot14p dprot14p ;

array gram_f  tfat14n afat14n vfat14n dfat14n sat14n mon14n mufa_pl14n mufa_a14n poly14n omega14n ttrn1114n;
array perc_f  tfat14p afat14p vfat14p dfat14p sat14p mon14p mufa_pl14p mufa_a14p poly14p omega14p trnss14p;

do over gram_cp;
perc_cp=(4*gram_cp)/calor14n*100;
end;

do over gram_f;
perc_f=(9*gram_f)/calor14n*100;
end;

hfat14p=sum(mufa_pl14p,poly14p);
ufat14p=sum(sat14p,mufa_a14p,trnss14p);

hcarbo14p=sum(leg_carb14p,frtcb_njs14p,vegcb_npl14p,sum_whgcarb14p);

run;

proc rank data=h14_nts groups=11 out=h14_nts;
var    prot14p  aprot14p  vprot14p  tfat14p  afat14p  vfat14p  hfat14p  ufat14p;
ranks  prot14pd aprot14pd vprot14pd tfat14pd afat14pd vfat14pd hfat14pd ufat14pd;
run;

proc rank data=h14_nts descending groups=11 out=h14_nts;
var    carbo14p  sumref_carb14p  hcarbo14p;
ranks  carbo14pd sumref_carb14pd hcarbo14pd;
run;

data h14_nts;
set h14_nts;

lcs14=sum(carbo14pd,prot14pd,tfat14pd);
alcs14=sum(carbo14pd,aprot14pd,afat14pd);
vlcs14=sum(carbo14pd,vprot14pd,vfat14pd);

hlcs14=sum(sumref_carb14pd,vprot14pd,vfat14pd);
ulcs14=sum(hcarbo14pd     ,aprot14pd,afat14pd);

run;

proc rank data=h14_nts groups=11 out=h14_nts;
var    carbo14p  sumref_carb14p  hcarbo14p   prot14p  aprot14p  vprot14p   ;
ranks  carbo14pdd sumref_carb14pdd hcarbo14pdd  prot14pdd aprot14pdd vprot14pdd  ;
run;

proc rank data=h14_nts descending groups=11 out=h14_nts ;
var    tfat14p  afat14p  vfat14p  ;
ranks  tfat14pdd afat14pdd vfat14pdd ;
run;

data h14_nts;
set h14_nts;

lfs14=sum(carbo14pdd,prot14pdd,tfat14pdd);
alfs14=sum(carbo14pdd,aprot14pdd,vfat14pdd);
vlfs14=sum(carbo14pdd,vprot14pdd,afat14pdd);
hlfs14=sum(hcarbo14pdd, vprot14pdd, afat14pdd);
ulfs14=sum(sumref_carb14pdd,aprot14pdd,vfat14pdd);
run;

/* nts18 */
%h18_nts(keep=id alco18n calor18n na18n

  carbo18n
  pot_carb18n/*Carbohydrate from potato*/
  sumref_carb18n/********Sum of carbohydrate from potato, added sugar in foods, and refined grain, foods only********/
  sumrf_nfr18n/*Sum of added sugar minus fruit carb, plus carb from refined grain, cereals only*/
  adds_cer18n/********Added sugar minus cereal fruit for Carb Quality Variables only********/
  adds_fds18n/********Added Sugar From Foods for Carb Quality Variables only********/
  ref_grn18n/*Carbohydrate from refined grain*/

  leg_carb18n/*Carbohdyrate from legumes*/
  frtcb_njs18n/*Fruit carbohydrate minus juice sugar*/
  vegcb_npl18n/*Vegetable Carbohydrate minus potato and legume carbohydrate*/
  sum_whgcarb18n/********Sum Carbohydrate From Wholegrain********/
  sumoth_crb_njs18n/********Sum of carbohydrate from fruit minus juice sugar, plus carb from vegetable, legume, and lact********/
  mill_grn18n/*Carbohydrate from milled whole grain*/
  int_grn18n/*Carbohydrate from intact whole grain*/

  tfat18n
  afat18n/*Animal Fat*/
  vfat18n/*Vegetable Fat*/
  dfat18n/*Dairy Fat*/
  sat18n/*Total Saturated Fat*/
  mon18n/*Total Monounsaturated Fat*/
  mufa_pl18n/*Plant Monounsaturated Fat*/
  mufa_a18n/*Animal Monounsaturated Fat*/
  poly18n/*Total Polyunsaturated Fat*/
  omega18n/*Omega 3 (f205+f226) gm, no alpha 18:3*/
  /*ttrn1118n*/ /*Total Trans Fat*/ 

  prot18n
  aprot18n/*Animal Protein*/
  vprot18n/*Vegetable Protein*/
  dprot18n/*Dairy Protein*/
  );

data h18_nts;
set h18_nts;

array gram_cp carbo18n
              pot_carb18n sumref_carb18n sumrf_nfr18n adds_cer18n adds_fds18n ref_grn18n
              leg_carb18n frtcb_njs18n vegcb_npl18n sum_whgcarb18n sumoth_crb_njs18n mill_grn18n int_grn18n
              prot18n aprot18n vprot18n dprot18n ;
array perc_cp carbo18p
              pot_carb18p sumref_carb18p sumrf_nfr18p adds_cer18p adds_fds18p ref_grn18p
              leg_carb18p frtcb_njs18p vegcb_npl18p sum_whgcarb18p sumoth_crb_njs18p mill_grn18p int_grn18p
              prot18p aprot18p vprot18p dprot18p ;

array gram_f  tfat18n afat18n vfat18n dfat18n sat18n mon18n mufa_pl18n mufa_a18n poly18n omega18n /*ttrn1118n*/;
array perc_f  tfat18p afat18p vfat18p dfat18p sat18p mon18p mufa_pl18p mufa_a18p poly18p omega18p /*trnss18p*/;

do over gram_cp;
perc_cp=(4*gram_cp)/calor18n*100;
end;

do over gram_f;
perc_f=(9*gram_f)/calor18n*100;
end;

hfat18p=sum(mufa_pl18p,poly18p);
ufat18p=sum(sat18p,mufa_a18p /*trnss18p*/);

hcarbo18p=sum(leg_carb18p,frtcb_njs18p,vegcb_npl18p,sum_whgcarb18p);

run;

proc rank data=h18_nts groups=11 out=h18_nts;
var    prot18p  aprot18p  vprot18p  tfat18p  afat18p  vfat18p  leg_carb18p  frtcb_njs18p  vegcb_npl18p  sum_whgcarb18p  sumoth_crb_njs18p  mill_grn18p  int_grn18p  hfat18p  ufat18p;
ranks  prot18pd aprot18pd vprot18pd tfat18pd afat18pd vfat18pd leg_carb18pd frtcb_njs18pd vegcb_npl18pd sum_whgcarb18pd sumoth_crb_njs18pd mill_grn18pd int_grn18pd hfat18pd ufat18pd;
run;

proc rank data=h18_nts descending groups=11 out=h18_nts;
var    carbo18p  pot_carb18p  sumref_carb18p  sumrf_nfr18p  adds_cer18p  adds_fds18p  ref_grn18p  hcarbo18p;
ranks  carbo18pd pot_carb18pd sumref_carb18pd sumrf_nfr18pd adds_cer18pd adds_fds18pd ref_grn18pd hcarbo18pd;
run;

data h18_nts;
  set h18_nts;
  lcs18  = sum(carbo18pd, prot18pd,  tfat18pd);
  alcs18 = sum(carbo18pd, aprot18pd, afat18pd);
  vlcs18 = sum(carbo18pd, vprot18pd, vfat18pd);
  hlcs18 = sum(sumref_carb18pd, vprot18pd, vfat18pd);
  ulcs18 = sum(hcarbo18pd,      aprot18pd, afat18pd);
run;

proc rank data=h18_nts groups=11 out=h18_nts;
  var    carbo18p  sumref_carb18p  hcarbo18p   prot18p  aprot18p  vprot18p;
  ranks  carbo18pdd sumref_carb18pdd hcarbo18pdd  prot18pdd aprot18pdd vprot18pdd;
run;

proc rank data=h18_nts descending groups=11 out=h18_nts ;
  var    tfat18p  afat18p  vfat18p;
  ranks  tfat18pdd afat18pdd vfat18pdd;
run;

data h18_nts;
  set h18_nts;

  lfs18  = sum(carbo18pdd, prot18pdd,  tfat18pdd);
  alfs18 = sum(carbo18pdd, aprot18pdd, vfat18pdd);
  vlfs18 = sum(carbo18pdd, vprot18pdd, afat18pdd);
  hlfs18 = sum(hcarbo18pdd,  vprot18pdd, afat18pdd);
  ulfs18 = sum(sumref_carb18pdd, aprot18pdd, vfat18pdd);
run;



                
*----------------3.Covariables-----------------------------;
** Endoscopy **;

/** Dietary Scores **/   
libname AHEIdat 'DIETSCORE_HPFS/AHEI_HPFS/';
data ahei86; set AHEIdat.hnahei86l; ahei86=nAHEI86_noal; run; proc sort nodupkey; by id; run; 
data ahei90; set AHEIdat.hnahei90l; ahei90=nAHEI90_noal; run; proc sort nodupkey; by id; run; 
data ahei94; set AHEIdat.hnahei94l; ahei94=nAHEI94_noal; run; proc sort nodupkey; by id; run; 
data ahei98; set AHEIdat.hnahei98l; ahei98=nAHEI98_noal; run; proc sort nodupkey; by id; run; 
data ahei02; set AHEIdat.hnahei02l; ahei02=nAHEI02_noal; run; proc sort nodupkey; by id; run; 
data ahei06; set AHEIdat.hnahei06l; ahei06=nAHEI06_noal; run; proc sort nodupkey; by id; run; 
data ahei10; set AHEIdat.hnahei10l; ahei10=nAHEI10_noal; run; proc sort nodupkey; by id; run; 
data ahei14; set AHEIdat.hnahei14l; ahei14=nAHEI14_noal; run; proc sort nodupkey; by id; run; 
data ahei18; set AHEIdat.hnahei18l; ahei18=nAHEI18_noal; run; proc sort nodupkey; by id; run; 
    

/** Physical activity **/ 
%hmet8620(keep= act86 act88 act90 act92 act94 act96 act98 act00 act02 act04 act06 act08 act10 act12 act14 act16);
data hp_mets;  set hmet8620;
array acts{16} act86 act88 act90 act92 act94 act96 act98 act00 act02 act04 act06 act08 act10 act12 act14 act16;
    do i=1 to 16; if acts{i}>900 then acts{i}=.; end;
    /* carry forward missing physical activity (following Ambika's plant based - wtchg program) '*/
    do i=2 to 16; if acts{i}=. then acts{i}=acts{i-1}; end;
    proc sort; by id;
run;


/** master and FFQ **/
%hp_der (keep= dbmy86 dbmy09 
               rtmnyr86 rtmnyr88 rtmnyr90 rtmnyr92 rtmnyr94 rtmnyr96 rtmnyr98 rtmnyr00 rtmnyr02 rtmnyr04 rtmnyr06 
               height bmi2186 wt86 wt88 wt90 wt92 wt94 wt96 wt98 wt00 wt02 wt04 wt06 
               bmi86  bmi88  bmi90  bmi92 bmi94 bmi96 bmi98  bmi00  bmi02  bmi04 bmi06
               smoke86 smoke88 smoke90 smoke92 smoke94 smoke96 smoke98 smoke00 smoke02 smoke04 smoke06 
               /* 1.Never 2.Past; 3.No, unknown past history; 4.Current; 5.PASSTHRU */
               cgnm86 cgnm88 cgnm90 cgnm92 cgnm94 cgnm96 cgnm98 cgnm00 cgnm02 cgnm04 cgnm06 /* no. of cigarettes per day */
               /*pckyr86 pckyr88 pckyr90 pckyr92 pckyr94 pckyr96 pckyr98 pckyr00 pckyr02 pckyr04 pckyr06*/
               fdb3087 mdb3087 sdb3087 bdb3087 dbfh87
               fclc86 afclc86 mclc86 amclc86 cafh86
               fclc90 fpro90  mclc90 sclc90 spro90 cafh90 
               stpyr86);
               
               /* family history of db */
               dbfh87=0;    * onset after age 30 presumed to be NIDDM;                  
               if fdb3087=1 or mdb3087=1 or sdb3087=1 or bdb3087=1 then dbfh87=1;
               /* family history of Cancer */
               cafh86=0; if fclc86=1 or mclc86=1 then cafh86=1;
               cafh90=0; if fclc90=1 or fpro90=1 or  mclc90=1 or sclc90=1 or spro90=1 then cafh90=1;run;


%hp_der_2 (keep=rtmnyr08 rtmnyr10 rtmnyr12 rtmnyr14 rtmnyr16 rtmnyr18 
                smoke08  smoke10 smoke12 smoke14 smoke16 smoke18 
                cgnm08   cgnm10 cgnm12 cgnm14 cgnm16 cgnm18 
                wt08     wt10   wt12   wt14   wt16   wt18   
                bmi08    bmi10  bmi12  bmi14  bmi16  bmi18);


%hp86(keep=hbp86 db86 chol86 trig86 mi86 ang86 str86 cabg86 pe86
           colc86  pros86  lymp86  ocan86  mel86
           seuro86 scand86 ocauc86 afric86 asian86 oanc86  /* races */
           mmi86 fmi86
           asp86 mar86 livng86
           betab86 lasix86 diur86 calcb86 ald86 /* anti-hypertensive */
           antch86); 
           if hbp86 ne 1 then hbp86=0; if chol86 ne 1 then chol86=0;
           if db86 ne 1 then db86=0; if mi86 ne 1 then mi86=0; if str86 ne 1 then str86=0;
           if colc86 ne 1 then colc86=0; if mel86 ne 1 then mel86=0; if pros86 ne 1 then pros86=0;
           if lymp86 ne 1 then lymp86=0; if ocan86 ne 1 then ocan86=0;
           if asp86=1 then asp86=1;else asp86=0;run;


%h86_dt(keep=nblnk86d mvt86d mvt86  /* number of unanswered food, multivitamin intakes 0=never, 1=past, 2=current;*/
             rais86d   avoc86d    ban86d    cant86d   h2om86d   appl86d    aj86d   peach86d 
             straw86d  blueb86d   oran86d      oj86d   grfr86d    grj86d  othj86d                      
             sbean86d  brocc86d  kraut86d    slaw86d  ccabb86d  cauli86d  bruss86d   car86d     corn86d   peas86d   mixv86d    
             bean86d  alfsp86d    cel86d    mush86d  ysqua86d  eggpl86d    yam86d cspin86d    rspin86d   kale86d  ilett86d  rlett86d    
             grpep86d   garl86d    tom86d     toj86d  tosau86d  chsau86d   tofu86d   
             bacon86d    dog86d   hamb86d    bmix86d  bmain86d  procm86d     
             ctuna86d  dkfsh86d  ofish86d   shrim86d     
             fries86d    pot86d        
             pchip86d  crack86d    popc86d       
             choco86d   cdyw86d   cdywo86d   cokh86d   cokr86d   brwn86d   donut86d  cakh86d    cakr86d     pieh86d   pier86d   srol86d
             cerbr86d  cer86d    oat86d   ckcer86d   dkbr86d  brice86d  otgrn86d    bran86d wgerm86d     
             whbr86d   engl86d    muff86d  wrice86d  pasta86d  pcake86d    
             lccaf86d   lcnoc86d   lcoth86d   cola86d    cnoc86d    otsug86d   
             skim86d whole86d sherb86d icecr86d yog86d cream86d sour86d cotch86d crmch86d otch86d    egg86d
             pnut86d   onut86d 
             beer86d  rwine86d   wwine86d   liq86d   h2o86d tea86d cof86d dcaf86d punch86d                 
             chwi86d   chwo86d
             q2286d    q2386d  ); if mvt86d=0 or mvt86d=1 then mvt86=1; else if mvt86d=2 then mvt86=2; else mvt86=3; run;

data h86_dt;
set  h86_dt;

if cerbr86d in ( 0, 4, 11, 12, 13, 15, 16, 17, 19, 20, 21, 22, 23, 25, 29, 30, 31, 33, 35,37,38, 42, 43, 45, 46, 47,
59, 60,61,62,63,64,65,66,68,73,74,75,76,79,80,81,82,83, 84,85,86,88,95,96,97,98,
101,109,111, 112,113,114, 115,119, 123,124, 127, 128, 129,130,131, 132, 133, 134, 146, 147, 148,149,150,
151, 152,153, 154, 155, 157,158, 164,165, 166, 167,170,173, 176, 177,182, 184, 185, 186, 188 /*,190,197,198,997*/ )
   then rcer86d=cer86d;
   else wcer86d=cer86d;
   if   rcer86d=. then rcer86d=0;
   if   wcer86d=. then wcer86d=0;

array foodd 
      rais86d   avoc86d    ban86d    cant86d   h2om86d   appl86d    aj86d   peach86d 
      straw86d  blueb86d   oran86d      oj86d   grfr86d    grj86d  othj86d                      
      sbean86d  brocc86d  kraut86d    slaw86d  ccabb86d  cauli86d  bruss86d   car86d     corn86d   peas86d   mixv86d    
      bean86d  alfsp86d    cel86d    mush86d  ysqua86d  eggpl86d    yam86d cspin86d    rspin86d   kale86d  ilett86d  rlett86d    
      grpep86d   garl86d    tom86d     toj86d  tosau86d  chsau86d   tofu86d   
      bacon86d    dog86d   hamb86d    bmix86d  bmain86d  procm86d     
      ctuna86d  dkfsh86d  ofish86d   shrim86d     
      fries86d    pot86d        
      pchip86d  crack86d    popc86d       
      choco86d   cdyw86d   cdywo86d   cokh86d   cokr86d   brwn86d   donut86d  cakh86d    cakr86d     pieh86d   pier86d   srol86d
      wcer86d   rcer86d    oat86d   ckcer86d   dkbr86d  brice86d  otgrn86d     bran86d wgerm86d     
      whbr86d   engl86d    muff86d  wrice86d  pasta86d  pcake86d    
      lccaf86d   lcnoc86d   lcoth86d   cola86d    cnoc86d    otsug86d   
      skim86d whole86d sherb86d icecr86d yog86d cream86d sour86d cotch86d crmch86d otch86d   egg86d
      pnut86d   onut86d 
      beer86d  rwine86d   wwine86d   liq86d  h2o86d tea86d cof86d dcaf86d punch86d                    
      chwi86d   chwo86d
;
   
array food 
      rais86   avoc86    ban86    cant86   h2om86   appl86    aj86   peach86 
      straw86  blueb86   oran86      oj86   grfr86    grj86  othj86                       
      sbean86  brocc86  kraut86    slaw86  ccabb86  cauli86  bruss86   car86     corn86   peas86   mixv86     
      bean86  alfsp86    cel86    mush86  ysqua86  eggpl86    yam86 cspin86    rspin86   kale86  ilett86  rlett86    
      grpep86   garl86    tom86     toj86  tosau86  chsau86   tofu86   
      bacon86    dog86   hamb86    bmix86  bmain86  procm86     
      ctuna86  dkfsh86  ofish86   shrim86     
      fries86    pot86        
      pchip86  crack86    popc86      
      choco86   cdyw86   cdywo86   cokh86   cokr86   brwn86   donut86  cakh86    cakr86      pieh86   pier86   srol86
      wcer86   rcer86    oat86   ckcer86   dkbr86  brice86  otgrn86    bran86 wgerm86     
      whbr86   engl86    muff86  wrice86  pasta86  pcake86    
      lccaf86   lcnoc86   lcoth86   cola86    cnoc86    otsug86   
      skim86 whole86 sherb86 icecr86 yog86 cream86 sour86 cotch86 crmch86 otch86    egg86
      pnut86   onut86 
      beer86  rwine86   wwine86   liq86  h2o86 tea86 coff86 dcaf86  punch86                  
      chwi86   chwo86
;
   
do over foodd; 
        if foodd=8  then food=6;     /*6+/day           */
   else if foodd=7  then food=4.5;   /*4-5/day          */
   else if foodd=6  then food=2.5;   /*2-3/day          */
   else if foodd=5  then food=1;     /*1/day            */
   else if foodd=4  then food=0.79;  /*5-6/wk           */
   else if foodd=3  then food=0.43;  /*2-4/wk           */
   else if foodd=2  then food=0.14;  /*1/wk             */
   else if foodd=1  then food=0.07;  /*1-3/month        */
   else if foodd=0  then food=0;     /*Never or <1/Month*/
   else if foodd=9  then food=.;     /*Pass Through     */
end; 

array friedd  q2386d   q2286d;
array fried   fhome86  faway86;

do over friedd; 
        if friedd=4  then fried=1;      /*Daily         */
   else if friedd=3  then fried=0.71;   /*4-6           */
   else if friedd=2  then fried=0.29;   /*1-3/wk        */
   else if friedd=1  then fried=0.07;   /*<1/wk         */
   else if friedd=9  then fried=.;      /*Pass Through  */
end; 

/* take sum of caffeinated and decaffeinated tea & coffee */
tea86d  = SUM(tea86);
coff86d  = sum(coff86, dcaf86);
/* whole grains */
whgrn86d=sum(wcer86,    oat86,  ckcer86,   dkbr86,  brice86,  otgrn86,    bran86, wgerm86,popc86);
/* refined grains */
rgrain86d  = sum(rcer86, whbr86 ,  wrice86,  pasta86,  engl86 ,  muff86 ,  pcake86 );
/* fruit and vegetable */
frui86d=sum(rais86, ban86, cant86, h2om86, avoc86, appl86, peach86, oran86, grfr86, straw86, blueb86);
vege86d=sum(sbean86,brocc86,  kraut86,    slaw86,  ccabb86,  cauli86,  bruss86, car86,    /*corn86,*/   peas86,   mixv86,     
                 bean86,  alfsp86,    cel86,    mush86,  ysqua86,  eggpl86,    yam86, cspin86,    rspin86,    kale86,  ilett86,  rlett86,   
                grpep86,   tom86,     toj86,  tosau86,  chsau86,   tofu86  );
fv86d=sum(frui86d, vege86d);
/* potatoes */
pota86d=sum(fries86, pot86, pchip86);
/* legumes */
leg86d=sum(sbean86, tofu86, bean86, peas86);
/* fruit juices */
fj86d=sum(aj86, oj86, grj86, othj86  );
/* vegetables minus potatoes and legumes */
veg_noleg86d=vege86d-leg86d;

/* dairy and fish */
dairy86d=sum(skim86, whole86, sherb86, icecr86, yog86, cream86, sour86d ,cotch86d ,crmch86, otch86);
fish86d=sum(ctuna86, dkfsh86, ofish86, shrim86);
egg86d=egg86;

/* poultry */
poul86d=sum(chwi86,   chwo86);

/* red and processed meat */
rpmeats86d = sum(bacon86, dog86, procm86, hamb86, bmix86, bmain86*1.67); 
/* ssb */
ssb86d = sum(cola86, cnoc86, otsug86, punch86); 
/* SWEETS & DESERTS */
swtsvg86d= sum(choco86, cdyw86, cdywo86, cokh86, cokr86, brwn86, donut86, cakh86, cakr86, pieh86, pier86, srol86); 
swtsnochoc86d=sum(cdyw86, cdywo86, cokh86, cokr86, brwn86, donut86, cakh86, cakr86, pieh86, pier86, srol86);
run;


%hp88(keep=hbp88 db88 chol88 trig88 mi88 ang88 str88 cabg88 pe88
              colc88  pros88  lymp88  ocan88  mel88  
              mvt88 mvtr88 asp88 mar88  livng88
              betab88 lasix88 diur88 ccblo88 ald88
              chrx88
              flag88 sigsy88 sigsc88 sigx88 sigpt88 colsy88 colsc88 colx88 colpt88 q25pt88 ); 
              if asp88=1 then asp88=1;else asp88=0;
              if mvt88=1 then mvtr88=2; else if mvt88=2 then mvtr88=1; else mvtr88=3; 
              mvt88=mvtr88;
              run;

%hp90(keep=hbp90 db90 chol90 trig90 mi90 ang90 str90 cabg90 pe90
              colc90  pros90  lymp90  ocan90  mel90  
              asp90 mar90 livng90 
              mdb90 fdb90 sdb90 dbfh90 
              betab90 lasix90 diur90 ccblo90 ald90
              chrx90
              flag90 sigm90 sigm190 sigmr90 vbld90 fecal90 abdpn90 diarr90 hxclc90 rout90 q13bp90 ); 
              if asp90=1 then asp90=1;else asp90=0;
              if fdb90=1 or mdb90=1 or sdb90=1 then dbfh90=1; else dbfh90=0;run;


%h90_dt(keep=mvt90d mvt90   
          rais90d    ban90d    cant90d   h2om90d   appl90d    aj90d   peach90d asau90d
          straw90d  blueb90d   oran90d      oj90d   grfr90d    grj90d  othj90d  prune90d                      
          sbean90d  brocc90d   slaw90d   cauli90d  bruss90d   rcar90d  ccar90d   corn90d   peas90d   mixv90d    
           bean90d    cel90d  osqua90d   eggpl90d    yam90d cspin90d    rspin90d   kale90d  ilett90d  rlett90d    
           beet90d    tom90d     toj90d  tosau90d  chsau90d   tofu90d    oniog90d oniov90d
          bacon90d    dog90d   hamb90d    bmix90d  bmain90d  procm90d   pmain90d  
          ctuna90d  dkfsh90d  ofsh90d /**/  shrim90d      
          fries90d    potb90d /**/       
          pchip90d  crack90d    popc90d       
          choco90d   cdyw90d   cdywo90d   cokh90d   cokr90d   brwn90d   donut90d  cakh90d    cakr90d     pieh90d   pier90d   srolh90d srolr90d/**/
          cerbr90d  cer90d    oat90d   ckcer90d   drbr90d /**/  brice90d  otgrn90d     bran90d wgerm90d  oatbr90d     
           whbr90d   engl90d    muff90d  wrice90d  pasta90d  pcake90d    
           lccaf90d   lcnoc90d   lcoth90d   cola90d    cnoc90d    otsug90d   
           skim90d whole90d sherb90d icecr90d yog90d cream90d sour90d cotch90d crmch90d otch90d   egg90d
           pnut90d   onut90d 
            q3690d    q3790d
           beer90d  rwine90d   wwine90d   liq90d  h2o90d tea90d cof90d decaf90d punch90d                     
           chwi90d chwo90d ); mvt90=mvt90d; run;

data h90_dt;
  set h90_dt;

if cerbr90d in ( 0, 4, 11, 12, 13, 15, 16, 17, 19, 20, 21, 22, 23, 25, 29, 30, 31, 33, 35,37,38, 42, 43, 45, 46, 47,
59, 60,61,62,63,64,65,66,68,73,74,75,76,79,80,81,82,83, 84,85,86,88,95,96,97,98,
101,109,111, 112,113,114, 115,119, 123,124, 127, 128, 129,130,131, 132, 133, 134, 146, 147, 148,149,150,
151, 152,153, 154, 155, 157,158, 164,165, 166, 167,170,173, 176, 177,182, 184, 185, 186, 188/*,
190,197,198,997*/)
   then rcer90d=cer90d;
   else wcer90d=cer90d;
   if rcer90d=. then rcer90d=0;
   if wcer90d=. then wcer90d=0;   
array foodd  
           rais90d    ban90d    cant90d   h2om90d   appl90d    aj90d   peach90d asau90d
          straw90d  blueb90d   oran90d      oj90d   grfr90d    grj90d  othj90d  prune90d                    
          sbean90d  brocc90d   slaw90d   cauli90d  bruss90d   rcar90d  ccar90d   corn90d   peas90d   mixv90d    
           bean90d    cel90d  osqua90d   eggpl90d    yam90d cspin90d    rspin90d   kale90d  ilett90d  rlett90d    
           beet90d    tom90d     toj90d  tosau90d  chsau90d   tofu90d    oniog90d oniov90d
          bacon90d    dog90d   hamb90d    bmix90d  bmain90d  procm90d  pmain90d   
          ctuna90d  dkfsh90d  ofsh90d /**/  shrim90d      
          fries90d    potb90d/**/       
          pchip90d  crack90d    popc90d       
          choco90d   cdyw90d   cdywo90d   cokh90d   cokr90d   brwn90d   donut90d  cakh90d    cakr90d     pieh90d   pier90d   srolh90d srolr90d/**/
          rcer90d  wcer90d    oat90d   ckcer90d   drbr90d /**/  brice90d  otgrn90d     bran90d wgerm90d oatbr90d      
           whbr90d   engl90d    muff90d  wrice90d  pasta90d  pcake90d    
           lccaf90d   lcnoc90d   lcoth90d   cola90d    cnoc90d    otsug90d   
           skim90d whole90d sherb90d icecr90d yog90d cream90d sour90d cotch90d crmch90d otch90d   egg90d
           pnut90d   onut90d 

           beer90d  rwine90d   wwine90d   liq90d  h2o90d tea90d cof90d decaf90d punch90d                     
           chwi90d chwo90d
;
array food 
           rais90    ban90    cant90   h2om90   appl90    aj90   peach90 asau90
          straw90  blueb90   oran90      oj90   grfr90    grj90  othj90  prune90                    
          sbean90  brocc90   slaw90   cauli90  bruss90   rcar90  ccar90   corn90   peas90   mixv90    
           bean90    cel90  osqua90   eggpl90    yam90 cspin90    rspin90   kale90  ilett90  rlett90    
           beet90    tom90     toj90  tosau90  chsau90    tofu90    oniog90 oniov90
          bacon90    dog90   hamb90    bmix90  bmain90  procm90  pmain90    
          ctuna90  dkfsh90  ofish90 /**/  shrim90     
          fries90    pot90 /**/        
          pchip90  crack90    popc90      
          choco90   cdyw90   cdywo90   cokh90   cokr90   brwn90   donut90  cakh90    cakr90      pieh90   pier90   srolh90 srolr90/**/
          rcer90  wcer90    oat90   ckcer90   dkbr90 /**/  brice90  otgrn90    bran90 wgerm90 oatbr90     
           whbr90   engl90    muff90  wrice90  pasta90  pcake90    
           lccaf90   lcnoc90   lcoth90   cola90    cnoc90    otsug90   
           skim90 whole90 sherb90 icecr90 yog90 cream90 sour90 cotch90 crmch90 otch90   egg90
           pnut90   onut90 
           beer90  rwine90   wwine90   liq90    h2o90 tea90 coff90 decaf90 punch90                   
           chwi90 chwo90
;
     
do over foodd; 
   if foodd=8  then food=6;          /*6+/day           */
   else if foodd=7  then food=4.5;   /*4-5/day          */
   else if foodd=6  then food=2.5;   /*2-3/day          */
   else if foodd=5  then food=1;     /*1/day            */
   else if foodd=4  then food=0.79;  /*5-6/wk           */
   else if foodd=3  then food=0.43;  /*2-4/wk           */
   else if foodd=2  then food=0.14;  /*1/wk             */
   else if foodd=1  then food=0.07;  /*1-3/month        */
   else if foodd=0  then food=0;     /*Never or <1/Month*/
   else if foodd=9  then food=.;     /*Pass Through     */
end; 

array friedd  q3690d   q3790d;
array fried   fhome90  faway90;

do over friedd; 
        if friedd=4  then fried=1;      /*Daily         */
   else if friedd=3  then fried=0.71;   /*4-6           */
   else if friedd=2  then fried=0.29;   /*1-3/wk        */
   else if friedd=1  then fried=0.07;   /*<1/wk         */
   else if friedd=5  then fried=.;      /*Pass Through  */
end; 


/* take sum of caffeinated and decaffeinated tea & coffee */
tea90d  = SUM(tea90);
coff90d  = sum(coff90, decaf90);
/** whole grains **/
whgrn90d=sum( wcer90,  oat90,ckcer90, dkbr90,brice90,otgrn90,oatbr90,bran90,wgerm90);
/** refined grains **/
rgrain90d  = sum(rcer90, whbr90 ,  wrice90,  pasta90,  engl90 ,  muff90 ,  pcake90);

/* fruit and vegetable*/
frui90d=sum(rais90, prune90, ban90, cant90, h2om90, appl90, asau90, peach90, oran90, grfr90, straw90, blueb90);
vege90d=sum(tom90, toj90, tosau90, chsau90, tofu90,sbean90,brocc90, cauli90,bruss90, rcar90,ccar90, /*corn90,*/
                 peas90, mixv90, bean90,osqua90,eggpl90,  yam90,cspin90,rspin90,kale90,ilett90,  
                rlett90,  cel90, beet90,oniog90,oniov90);
fv90d=sum(frui90d, vege90d);


/* potatoes */
pota90d=sum(fries90, pot90, pchip90);
/*legumes */
leg90d=sum(sbean90, tofu90, bean90, peas90);
/* fruit juices */
fj90d=sum(aj90, oj90, grj90, othj90  );

/* vegetables minus potatoes and legumes */
veg_noleg90d=vege90d-leg90d;

/* dairy and fish */
dairy90d=sum(skim90, whole90, sherb90, icecr90, yog90, cream90, sour90d ,cotch90d ,crmch90, otch90);
fish90d=sum(ctuna90, dkfsh90, ofish90, shrim90);
egg90d=egg90;

/* poultry */
poul90d=sum(chwi90,   chwo90);

/* red and processed meat */
rpmeats90d = sum(bacon90, dog90, procm90, hamb90, bmix90, pmain90*1.67, bmain90*1.67); 
/* ssb*/
ssb90d = sum(cola90, cnoc90, otsug90, punch90); 
/*** SWEETS & DESERTS  ***/
swtsvg90d= sum(choco90, cdyw90, cdywo90, cokh90, cokr90, brwn90, donut90, cakh90, cakr90, pieh90, pier90, srolr90, srolh90); 
swtsnochoc90d=sum( cdyw90, cdywo90, cokh90, cokr90, brwn90, donut90, cakh90, cakr90, pieh90, pier90, srolr90, srolh90);
run;

%hp92(keep=hbp92 db92 chol92 trig92 mi92 ang92 str92 cabg92 
              colc92  pros92  lymp92  ocan92  mel92  
              mvt92 asp92 mar92 livng92
              betab92 lasix92 thiaz92 ccblo92 ald92  
              chrx92
              fdb92 mdb92 sdb92 dbfh92
              mlng92 mclc92 flng92 fclc92 fpro92 slng92 sclc92 spro92 cafh92
              flag92 sigm92 sig9092 vbld92 fecal92 abdpn92 diarr92 hxclc92 rout92 q13bp92 );
              if fdb92=1 or mdb92=1 or sdb92=1 then dbfh92=1; else dbfh92=0; 
              if mlng92=1 or mclc92=1 or flng92=1 or fclc92=1 or fpro92=1 or slng92=1 or sclc92=1 or spro92=1 then cafh92=1; 
              else cafh92=0;
              if asp92=1 then asp92=1;else asp92=0;run;

%hp94(keep=hbp94 db94 chol94 trig94 mi94 ang94 str94 cabg94 pe94
              colc94  pros94  lymp94  ocan94  mel94  
              asp94 mar94 livng94
              betab94 lasix94 thiaz94 ccblo94 ald94 
              chrx94
              flag94 sigm94 vbld94 fecal94 abdpn94 diarr94 hxclc94 rout94 q9apt94 );
              if asp94=1 then asp94=1;else asp94=0;run;


%h94_dt(keep=  mvt94d mvt94 
rais94d    prune94d    ban94d      cant94d     avo94d      appl94d     aj94d       
oran94d    oj94d       grfr94d     grj94d      othj94d     straw94d    blueb94d    
peach94d     
tom94d     toj94d      tosau94d    salsa94d    tofu94d     sbean94d    brocc94d    
cabb94d    cauli94d    bruss94d    rcar94d     ccar94d     corn94d     peas94d     
mixv94d    bean94d     osqua94d    eggpl94d    yam94d      cspin94d    rspin94d    
kale94d    ilett94d    rlett94d    cel94d      grpep94d    oniog94d    oniov94d      
bacon94d   dog94d      ctdog94d    procm94d    hamb94d     hambl94d    bmix94d     
pmain94d   bmain94d      
ctuna94d   dkfsh94d    ofish94d    shrim94d      
fries94d   pot94d        
pchip94d   crack94d    popc94d 
choco94d   cdyw94d     cdywo94d    cokh94d     cokr94d     brwn94d     donut94d    
cakh94d    cakr94d     pieh94d     pier94d     srolh94d    srolr94d        
cerbr94d   cer94d      oat94d      ckcer94d    dkbr94d     brice94d    otgrn94d    oatbr94d    
bran94d    wgerm94d        
whbr94d    engl94d     muff94d     wrice94d    pasta94d    pcake94d      
cola94d    cnoc94d     otsug94d 
lccaf94d        lcnoc94d         lcoth94d 
skim194d milk294d whole94d sherb94d icecr94d flyog94d plyog94d cream94d cotch94d crmch94d otch94d but94d 
egg94d
pnut94d onut94d 
q4094d       q4194d
chtyp94d 
beer94d    lbeer94d   rwine94d   wwine94d   liq94d   h2o94d tea94d coff94d dcaf94d punch94d    
chwi94d chwo94d ctdog94d );
mvt94=mvt94d;
run;

data h94_dt;
  set h94_dt;

if cerbr94d in ( 0, 4, 11, 12, 13, 15, 16, 17, 19, 20, 21, 22, 23, 25, 29, 30, 31, 33, 35,37,38, 42, 43, 45, 46, 47,
59, 60,61,62,63,64,65,66,68,73,74,75,76,79,80,81,82,83, 84,85,86,88,95,96,97,98,
101,109,111, 112,113,114, 115,119, 123,124, 127, 128, 129,130,131, 132, 133, 134, 146, 147, 148,149,150,
151, 152,153, 154, 155, 157,158, 164,165, 166, 167,170,173, 176, 177,182, 184, 185, 186, 188/*,
190,197,198,997*/)
   then rcer94d=cer94d;
   else wcer94d=cer94d;
   if rcer94d=. then rcer94d=0;
   if wcer94d=. then wcer94d=0;   

   array foodd 
rais94d    prune94d    ban94d      cant94d     avo94d      appl94d     aj94d       
oran94d    oj94d       grfr94d     grj94d      othj94d     straw94d    blueb94d    
peach94d     
tom94d     toj94d      tosau94d    salsa94d    tofu94d     sbean94d    brocc94d    
cabb94d    cauli94d    bruss94d    rcar94d     ccar94d     corn94d     peas94d     
mixv94d    bean94d     osqua94d    eggpl94d    yam94d      cspin94d    rspin94d    
kale94d    ilett94d    rlett94d    cel94d      grpep94d    oniog94d    oniov94d      
bacon94d   dog94d      ctdog94d    procm94d    hamb94d     hambl94d    bmix94d     
pmain94d   bmain94d      
ctuna94d   dkfsh94d    ofish94d    shrim94d      
fries94d   pot94d        
pchip94d   crack94d    popc94d 
choco94d   cdyw94d     cdywo94d    cokh94d     cokr94d     brwn94d     donut94d    
cakh94d    cakr94d     pieh94d     pier94d     srolh94d    srolr94d        
wcer94d    rcer94d     oat94d      ckcer94d    dkbr94d     brice94d    otgrn94d    oatbr94d    
bran94d    wgerm94d        
whbr94d    engl94d     muff94d     wrice94d    pasta94d    pcake94d      
cola94d    cnoc94d     otsug94d 
lccaf94d        lcnoc94d         lcoth94d 
skim194d milk294d whole94d sherb94d icecr94d flyog94d plyog94d cream94d cotch94d crmch94d otch94d but94d    egg94d
pnut94d onut94d 
beer94d    lbeer94d   rwine94d   wwine94d   liq94d    h2o94d tea94d coff94d dcaf94d punch94d    
chwi94d chwo94d ctdog94d
;
   
array food 
rais94     prune94   ban94     cant94    avo94     appl94    aj94      
oran94     oj94      grfr94    grj94     othj94    straw94   blueb94   
peach94      
tom94      toj94     tosau94   salsa94   tofu94    sbean94   brocc94   
cabb94     cauli94   bruss94   rcar94    ccar94    corn94    peas94    
mixv94     bean94    osqua94   eggpl94   yam94     cspin94   rspin94   
kale94     ilett94   rlett94   cel94     grpep94   oniog94   oniov94     
bacon94    dog94     ctdog94   procm94   hamb94    hambl94   bmix94    
pmain94    bmain94     
ctuna94    dkfsh94   ofish94   shrim94     
fries94    pot94       
pchip94    crack94   popc94 
choco94    cdyw94    cdywo94   cokh94    cokr94    brwn94    donut94   
cakh94     cakr94    pieh94    pier94    srolh94   srolr94       
wcer94     rcer94      oat94     ckcer94   dkbr94    brice94   otgrn94   oatbr94   
bran94     wgerm94       
whbr94     engl94    muff94    wrice94   pasta94   pcake94     
cola94     cnoc94    otsug94 
lccaf94        lcnoc94         lcoth94 
skim194 milk294 whole94 sherb94 icecr94 flyog94 plyog94 cream94 cotch94 crmch94 otch94  but94    egg94
pnut94 onut94 
beer94    lbeer94   rwine94   wwine94   liq94    h2o94 tea94 coff94 dcaf94 punch94    
chwi94 chwo94 ctdog94
; 
     
do over foodd; 
   if foodd=8  then food=6;          /*6+/day           */
   else if foodd=7  then food=4.5;   /*4-5/day          */
   else if foodd=6  then food=2.5;   /*2-3/day          */
   else if foodd=5  then food=1;     /*1/day            */
   else if foodd=4  then food=0.79;  /*5-6/wk           */
   else if foodd=3  then food=0.43;  /*2-4/wk           */
   else if foodd=2  then food=0.14;  /*1/wk             */
   else if foodd=1  then food=0.07;  /*1-3/month        */
   else if foodd=0  then food=0;     /*Never or <1/Month*/
   else if foodd=9  then food=.;     /*Pass Through     */
end; 

array friedd  q4094d   q4194d;
array fried   fhome94  faway94;

do over friedd; 
        if friedd=4  then fried=1;      /*Daily          */
   else if friedd=3  then fried=0.71;   /*4-6/wk         */
   else if friedd=2  then fried=0.29 ;  /*1-3/wk         */
   else if friedd=1  then fried=0.07;   /*never or <1/wk */
   else if friedd=5  then fried=.;      /*Pass Through   */
end; 


/* take sum of caffeinated and decaffeinated tea & coffee */
tea94d  = SUM(tea94);
coff94d  = sum(coff94, dcaf94);
/** whole grains **/
whgrn94d=sum( wcer94,  oat94,ckcer94, dkbr94,brice94,otgrn94,oatbr94,bran94,wgerm94);
/** refined grains **/
rgrain94d  = sum(rcer94, whbr94 ,  wrice94,  pasta94,  engl94 ,  muff94 ,  pcake94);
/* fruit and vegetable*/
frui94d=sum(rais94, prune94, ban94, cant94, avo94, appl94, peach94, oran94, grfr94, straw94, blueb94);
vege94d=sum(tom94 ,   toj94  ,   tosau94,   salsa94,   tofu94 ,   sbean94,   brocc94,   
               cabb94 ,   cauli94,   bruss94,   rcar94 ,   ccar94 ,   /*corn94 ,*/   peas94 ,   
               mixv94 ,   bean94 ,   osqua94,   eggpl94,   yam94  ,   cspin94,   rspin94,   
               kale94 ,   ilett94,   rlett94,   cel94  ,   grpep94,   oniog94,   oniov94 );
fv94d=sum(frui94d, vege94d);

/* potatoes */
pota94d=sum(fries94, pot94, pchip94);
/*legumes */
leg94d=sum(sbean94, tofu94, bean94, peas94);
/* fruit juices */
fj94d=sum(aj94, oj94, grj94, othj94  );

/* vegetables minus potatoes and legumes */
veg_noleg94d=vege94d-leg94d;

/* dairy and fish */
dairy94d=sum(skim194, milk294, whole94, sherb94, icecr94, flyog94, plyog94, cream94, cotch94, crmch94, otch94);
fish94d=sum(ctuna94, dkfsh94, ofish94, shrim94);
egg94d=egg94;

/* poultry */
poul94d=sum(chwi94, chwo94, ctdog94);

/* red and processed meat */
rpmeats94d = sum(bacon94, dog94, procm94, hamb94,hambl94, bmix94, pmain94*1.67, bmain94*1.67); 
/* ssb*/
ssb94d = sum(cola94, cnoc94, otsug94, punch94); 
/*** SWEETS & DESERTS  ***/
swtsvg94d= sum(choco94, cdyw94,cdywo94, cokh94, cokr94, brwn94,donut94, cakh94,cakr94, pieh94, pier94,srolh94,srolr94, icecr94); 
swtsnochoc94d=sum(cdyw94,cdywo94, cokh94, cokr94, brwn94,donut94, cakh94,cakr94, pieh94, pier94,srolh94,srolr94, icecr94);
run;


%hp96(keep=hbp96 db96 chol96 trig96 mi96 ang96 str96 cabg96 pe96 
              colc96  pros96  lymp96  ocan96  mel96 
              asp96 aspd96 mvt96 mar96  livar96 
              betab96 lasix96 thiaz96 ccblo96 ald96 
              chrx96
              pclc96 sclc196 sclc296 fpro96 bpro196 bpro296 sbrc196 mbrcn96 cafh96
              flag96 sig9496 vbld96 fecal96 abdpn96 diarr96 hxclc96 rout96 q15pt96 ); 
              if aspd96 in (2,3,4,5) then asp96=1; else asp96=0;
              if pclc96=1 or sclc196=1 or sclc296=1 or fpro96=1 or bpro196=1 or bpro296=1 or sbrc196=1 or mbrcn96=1 then cafh96=1; 
              else cafh96=0;
              run;

%hp98(keep=hbp98 db98 chol98 trig98 mi98 ang98 str98 cabg98 pe98 
              colc98  pros98  lymp98  ocan98  mel98 
              mar98  livng98  asp98 aspfr98
              betab98 lasix98 thiaz98 ccblo98 ald98 
              chrx98
              flag98 sigm98 vbld98 fecal98 abdpn98 diar98 hxclc98 rout98 sigmp98 );
              if aspfr98 in (2,3,4,5,6) then asp98=1; else asp98=0;run;



%h98_dt(keep=  mvt98d mvt98 
rais98d       prune98d      ban98d        cant98d       avo98d        appl98d       aj98d         oran98d       
oj98d        ojca98d       grfr98d       grj98d        othj98d       straw98d      blueb98d      peach98d        
tom98d       toj98d        tosau98d      salsa98d      tofu98d       sbean98d      brocc98d      cabb98d       
cauli98d     bruss98d      rcar98d       ccar98d       corn98d       peas98d       mixv98d       bean98d       
yam98d       osqua98d      eggpl98d      kale98d       cspin98d      rspin98d      ilett98d      rlett98d      
cel98d       grpep98d      oniog98d      oniov98d          
dog98d  ctdog98d  ctsan98d  chwi98d   chwo98d   bacon98d  sbol98d      oprom98d    hamb98d   hambl98d   bmix98d  pmain98d   bmain98d       
ctuna98d     bfsh98d       shrim98d      dkfsh98d      ofish98d        
fries98d     pot98d          
pchip98d     crack98d      popc98d    
choco98d     cdyw98d       cdywo98d      cokff98d      cokr98d       cokh98d       brwn98d       donut98d      
cakr98d      cakh98d       piehr98d      srolh98d      srolr98d      srolh98d        
cerbr98d   cer98d        oat98d        ckcer98d      dkbr98d       brice98d      otgrn98d      oatbr98d      bran98d       wgerm98d       
ooil98d
whbr98d      engl98d       muff98d       pcake98d      wrice98d      pasta98d        egg98d eggbt98d
cola98d      otsug98d       cnoc98d 
lccaf98d         lcoth98d        lcnoc98d 
skim198d milk298d whole98d sherb98d icecr98d flyog98d plyog98d cream98d cotch98d crmch98d otch98d but98d
pnut98d wnut98d onut98d 
q3498d       q3798d
chreg98d chlf98d chnf98d noch98d
beer98d      lbeer98d     rwine98d     wwine98d     liq98d    h2o98d htea98d tea98d coff98d dcaf98d punch98d      
chwi98d chwo98d ctdog98d ctsan98d);
mvt98=mvt98d;

data h98_dt;
  set h98_dt;

if cerbr98d in ( 0, 4, 11, 12, 13, 15, 16, 17, 19, 20, 21, 22, 23, 25, 29, 30, 31, 33, 35,37,38, 42, 43, 45, 46, 47,
59, 60,61,62,63,64,65,66,68,73,74,75,76,79,80,81,82,83, 84,85,86,88,95,96,97,98,
101,109,111, 112,113,114, 115,119, 123,124, 127, 128, 129,130,131, 132, 133, 134, 146, 147, 148,149,150,
151, 152,153, 154, 155, 157,158, 164,165, 166, 167,170,173, 176, 177,182, 184, 185, 186, 188/*,
190,197,198,997*/)
   then rcer98d=cer98d;
   else wcer98d=cer98d;
   if rcer98d=. then rcer98d=0;
   if wcer98d=. then wcer98d=0;   
array foodd  
rais98d      prune98d      ban98d        cant98d       avo98d        appl98d       aj98d         oran98d       
oj98d        ojca98d       grfr98d       grj98d        othj98d       straw98d      blueb98d      peach98d        
tom98d       toj98d        tosau98d      salsa98d      tofu98d       sbean98d      brocc98d      cabb98d       
cauli98d     bruss98d      rcar98d       ccar98d       corn98d       peas98d       mixv98d       bean98d       
yam98d       osqua98d      eggpl98d      kale98d       cspin98d      rspin98d      ilett98d      rlett98d      
cel98d       grpep98d      oniog98d      oniov98d          
dog98d       ctdog98d      ctsan98d      chwi98d   chwo98d   bacon98d  sbol98d   oprom98d  hamb98d   hambl98d bmix98d  pmain98d  bmain98d       
ctuna98d     bfsh98d       shrim98d      dkfsh98d      ofish98d        egg98d eggbt98d
fries98d     pot98d          
pchip98d     crack98d      popc98d    
choco98d     cdyw98d       cdywo98d      cokff98d      cokr98d       cokh98d       brwn98d       donut98d      
cakr98d      cakh98d       piehr98d      srolh98d      srolr98d      srolh98d        
rcer98d       wcer98d        oat98d        ckcer98d      dkbr98d       brice98d      otgrn98d      oatbr98d      bran98d       wgerm98d       
ooil98d
whbr98d      engl98d       muff98d       pcake98d      wrice98d      pasta98d        
cola98d      otsug98d      cnoc98d 
lccaf98d         lcoth98d        lcnoc98d 
skim198d milk298d whole98d sherb98d icecr98d flyog98d plyog98d cream98d cotch98d crmch98d otch98d but98d
pnut98d wnut98d onut98d 
beer98d      lbeer98d     rwine98d     wwine98d     liq98d   h2o98d htea98d tea98d coff98d dcaf98d punch98d     
chwi98d chwo98d ctdog98d ctsan98d
;
   
array food 
rais98       prune98     ban98       cant98      avo98       appl98      aj98        oran98      
oj98         ojca98      grfr98      grj98       othj98      straw98     blueb98     peach98       
tom98        toj98       tosau98     salsa98     tofu98      sbean98     brocc98     cabb98      
cauli98      bruss98     rcar98      ccar98      corn98      peas98      mixv98      bean98      
yam98       osqua98     eggpl98     kale98      cspin98     rspin98     ilett98     rlett98     
cel98        grpep98     oniog98     oniov98         
dog98          ctdog98          ctsan98   chwi98          chwo98         bacon98         sbol98      oprocm98      
hamb98       hambl98     bmix98      pmain98     bmain98      
ctuna98      bfsh98      shrim98     dkfsh98     ofish98       egg98 eggbt98
fries98      pot98         
pchip98      crack98     popc98    
choco98      cdyw98      cdywo98     coknf98     cokr98      cokh98      brwn98      donut98     
cakr98       cakh98      pieh98      srolf98     srolr98     srolh98       
rcer98       wcer98        oat98       ckcer98     dkbr98      brice98     otgrn98     oatbr98     bran98      wgerm98      
ooil98 
whbr98       engl98      muff98      pcake98     wrice98     pasta98       
cola98       otsug98    cnoc98  
lccaf98         lcoth98        lcnoc98  
skim198 milk298 whole98 sherb98 icecr98 flyog98 plyog98 cream98 cotch98 crmch98 otch98 but98d
pnut98 wnut98 onut98 
beer98      lbeer98     rwine98     wwine98     liq98  h2o98 htea98 tea98 coff98 dcaf98 punch98      
chwi98 chwo98  ctdog98  ctsan98
; 
     
do over foodd; 
   if foodd=8  then food=6;          /*6+/day           */
   else if foodd=7  then food=4.5;   /*4-5/day          */
   else if foodd=6  then food=2.5;   /*2-3/day          */
   else if foodd=5  then food=1;     /*1/day            */
   else if foodd=4  then food=0.79;  /*5-6/wk           */
   else if foodd=3  then food=0.43;  /*2-4/wk           */
   else if foodd=2  then food=0.14;  /*1/wk             */
   else if foodd=1  then food=0.07;  /*1-3/month        */
   else if foodd=0  then food=0;     /*Never or <1/Month*/
   else if foodd=9  then food=.;     /*Pass Through     */
end; 

array friedd  q3498d   q3798d;
array fried   fhome98  faway98;

do over friedd; 
        if friedd=4  then fried=1;      /*Daily          */
   else if friedd=3  then fried=0.71;   /*4-6/wk         */
   else if friedd=2  then fried=0.29;   /*1-3/wk         */
   else if friedd=1  then fried=0.07;   /*never or <1/wk */
   else if friedd=5  then fried=.;      /*Pass Through   */
end; 


/* take sum of caffeinated and decaffeinated tea & coffee */
tea98d  = SUM(tea98, htea98);
coff98d  = sum(coff98, dcaf98);
/** whole grains **/
whgrn98d=sum(wcer98, oat98, ckcer98, dkbr98, brice98, otgrn98, oatbr98, bran98, wgerm98 );
/** refined grains **/
rgrain98d  = sum(rcer98, whbr98 ,  wrice98,  pasta98,  engl98 ,  muff98 ,  pcake98);

/* fruit and vegetable*/
frui98d=sum(rais98, prune98, ban98, cant98, avo98, appl98, peach98, oran98, grfr98, straw98, blueb98);
vege98d=sum(tom98   ,  toj98   ,  tosau98 ,  salsa98 ,  tofu98  ,  sbean98 ,  brocc98 ,  cabb98  ,  
               cauli98 ,  bruss98 ,  rcar98  ,  ccar98  ,  /*corn98,*/  peas98  ,  mixv98  ,  bean98  ,  
                yam98  ,  osqua98 ,  eggpl98 ,  kale98  ,  cspin98 ,  rspin98 ,  ilett98 ,  rlett98 ,  
                 cel98   ,  grpep98 ,  oniog98 ,  oniov98);
fv98d=sum(frui98d, vege98d);

/* potatoes */
pota98d=sum(fries98, pot98, pchip98);
/*legumes */
leg98d=sum(sbean98, tofu98, bean98, peas98);
/* fruit juices */
fj98d=sum(aj98, oj98, grj98, othj98  );

/* vegetables minus potatoes and legumes */
veg_noleg98d=vege98d-leg98d;

/* dairy and fish */
dairy98d=sum(skim198, milk298, whole98, sherb98, icecr98, flyog98, plyog98, cream98, cotch98, crmch98, otch98);
fish98d=sum(ctuna98, dkfsh98, ofish98, shrim98);
egg98d=sum(egg98, eggbt98);

/* poultry */
poul98d=sum(chwi98, chwo98,  ctdog98,  ctsan98);
/* red and processed meat */
rpmeats98d = sum(bacon98, dog98, sbol98, oprocm98, hamb98,hambl98, bmix98, pmain98*1.67, bmain98*1.67); 
/* ssb*/
ssb98d = sum(cola98, cnoc98, otsug98, punch98d); 
/*** SWEETS & DESERTS  ***/
swtsvg98d= sum(choco98, cdyw98, cdywo98, coknf98, cokr98, cokh98, brwn98, donut98, cakr98, cakh98, pieh98, srolf98, srolr98, srolh98, icecr98); 
swtsnochoc98d=sum(cdyw98, cdywo98, coknf98, cokr98, cokh98, brwn98, donut98, cakr98, cakh98, pieh98, srolf98, srolr98, srolh98, icecr98);
run;


%hp00(keep=hbp00 db00 chol00 trig00 mi00 ang00 strk00 cabg00 pe00 /* str=strk */
            colc00  pros00  lymp00  ocan00  mel00 
            asp00 mvt00 mar00 livng00
            betab00 lasix00 thiaz00 ccblo00 oanth00
            chrx00 ochrx00
            flag00 sigm00 vbld00 fecal00 abdpn00 diar00 hxclc00 rout00 sigmp00 );


%hp02(keep=hbp02 db02 chol02 trig02 mi02 ang02 strk02 cabg02 pe02 /* str=strk */
              colc02  pros02  lymp02  ocan02  mel02 
              betab02 lasix02 thiaz02 calcb02 anthp02  
              asp02 mvt02 mar02 lalon02
              stat02 ochrx02
              flag02 sigmd02 vbld02 fecal02 abdpn02 diarr02 hxclc02 rout02 scopt02 

              rais02d       prune02d      prunj02d      ban02d        cant02d       appl02d       
aj02d        oran02d       ojca02d       ojreg02d        grfrj02d      othj02d     straw02d      blueb02d      peach02d        
tom02d       toj02d        tosau02d      salsa02d      tofu02d       sbean02d      bean02d       peas02d       soym02d
brocc02d     cauli02d      cabb02d       bruss02d      rcar02d       ccar02d       corn02d       mixv02d       
yam02d       osqua02d      eggpl02d      kale02d       cspin02d      rspin02d      ilett02d      rlett02d      
cel02d       grpep02d      oniog02d      oniov02d        
dog02d  ctdog02d ctsan02d chwi02d  chwo02d  bacon02d     sbol02d     oprom02d      hambl02d      hamb02d       
bmix02d       pmain02d      bmain02d        
ctuna02d     bfsh02d       shrim02d      dkfsh02d      ofish02d        
fries02d     pot02d          
pchip02d     crklt02d      crack02d      popcl02d    popcr02d          
choco02d     cdyw02d       cdywo02d      cokff02d      cokr02d       cokh02d       brwn02d       donut02d      
cakr02d      cakh02d       piehr02d      srolh02d      srolr02d  srlff02d      
cerbr02d cer02d        oat02d        ckcer02d      dkbr02d       brice02d      oatbr02d      bran02d       wgerm02d        
whbr02d      wrice02d      pasta02d      engl02d       muff02d       pcake02d        
colas02d       otsug02d   obvcs02d 
lccaf02d         lcbnc02d   lcbcf02d
skim02d mlk1202d whole02d sherb02d icecr02d plyog02d flyog02d cream02d cotch02d crmch02d otch02d  eggbt02d egg02d
pnut02d wnut02d onut02d 
ffh02d ffa02d
chreg02d chlf02d 
chnf02d          
beer02d     lbeer02d    rwine02d    wwine02d    liq02d   h2o02d htea02d tea02d coff02d dcaf02d punch02d     
chwi02d chwo02d ctdog02d ctsan02d ); 

data hp02;
  set hp02;


if cerbr02d in ( 0, 4, 11, 12, 13, 15, 16, 17, 19, 20, 21, 22, 23, 25, 29, 30, 31, 33, 35,37,38, 42, 43, 45, 46, 47,
59, 60,61,62,63,64,65,66,68,73,74,75,76,79,80,81,82,83, 84,85,86,88,95,96,97,98,
101,109,111, 112,113,114, 115,119, 123,124, 127, 128, 129,130,131, 132, 133, 134, 146, 147, 148,149,150,
151, 152,153, 154, 155, 157,158, 164,165, 166, 167,170,173, 176, 177,182, 184, 185, 186, 188/*,
190,197,198,997*/)
   then rcer02d=cer02d;
   else wcer02d=cer02d;
   if rcer02d=. then rcer02d=0;
   if wcer02d=. then wcer02d=0;

array foodd  
rais02d      prune02d      prunj02d      ban02d        cant02d       appl02d       
aj02d        oran02d       ojca02d       ojreg02d        grfrj02d      othj02d       straw02d      blueb02d      peach02d        
tom02d       toj02d        tosau02d      salsa02d      tofu02d       sbean02d      bean02d       peas02d       soym02d
brocc02d     cauli02d      cabb02d       bruss02d      rcar02d       ccar02d       corn02d       mixv02d       
yam02d       osqua02d      eggpl02d      kale02d       cspin02d      rspin02d      ilett02d      rlett02d      
cel02d       grpep02d      oniog02d      oniov02d        
dog02d  ctdog02d ctsan02d chwi02d  chwo02d  bacon02d     sbol02d     oprom02d      hambl02d      
hamb02d       bmix02d       pmain02d      bmain02d        
ctuna02d     bfsh02d       shrim02d      dkfsh02d      ofish02d        
fries02d     pot02d          
pchip02d     crklt02d      crack02d      popcl02d    popcr02d        
choco02d     cdyw02d       cdywo02d      cokff02d      cokr02d       cokh02d       brwn02d       donut02d      
cakr02d      cakh02d       piehr02d      srolh02d      srolr02d      srlff02d            
wcer02d     rcer02d        oat02d        ckcer02d      dkbr02d       brice02d      oatbr02d      bran02d       wgerm02d        
whbr02d      wrice02d      pasta02d      engl02d       muff02d       pcake02d        
colas02d       otsug02d   obvcs02d 
lccaf02d         lcbnc02d   lcbcf02d
skim02d mlk1202d whole02d sherb02d icecr02d plyog02d flyog02d cream02d cotch02d crmch02d otch02d  eggbt02d egg02d
pnut02d wnut02d onut02d 
beer02d     lbeer02d    rwine02d    wwine02d    liq02d  h2o02d htea02d tea02d coff02d dcaf02d punch02d     
chwi02d chwo02d ctdog02d ctsan02d
;
   
array food  

rais02       prune02     prunj02     ban02       cant02      appl02      
aj02         oran02      ojca02      ojreg02         grfrj02       othj02      straw02     blueb02     peach02       
tom02        toj02       tosau02     salsa02     tofu02      sbean02     bean02      peas02      soym02
brocc02      cauli02     cabb02      bruss02     rcar02      ccar02      corn02      mixv02      
yams02       osqua02     eggpl02     kale02      cspin02     rspin02     ilett02     rlett02     
cel02        grpep02     oniog02     oniov02       
dog02  ctdog02 ctsan02 chwi02  chwo02  bacon02     sbol02      oprom02     hambl02     hamb02      bmix02      pmain02     bmain02       
ctuna02      bfsh02      shrim02     dkfsh02     ofish02       
fries02      pot02         
pchip02      crklt02     crack02     popcl02    popcr02            
choco02      cdyw02      cdywo02     cokff02     cokr02      cokh02      brwn02      donut02     
cakr02       cakh02      pieh02      srolh02     srolr02     srlff02           
wcer02      rcer02         oat02       ckcer02     dkbr02      brice02     oatbr02     bran02      wgerm02       
whbr02       wrice02     pasta02     engl02      muff02      pcake02       
colas02      otsug02 obvcs02
lccaf02         lcbnc02  lcbcf02
skim02 mlk1202 whole02 sherb02 icecr02 plyog02 flyog02 cream02 cotch02 crmch02 otch02  eggbt02 egg02
pnut02 wnut02 onut02 
beer02     lbeer02    rwine02    wwine02    liq02  h2o02 htea02 tea02 coff02 dcaf02 punch02     
chwi02  chwo02  ctdog02  ctsan02 
; 
     
do over foodd; 
   if foodd=8  then food=6;          /*6+/day           */
   else if foodd=7  then food=4.5;   /*4-5/day          */
   else if foodd=6  then food=2.5;   /*2-3/day          */
   else if foodd=5  then food=1;     /*1/day            */
   else if foodd=4  then food=0.79;  /*5-6/wk           */
   else if foodd=3  then food=0.43;  /*2-4/wk           */
   else if foodd=2  then food=0.14;  /*1/wk             */
   else if foodd=1  then food=0.07;  /*1-3/month        */
   else if foodd=0  then food=0;     /*Never or <1/Month*/
   else if foodd=9  then food=.;     /*Pass Through     */
end; 

array friedd  ffh02d    ffa02d;
array fried   fhome02  faway02;

do over friedd; 
        if friedd=4  then fried=1;      /*Daily          */
   else if friedd=3  then fried=0.71;   /*4-6/wk         */
   else if friedd=2  then fried=0.29 ;  /*1-3/wk         */
   else if friedd=1  then fried=0.07;   /*never or <1/wk */
   else if friedd=5  then fried=.;      /*Pass Through   */
end; 

/* take sum of caffeinated and decaffeinated tea & coffee */
tea02d  = SUM(tea02, htea02);
coff02d  = sum(coff02, dcaf02);
/** whole grains **/
whgrn02d=sum(wcer02  ,  oat02  ,  ckcer02,  dkbr02 ,  brice02,  oatbr02,  bran02 ,  wgerm02);
/** refined grains **/
rgrain02d  = sum(rcer02, whbr02 ,  wrice02,  pasta02,  engl02 ,  muff02 ,  pcake02);

/* fruit and vegetable*/
frui02d=sum(rais02, prune02, ban02, cant02, appl02, peach02, oran02, grfrj02, straw02, blueb02);
vege02d=sum(tom02,  toj02,   tosau02,  salsa02,   tofu02, sbean02,  bean02,    peas02,  soym02,
                brocc02,cauli02,    cabb02,  bruss02,   rcar02,  ccar02,  /*corn02,*/    mixv02,  
                 yams02,osqua02,   eggpl02,   kale02,  cspin02, rspin02, ilett02,   rlett02,  
                  cel02,grpep02,   oniog02,  oniov02 );
fv02d=sum(frui02d, vege02d);

/* potatoes */
pota02d=sum(fries02, pot02, pchip02);
/*legumes */
leg02d=sum(sbean02, tofu02, bean02, peas02, soym02);
/* fruit juices */
fj02d=sum(aj02, ojca02, ojreg02, prunj02, othj02);

/* vegetables minus potatoes and legumes */
veg_noleg02d=vege02d-leg02d;

/* dairy and fish */
dairy02d=sum(skim02, mlk1202, whole02, sherb02, icecr02, plyog02, flyog02, cream02, cotch02, crmch02, otch02);
fish02d=sum(ctuna02, dkfsh02, ofish02, shrim02);
egg02d=sum(egg02, eggbt02);

/* poultry */
poul02d=sum(chwi02, chwo02,  ctdog02,  ctsan02);

/* red and processed meat */
rpmeats02d = sum(bacon02, dog02, sbol02, oprom02, hamb02,hambl02, bmix02, pmain02*1.67, bmain02*1.67); 
/* ssb*/
ssb02d = sum(colas02, otsug02,obvcs02, punch02); 
/*** SWEETS & DESERTS  ***/
swtsvg02d= sum(choco02, cdyw02, cdywo02, cokff02, cokr02, cokh02, brwn02, donut02, cakr02, cakh02, pieh02,  srolr02, srolh02, icecr02); 
swtsnochoc02d=sum(cdyw02, cdywo02, cokff02, cokr02, cokh02, brwn02, donut02, cakr02, cakh02, pieh02,  srolr02, srolh02, icecr02);
run;


%hp04(keep=hbp04 db04 chol04 trig04 mi04 mi1d04 ang04 strk04 cabg04 pe04
              colc04  pros04  lymp04  ocan04  mel04 
              mar04 lalon04 asp04 mvt04 
              betab04 ace04 lasix04 thiaz04 calcb04 anthp04 
              stat04 mev04 zoc04 crest04 prav04 lip04 ostat04 ochrx04
              flag04 illpt04 sigmd04 colsc04 scopt04 vbld04 hxclc04 virt04 fecal04 cpolp04 bari04 diarr04 abdpn04 rout04
              sighx04 sig8604 sig8704 sig8904 sig9104 sig9304 sig9504 sig9704 sig9904 sig0104 sig0204 sig0304 sig0404 sigpt04
              colhx04 col8604 col8704 col8904 col9104 col9304 col9504 col9704 col9904 col0104 col0204 col0304 col0404 colpt04
           );

%hp06(keep=hbp06 hbpd06 db06 chol06 trig06 mi06 ang06 strk06 cabg06  pe06
              colc06  pros06  lymp06  ocan06  mel06 
              betab06 ace06 lasix06 thiaz06 calcb06 anthp06 
              asp06 mvt06 mar06 lalon06 
              stat06 mev06 zoc06 crest06 prav06 lip06 ostat06 ochrx06
              flag06 id colsc06 sigmd06 vbld06 hxclc06 virt06 cpol06 fecal06 diarr06 cpolp06 abdpn06 bari06 rout06 
mchoc06d /*milk chocolate*/
dchoc06d /*dark chocolate*/
cdyw06d /*Candy bar*/
cdywo06d/*Candy wo Chocolate*/
cokff06d/*Cookies, Fat-free/reduced fat*/
cokr06d/*Cookies, Commercial*/
cokh06d/*Cookies, Home Baked*/
brwn06d /*Brownies*/
donut06d/*Doughnuts*/
cake06d/*Cake, Homemade*/
piehr06d/*Pie, Ready-Made*/
jam06d  /*Jams*/
pbut06d /*Peanut butter*/
popcl06d /*Popcorn: Fat Free or light*/
popcr06d/*Popcorn: Regular*/
srlff06d/*Sweetroll Fat Free*/
srolr06d/*Sweetroll, Coffee Cake-Ready Made*/
srolh06d/*Sweetroll, Coffee Cake-Home Made*/
pretz06d/*Pretzels*/
pnut06d/*Peanuts*/
wnut06d/*Walnuts*/
onut06d/*Other Nuts*/
fries06d pot06d 
pchip06d /*chips*/
icecr06d /*Ice cream*/
plyog06d /*DAIRY - Yogurt-plain, lo-carb*/
flyog06d /*Yogurt*/

rpmeats06d whgrn06d ssb06d 
rais06d  prune06d  prunj06d  ban06d  cant06d   appl06d  aj06d  oran06d  ojca06d  ojreg06d  grfrj06d  othj06d  straw06d  blueb06d  peach06d
tom06d  toj06d  tosau06d  salsa06d  tofu06d  sbean06d  bean06d  peas06d soym06d brocc06d  cauli06d  cabb06d  bruss06d  rcar06d  ccar06d  corn06d  
mixv06d  yam06d  osqua06d  eggpl06d kale06d  cspin06d  rspin06d  ilett06d  rlett06d  cel06d  grpep06d oniog06d oniov06d
rwine06d wwine06d

cola06d  otsug06d  punch06d
bmix06d  bmain06d  pmain06d  oprom06d  bacon06d  dog06d  hamb06d  hambl06d  sbol06d
wcer06d rcer06d cerbr06d cer06d   ckcer06d  oat06d  dkbr06d  brice06d  otgrn06d  wgerm06d  popc06d  bran06d rye06d oatbr06d
whbr06d  wrice06d  pasta06d engl06d muff06d  pcake06d

tea06d htea06d  coff06d dcaf06d
swtsvg06d swtsnochoc06d frui06d vege06d fv06d
/* dairy */ skim06d mlk1206d whole06d sherb06d icecr06d plyog06d flyog06d cream06d cotch06d crmch06d otch06d
/* eggs */ eggom06d egg06d
/* fish */ ctuna06d bfsh06d dkfsh06d ofish06d shrim06d
/* fruit juices */ aj06d ojreg06d ojca06d prunj06d othj06d   
chwi06d chwo06d ctsan06d ctdog06d);

if cerbr06d in ( 0, 4, 11, 12, 13, 15, 16, 17, 19, 20, 21, 22, 23, 25, 29, 30, 31, 33, 35,37,38, 42, 43, 45, 46, 47,
59, 60,61,62,63,64,65,66,68,73,74,75,76,79,80,81,82,83, 84,85,86,88,95,96,97,98,
101,109,111, 112,113,114, 115,119, 123,124, 127, 128, 129,130,131, 132, 133, 134, 146, 147, 148,149,150,
151, 152,153, 154, 155, 157,158, 164,165, 166, 167,170,173, 176, 177,182, 184, 185, 186, 188)

   then rcer06d=cer06d;
   else wcer06d=cer06d;
   if rcer06d=. then rcer06d=0;
   if wcer06d=. then wcer06d=0;

array food mchoc06d dchoc06d cdyw06d cdywo06d cokff06d cokr06d cokh06d brwn06d donut06d cake06d piehr06d srlff06d srolr06d srolh06d jam06d pbut06d  
           popcl06d popcr06d srlff06d srolr06d srolh06d pretz06d pnut06d wnut06d onut06d fries06d pot06d  pchip06d icecr06d plyog06d flyog06d
           bmix06d  bmain06d  pmain06d  oprom06d  bacon06d  dog06d  hamb06d  hambl06d  sbol06d
           wcer06d  rcer06d ckcer06d  oat06d  dkbr06d  brice06d  otgrn06d  wgerm06d  popc06d  bran06d rye06d oatbr06d
           whbr06d  wrice06d  pasta06d engl06d muff06d  pcake06d
           rais06d  prune06d  prunj06d  ban06d  cant06d   appl06d  aj06d  oran06d  ojca06d  ojreg06d  grfrj06d  othj06d  straw06d  blueb06d  peach06d
           tom06d  toj06d  tosau06d  salsa06d  tofu06d  sbean06d  bean06d  peas06d soym06d brocc06d  cauli06d  cabb06d  bruss06d  rcar06d  ccar06d  corn06d  
           mixv06d  yam06d  osqua06d  eggpl06d kale06d  cspin06d  rspin06d  ilett06d  rlett06d  cel06d  grpep06d oniog06d oniov06d
           rwine06d wwine06d
           cola06d  otsug06d  punch06d tea06d htea06d coff06d dcaf06d
           skim06d mlk1206d whole06d sherb06d icecr06d plyog06d flyog06d cream06d cotch06d crmch06d otch06d
            eggom06d egg06d
            ctuna06d bfsh06d dkfsh06d ofish06d shrim06d
            aj06d ojreg06d ojca06d prunj06d othj06d   
chwi06d chwo06d ctsan06d ctdog06d;

    do over food; /* follow the same definition about continuous servings in NHS and NHS2 */
        if      food=0 then food=0;       /*never or <1/month*/
        else if food=1 then food=0.07;    /*1-3/month*/
        else if food=2 then food=0.14;     /*1/week*/
        else if food=3 then food=0.43;     /*2-4/week*/
        else if food=4 then food=0.79;   /*5-6/week*/
        else if food=5 then food=1;       /*1/day*/
        else if food=6 then food=2.5;     /*2-3/day*/
        else if food=7 then food=4.5;     /*4-5/day*/
        else if food=8 then food=6;     /*6+ per day*/
        else if food=9 then food=.;       /*missing*/
    end;



/* take sum of caffeinated and decaffeinated tea & coffee */
tea06d  = SUM(tea06d, htea06d);
coff06d  = sum(coff06d, dcaf06d);
/** whole grains **/
whgrn06d=sum(wcer06d  , dkbr06d, oat06d ,  ckcer06d,  rye06d ,  brice06d,  oatbr06d, popc06d, bran06d, wgerm06d, otgrn06d);
/** refined grains **/
rgrain06d  = sum(rcer06d, whbr06d ,  wrice06d,  pasta06d,  engl06d ,  muff06d ,  pcake06d);

/* fruit and vegetable*/
frui06d=sum(rais06d,  prune06d, ban06d , cant06d , appl06d, oran06d , straw06d, blueb06d, peach06d);
vege06d=sum( tom06d,  toj06d,   tosau06d,  salsa06d,   tofu06d, sbean06d,  bean06d,    peas06d,  soym06d,
                brocc06d,cauli06d,    cabb06d,  bruss06d,   rcar06d,  ccar06d,  /*corn06d,*/    mixv06d,  
                 yam06d,osqua06d,   eggpl06d,   kale06d,  cspin06d, rspin06d, ilett06d,   rlett06d,  
                  cel06d,grpep06d,   oniog06d,  oniov06d    );
fv06d=sum(frui06d, vege06d);

/* potatoes */
pota06d=sum(fries06d, pot06d, pchip06d);
/*legumes */
leg06d=sum(sbean06d, tofu06d, bean06d, peas06d, soym06d);
/* fruit juices */
fj06d=sum(aj06d, ojreg06d, ojca06d, prunj06d, othj06d  );

/* vegetables minus potatoes and legumes */
veg_noleg06d=vege06d-leg06d;

/* dairy and fish */
dairy06d=sum(skim06d, mlk1206d, whole06d, sherb06d, icecr06d, plyog06d, flyog06d, cream06d, cotch06d, crmch06d, otch06d);
fish06d=sum(ctuna06d, bfsh06d, dkfsh06d, ofish06d, shrim06d);
egg06d=sum(egg06d, eggom06d);

/* poultry */
poul06d=sum(chwi06d, chwo06d, ctsan06d, ctdog06d);

/* red and processed meat */
rpmeats06d=sum(bmix06d, bmain06d*1.67, pmain06d*1.67, oprom06d, bacon06d, dog06d, hamb06d, hambl06d, sbol06d);
/* ssb*/
ssb06d= sum(cola06d, otsug06d, punch06d);

/*** SWEETS & DESERTS  ***/
swtsvg06d= sum(mchoc06d, dchoc06d, cdyw06d, cdywo06d, cokff06d, cokh06d, cokr06d, brwn06d, donut06d, cake06d, piehr06d, srlff06d, srolr06d, srolh06d, jam06d);
swtsnochoc06d= sum(cdyw06d, cdywo06d, cokff06d, cokh06d, cokr06d, brwn06d, donut06d, cake06d, piehr06d, srlff06d, srolr06d, srolh06d, jam06d);
run;

%hp08(keep=hbp08 db08 chol08 trig08 mi08 ang08 strk08 cabg08 pe08
              colc08  pros08  lymp08  ocan08  mel08
              sclc108 pmel08 smel08 mpnc08 fpnc08 spnc08 bpnc08 cafh08
              asp08 mvt08  mar08 lalon08 
              betab08 ace08 lasix08 thiaz08 calcb08 arb08 anthp08 
              stat08 mev08 zoc08 crest08 prav08 lip08 ostat08 ochrx08
              flag08 id colsc08 sigmd08 vbld08 hxclc08 virt08 cpol08 fecal08 diarr08 cpolp08 abdpn08 bari08 rout08 ctcol08 );
              if sclc108=1 or pmel08=1 or smel08=1 or mpnc08=1 or fpnc08=1 or spnc08=1 or bpnc08=1 then cafh08=1; 
              else cafh08=0;
              run;


%hp10(keep=hbp10 db10  chol10  trig10  cabg10 ang10 mi10 strk10 pe10 
              colc10  pros10  blad10  mel10   lymp10  ocan10  /* bladder cancer new */
              betab10 calcb10 ace10   arb10   thiaz10 lasix10 anthp10 
              stat10  mev10   zoc10   crest10 prav10  lip10   ostat10 ochrx10 
              mar10  lalon10   mvt10 asp10
              flag10 id colsc10 sigmd10 vbld10 hxclc10 virt10 cpol10 fecal10 diarr10 cpolp10 abdpn10 bari10 rout10 ctcol10 
mchoc10d /*milk chocolate*/
dchoc10d /*dark chocolate*/
cdyw10d /*Candy bar*/
cdywo10d/*Candy wo Chocolate*/
cokff10d/*Cookies, Fat-free/reduced fat*/
cokr10d/*Cookies, Commercial*/
cokh10d/*Cookies, Home Baked*/
donut10d/*Doughnuts*/
cake10d /*Cake, Homemade*/
piehr10d/*Pie, Ready-Made*/
engl10d /* Bagels/English muffins/rolls */
muff10d /* muffins or biscuits */
jam10d  /*Jams, Jellies*/
pbut10d /*P.bu*/
popcl10d/*Popcorn, fat free or light*/
popcr10d/*Popcorn: Regular*/
sroll10d/*Sweetroll Fat Free*/
bbars10d/*Breakfast Bars*/
ebars10d/*Energy Bars*/
hpbars10d/*Protein Bars*/
pretz10d/*Pretzels*/
fd2slm10d/*Slimfast*/
fd3en10d/*Ensure*/
fd3enp10d/*Ensure plus*/
fd3mdf10d /*Mixed Dried Fruit*/
pnut10d/*Peanuts*/
wnut10d/*Walnuts*/
onut10d/*Other Nuts*/
fries10d pot10d 
pchip10d /*chips*/
icecr10d /*Ice cream*/
plyog10d /*DAIRY - Yogurt-plain, lo-carb*/
artyog10d /*Yogurt, artificially sweetened*/
sweyog10d /*Yogurt, sweetened*/

rpmeats10d frui10d vege10d fv10d whgrn10d ssb10d  
bmix10d  bmain10d  pmain10d  oprom10d  bacon10d  dog10d  hamb10d  hambl10d  sbol10d
rais10d  prune10d  prunj10d  ban10d  cant10d   appl10d  aj10d  oran10d  ojca10d  oj10d  grfr10d  othj10d  straw10d  blueb10d  peach10d
tom10d  toj10d  tosau10d  salsa10d  tofu10d  sbean10d  bean10d  peas10d  brocc10d  cauli10d  cabb10d  bruss10d  rcar10d  ccar10d  corn10d  
mixv10d  yam10d  osqua10d  eggpl10d kale10d  cspin10d  rspin10d  ilett10d  rlett10d  cel10d  grpep10d  oniov10d
cerbr10d rcer10d wcer10d cer10d  ckcer10d  oat10d  dkbr10d  brice10d  otgrn10d  wgerm10d  popc10d  bran10d
whbr10d  wrice10d  pasta10d  engl10d  muff10d  pcake10d
olcnc10d  cbvcfs10d   otsug10d   punch10d
tea10d htea10d  coff10d dcaf10d
rwine10d wwine10d

swtsvg10d swtsnochoc10d 
soym10d
/* dairy */ skim10d mlk1210d whole10d sherb10d icecr10d plyog10d artyog10d sweyog10d cream10d cotch10d crmch10d otch10d dcofdk10d /*yoghurt has three categories on FFQ 2010, dairy coffee drink is a new item on FFQ 2010 so it is added here (not to tea-coffee group)*/
/* eggs */ eggom10d egg10d
/* fish */ ctuna10d bfsh10d dkfsh10d ofish10d shrim10d
/* fruit juices */ aj10d ojreg10d ojca10d prunj10d othj10d   
chwi10d chwo10d ctsan10d ctdog10d);       

if cerbr10d in ( 0, 4, 11, 12, 13, 15, 16, 17, 19, 20, 21, 22, 23, 25, 29, 30, 31, 33, 35,37,38, 42, 43, 45, 46, 47,
59, 60,61,62,63,64,65,66,68,73,74,75,76,79,80,81,82,83, 84,85,86,88,95,96,97,98,
101,109,111, 112,113,114, 115,119, 123,124, 127, 128, 129,130,131, 132, 133, 134, 146, 147, 148,149,150,
151, 152,153, 154, 155, 157,158, 164,165, 166, 167,170,173, 176, 177,182, 184, 185, 186, 188)

   then rcer10d=cer10d;
   else wcer10d=cer10d;
   if rcer10d=. then rcer10d=0;
   if wcer10d=. then wcer10d=0;

array food mchoc10d dchoc10d cdyw10d cdywo10d cokff10d cokr10d cokh10d donut10d cake10d piehr10d jam10d pbut10d
               popcl10d popcr10d sroll10d bbars10d ebars10d hpbars10d pretz10d fd2slm10d fd3en10d fd3enp10d fd3mdf10d
               pnut10d  wnut10d  onut10d  fries10d pot10d pchip10d icecr10d plyog10d  artyog10d sweyog10d 
               bmix10d  bmain10d  pmain10d  oprom10d  bacon10d  dog10d  hamb10d  hambl10d  sbol10d
               rais10d  prune10d  prunj10d  ban10d  cant10d   appl10d  aj10d  oran10d  ojca10d  oj10d  grfr10d  othj10d  straw10d  blueb10d  peach10d
               tom10d  toj10d  tosau10d  salsa10d  tofu10d  sbean10d  bean10d  peas10d  brocc10d  cauli10d  cabb10d  bruss10d  rcar10d  ccar10d  corn10d  
               mixv10d  yam10d  osqua10d  eggpl10d kale10d  cspin10d  rspin10d  ilett10d  rlett10d  cel10d  grpep10d  oniov10d
               cer10d  rcer10d wcer10d ckcer10d  oat10d  dkbr10d  brice10d  otgrn10d  wgerm10d  popc10d  bran10d
               whbr10d  wrice10d  pasta10d  engl10d  muff10d  pcake10d
               ola10d   otsug10d   punch10d  tea10d htea10d  coff10d dcaf10d rwine10d wwine10d
               soym10d
/* dairy */ skim10d mlk1210d whole10d sherb10d icecr10d flyog10d plyog10d artyog10d sweyog10d cream10d cotch10d crmch10d otch10d dcofdk10d /*yoghurt has three categories on FFQ 2010, dairy coffee drink is a new item on FFQ 2010 so it is added here (not to tea-coffee group)*/
/* eggs */ eggom10d egg10d
/* fish */ ctuna10d bfsh10d dkfsh10d ofish10d shrim10d
/* fruit juices */ aj10d ojreg10d ojca10d prunj10d othj10d   
chwi10d chwo10d ctsan10d ctdog10d

;

    do over food; /* follow the same definition about continuous servings in NHS and NHS2 */
        if      food=0 then food=0;       /*never or <1/month*/
        else if food=1 then food=0.07;    /*1-3/month*/
        else if food=2 then food=0.14;     /*1/week*/
        else if food=3 then food=0.43;     /*2-4/week*/
        else if food=4 then food=0.79;   /*5-6/week*/
        else if food=5 then food=1;       /*1/day*/
        else if food=6 then food=2.5;     /*2-3/day*/
        else if food=7 then food=4.5;     /*4-5/day*/
        else if food=8 then food=6;     /*6+ per day*/
        else if food=9 then food=.;       /*missing*/
    end;


/* take sum of caffeinated and decaffeinated tea & coffee */
tea10d  = SUM(tea10d, htea10d);
coff10d  = sum(coff10d, dcaf10d);

/* red and processed meat */
rpmeats10d=sum(bmix10d, bmain10d*1.67, pmain10d*1.67, oprom10d, bacon10d, dog10d, hamb10d, hambl10d, sbol10d);
/* fruit and vegetable*/
frui10d=sum(rais10d, prune10d, prunj10d, ban10d, cant10d , appl10d, aj10d, oran10d, ojca10d, oj10d, grfr10d, othj10d, straw10d, blueb10d, peach10d);
vege10d=sum(tom10d, toj10d, tosau10d, salsa10d, tofu10d, sbean10d, bean10d, peas10d, soym10d, brocc10d, cauli10d, cabb10d, bruss10d, rcar10d, ccar10d, corn10d, mixv10d, yam10d, osqua10d, eggpl10d, 
        kale10d, cspin10d, rspin10d, ilett10d, rlett10d, cel10d, grpep10d, oniov10d);
fv10d=sum(frui10d, vege10d);

/* potatoes */
pota10d=sum(fries10d, pot10d, pchip10d);
/*legumes */
leg10d=sum(sbean10d, tofu10d, bean10d, peas10d, soym10d);
/* fruit juices */
fj10d=sum(aj10d, ojreg10d, ojca10d, prunj10d, othj10d  );

/* vegetables minus potatoes and legumes */
veg_noleg10d=vege10d-leg10d;

/* dairy and fish */
dairy10d=sum(skim10d, mlk1210d, whole10d, sherb10d, icecr10d, plyog10d, flyog10d, cream10d, cotch10d, crmch10d, otch10d);
fish10d=sum(ctuna10d, bfsh10d, dkfsh10d, ofish10d, shrim10d);
egg10d=sum(egg10d, eggom10d);

/* poultry */
poul10d=sum(chwi10d, chwo10d, ctsan10d, ctdog10d);


/** whole grains **/
whgrn10d=sum(wcer10d, ckcer10d, oat10d, dkbr10d, brice10d, otgrn10d, popc10d, wgerm10d, bran10d);
/** refined grains **/
rgrain10d  = sum(rcer10d, whbr10d ,  wrice10d,  pasta10d,  engl10d ,  muff10d ,  pcake10d);

/* ssb*/
ssb10d=sum(olcnc10d,  cbvcfs10d, otsug10d,  punch10d);

/*** SWEETS & DESERTS  ***/
swtsvg10d= sum(mchoc10d, dchoc10d, cdyw10d, cdywo10d, cokff10d, cokh10d, cokr10d, donut10d, cake10d, piehr10d, sroll10d, engl10d, muff10d, jam10d);
swtsnochoc10d= sum(cdyw10d, cdywo10d, cokff10d, cokh10d, cokr10d, donut10d, cake10d, piehr10d, sroll10d, engl10d, muff10d, jam10d);
run;

%hp12(keep= db12  hbp12 chol12 mi12  ang12 cabg12  strk12  
colc12  pros12   blad12  mel12 lymp12  ocan12    
betab12 calcb12 ace12   arb12   thiaz12 lasix12 anthp12
stat12 mev12  zoc12 crest12  prav12  lip12  ostat12 ochrx12
mvt12 asp12 mar12 lalon12
flag12 id colsc12 sigmd12 vbld12 hxclc12 virt12 cpol12 fecal12 diarr12 cpolp12 abdpn12 bari12 rout12 ctcol12
);   *renal failure no longer asked asfter 2010;



%hp14(keep=  db14  hbp14  chol14 mi14  ang14  cabg14  strk14 
colc14  pros14   blad14  mel14 lymp14  ocan14  
betab14 calcb14 ace14   arb14   thiaz14 lasix14 anthp14
stat14 mev14 prav14 zoc14 crest14  lip14  ostat14 ochrx14
mvt14 asp14 mar14 lalon14
flag14 id colsc14 sigmd14 vbld14 hxclc14 ctcol14 cpol14 fecal14 diarr14 cpolp14 abdpn14         rout14 ctcol14 
mchoc14d /*milk chocolate*/
dchoc14d /*dark chocolate*/
cdyw14d /*Candy bar*/
cdywo14d/*Candy wo Chocolate*/
cokff14d/*Cookies, Fat-free/reduced fat*/
cokr14d/*Cookies, Commercial*/
cokh14d/*Cookies, Home Baked*/
donut14d/*Doughnuts*/
cake14d /*Cake, Homemade*/
piehr14d/*Pie, Ready-Made*/
engl14d /* Bagels/English muffins/rolls */
muff14d /* muffins or biscuits */
jam14d  /*Jams, Jellies*/
pbut14d /*P.bu*/
popcl14d/*Popcorn, fat free or light*/
popcr14d/*Popcorn: Regular*/
sroll14d/*Sweetroll Fat Free*/
bbars14d/*Breakfast Bars*/
ebars14d/*Energy Bars*/
hpbars14d/*Protein Bars*/
pretz14d/*Pretzels*/
fd2slm14d/*Slimfast*/
fd3en14d/*Ensure*/
fd3enp14d/*Ensure plus*/
fd3mdf14d /*Mixed Dried Fruit*/
pnut14d/*Peanuts*/
wnut14d/*Walnuts*/
onut14d/*Other Nuts*/
fries14d pot14d 
pchip14d /*chips*/
icecr14d /*Ice cream*/
plyog14d /*DAIRY - Yogurt-plain, lo-carb*/
artyog14d /*Yogurt, artificially sweetened*/
sweyog14d /*Yogurt, sweetened*/

rpmeats14d frui14d vege14d fv14d whgrn14d ssb14d cerbr14d cer14d wcer14d rcer14d
bmix14d  bmain14d  pmain14d  oprom14d  bacon14d  dog14d  hamb14d  hambl14d  sbol14d
rais14d  prune14d  prunj14d  ban14d  cant14d   appl14d  aj14d  oran14d  ojca14d  oj14d  grfr14d  othj14d  straw14d  blueb14d  peach14d
tom14d  toj14d  tosau14d  salsa14d  tofu14d  sbean14d  bean14d  peas14d  brocc14d  cauli14d  cabb14d  bruss14d  rcar14d  ccar14d  corn14d  
mixv14d  yam14d  osqua14d  eggpl14d kale14d  cspin14d  rspin14d  ilett14d  rlett14d  cel14d  grpep14d  oniov14d
cer14d wcer14d rcer14d ckcer14d  oat14d  dkbr14d  brice14d  otgrn14d  wgerm14d  popc14d  bran14d
whbr14d  wrice14d  pasta14d  engl14d  muff14d  pcake14d
olcnc14d  cbvcfs14d   otsug14d   punch14d
tea14d htea14d  coff14d dcaf14d
rwine14d wwine14d
swtsvg14d swtsnochoc14d 
               soym14d
/* dairy */ skim14d mlk1214d whole14d sherb14d icecr14d  plyog14d artyog14d sweyog14d cream14d cotch14d crmch14d otch14d dcofdk14d /*yoghurt has three categories on FFQ 2014, dairy coffee drink is a new item on FFQ 2014 so it is added here (not to tea-coffee group)*/
/* eggs */ eggom14d egg14d
/* fish */ ctuna14d bfsh14d dkfsh14d ofish14d shrim14d
/* fruit juices */ aj14d ojreg14d ojca14d prunj14d othj14d   
chwi14d chwo14d ctsan14d ctdog14d);
if cerbr14d in ( 0, 4, 11, 12, 13, 15, 16, 17, 19, 20, 21, 22, 23, 25, 29, 30, 31, 33, 35,37,38, 42, 43, 45, 46, 47,
59, 60,61,62,63,64,65,66,68,73,74,75,76,79,80,81,82,83, 84,85,86,88,95,96,97,98,
101,109,111, 112,113,114, 115,119, 123,124, 127, 128, 129,130,131, 132, 133, 134, 146, 147, 148,149,150,
151, 152,153, 154, 155, 157,158, 164,165, 166, 167,170,173, 176, 177,182, 184, 185, 186, 188)

   then rcer14d=cer14d;
   else wcer14d=cer14d;
   if rcer14d=. then rcer14d=0;
   if wcer14d=. then wcer14d=0;


    array food mchoc14d dchoc14d cdyw14d cdywo14d cokff14d cokr14d cokh14d donut14d cake14d piehr14d engl14d muff14d jam14d pbut14d
               popcl14d popcr14d sroll14d bbars14d ebars14d hpbars14d pretz14d fd2slm14d fd3en14d fd3enp14d fd3mdf14d
               pnut14d  wnut14d  onut14d fries14d pot14d  pchip14d icecr14d plyog14d  artyog14d sweyog14d 
               bmix14d  bmain14d  pmain14d  oprom14d  bacon14d  dog14d  hamb14d  hambl14d  sbol14d
               rais14d  prune14d  prunj14d  ban14d  cant14d   appl14d  aj14d  oran14d  ojca14d  oj14d  grfr14d  othj14d  straw14d  blueb14d  peach14d
               tom14d  toj14d  tosau14d  salsa14d  tofu14d  sbean14d  bean14d  peas14d  brocc14d  cauli14d  cabb14d  bruss14d  rcar14d  ccar14d  corn14d  
               mixv14d  yam14d  osqua14d  eggpl14d kale14d  cspin14d  rspin14d  ilett14d  rlett14d  cel14d  grpep14d  oniov14d
               cer14d  wcer14d rcer14d ckcer14d  oat14d  dkbr14d  brice14d  otgrn14d  wgerm14d  popc14d  bran14d
               whbr14d  wrice14d  pasta14d  engl14d  muff14d  pcake14d
               ola14d   otsug14d   punch14d tea14d htea14d  coff14d dcaf14d rwine14d wwine14d
               soym14d
/* dairy */ skim14d mlk1214d whole14d sherb14d icecr14d  plyog14d artyog14d sweyog14d cream14d cotch14d crmch14d otch14d dcofdk14d /*yoghurt has three categories on FFQ 2014, dairy coffee drink is a new item on FFQ 2014 so it is added here (not to tea-coffee group)*/
/* eggs */ eggom14d egg14d
/* fish */ ctuna14d bfsh14d dkfsh14d ofish14d shrim14d
/* fruit juices */ aj14d ojreg14d ojca14d prunj14d othj14d   
chwi14d chwo14d ctsan14d ctdog14d

;

    do over food; /* follow the same definition about continuous servings in NHS and NHS2 */
        if      food=0 then food=0;       /*never or <1/month*/
        else if food=1 then food=0.07;    /*1-3/month*/
        else if food=2 then food=0.14;     /*1/week*/
        else if food=3 then food=0.43;     /*2-4/week*/
        else if food=4 then food=0.79;   /*5-6/week*/
        else if food=5 then food=1;       /*1/day*/
        else if food=6 then food=2.5;     /*2-3/day*/
        else if food=7 then food=4.5;     /*4-5/day*/
        else if food=8 then food=6;     /*6+ per day*/
        else if food=9 then food=.;       /*missing*/
    end;


/* take sum of caffeinated and decaffeinated tea & coffee */
tea14d  = SUM(tea14d, htea14d);
coff14d  = sum(coff14d, dcaf14d);

/* red and processed meat */
rpmeats14d=sum(bmix14d, bmain14d*1.67, pmain14d*1.67, oprom14d, bacon14d, dog14d, hamb14d, hambl14d, sbol14d);
/* fruit and vegetable*/
frui14d=sum(rais14d, prune14d, prunj14d, ban14d, cant14d , appl14d, aj14d, oran14d, ojca14d, oj14d, grfr14d, othj14d, straw14d, blueb14d, peach14d);
vege14d=sum(tom14d, toj14d, tosau14d, salsa14d, tofu14d, sbean14d, bean14d, peas14d, brocc14d, cauli14d, cabb14d, bruss14d, rcar14d, ccar14d, corn14d, mixv14d, yam14d, osqua14d, eggpl14d, 
        kale14d, cspin14d, rspin14d, ilett14d, rlett14d, cel14d, grpep14d, oniov14d);
fv14d=sum(frui14d, vege14d);

/* potatoes */
pota14d=sum(fries14d, pot14d, pchip14d);
/*legumes */
leg14d=sum(sbean14d, tofu14d, bean14d, peas14d, soym14d);
/* fruit juices */
fj14d=sum(aj14d, ojreg14d, ojca14d, prunj14d, othj14d  );

/* vegetables minus potatoes and legumes */
veg_noleg14d=vege14d-leg14d;

/* dairy and fish */
dairy14d=sum(skim14d, mlk1214d, whole14d, sherb14d, icecr14d, plyog14d, artyog14d, sweyog14d , cream14d, cotch14d, crmch14d, otch14d);
fish14d=sum(ctuna14d, bfsh14d, dkfsh14d, ofish14d, shrim14d);
egg14d=sum(egg14d, eggom14d);

/* poultry */
poul14d=sum(chwi14d, chwo14d, ctsan14d, ctdog14d);

/** whole grains **/
whgrn14d=sum(wcer14d, ckcer14d, oat14d, dkbr14d, brice14d, otgrn14d, wgerm14d, popc14d, bran14d);
/** refined grains **/
rgrain14d  = sum(rcer14d, whbr14d ,  wrice14d,  pasta14d,  engl14d ,  muff14d ,  pcake14d);

/* ssb*/
ssb14d=sum(olcnc14d,  cbvcfs14d, otsug14d,  punch14d);
/*** SWEETS & DESERTS  ***/
swtsvg14d = sum(mchoc14d, dchoc14d, cdyw14d, cdywo14d, cokff14d, cokr14d, cokh14d, donut14d, cake14d, piehr14d, engl14d, muff14d, jam14d);
swtsnochoc14d= sum(cdyw14d, cdywo14d, cokff14d, cokr14d, cokh14d, donut14d, cake14d, piehr14d, engl14d, muff14d, jam14d);
run;

%hp16(keep= db16    hbp16    chol16 mi16    ang16    cabg16  stk16 /* change */
colc16  pros16   blad16  mel16 lymp16  ocan16  
/* multivitamin emphysema/COPD NOT asked in 2016 */ 
bblock16 cblock16 aceinhb16  angio16  thiaz16 lasix16 antihy16
stat16 prava16 zocor16 crestor16  lipit16  othmed16 ochrx16
asp16 mar16 lalon16
flag16 id colsc16 sigmd16 vbld16 hxclc16 ctcol16 cpol16 fecal16 diarr16 cpolp16 abdpn16 baren16 routscr16 ctcol16 
); 


%hp18(keep= db18    hbp18    chol18 mi18    ang18    cabg18  strk18 
colc18  pros18   blad18  mel18 lymp18  ocan18   
mvt18  asp18  mar18 lalon18
bblock18 cblock18 aceinhb18  angio18  thiaz18 lasix18 antihy18
stat18 prava18 zocor18 crestor18  lipit18  othmed18 ochrx18
flag18 id colsc18 sigmd18 vbld18 hxclc18 ctcol18 cpol18 fecal18 diarr18 cpolp18 abdpn18 baren18 routscr18 ctcol18
mchoc18d /*milk chocolate*/
dchoc18d /*dark chocolate*/
cdyw18d /*Candy bar*/
cdywo18d/*Candy wo Chocolate*/
cokff18d/*Cookies, Fat-free/reduced fat*/
cokr18d/*Cookies, Commercial*/
cokh18d/*Cookies, Home Baked*/
donut18d/*Doughnuts*/
cake18d /*Cake, Homemade*/
piehr18d/*Pie, Ready-Made*/
engl18d /* Bagels/English muffins/rolls */
muff18d /* muffins or biscuits */
jam18d  /*Jams, Jellies*/
pbut18d /*P.bu*/
popcl18d/*Popcorn, fat free or light*/
popcr18d/*Popcorn: Regular*/
sroll18d/*Sweetroll Fat Free*/
bbars18d/*Breakfast Bars*/
ebars18d/*Energy Bars*/
hpbars18d/*Protein Bars*/
pretz18d/*Pretzels*/
fd2slm18d/*Slimfast*/
fd3en18d/*Ensure*/
fd3enp18d/*Ensure plus*/
fd3mdf18d /*Mixed Dried Fruit*/
pnut18d/*Peanuts*/
wnut18d/*Walnuts*/
onut18d/*Other Nuts*/

pchip18d /*chips*/
icecr18d /*Ice cream*/
plyog18d /*DAIRY - Yogurt-plain, lo-carb*/
artyog18d /*Yogurt, artificially sweetened*/
sweyog18d /*Yogurt, sweetened*/

rpmeats18d frui18d vege18d fv18d whgrn18d ssb18d cerbr18d cer18d wcer18d rcer18d
bmix18d  bmain18d  pmain18d  oprom18d  bacon18d  dog18d  hamb18d  hambl18d  sbol18d
rais18d  prune18d  prunj18d  ban18d  cant18d   appl18d  aj18d  oran18d  ojca18d  oj18d  grfr18d  othj18d  straw18d  blueb18d  peach18d
tom18d  toj18d  tosau18d  salsa18d  tofu18d  sbean18d  bean18d  peas18d  brocc18d  cauli18d  cabb18d  bruss18d  rcar18d  ccar18d  corn18d  
mixv18d  yam18d  osqua18d  eggpl18d kale18d  cspin18d  rspin18d  ilett18d  rlett18d  cel18d  grpep18d  oniov18d
cer18d  ckcer18d  oat18d  dkbr18d  brice18d  otgrn18d  wgerm18d  popc18d  bran18d
olcnc18d  cbvcfs18d   otsug18d   punch18d
tea18d htea18d  coff18d dcaf18d
rwine18d wwine18d
swtsvg18d
chwi18d chwo18d ctsan18d ctdog18d);

if cerbr18d in ( 0, 4, 11, 12, 13, 15, 16, 17, 19, 20, 21, 22, 23, 25, 29, 30, 31, 33, 35,37,38, 42, 43, 45, 46, 47,
59, 60,61,62,63,64,65,66,68,73,74,75,76,79,80,81,82,83, 84,85,86,88,95,96,97,98,
101,109,111, 112,113,114, 115,119, 123,124, 127, 128, 129,130,131, 132, 133, 134, 146, 147, 148,149,150,
151, 152,153, 154, 155, 157,158, 164,165, 166, 167,170,173, 176, 177,182, 184, 185, 186, 188)

   then rcer18d=cer18d;
   else wcer18d=cer18d;
   if rcer18d=. then rcer18d=0;
   if wcer18d=. then wcer18d=0;


    array food mchoc18d dchoc18d cdyw18d cdywo18d cokff18d cokr18d cokh18d donut18d cake18d piehr18d engl18d muff18d jam18d pbut18d
               popcl18d popcr18d sroll18d bbars18d ebars18d hpbars18d pretz18d fd2slm18d fd3en18d fd3enp18d fd3mdf18d
               pnut18d  wnut18d  onut18d  pchip18d icecr18d plyog18d  artyog18d sweyog18d 
               bmix18d  bmain18d  pmain18d  oprom18d  bacon18d  dog18d  hamb18d  hambl18d  sbol18d
               rais18d  prune18d  prunj18d  ban18d  cant18d   appl18d  aj18d  oran18d  ojca18d  oj18d  grfr18d  othj18d  straw18d  blueb18d  peach18d
               tom18d  toj18d  tosau18d  salsa18d  tofu18d  sbean18d  bean18d  peas18d  brocc18d  cauli18d  cabb18d  bruss18d  rcar18d  ccar18d  corn18d  
               mixv18d  yam18d  osqua18d  eggpl18d kale18d  cspin18d  rspin18d  ilett18d  rlett18d  cel18d  grpep18d  oniov18d
               cer18d  ckcer18d  oat18d  dkbr18d  brice18d  otgrn18d  wgerm18d  popc18d  bran18d
               ola18d   otsug18d   punch18d  tea18d htea18d  coff18d dcaf18d rwine18d wwine18d
               chwi18d chwo18d ctsan18d ctdog18d

;

    do over food; /* follow the same definition about continuous servings in NHS and NHS2 */
        if      food=0 then food=0;       /*never or <1/month*/
        else if food=1 then food=0.07;    /*1-3/month*/
        else if food=2 then food=0.14;     /*1/week*/
        else if food=3 then food=0.43;     /*2-4/week*/
        else if food=4 then food=0.79;   /*5-6/week*/
        else if food=5 then food=1;       /*1/day*/
        else if food=6 then food=2.5;     /*2-3/day*/
        else if food=7 then food=4.5;     /*4-5/day*/
        else if food=8 then food=6;     /*6+ per day*/
        else if food=9 then food=.;       /*missing*/
    end;

choco18=sum(mchoc18d, dchoc18d);

/* take sum of caffeinated and decaffeinated tea & coffee */
tea18d  = SUM(tea18d, htea18d);
coff18d  = sum(coff18d, dcaf18d);

/* red and processed meat */
rpmeats18d=sum(bmix18d, bmain18d*1.67, pmain18d*1.67, oprom18d, bacon18d, dog18d, hamb18d, hambl18d, sbol18d);
/* fruit and vegetable*/
frui18d=sum(rais18d, prune18d, prunj18d, ban18d, cant18d , appl18d, aj18d, oran18d, ojca18d, oj18d, grfr18d, othj18d, straw18d, blueb18d, peach18d);
vege18d=sum(tom18d, toj18d, tosau18d, salsa18d, tofu18d, sbean18d, bean18d, peas18d, brocc18d, cauli18d, cabb18d, bruss18d, rcar18d, ccar18d, corn18d, mixv18d, yam18d, osqua18d, eggpl18d, 
        kale18d, cspin18d, rspin18d, ilett18d, rlett18d, cel18d, grpep18d, oniov18d);
fv18d=sum(frui18d, vege18d);

/* poultry */
poul18d=sum(chwi18d, chwo18d, ctsan18d, ctdog18d);

/** whole grains **/
whgrn18d=sum(cer18d, ckcer18d, oat18d, dkbr18d, brice18d, otgrn18d, wgerm18d, popc18d, bran18d);
/* ssb*/
ssb18d=sum(olcnc18d,  cbvcfs18d, otsug18d,  punch18d);
/*** SWEETS & DESERTS  ***/
swtsvg18d = sum(mchoc18d, dchoc18d, cdyw18d, cdywo18d, cokff18d, cokr18d, cokh18d, donut18d, cake18d, piehr18d, engl18d, muff18d, jam18d);
swtsnochoc18d= sum(cdyw18d, cdywo18d, cokff18d, cokr18d, cokh18d, donut18d, cake18d, piehr18d, engl18d, muff18d, jam18d);
run;



%nses8818; 
run;

data nses8818;
set nses8818;

    *invert to make increases good*;
    array raw {*} ppov90_88 ppov90_90 ppov90_92 ppov90_94 ppov00_96 ppov00_98 ppov00_00 ppov00_02 
                  ppov00_04 ppov10_06 ppov10_08 ppov10_10 ppov10_12 ppov10_14 ppov10_16 ppov10_18;

    array der {*} dppov90_88 dppov90_90 dppov90_92 dppov90_94 dppov00_96 dppov00_98 dppov00_00 dppov00_02 
                  dppov00_04 dppov10_06 dppov10_08 dppov10_10 dppov10_12 dppov10_14 dppov10_16 dppov10_18;

    do j=1 to dim(raw);
        der{j}=1-raw{j};
    end;

   array inc {16} mdinc90_88 mdinc90_90 mdinc90_92 mdinc90_94 mdinc00_96 mdinc00_98 mdinc00_00 mdinc00_02 
                  mdinc00_04 mdinc10_06 mdinc10_08 mdinc10_10 mdinc10_12 mdinc10_14 mdinc10_16 mdinc10_18;
   array vhs {16} mdvhs90_88 mdvhs90_90 mdvhs90_92 mdvhs90_94 mdvhs00_96 mdvhs00_98 mdvhs00_00 mdvhs00_02 
                  mdvhs00_04 mdvhs10_06 mdvhs10_08 mdvhs10_10 mdvhs10_12 mdvhs10_14 mdvhs10_16 mdvhs10_18;
   array pcoll {16} pcolled90_88 pcolled90_90 pcolled90_92 pcolled90_94 pcolled00_96 pcolled00_98 pcolled00_00 pcolled00_02 
                    pcolled00_04 pcolled10_06 pcolled10_08 pcolled10_10 pcolled10_12 pcolled10_14 pcolled10_16 pcolled10_18;
   array pfint {16} pfaminterest90_88 pfaminterest90_90 pfaminterest90_92 pfaminterest90_94 pfaminterest00_96 pfaminterest00_98 pfaminterest00_00 pfaminterest00_02 
                    pfaminterest00_04 pfaminterest10_06 pfaminterest10_08 pfaminterest10_10 pfaminterest10_12 pfaminterest10_14 pfaminterest10_16 pfaminterest10_18;
   array pohse {16} pohse90_88 pohse90_90 pohse90_92 pohse90_94 pohse00_96 pohse00_98 pohse00_00 pohse00_02 
                    pohse00_04 pohse10_06 pohse10_08 pohse10_10 pohse10_12 pohse10_14 pohse10_16 pohse10_18;
   array pw {16}   pwht90_88 pwht90_90 pwht90_92 pwht90_94 pwht00_96 pwht00_98 pwht00_00 pwht00_02 
                   pwht00_04 pwht10_06 pwht10_08 pwht10_10 pwht10_12 pwht10_14 pwht10_16 pwht10_18;

   array zinc {16} zmdinc90_88 zmdinc90_90 zmdinc90_92 zmdinc90_94 zmdinc00_96 zmdinc00_98 zmdinc00_00 zmdinc00_02 
                   zmdinc00_04 zmdinc10_06 zmdinc10_08 zmdinc10_10 zmdinc10_12 zmdinc10_14 zmdinc10_16 zmdinc10_18;
   array zvhs {16} zmdvhs90_88 zmdvhs90_90 zmdvhs90_92 zmdvhs90_94 zmdvhs00_96 zmdvhs00_98 zmdvhs00_00 zmdvhs00_02 
                   zmdvhs00_04 zmdvhs10_06 zmdvhs10_08 zmdvhs10_10 zmdvhs10_12 zmdvhs10_14 zmdvhs10_16 zmdvhs10_18;
   array zpcoll {16} zpcolled90_88 zpcolled90_90 zpcolled90_92 zpcolled90_94 zpcolled00_96 zpcolled00_98 zpcolled00_00 zpcolled00_02 
                     zpcolled00_04 zpcolled10_06 zpcolled10_08 zpcolled10_10 zpcolled10_12 zpcolled10_14 zpcolled10_16 zpcolled10_18;
   array zpfint {16} zpfaminterest90_88 zpfaminterest90_90 zpfaminterest90_92 zpfaminterest90_94 zpfaminterest00_96 zpfaminterest00_98 zpfaminterest00_00 zpfaminterest00_02 
                     zpfaminterest00_04 zpfaminterest10_06 zpfaminterest10_08 zpfaminterest10_10 zpfaminterest10_12 zpfaminterest10_14 zpfaminterest10_16 zpfaminterest10_18;
   array zpohse {16} zpohse90_88 zpohse90_90 zpohse90_92 zpohse90_94 zpohse00_96 zpohse00_98 zpohse00_00 zpohse00_02 
                     zpohse00_04 zpohse10_06 zpohse10_08 zpohse10_10 zpohse10_12 zpohse10_14 zpohse10_16 zpohse10_18;
   array zpw {16}   zpwht90_88 zpwht90_90 zpwht90_92 zpwht90_94 zpwht00_96 zpwht00_98 zpwht00_00 zpwht00_02 
                    zpwht00_04 zpwht10_06 zpwht10_08 zpwht10_10 zpwht10_12 zpwht10_14 zpwht10_16 zpwht10_18;
   array zdppov {16} zdppov90_88 zdppov90_90 zdppov90_92 zdppov90_94 zdppov00_96 zdppov00_98 zdppov00_00 zdppov00_02 
                     zdppov00_04 zdppov10_06 zdppov10_08 zdppov10_10 zdppov10_12 zdppov10_14 zdppov10_16 zdppov10_18;

    do j=1 to dim(inc);
        zinc{j}  =inc{j};
        zvhs{j}  =vhs{j};
        zpcoll{j}=pcoll{j};
        zpfint{j}=pfint{j};
        zpohse{j}=pohse{j};
        zpw{j}   =pw{j};
        zdppov{j}=der{j};
    end;

        keep id  region: mdinc:  mdvhs:  pcolled:  pfaminterest:  pohse:  pwht:  dppov:
                         zmdinc: zmdvhs: zpcolled: zpfaminterest: zpohse: zpwht: zdppov:;
        run;

  /*
  proc means n nmiss min max mean median std; var mdinc:  mdvhs:  pcolled:  pfaminterest:  pohse:  pwht:  dppov:; run;
  proc freq;
    tables region:;
  run;
  */
        proc sort data=nses8818; by id; run;

* generate z-scores of variables and overwrites the values of the existing variable with z-scores*;
proc standard data=nses8818 mean=0 std=1 out=nses8818;
        var       zmdinc90_88 zmdinc90_90 zmdinc90_92 zmdinc90_94 zmdinc00_96 zmdinc00_98 zmdinc00_00 zmdinc00_02 
                  zmdinc00_04 zmdinc10_06 zmdinc10_08 zmdinc10_10 zmdinc10_12 zmdinc10_14 zmdinc10_16 zmdinc10_18
                  zmdvhs90_88 zmdvhs90_90 zmdvhs90_92 zmdvhs90_94 zmdvhs00_96 zmdvhs00_98 zmdvhs00_00 zmdvhs00_02 
                  zmdvhs00_04 zmdvhs10_06 zmdvhs10_08 zmdvhs10_10 zmdvhs10_12 zmdvhs10_14 zmdvhs10_16 zmdvhs10_18
                  zpcolled90_88 zpcolled90_90 zpcolled90_92 zpcolled90_94 zpcolled00_96 zpcolled00_98 zpcolled00_00 zpcolled00_02 
                  zpcolled00_04 zpcolled10_06 zpcolled10_08 zpcolled10_10 zpcolled10_12 zpcolled10_14 zpcolled10_16 zpcolled10_18
                  zpfaminterest90_88 zpfaminterest90_90 zpfaminterest90_92 zpfaminterest90_94 zpfaminterest00_96 zpfaminterest00_98 zpfaminterest00_00 zpfaminterest00_02 
                  zpfaminterest00_04 zpfaminterest10_06 zpfaminterest10_08 zpfaminterest10_10 zpfaminterest10_12 zpfaminterest10_14 zpfaminterest10_16 zpfaminterest10_18
                  zpohse90_88 zpohse90_90 zpohse90_92 zpohse90_94 zpohse00_96 zpohse00_98 zpohse00_00 zpohse00_02 
                  zpohse00_04 zpohse10_06 zpohse10_08 zpohse10_10 zpohse10_12 zpohse10_14 zpohse10_16 zpohse10_18
                  zpwht90_88 zpwht90_90 zpwht90_92 zpwht90_94 zpwht00_96 zpwht00_98 zpwht00_00 zpwht00_02 
                  zpwht00_04 zpwht10_06 zpwht10_08 zpwht10_10 zpwht10_12 zpwht10_14 zpwht10_16 zpwht10_18
                  zdppov90_88 zdppov90_90 zdppov90_92 zdppov90_94 zdppov00_96 zdppov00_98 zdppov00_00 zdppov00_02 
                  zdppov00_04 zdppov10_06 zdppov10_08 zdppov10_10 zdppov10_12 zdppov10_14 zdppov10_16 zdppov10_18;
run;

data nses8818;
   set nses8818;

   array zinc {16} zmdinc90_88 zmdinc90_90 zmdinc90_92 zmdinc90_94 zmdinc00_96 zmdinc00_98 zmdinc00_00 zmdinc00_02 
                   zmdinc00_04 zmdinc10_06 zmdinc10_08 zmdinc10_10 zmdinc10_12 zmdinc10_14 zmdinc10_16 zmdinc10_18;
   array zvhs {16} zmdvhs90_88 zmdvhs90_90 zmdvhs90_92 zmdvhs90_94 zmdvhs00_96 zmdvhs00_98 zmdvhs00_00 zmdvhs00_02 
                   zmdvhs00_04 zmdvhs10_06 zmdvhs10_08 zmdvhs10_10 zmdvhs10_12 zmdvhs10_14 zmdvhs10_16 zmdvhs10_18;
   array zpcoll {16} zpcolled90_88 zpcolled90_90 zpcolled90_92 zpcolled90_94 zpcolled00_96 zpcolled00_98 zpcolled00_00 zpcolled00_02 
                     zpcolled00_04 zpcolled10_06 zpcolled10_08 zpcolled10_10 zpcolled10_12 zpcolled10_14 zpcolled10_16 zpcolled10_18;
   array zpfint {16} zpfaminterest90_88 zpfaminterest90_90 zpfaminterest90_92 zpfaminterest90_94 zpfaminterest00_96 zpfaminterest00_98 zpfaminterest00_00 zpfaminterest00_02 
                     zpfaminterest00_04 zpfaminterest10_06 zpfaminterest10_08 zpfaminterest10_10 zpfaminterest10_12 zpfaminterest10_14 zpfaminterest10_16 zpfaminterest10_18;
   array zpohse {16} zpohse90_88 zpohse90_90 zpohse90_92 zpohse90_94 zpohse00_96 zpohse00_98 zpohse00_00 zpohse00_02 
                     zpohse00_04 zpohse10_06 zpohse10_08 zpohse10_10 zpohse10_12 zpohse10_14 zpohse10_16 zpohse10_18;
   array zpw {16}   zpwht90_88 zpwht90_90 zpwht90_92 zpwht90_94 zpwht00_96 zpwht00_98 zpwht00_00 zpwht00_02 
                    zpwht00_04 zpwht10_06 zpwht10_08 zpwht10_10 zpwht10_12 zpwht10_14 zpwht10_16 zpwht10_18;
   array zdppov {16} zdppov90_88 zdppov90_90 zdppov90_92 zdppov90_94 zdppov00_96 zdppov00_98 zdppov00_00 zdppov00_02 
                     zdppov00_04 zdppov10_06 zdppov10_08 zdppov10_10 zdppov10_12 zdppov10_14 zdppov10_16 zdppov10_18;

   array nses {16} nSES_88 nSES_90 nSES_92 nSES_94 nSES_96 nSES_98 nSES_00 nSES_02 
                   nSES_04 nSES_06 nSES_08 nSES_10 nSES_12 nSES_14 nSES_16 nSES_18;

   *generate ses score using the sum of the z-scores of the 7 census-tract factors*;
   *not doing cumulative average as where you live close to diagnosis may impact short-term diagnosis*;
   do j=1 to 16;
        nses{j}=zinc{j}+zvhs{j}+zpcoll{j}+zpfint{j}+zpohse{j}+zpw{j}+zdppov{j};
   end;
run;

data nses8818;
        set nses8818;
keep id nses:;
run;

proc sort; by id;


*------------------------------------4.Merge-----------------------------;
data base;
   merge hp_der(in=mst) hp_der_2 hp_dead cancer cvd16 diabetes hp_mets  
                        hp86 hp88 hp90 hp92 hp94 hp96 hp98 hp00 hp02 hp04 hp06 hp08 hp10 hp12 hp14 hp16 hp18 
                        h86_dt h90_dt h94_dt h98_dt
                        h86_nts(in=ffq) h90_nts h94_nts h98_nts h02_nts h06_nts h10_nts h14_nts h18_nts 
                        ahei86 ahei90 ahei94 ahei98 ahei02 ahei06 ahei10 ahei14 ahei18 
                        nses8818 end=_end_;
    by id;
    exrec=1; if  first.id and mst then exrec=0; 
    baseffq=1; if  ffq then baseffq=0; 


/*********************** covariates ***********************/

/* ethnicity data */
     if seuro86=1 then ethnic =1;
     else if scand86=1 then ethnic =1;
     else if ocauc86=1 then ethnic =1;
     else if afric86=1 then ethnic =2;
     else if asian86=1 then ethnic =3;
     else if oanc86=1  then ethnic =4;
     else ethnic=.;

     if ethnic>3 then eth3g=3;
     else             eth3g=ethnic;
     
     /*if eth3g=1 then white=1; else white=0; if ethnic=. then white=.; */
     if eth3g=1 then white=1; else white=0; 

/* antihypertensive */
   if (betab86 =1 or lasix86=1 or diur86 =1 or calcb86=1 or ald86=1) then htnrx86=1;
       else htnrx86=0;
   if (betab88 =1 or lasix88=1 or diur88 =1 or ccblo88=1 or ald88=1) then htnrx88=1;
       else htnrx88=0;
   if (betab90 =1 or lasix90=1 or diur90 =1 or ccblo90=1 or ald90=1) then htnrx90=1;
       else htnrx90=0;
   if (betab92 =1 or lasix92=1 or thiaz92=1 or ccblo92=1 or ald92=1) then htnrx92=1;
       else htnrx92=0;
   if (betab94 =1 or lasix94=1 or thiaz94=1 or ccblo94=1 or ald94=1) then htnrx94=1;
       else htnrx94=0;
   if (betab96 =1 or lasix96=1 or thiaz96=1 or ccblo96=1 or ald96=1) then htnrx96=1;
       else htnrx96=0;
   if (betab98 =1 or lasix98=1 or thiaz98=1 or ccblo98=1 or ald98=1) then htnrx98=1;
       else htnrx98=0;
   if (betab00 =1 or lasix00=1 or thiaz00=1 or ccblo00=1 or oanth00=1) then htnrx00=1;
       else htnrx00=0;
   if (betab02 =1 or lasix02=1 or thiaz02 =1 or calcb02=1 or anthp02=1) then htnrx02=1;
       else htnrx02=0;
   if (betab04 =1 or ace04=1 or lasix04=1 or thiaz04 =1 or calcb04=1 or anthp04=1) then htnrx04=1;
       else htnrx04=0; 
   if (betab06 =1 or ace06=1 or lasix06=1 or thiaz06 =1 or calcb06=1 or anthp06=1) then htnrx06=1;
       else htnrx06=0; 
   if (betab08 =1 or ace08=1 or lasix08=1 or thiaz08 =1 or calcb08=1 or arb08=1 or anthp08=1) then htnrx08=1;
       else htnrx08=0; 
   if (betab10 =1 or ace10=1 or lasix10=1 or thiaz10 =1 or calcb10=1 or arb10=1 or anthp10=1) then htnrx10=1;
       else htnrx10=0;
   if (betab12 =1 or ace12=1 or lasix12=1 or thiaz12 =1 or calcb12=1 or arb12=1 or anthp12=1) then htnrx12=1;
       else htnrx12=0;
   if (betab14 =1 or ace14=1 or lasix14=1 or thiaz14 =1 or calcb14=1 or arb14=1 or anthp14=1) then htnrx14=1;
       else htnrx14=0;
   
/* cholesterol lowering */
if antch86=1 then hchtx86=1;else hchtx86=0;
if chrx88=1 then hchtx88=1;else hchtx88=0;
if chrx90=1 then hchtx90=1;else hchtx90=0;
if chrx92=1 then hchtx92=1;else hchtx92=0;
if chrx94=1 then hchtx94=1;else hchtx94=0;
if chrx96=1 then hchtx96=1;else hchtx96=0;
if chrx98=1 then hchtx98=1;else hchtx98=0;
if chrx00=1 or ochrx00=1 then hchtx00=1;else hchtx00=0;
if stat02=1 or ochrx02=1 then hchtx02=1;else hchtx02=0;
if stat04=1 or mev04=1 or zoc04=1 or crest04=1 or prav04=1 or lip04=1 or ostat04=1 or ochrx04=1 then hchtx04=1; else hchtx04=0;
if stat06=1 or mev06=1 or zoc06=1 or crest06=1 or prav06=1 or lip06=1 or ostat06=1 or ochrx06=1 then hchtx06=1; else hchtx06=0;
if stat08=1 or mev08=1 or zoc08=1 or crest08=1 or prav08=1 or lip08=1 or ostat08=1 or ochrx08=1 then hchtx08=1; else hchtx08=0;
if stat10=1 or mev10=1 or zoc10=1 or crest10=1 or prav10=1 or lip10=1 or ostat10=1 or ochrx10=1 then hchtx10=1; else hchtx10=0;
if stat12=1 or mev12=1 or zoc12=1 or crest12=1 or prav12=1 or lip12=1 or ostat12=1 or ochrx12=1 then hchtx12=1; else hchtx12=0;
if stat14=1 or mev14=1 or zoc14=1 or crest14=1 or prav14=1 or lip14=1 or ostat14=1 or ochrx14=1 then hchtx14=1; else hchtx14=0;

/* Family history of diabetes, MI, stroke, and cancer */
   mifh=0; if mmi86=1 or fmi86=1 then mifh=1;
   dbfh=0; if dbfh87=1 or dbfh90=1 or dbfh92=1 then dbfh=1;
   parentaldb=0; if fdb3087=1 or mdb3087=1 or fdb90=1 or mdb90=1 or fdb92=1 or mdb92=1 then parentaldb=1;
   mdm=0; if mdb3087=1 or mdb90=1 or mdb92=1 then mdm=1;
   fdm=0; if fdb3087=1 or fdb90=1 or fdb92=1 then fdm=1;  
   cafh=0; if cafh86=1 or cafh90=1 or cafh92=1 or cafh96=1 or cafh08=1 then cafh=1;

/************** cancer: any confirmed cancer (from registry and self-reported) *****************/

  array canc_prev {*} can86  can88  can90  can92  can94  can96  can98  can00  can02  can04  can06  can08  can10  can12  can14  can16 can18;
  array canc_new  {*} cann86 cann88 cann90 cann92 cann94 cann96 cann98 cann00 cann02 cann04 cann06 cann08 cann10 cann12 cann14 can16 can18;

/* cancer: self-reported */
   can86=0;
   if pros86=1 or lymp86=1 or ocan86=1 or colc86=1 or mel86=1 then can86=1;
   can88=0;
   if pros88=1 or lymp88=1 or ocan88=1 or colc88=1 or mel88=1 then can88=1;
   can90=0;
   if pros90=1 or lymp90=1 or ocan90=1 or colc90=1 or mel90=1 then can90=1;
   can92=0;
   if pros92=1 or lymp92=1 or ocan92=1 or colc92=1 or mel92=1 then can92=1;
   can94=0;
   if pros94=1 or lymp94=1 or ocan94=1 or colc94=1 or mel94=1 then can94=1;
   can96=0;
   if pros96=1 or lymp96=1 or ocan96=1 or colc96=1 or mel96=1 then can96=1;
   can98=0;
   if pros98=1 or lymp98=1 or ocan98=1 or colc98=1 or mel98=1 then can98=1;
   can00=0;
   if pros00=1 or lymp00=1 or ocan00=1 or colc00=1 or mel00=1 then can00=1;
   can02=0;
   if pros02=1 or lymp02=1 or ocan02=1 or colc02=1 or mel02=1 then can02=1;
   can04=0;
   if pros04=1 or lymp04=1 or ocan04=1 or colc04=1 or mel04=1 then can04=1;
   can06=0;
   if pros06=1 or lymp06=1 or ocan06=1 or colc06=1 or mel06=1 then can06=1;
   can08=0;
   if pros08=1 or lymp08=1 or ocan08=1 or colc08=1 or mel08=1 then can08=1;
   can10=0;
   if pros10=1 or lymp10=1 or ocan10=1 or colc10=1 or mel10=1 or blad10=1 then can10=1;
   can12=0;
   if pros12=1 or lymp12=1 or ocan12=1 or colc12=1 or mel12=1 or blad12=1 then can12=1;
   can14=0;
   if pros14=1 or lymp14=1 or ocan14=1 or colc14=1 or mel14=1 or blad14=1 then can14=1;
   can16=0;
   if pros16=1 or lymp16=1 or ocan16=1 or colc16=1 or mel16=1 or blad16=1 then can16=1;
   can18=0;
   if pros18=1 or lymp18=1 or ocan18=1 or colc18=1 or mel18=1 or blad18=1 then can18=1;


   if dtdx_allca > 0 then first_dtdx_ca = dtdx_allca; /* ANY confirmed cancer */
   else if dtdxca > 0 then first_dtdx_ca = dtdxca;
   else first_dtdx_ca = .;

/* set first cancer */
    array qqyear {*} rtmnyr86 rtmnyr88 rtmnyr90 rtmnyr92 rtmnyr94 rtmnyr96 rtmnyr98 
                     rtmnyr00 rtmnyr02 rtmnyr04 rtmnyr06 rtmnyr08 rtmnyr10 rtmnyr12 rtmnyr14 rtmnyr16 rtmnyr18;

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
   array hbp     {*}  hbp86   hbp88   hbp90   hbp92   hbp94   hbp96   hbp98   hbp00    hbp02    hbp04    hbp06     hbp08   hbp10  hbp12  hbp14  hbp16  hbp18; 
   array chol    {*}  chol86  chol88  chol90  chol92  chol94  chol96  chol98  chol00   chol02   chol04   chol06    chol08  chol10 chol12 chol14 chol16 chol18;
   array diab    {*}  db86    db88    db90    db92    db94    db96    db98    db00     db02     db04     db06      db08    db10   db12   db14   db16   db18;
   array mi      {*}  mi86    mi88    mi90    mi92    mi94    mi96    mi98    mi00     mi02     mi04     mi06      mi08    mi10   mi12   mi14   mi16   mi18;
   array cabg    {*}  cabg86  cabg88  cabg90  cabg92  cabg94  cabg96  cabg98  cabg00   cabg02   cabg04   cabg06    cabg08  cabg10 cabg12 cabg14 cabg16 cabg18;
   array strk    {*}  str86   str88   str90   str92   str94   str96   str98   strk00   strk02   strk04   strk06    strk08  strk10 strk12 strk14 strk16 strk18;
   /*array canc_prev {*} can86  can88   can90   can92   can94   can96   can98   can00    can02    can04    can06     can08   can10  can12  can14 can16  can18;*/
   array chrodis {*}  chrodis86 chrodis88 chrodis90 chrodis92 chrodis94 chrodis96 chrodis98 chrodis00 chrodis02 
                      chrodis04 chrodis06 chrodis08 chrodis10 chrodis12 chrodis14 chrodis16 chrodis18; /* for stop update */
   
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
        if canc_prev {i} ne 1  then canc_prev {i} = 0;
        if hbp  {i} ne 1    then hbp  {i} = 0;
        if chol {i} ne 1    then chol {i} = 0;
        if diab {i} ne 1    then diab {i} = 0;
        if mi   {i} ne 1    then mi   {i} = 0;
        if cabg {i} ne 1    then cabg {i} = 0;
        if strk {i} ne 1    then strk {i} = 0;
        if mi(i) =1 or strk(i)=1 or cabg(i)=1 or diab(i)=1 or canc_prev(i)=1 then chrodis(i)=1; else chrodis(i)=0; 
    end;

    do i = 2 to dim(chrodis); if chrodis{i-1} = 1 then chrodis{i} = 1; 
    end; /* Make chrodis cumulative */


if nSES_88 = . then nSES_88m = 0.0296195; else nSES_88m = nSES_88;


/***************** Diet scores: carry forward + STOP updating after cancer/CVD *******************/

/* simple updated scores */   
array dis   {*}  chrodis86 chrodis90  chrodis94   chrodis98   chrodis02   chrodis06  chrodis10 chrodis14 chrodis18;
array calor {*}  calor86n  calor90n   calor94n    calor98n    calor02n    calor06n   calor10n  calor14n calor18n;   
array aheim {*}  ahei86   ahei90   ahei94   ahei98   ahei02   ahei06   ahei10  ahei14 ahei18;
array lcs   {*}  lcs86    lcs90    lcs94    lcs98    lcs02    lcs06    lcs10   lcs14  lcs18;
array alcs  {*}  alcs86   alcs90   alcs94   alcs98   alcs02   alcs06   alcs10  alcs14 alcs18;
array vlcs  {*}  vlcs86   vlcs90   vlcs94   vlcs98   vlcs02   vlcs06   vlcs10  vlcs14 vlcs18;
array ulcs  {*}  ulcs86   ulcs90   ulcs94   ulcs98   ulcs02   ulcs06   ulcs10  ulcs14 ulcs18;
array hlcs  {*}  hlcs86   hlcs90   hlcs94   hlcs98   hlcs02   hlcs06   hlcs10  hlcs14 hlcs18;
array lfs   {*}  lfs86    lfs90    lfs94    lfs98    lfs02    lfs06    lfs10   lfs14  lfs18;
array ulfs  {*}  ulfs86   ulfs90   ulfs94   ulfs98   ulfs02   ulfs06   ulfs10  ulfs14 ulfs18;
array hlfs  {*}  hlfs86   hlfs90   hlfs94   hlfs98   hlfs02   hlfs06   hlfs10  hlfs14 hlfs18;
array alfs  {*}  alfs86   alfs90   alfs94   alfs98   alfs02   alfs06   alfs10  alfs14 alfs18;
array vlfs  {*}  vlfs86   vlfs90   vlfs94   vlfs98   vlfs02   vlfs06   vlfs10  vlfs14 vlfs18;

array    carbp  {*}   carbo86p         carbo90p        carbo94p         carbo98p         carbo02p        carbo06p         carbo10p        carbo14p  carbo18p;
array    hcarbp  {*}   hcarbo86p        hcarbo90p       hcarbo94p        hcarbo98p        hcarbo02p       hcarbo06p        hcarbo10p       hcarbo14p hcarbo18p ;
array sumref_carbp {*}  sumref_carb86p   sumref_carb90p  sumref_carb94p   sumref_carb98p   sumref_carb02p  sumref_carb06p   sumref_carb10p  sumref_carb14p  sumref_carb18p;
array     protp  {*}   prot86p          prot90p         prot94p          prot98p          prot02p         prot06p          prot10p         prot14p  prot18p;
array    aprotp  {*}   aprot86p         aprot90p        aprot94p         aprot98p         aprot02p        aprot06p         aprot10p        aprot14p aprot18p ;
array    vprotp  {*}   vprot86p         vprot90p        vprot94p         vprot98p         vprot02p        vprot06p         vprot10p        vprot14p  vprot18p ;
array    tfatp   {*}   tfat86p          tfat90p         tfat94p          tfat98p          tfat02p         tfat06p          tfat10p         tfat14p  tfat18p;
array    afatp   {*}   afat86p          afat90p         afat94p          afat98p          afat02p         afat06p          afat10p         afat14p   afat18p;
array    vfatp   {*}   vfat86p          vfat90p         vfat94p          vfat98p          vfat02p         vfat06p          vfat10p         vfat14p   vfat18p;

array    carbs  {*}   carbo86s         carbo90s        carbo94s         carbo98s         carbo02s        carbo06s         carbo10s        carbo14s  carbo18s ;
array    hcarbs  {*}   hcarbo86s        hcarbo90s       hcarbo94s        hcarbo98s        hcarbo02s       hcarbo06s        hcarbo10s       hcarbo14s hcarbo18s ;
array sumref_carbs {*}  sumref_carb86s   sumref_carb90s  sumref_carb94s   sumref_carb98s   sumref_carb02s  sumref_carb06s   sumref_carb10s  sumref_carb14s  sumref_carb18s ;
array     prots  {*}   prot86s          prot90s         prot94s          prot98s          prot02s         prot06s          prot10s         prot14s  prot18s ;
array    aprots  {*}   aprot86s         aprot90s        aprot94s         aprot98s         aprot02s        aprot06s         aprot10s        aprot14s aprot18s ;
array    vprots  {*}   vprot86s         vprot90s        vprot94s         vprot98s         vprot02s        vprot06s         vprot10s        vprot14s  vprot18s ;
array    tfats   {*}   tfat86s          tfat90s         tfat94s          tfat98s          tfat02s         tfat06s          tfat10s         tfat14s  tfat18s ;
array    afats   {*}   afat86s          afat90s         afat94s          afat98s          afat02s         afat06s          afat10s         afat14s   afat18s ;
array    vfats   {*}   vfat86s          vfat90s         vfat94s          vfat98s          vfat02s         vfat06s          vfat10s         vfat14s   vfat18s ;

array    whgrn_o  {*}   whgrn86d       whgrn90d      whgrn94d       whgrn98d       whgrn02d      whgrn06d    whgrn10d  whgrn14d whgrn18d;
array    rgrain_o {*}   rgrain86d      rgrain90d     rgrain94d      rgrain98d      rgrain02d     rgrain06d   rgrain10d rgrain14d rgrain18d ;
array    frui_o   {*}   frui86d    frui90d    frui94d frui98d frui02d frui06d frui10d  frui14d frui18d ;
array     leg_o   {*}   leg86d        leg90d         leg94d        leg98d        leg02d       leg06d     leg10d  leg14d leg18d ;
array    pota_o   {*}   pota86d       pota90d        pota94d       pota98d       pota02d      pota06d    pota10d  pota14d pota18d ;
array  veg_noleg_o{*}   veg_noleg86d   veg_noleg90d  veg_noleg94d       veg_noleg98d       veg_noleg02d      veg_noleg06d    veg_noleg10d  veg_noleg14d veg_noleg18d ;
array    fish_o   {*}   fish86d        fish90d       fish94d        fish98d        fish02d       fish06d     fish10d fish14d fish18d ;
array   rpmeats_o {*}   rpmeats86d     rpmeats90d    rpmeats94d        rpmeats98d        rpmeats02d       rpmeats06d     rpmeats10d rpmeats14d rpmeats18d ;
array    ssb_o    {*}   ssb86d         ssb90d        ssb94d        ssb98d        ssb02d       ssb06d     ssb10d  ssb14d ssb18d ;
array    dairy_o   {*}   dairy86d        dairy90d       dairy94d        dairy98d        dairy02d       dairy06d     dairy10d dairy14d dairy18d;
array swtsnochoc_o {*}   swtsnochoc86d     swtsnochoc90d    swtsnochoc94d        swtsnochoc98d        swtsnochoc02d       swtsnochoc06d     swtsnochoc10d swtsnochoc14d swtsnochoc18d ;
array    poul_o    {*}   poul86d         poul90d        poul94d        poul98d        poul02d       poul06d     poul10d  poul14d poul18d;

array    whgrn_os  {*}   whgrn86ds       whgrn90ds      whgrn94ds       whgrn98ds       whgrn02ds      whgrn06ds    whgrn10ds whgrn14ds whgrn18ds;
array    rgrain_os {*}   rgrain86ds      rgrain90ds     rgrain94ds      rgrain98ds      rgrain02ds     rgrain06ds   rgrain10ds rgrain14ds rgrain18ds ;
array    frui_os   {*}   frui86ds    frui90ds    frui94ds frui98ds frui02ds frui06ds frui10ds frui14ds frui18ds ;
array     leg_os   {*}   leg86ds        leg90ds         leg94ds        leg98ds        leg02ds       leg06ds     leg10ds  leg14ds leg18ds ;
array    pota_os   {*}   pota86ds       pota90ds        pota94ds       pota98ds       pota02ds      pota06ds    pota10ds pota14ds pota18ds ;
array  veg_noleg_os{*}   veg_noleg86ds   veg_noleg90ds  veg_noleg94ds       veg_noleg98ds       veg_noleg02ds      veg_noleg06ds    veg_noleg10ds  veg_noleg14ds veg_noleg18ds ;
array    fish_os   {*}   fish86ds        fish90ds       fish94ds        fish98ds        fish02ds       fish06ds     fish10ds fish14ds fish18ds ;
array   rpmeats_os {*}   rpmeats86ds     rpmeats90ds    rpmeats94ds        rpmeats98ds        rpmeats02ds       rpmeats06ds     rpmeats10ds rpmeats14ds rpmeats18ds ;
array    ssb_os    {*}   ssb86ds         ssb90ds        ssb94ds        ssb98ds        ssb02ds       ssb06ds     ssb10ds ssb14ds ssb18ds ;
array    dairy_os   {*}   dairy86ds        dairy90ds       dairy94ds        dairy98ds        dairy02ds       dairy06ds     dairy10ds dairy14ds dairy18ds ;
array swtsnochoc_os {*}   swtsnochoc86ds     swtsnochoc90ds    swtsnochoc94ds        swtsnochoc98ds        swtsnochoc02ds       swtsnochoc06ds     swtsnochoc10ds  swtsnochoc14ds swtsnochoc18ds ;
array    poul_os    {*}   poul86ds         poul90ds        poul94ds        poul98ds        poul02ds       poul06ds     poul10ds poul14ds poul18ds ;

/* stop updated scores */
array calors {*}  calor86ns  calor90ns  calor94ns  calor98ns  calor02ns  calor06ns   calor10ns  calor14ns calor18ns;   
array aheims {*}  ahei86s  ahei90s  ahei94s  ahei98s  ahei02s  ahei06s ahei10s ahei14s ahei18s;
array lcss   {*}  lcs86s    lcs90s    lcs94s    lcs98s    lcs02s    lcs06s    lcs10s  lcs14s lcs18s;
array alcss  {*}  alcs86s   alcs90s   alcs94s   alcs98s   alcs02s   alcs06s   alcs10s  alcs14s alcs18s;
array vlcss  {*}  vlcs86s   vlcs90s   vlcs94s   vlcs98s   vlcs02s   vlcs06s   vlcs10s  vlcs14s vlcs18s;
array ulcss  {*}  ulcs86s   ulcs90s   ulcs94s   ulcs98s   ulcs02s   ulcs06s   ulcs10s  ulcs14s ulcs18s;
array hlcss  {*}  hlcs86s   hlcs90s   hlcs94s   hlcs98s   hlcs02s   hlcs06s   hlcs10s  hlcs14s hlcs18s;
array lfss   {*}  lfs86s    lfs90s    lfs94s    lfs98s    lfs02s    lfs06s    lfs10s  lfs14s lfs18s;
array ulfss  {*}  ulfs86s   ulfs90s   ulfs94s   ulfs98s   ulfs02s   ulfs06s   ulfs10s  ulfs14s ulfs18s;
array hlfss  {*}  hlfs86s   hlfs90s   hlfs94s   hlfs98s   hlfs02s   hlfs06s   hlfs10s  hlfs14s hlfs18s;
array alfss {*} alfs86s alfs90s alfs94s alfs98s alfs02s alfs06s alfs10s alfs14s alfs18s;
array vlfss {*} vlfs86s vlfs90s vlfs94s vlfs98s vlfs02s vlfs06s vlfs10s vlfs14s vlfs18s;

/* carry forward */
do i=2 to DIM(calor);
if calor{i}=.  then calor{i}=calor{i-1};
if aheim{i}=.  then aheim{i}=aheim{i-1};
if lcs{i}=.    then lcs{i}=lcs{i-1};
if alcs{i}=.   then alcs{i}=alcs{i-1};
if vlcs{i}=.   then vlcs{i}=vlcs{i-1};
if ulcs{i}=.   then ulcs{i}=ulcs{i-1};
if hlcs{i}=.   then hlcs{i}=hlcs{i-1};
if lfs{i}=.    then lfs{i}=lfs{i-1};
if ulfs{i}=.   then ulfs{i}=ulfs{i-1};
if hlfs{i}=.   then hlfs{i}=hlfs{i-1};
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
do i=1 to DIM(calor); /* initialize stop-updated scores */
  calors{i}=calor{i};  
  aheims{i}=aheim{i};
  lcss{i}=lcs{i};  
  alcss{i}=alcs{i};
  vlcss{i}=vlcs{i};
  ulcss{i}=ulcs{i};
  hlcss{i}=hlcs{i};
  lfss{i} =lfs{i};
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

do i=2 to DIM(calor); /* STOP updating scores after cancer OR CVD */
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

do i=2 to DIM(carbs); /* STOP updating scores after cancer OR CVD */
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


/* calculate cumulative average */
array nut{9,12}
calor86ns ahei86s lcs86s  alcs86s vlcs86s  ulcs86s  hlcs86s  lfs86s  ulfs86s  hlfs86s alfs86s  vlfs86s
calor90ns ahei90s lcs90s  alcs90s vlcs90s  ulcs90s  hlcs90s  lfs90s  ulfs90s  hlfs90s alfs90s  vlfs90s
calor94ns ahei94s lcs94s  alcs94s vlcs94s  ulcs94s  hlcs94s  lfs94s  ulfs94s  hlfs94s alfs94s  vlfs94s
calor98ns ahei98s lcs98s  alcs98s vlcs98s  ulcs98s  hlcs98s  lfs98s  ulfs98s  hlfs98s alfs98s  vlfs98s
calor02ns ahei02s lcs02s  alcs02s vlcs02s  ulcs02s  hlcs02s  lfs02s  ulfs02s  hlfs02s alfs02s  vlfs02s
calor06ns ahei06s lcs06s  alcs06s vlcs06s  ulcs06s  hlcs06s  lfs06s  ulfs06s  hlfs06s alfs06s  vlfs06s
calor10ns ahei10s lcs10s  alcs10s vlcs10s  ulcs10s  hlcs10s  lfs10s  ulfs10s  hlfs10s alfs10s  vlfs10s
calor14ns ahei14s lcs14s  alcs14s vlcs14s  ulcs14s  hlcs14s  lfs14s  ulfs14s  hlfs14s alfs14s  vlfs14s
calor18ns ahei18s lcs18s  alcs18s vlcs18s  ulcs18s  hlcs18s  lfs18s  ulfs18s  hlfs18s alfs18s  vlfs18s
;
array nutv{9,12}
calor86nsv  ahei86sv  lcs86sv alcs86sv vlcs86sv  ulcs86sv  hlcs86sv  lfs86sv  ulfs86sv  hlfs86sv  alfs86sv  vlfs86sv
calor90nsv  ahei90sv  lcs90sv alcs90sv vlcs90sv  ulcs90sv  hlcs90sv  lfs90sv  ulfs90sv  hlfs90sv  alfs90sv  vlfs90sv
calor94nsv  ahei94sv  lcs94sv alcs94sv vlcs94sv  ulcs94sv  hlcs94sv  lfs94sv  ulfs94sv  hlfs94sv  alfs94sv  vlfs94sv
calor98nsv  ahei98sv  lcs98sv alcs98sv vlcs98sv  ulcs98sv  hlcs98sv  lfs98sv  ulfs98sv  hlfs98sv  alfs98sv  vlfs98sv
calor02nsv  ahei02sv  lcs02sv alcs02sv vlcs02sv  ulcs02sv  hlcs02sv  lfs02sv  ulfs02sv  hlfs02sv  alfs02sv  vlfs02sv
calor06nsv  ahei06sv  lcs06sv alcs06sv vlcs06sv  ulcs06sv  hlcs06sv  lfs06sv  ulfs06sv  hlfs06sv  alfs06sv  vlfs06sv
calor10nsv  ahei10sv  lcs10sv alcs10sv vlcs10sv  ulcs10sv  hlcs10sv  lfs10sv  ulfs10sv  hlfs10sv  alfs10sv  vlfs10sv
calor14nsv  ahei14sv  lcs14sv alcs14sv vlcs14sv  ulcs14sv  hlcs14sv  lfs14sv  ulfs14sv  hlfs14sv  alfs14sv  vlfs14sv
calor18nsv  ahei18sv  lcs18sv alcs18sv vlcs18sv  ulcs18sv  hlcs18sv  lfs18sv  ulfs18sv  hlfs18sv  alfs18sv  vlfs18sv
;

/* calculate cumulative average intake */
do j=1 to 12;
  nutv{1,j}=nut{1,j};
 do i=2 to 9;
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
array foodlist{9,12}
whgrn86ds  rgrain86ds  frui86ds  leg86ds pota86ds  veg_noleg86ds  fish86ds  rpmeats86ds  ssb86ds  dairy86ds swtsnochoc86ds  poul86ds 
whgrn90ds  rgrain90ds  frui90ds  leg90ds pota90ds  veg_noleg90ds  fish90ds  rpmeats90ds  ssb90ds  dairy90ds swtsnochoc90ds  poul90ds 
whgrn94ds  rgrain94ds  frui94ds  leg94ds pota94ds  veg_noleg94ds  fish94ds  rpmeats94ds  ssb94ds  dairy94ds swtsnochoc94ds  poul94ds 
whgrn98ds  rgrain98ds  frui98ds  leg98ds pota98ds  veg_noleg98ds  fish98ds  rpmeats98ds  ssb98ds  dairy98ds swtsnochoc98ds  poul98ds 
whgrn02ds  rgrain02ds  frui02ds  leg02ds pota02ds  veg_noleg02ds  fish02ds  rpmeats02ds  ssb02ds  dairy02ds swtsnochoc02ds  poul02ds 
whgrn06ds  rgrain06ds  frui06ds  leg06ds pota06ds  veg_noleg06ds  fish06ds  rpmeats06ds  ssb06ds  dairy06ds swtsnochoc06ds  poul06ds 
whgrn10ds  rgrain10ds  frui10ds  leg10ds pota10ds  veg_noleg10ds  fish10ds  rpmeats10ds  ssb10ds  dairy10ds swtsnochoc10ds  poul10ds 
whgrn14ds  rgrain14ds  frui14ds  leg14ds pota14ds  veg_noleg14ds  fish14ds  rpmeats14ds  ssb14ds  dairy14ds swtsnochoc14ds  poul14ds 
whgrn18ds  rgrain18ds  frui18ds  leg18ds pota18ds  veg_noleg18ds  fish18ds  rpmeats18ds  ssb18ds  dairy18ds swtsnochoc18ds  poul18ds 
;
array foodlistv{9,12}
whgrn_s86sv  rgrain_s86sv  frui_s86sv  leg_s86sv pota_s86sv  veg_noleg_s86sv  fish_s86sv  rpmeats_s86sv  ssb_s86sv  dairy_s86sv swtsnochoc_s86sv  poul_s86sv
whgrn_s90sv  rgrain_s90sv  frui_s90sv  leg_s90sv pota_s90sv  veg_noleg_s90sv  fish_s90sv  rpmeats_s90sv  ssb_s90sv  dairy_s90sv swtsnochoc_s90sv  poul_s90sv
whgrn_s94sv  rgrain_s94sv  frui_s94sv  leg_s94sv pota_s94sv  veg_noleg_s94sv  fish_s94sv  rpmeats_s94sv  ssb_s94sv  dairy_s94sv swtsnochoc_s94sv  poul_s94sv
whgrn_s98sv  rgrain_s98sv  frui_s98sv  leg_s98sv pota_s98sv  veg_noleg_s98sv  fish_s98sv  rpmeats_s98sv  ssb_s98sv  dairy_s98sv swtsnochoc_s98sv  poul_s98sv
whgrn_s02sv  rgrain_s02sv  frui_s02sv  leg_s02sv pota_s02sv  veg_noleg_s02sv  fish_s02sv  rpmeats_s02sv  ssb_s02sv  dairy_s02sv swtsnochoc_s02sv  poul_s02sv
whgrn_s06sv  rgrain_s06sv  frui_s06sv  leg_s06sv pota_s06sv  veg_noleg_s06sv  fish_s06sv  rpmeats_s06sv  ssb_s06sv  dairy_s06sv swtsnochoc_s06sv  poul_s06sv
whgrn_s10sv  rgrain_s10sv  frui_s10sv  leg_s10sv pota_s10sv  veg_noleg_s10sv  fish_s10sv  rpmeats_s10sv  ssb_s10sv  dairy_s10sv swtsnochoc_s10sv  poul_s10sv
whgrn_s14sv  rgrain_s14sv  frui_s14sv  leg_s14sv pota_s14sv  veg_noleg_s14sv  fish_s14sv  rpmeats_s14sv  ssb_s14sv  dairy_s14sv swtsnochoc_s14sv  poul_s14sv
whgrn_s18sv  rgrain_s18sv  frui_s18sv  leg_s18sv pota_s18sv  veg_noleg_s18sv  fish_s18sv  rpmeats_s18sv  ssb_s18sv  dairy_s18sv swtsnochoc_s18sv  poul_s18sv
;

do j=1 to 12;
  foodlistv{1,j}=foodlist{1,j};
 do i=2 to 9;
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
array macron{9,9}
carbo86s hcarbo86s sumref_carb86s  prot86s aprot86s  vprot86s  tfat86s  afat86s  vfat86s
carbo90s hcarbo90s sumref_carb90s  prot90s aprot90s  vprot90s  tfat90s  afat90s  vfat90s
carbo94s hcarbo94s sumref_carb94s  prot94s aprot94s  vprot94s  tfat94s  afat94s  vfat94s
carbo98s hcarbo98s sumref_carb98s  prot98s aprot98s  vprot98s  tfat98s  afat98s  vfat98s
carbo02s hcarbo02s sumref_carb02s  prot02s aprot02s  vprot02s  tfat02s  afat02s  vfat02s
carbo06s hcarbo06s sumref_carb06s  prot06s aprot06s  vprot06s  tfat06s  afat06s  vfat06s
carbo10s hcarbo10s sumref_carb10s  prot10s aprot10s  vprot10s  tfat10s  afat10s  vfat10s
carbo14s hcarbo14s sumref_carb14s  prot14s aprot14s  vprot14s  tfat14s  afat14s  vfat14s
carbo18s hcarbo18s sumref_carb18s  prot18s aprot18s  vprot18s  tfat18s  afat18s  vfat18s
;
array macronv{9,9}
carbo86sv  hcarbo86sv  sumref_carb86sv prot86sv aprot86sv  vprot86sv  tfat86sv  afat86sv  vfat86sv
carbo90sv  hcarbo90sv  sumref_carb90sv prot90sv aprot90sv  vprot90sv  tfat90sv  afat90sv  vfat90sv
carbo94sv  hcarbo94sv  sumref_carb94sv prot94sv aprot94sv  vprot94sv  tfat94sv  afat94sv  vfat94sv
carbo98sv  hcarbo98sv  sumref_carb98sv prot98sv aprot98sv  vprot98sv  tfat98sv  afat98sv  vfat98sv
carbo02sv  hcarbo02sv  sumref_carb02sv prot02sv aprot02sv  vprot02sv  tfat02sv  afat02sv  vfat02sv
carbo06sv  hcarbo06sv  sumref_carb06sv prot06sv aprot06sv  vprot06sv  tfat06sv  afat06sv  vfat06sv
carbo10sv  hcarbo10sv  sumref_carb10sv prot10sv aprot10sv  vprot10sv  tfat10sv  afat10sv  vfat10sv
carbo14sv  hcarbo14sv  sumref_carb14sv prot14sv aprot14sv  vprot14sv  tfat14sv  afat14sv  vfat14sv
carbo18sv  hcarbo18sv  sumref_carb18sv prot18sv aprot18sv  vprot18sv  tfat18sv  afat18sv  vfat18sv
;

/* calculate cumulative average intake */
do j=1 to 9;
  macronv{1,j}=macron{1,j};
 do i=2 to 9;
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


/* endoscopy */ 
endo86=.; endo88=.; endo90=.; endo92=.; endo94=.; endo96=.;
endo98=.; endo00=.; endo02=.; endo04=.; endo06=.; endo08=.; endo10=.;endo12=.; 
endo14=.;endo16=.; endo18=.;

array endo {*} endo86  endo88  endo90  endo92  endo94  endo96  endo98  
               endo00  endo02  endo04  endo06  endo08  endo10  endo12  endo14  endo16 endo18; 
if sigx88=1 and colx88=1 then endo88=1;  /*endo no*/
else if sigsy88=1 or sigsc88=1 or colsy88=1 or colsc88=1 then endo88=2;  /*endo yes*/
if sigm90=1 then endo90=1;  /*endo no*/
else if sigm90=2 then endo90=2;  /*endo yes*/
if sig9092=1 then endo92=1;  /*endo no*/
else if sig9092=2 then endo92=2;  /*endo yes*/
else if sigm92=1 then endo92=1; 
if sigm94=1  then endo94=1;  else if sigm94=2 then endo94=2;  
if sig9496=1 then endo96=1;  else if sig9496=2 then endo96=2; 
if sigm98=1  then endo98=1;  else if sigm98=2 then endo98=2; 
if sigm00=1  then endo00=1; else if sigm00=2 then endo00=2; 
if sigmd02=1 then endo02=1; else if sigmd02=2 then endo02=2; 
if colsc04=2 or col0204=1 or col0304=1 or sigmd04=2 or sig0204=1 or sig0304=1 then endo04=2;  /*endo yes*/
else if colsc04=1 and sigmd04=1 then endo04=1;  /*endo no*/
if colsc06=2 or  sigmd06=2 then endo06=2;  /*endo yes*/
else if colsc06=1 and  sigmd06=1 then endo06=1;  /*endo no*/
if colsc08=2 or  sigmd08=2  or  ctcol08=2 then endo08=2;  /*endo yes*/
else if colsc08=1 and  sigmd08=1 and ctcol08=1 then endo08=1;  /*endo no*/
if colsc10=2 or  sigmd10=2  or  ctcol10=2 then endo10=2;  /*endo yes*/
else if colsc10=1 and  sigmd10=1 and ctcol10=1 then endo10=1;  /*endo no*/
if colsc12=2 or  sigmd12=2  or  ctcol12=2 then endo12=2;  /*endo yes*/
else if colsc12=1 and  sigmd12=1 and ctcol12=1 then endo12=1;  /*endo no*/
if colsc14=2 or  sigmd14=2  or  ctcol14=2 then endo14=2;  /*endo yes*/
else if colsc14=1 and  sigmd14=1 and ctcol14=1 then endo14=1;  /*endo no*/
if colsc16=2 or  sigmd16=2  or  ctcol16=2 then endo16=2;  /*endo yes*/
else if colsc16=1 and  sigmd16=1 and ctcol16=1 then endo16=1;  /*endo no*/
if colsc18=2 or  sigmd18=2  or  ctcol18=2 then endo18=2;  /*endo yes*/
else if colsc18=1 and  sigmd18=1 and ctcol18=1 then endo18=1;  /*endo no*/

*-----------*
 Use 1990 Qx
*-----------*;
/*endo86 indicate experience of endoscopy before 1986 (1=BEFORE 1986)*/
if sigm190 in (1) or sigmr90 in (1) then endo86=2; 
else if (sigm190 in (2,3,4,5,6)) then endo86=1;  /*endo no*/

/*endo88-endo90 were missing, 1990 info was used to fill the answer.
 We don't impute endo88 if endo90=yes becase person shouldn't take
 endoscopy twice within 4 years. Same for endo90*/

if ((endo90 ne 2) and (sigm190 in (2,3) or sigmr90 in (2,3)) and endo88=. )
   then endo88=2;  /*yes endo*/
else if (endo88=. and sigm190 in (4,5,6)) then endo88=1;  /*when first sigmo is later than 1988, endo88 should be no*/

if ((endo88 ne 2) and (sigm190 in (4,5,6) or sigmr90 in (4,5,6)) and endo90=. )
   then endo90=2; 
else if (endo90=. and sigm190 in (2,3)) then endo90=1;  /*she won't take endo twice within 5 years*/
else if (endo90=. and sigmr90 in (1,2,3)) then endo90=2;  /*if most recent endo in 86 or 87, she won't take in 88, 89*/

*-----------*
 Use 2004 Qx
*-----------*;
/*For only missing variables of the past Qx, we use 2004 info*/
/*We don't impute endo[i] if person already endo[i-1]=yes or endo[i+1]=yes*/
/*never*/
/*pass through means "he never took" because there's no "never" option in 2004 Qx in HPFS*/
if endo86=2 or endo88=2 or endo90=2 or endo92=2 or endo94=2 or
   endo96=2 or endo98=2 or endo00=2 or endo02=2 or endo04=2
   then atleast=2;  /*yes edo at least once in some year. They cannot be "Never endo"*/

/*we fill from the latest information*/
/*04 can fill endo04 as we assume they took endo in the midepoint between qxs
  (i.e., midpoint between 04qx and 06qx is June 05. Thus, endo in 04 is included in 04qx cycle*/
if (endo04 ne 1 and (sig0404=1 or col0404=1)) then endo04=2;  /*yes*/   
/*03 endo can fill endo04*/
if (endo02 ne 2 and endo06 ne 2 and endo04=. and (sig0304=1 or col0304=1)) then endo04=2;  /*yes*/
else if endo04=. and flag04 in (1,2) then endo04=1;  /*participated endo q but no answer*/
/*02 endo can fill endo04*/
if (endo02 ne 2 and endo06 ne 2 and endo04=. and (sig0204=1 or col0204=1)) then endo04=2;  /*yes*/
else if endo04=. and flag04 in (1,2) then endo04=1;  /*participated endo q but no answer*/
/*2000-2001 endo can fill endo02*/
if (endo00 ne 2 and endo04 ne 2 and endo02=. and (sig0104=1 or col0104=1)) then endo02=2;  /*yes*/
else if endo02=. and flag04 in (1,2) then endo02=1;  /*participated endo q but no answer*/
/*1998-1999 endo can fill endo00*/
if (endo98 ne 2 and endo02 ne 2 and endo00=. and (sig9904=1 or col9904=1)) then endo00=2;  /*yes*/
else if endo00=. and flag04 in (1,2) then endo00=1;  /*participated endo q but no answer*/
/*1996-1997 endo can fill endo98*/
if (endo96 ne 2 and endo00 ne 2 and endo98=. and (sig9704=1 or col9704=1)) then endo98=2;  /*yes*/
else if endo98=. and flag04 in (1,2) then endo98=1;  /*participated endo q but no answer*/
/*1994-1995 endo can fill endo96*/
if (endo94 ne 2 and endo98 ne 2 and endo96=. and (sig9504=1 or col9504=1)) then endo96=2;  /*yes*/
else if endo96=. and flag04 in (1,2) then endo96=1;  /*participated endo q but no answer*/
/*92-93 endo can fill endo94*/
if (endo92 ne 2 and endo96 ne 2 and endo94=. and (sig9304=1 or col9304=1)) then endo94=2;  /*yes*/
else if endo94=. and flag04 in (1,2) then endo94=1;  /*participated endo q but no answer*/
/*90-91 endo can fill endo92*/
if (endo90 ne 2 and endo94 ne 2 and endo92=. and (sig9104=1 or col9104=1)) then endo92=2;  /*yes*/
else if endo92=. and flag04 in (1,2) then endo92=1;  /*participated endo q but no answer*/
/*88-89 endo can fill endo90*/
if (endo88 ne 2 and endo92 ne 2 and endo90=. and (sig8904=1 or col8904=1)) then endo90=2;  /*yes*/
else if endo90=. and flag04 in (1,2) then endo90=1;  /*participated endo q but no answer*/
/*86-87 endo can fill endo88*/
if (endo90 ne 2 and endo88=. and (sig8704=1 or col8704=1)) then endo88=2;  /*yes*/
else if endo88=. and flag04 in (1,2) then endo88=1;  /*participated endo q but no answer*/
/*before 1986, we fill endo86*/
if (endo86=. and (sig8604=1 or col8604=1)) then endo86=2;  /*yes*/                        /*illpt04=1 is all path through*/
else if endo86=. and flag04 in (1,2) then endo86=1;  /*participated endo q but no answer*/

*****fill in if missing*******;
if endo88=. and flag88=1 then endo88=1;  /*those who participated in long Qx endo no*/
if endo90=. and flag90=1 then endo90=1;  /*those who participated in long Qx endo no*/
if endo92=. and flag92=1 then endo92=1;  /*those who participated in long Qx endo no*/
if endo94=. and flag94=1 then endo94=1;  /*those who participated in long Qx endo no*/
if endo96=. and flag96 in (1,2) then endo96=1;  /*those who participated in long and short Qx endo no*/
if endo98=. and flag98 in (1,2) then endo98=1;  /*those who participated in long and short Qx endo no*/
if endo00=. and flag00 in (1,2) then endo00=1;  /*those who participated in long and short Qx endo no*/
if endo02=. and flag02 in (1,2) then endo02=1;  /*those who participated in long and short Qx endo no*/
if endo04=. and flag04 in (1,2) then endo04=1;  /*those who participated in long and short Qx endo no*/
if endo06=. and flag06 in (1,2) then endo06=1;  /*those who participated in long and short Qx endo no*/
if endo08=. and flag08 in (1,2) then endo08=1;  /*those who participated in long and short Qx endo no*/
if endo10=. and flag10 in (1,2) then endo10=1;  /*those who participated in long and short Qx endo no*/
if endo12=. and flag12 in (1,2) then endo12=1;  /*those who participated in long and short Qx endo no*/
if endo14=. and flag14 in (1,2) then endo14=1;  /*those who participated in long and short Qx endo no*/
if endo16=. and flag16 in (1,2) then endo16=1;  /*those who participated in long and short Qx endo no*/
if endo18=. and flag18 in (1,2) then endo18=1;  /*those who participated in long and short Qx endo no*/

*-----------*
 Use 1990 Qx
*-----------*;
if endo86=. and flag90=1 then endo86=1;  /*endo no*/
if endo88=. and flag90=1 then endo88=1;  /*long qx*/
if endo90=. and flag90=1 then endo90=1;  /*long qx*/

*----------------------------------------*
 reason for endoscopy, create dichotomous
*----------------------------------------*;
   if (vbld90=1 or fecal90=1 or abdpn90=1 or diarr90=1) then reasoned86=1; /*symptom*/
   else if hxclc90=1 or rout90=1 then reasoned86=2;                 /*screeting or family hx*/

   if sigsy88=1 or colsy88=1 then reasoned88=1; /*symptom. NOTE, Precise reasons were not asked in 1980 Qx*/
   if sigsc88=1 or colsc88=1 then reasoned88=2; /*screening*/

   if (vbld90=1 or fecal90=1 or abdpn90=1 or diarr90=1) then reasoned90=1; /*symptom*/
   else if hxclc90=1 or rout90=1 then reasoned90=2;                 /*screeting or family hx*/

   if (vbld92=1 or fecal92=1 or abdpn92=1 or diarr92=1) then reasoned92=1; /*symptom*/
   else if hxclc92=1 or rout92=1 then reasoned92=2;                 /*screeting or family hx*/

   if (vbld94=1 or fecal94=1 or abdpn94=1 or diarr94=1) then reasoned94=1; /*symptom*/
   else if hxclc94=1 or rout94=1 then reasoned94=2;                 /*screeting or family hx*/

   if (vbld96=1 or fecal96=1 or abdpn96=1 or diarr96=1) then reasoned96=1; /*symptom*/
   else if hxclc96=1 or rout96=1 then reasoned96=2;                 /*screeting or family hx*/

   if (vbld98=1 or fecal98=1 or abdpn98=1 or diar98=1) then reasoned98=1; /*symptom*/
   else if hxclc98=1 or rout98=1 then reasoned98=2;                 /*screeting or family hx*/

   if (vbld00=1 or fecal00=1 or abdpn00=1 or diar00=1) then reasoned00=1; /*symptom*/
   else if hxclc00=1 or rout00=1 then reasoned00=2;                 /*screeting or family hx*/

   if (vbld02=1 or fecal02=1 or abdpn02=1 or diarr02=1) then reasoned02=1; /*symptom*/
   else if hxclc02=1 or rout02=1 then reasoned02=2;                 /*screeting or family hx*/

   if (vbld04=1 or fecal04=1 or abdpn04=1 or diarr04=1 or bari04=1) then reasoned04=1; /*symptom*/
   else if (hxclc04=1 or rout04=1 or cpolp04=1 or virt04=1) then reasoned04=2;                 /*screeting or family hx*/

   if (vbld06=1 or fecal06=1 or abdpn06=1 or diarr06=1 or bari06=1) then reasoned06=1; /*symptom*/
   else if (hxclc06=1 or rout06=1 or cpolp06=1 or virt06=1) then reasoned06=2;                 /*screeting or family hx*/

   if (vbld08=1 or fecal08=1 or abdpn08=1 or diarr08=1 or bari08=1) then reasoned08=1; /*symptom*/
   else if (hxclc08=1 or rout08=1 or cpolp08=1 or virt08=1) then reasoned08=2;                 /*screeting or family hx*/

   if (vbld10=1 or fecal10=1 or abdpn10=1 or diarr10=1 or bari10=1) then reasoned10=1; /*symptom*/
   else if (hxclc10=1 or rout10=1 or cpolp10=1 or virt10=1) then reasoned10=2;                 /*screeting or family hx*/

   if (vbld12=1 or fecal12=1 or abdpn12=1 or diarr12=1 or bari12=1) then reasoned12=1; /*symptom*/
   else if (hxclc12=1 or rout12=1 or cpolp12=1 or virt12=1) then reasoned12=2;                 /*screeting or family hx*/

   if (vbld14=1 or fecal14=1 or abdpn14=1 or diarr14=1  ) then reasoned14=1; /*symptom*/
   else if (hxclc14=1 or rout14=1 or cpolp14=1  ) then reasoned14=2;                 /*screeting or family hx*/
   if (vbld16=1 or fecal16=1 or abdpn16=1 or diarr16=1 or baren16=1) then reasoned16=1; /*symptom*/
   else if ( routscr16=1 or cpolp16=1  ) then reasoned16=2;                 /*screeting or family hx*/
   if (vbld18=1 or fecal18=1 or abdpn18=1 or diarr18=1 or baren18=1) then reasoned18=1; /*symptom*/
   else if (hxclc18=1 or routscr18=1 or cpolp18=1 ) then reasoned18=2;                 /*screeting or family hx*/


*------------------*
 Reason for endo (updated and carry forward)
*------------------*;
/*If endo=yes and reason=symptom, endo-reason status is yes-endo-w/symptom later on until
 endo-reason status changes, i.e. endo=yes w/sympton in 90, endo=no in 92,
 endo=yes w/screening in 94, her/his endoreason is symptom in 90, sympton in 92,
 screening in 94*/
array endolist {*}    endo86    endo88    endo90    endo92    endo94    endo96    endo98    endo00    endo02    endo04    endo06    endo08    endo10  endo12 endo14 endo16 endo18 ;
array reasonxx[*] reasoned88 reasoned90 reasoned92 reasoned94 reasoned96
                  reasoned98 reasoned00 reasoned02 reasoned04 reasoned06 reasoned08 reasoned10 reasoned12 reasoned14 reasoned16 reasoned18 ;   /*1=symptom,2=screening*/
array endoreason[*] endoreason88 endoreason90 endoreason92 endoreason94 endoreason96
                    endoreason98 endoreason00 endoreason02 endoreason04 endoreason06 endoreason08 endoreason10 endoreason12 endoreason14 endoreason16 endoreason18 ;/*1=symptom,2=screening*/
do i=1 to dim(endoreason);
   if endolist[i]=2 and reasonxx[i]=. then do; /*missing*/
      do j=i to dim(endoreason);
         endoreason[j]=.;
      end;
   end;
   /*overwrite if status changes*/
   if endolist[i]=2 and reasonxx[i]=1 then do; /*symptom*/
      do j=i to dim(endoreason);
         endoreason[j]=1;
      end;
   end;
   /*overwrite if status changes*/
   if endolist[i]=2 and reasonxx[i]=2 then do; /*screening*/
      do j=i to dim(endoreason);
         endoreason[j]=2;
      end;
   end;
end;


/***1. time varying, updated every 2 years, 1=no endo, 2=endo for symptoms, 3=screening endo***/
array endo3c [*] endo3c88 endo3c90 endo3c92 endo3c94 endo3c96 endo3c98 endo3c00 endo3c02 endo3c04 endo3c06 endo3c08 endo3c10 endo3c12 endo3c14 endo3c16 endo3c18 ;

do i=1 to dim(endo3c);
   if endolist[i]=1 then endo3c[i]=1; *no endo;
   else if endolist[i]=2 and reasonxx[i]=1 then endo3c[i]=2; *endo for symptoms;
   else if endolist[i]=2 and reasonxx[i]=2 then endo3c[i]=3; *endo for screening;
   else endo3c[i]=.; *missing;
end;


/***2. ever endo3c, carry forward***/
array eendo3c {*}  eendo3c88    eendo3c90    eendo3c92    eendo3c94    eendo3c96    eendo3c98    eendo3c00    eendo3c02    eendo3c04    eendo3c06    eendo3c08  eendo3c10 eendo3c12 eendo3c14  eendo3c16 eendo3c18 ;

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
      if eendo3c[j]=. and endo3c[j]=1 then eendo3c[j]=1; *no history of endoscopy;
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
array eendo2c    {*} eendo2c88 eendo2c90 eendo2c92 eendo2c94 eendo2c96 eendo2c98 eendo2c00 eendo2c02 eendo2c04 eendo2c06 eendo2c08 eendo2c10 eendo2c12 eendo2c14 eendo2c16 eendo2c18  ;

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

proc datasets nolist;
   delete hp_der hp_der_2 cancer diabetes cvd16 hp_dead hp_mets
   hp86 hp88 hp90 hp92 hp94 hp96 hp98 hp00 hp02 hp04 hp06 hp08 hp12 hp14 hp16 hp18
   h86_nts h90_nts h94_nts h98_nts h02_nts h06_nts h10_nts h14_nts h18_nts 
   h86_dt h90_dt h94_dt h98_dt 
   ahei86 ahei90 ahei94 ahei98 ahei02 ahei06 ahei10 ahei18 ahei14
   nses8818 ;
   run;



data hpfs;
     set base end=_end_;

/*Siyu, revise the cutoff: NHS1 2020.6, NHS2 2021.6, HPFS 2020.1 */
cutoff=1441;    /*** Set cutoff at 2020.1 ***/

array irt   {*} rtmnyr86 rtmnyr88 rtmnyr90 rtmnyr92 rtmnyr94 rtmnyr96 rtmnyr98 rtmnyr00 rtmnyr02 rtmnyr04 rtmnyr06 rtmnyr08 rtmnyr10 rtmnyr12 rtmnyr14 rtmnyr16 rtmnyr18 cutoff;
array tvar  {*} t86 t88 t90 t92 t94 t96 t98 t00 t02 t04 t06 t08 t10 t12 t14 t16 t18;

array age   {*} age86   age88   age90   age92   age94   age96   age98   age00   age02   age04   age06   age08   age10 age12 age14 age16 age18;   
array smok  {*} smoke86 smoke88 smoke90 smoke92 smoke94 smoke96 smoke98 smoke00 smoke02 smoke04 smoke06 smoke08 smoke10 smoke12 smoke14 smoke16 smoke18;
array cignm {*} cgnm86  cgnm88  cgnm90  cgnm92  cgnm94  cgnm96  cgnm98  cgnm00   cgnm02   cgnm04   cgnm06    cgnm08 cgnm10 cgnm12 cgnm14 cgnm16 cgnm18;
array smm   {*} smm86 smm88 smm90 smm92 smm94 smm96 smm98 smm00 smm02 smm04 smm06 smm08 smm10 smm12 smm14 smm16 smm18;
array bmi   {*} bmi86 bmi88 bmi90 bmi92 bmi94 bmi96 bmi98 bmi00 bmi02 bmi04 bmi06 bmi08 bmi10 bmi12 bmi14 bmi16 bmi18;
array bmicq {*} bmic86 bmic88 bmic90 bmic92 bmic94 bmic96 bmic98 bmic00 bmic02 bmic04 bmic06 bmic08 bmic10 bmic12 bmic14 bmic16 bmic18;
array alcon {*} alco86n alco86n alco90n alco90n alco94n alco94n alco98n alco98n alco02n alco02n alco06n alco06n alco10n alco10n alco14n alco14n alco18n;
array alcq  {*} alc86 alc88 alc90 alc92 alc94 alc96 alc98 alc00 alc02 alc04 alc06 alc08 alc10 alc12 alc14 alc16 alc18; 
array act   {*} act86 act88 act90 act92 act94 act96 act98 act00 act02 act04 act06 act08 act10 act12 act14 act16 act18;
array mvit  {*} mvt86 mvt88 mvt90 mvt92 mvt94 mvt96 mvt98 mvt00 mvt02 mvt04 mvt06 mvt08 mvt10 mvt12 mvt14 mvt16 mvt18; 
array aspi  {*} asp86 asp88 asp90 asp92 asp94 asp96 asp98 asp00 asp02 asp04 asp06 asp08 asp10 asp12 asp14 asp16 asp18;
*array mar  {*} mar86 mar88  mar90  mar92  mar94  mar96  mar98  mar00  mar02  mar04  mar06  mar08 mar10 mar12 mar14 mar16 mar18;
*array alone{*} livng86 livng88 livng90 livng92 livng94 livng96 livar96 livng00  lalon02  lalon04  lalon06 lalon08 lalon10 lalon12 lalon14 lalon16 lalon18;
array nSESArr{*} nSES_88m nSES_88m nSES_90 nSES_92 nSES_94 nSES_96 nSES_98 nSES_00 nSES_02 nSES_04 nSES_06 nSES_08 nSES_10 nSES_12 nSES_14 nSES_16 nSES_18 ;
 
array hbp   {*}  hbp86   hbp88   hbp90   hbp92   hbp94   hbp96   hbp98   hbp00    hbp02    hbp04    hbp06     hbp08   hbp10  hbp12  hbp14  hbp16  hbp18; 
array chol  {*}  chol86  chol88  chol90  chol92  chol94  chol96  chol98  chol00   chol02   chol04   chol06    chol08  chol10 chol12 chol14 chol16 chol18;
array antihp {*} htnrx86 htnrx88 htnrx90 htnrx92 htnrx94 htnrx96 htnrx98 htnrx00 htnrx02 htnrx04 htnrx06 htnrx08 htnrx10 htnrx12 htnrx14 htnrx16 htnrx18;
array antitc {*} hchtx86 hchtx88 hchtx90 hchtx92 hchtx94 hchtx96 hchtx98 hchtx00 hchtx02 hchtx04 hchtx06 hchtx08 hchtx10 hchtx12 hchtx14 hchtx16 hchtx18;
array diab   {*} db86    db88    db90    db92    db94    db96    db98    db00     db02     db04     db06      db08    db10   db12   db14   db16   db18;
*array mi    {*} mi86    mi88    mi90    mi92    mi94    mi96    mi98    mi00     mi02     mi04     mi06      mi08    mi10   mi12   mi14   mi16   mi18;
*array cabg  {*} cabg86  cabg88  cabg90  cabg92  cabg94  cabg96  cabg98  cabg00   cabg02   cabg04   cabg06    cabg08  cabg10 cabg12 cabg14 cabg16 cabg18;
*array strk  {*} str86   str88   str90   str92   str94   str96   str98   strk00   strk02   strk04   strk06    strk08  strk10 strk12 strk14 strk16 strk18;
array canc   {*} can86  can88   can90   can92   can94   can96   can98   can00    can02    can04    can06     can08   can10  can12  can14 can16  can18;


array    carbp  {*}   carbo86p  carbo86p  carbo90p  carbo90p  carbo94p  carbo94p  carbo98p  carbo98p  carbo02p   carbo02p  carbo06p carbo06p  carbo10p carbo10p  carbo14p  carbo14p carbo18p;
array    hcarbp  {*}   hcarbo86p  hcarbo86p  hcarbo90p  hcarbo90p  hcarbo94p  hcarbo94p  hcarbo98p  hcarbo98p  hcarbo02p   hcarbo02p  hcarbo06p hcarbo06p  hcarbo10p hcarbo10p  hcarbo14p  hcarbo14p hcarbo18p;
array sumref_carbp {*}  sumref_carb86p  sumref_carb86p  sumref_carb90p  sumref_carb90p  sumref_carb94p  sumref_carb94p  sumref_carb98p  sumref_carb98p  sumref_carb02p   sumref_carb02p  sumref_carb06p sumref_carb06p  sumref_carb10p sumref_carb10p  sumref_carb14p  sumref_carb14p sumref_carb18p;
array     protp  {*}   prot86p  prot86p  prot90p  prot90p  prot94p  prot94p  prot98p  prot98p  prot02p   prot02p  prot06p prot06p  prot10p prot10p  prot14p  prot14p prot18p;
array    aprotp  {*}   aprot86p  aprot86p  aprot90p  aprot90p  aprot94p  aprot94p  aprot98p  aprot98p  aprot02p   aprot02p  aprot06p aprot06p  aprot10p aprot10p  aprot14p  aprot14p aprot18p;
array    vprotp  {*}   vprot86p  vprot86p  vprot90p  vprot90p  vprot94p  vprot94p  vprot98p  vprot98p  vprot02p   vprot02p  vprot06p vprot06p  vprot10p vprot10p  vprot14p  vprot14p vprot18p;
array    tfatp   {*}   tfat86p  tfat86p  tfat90p  tfat90p  tfat94p  tfat94p  tfat98p  tfat98p  tfat02p   tfat02p  tfat06p tfat06p  tfat10p tfat10p  tfat14p  tfat14p tfat18p;
array    afatp   {*}   afat86p  afat86p  afat90p  afat90p  afat94p  afat94p  afat98p  afat98p  afat02p   afat02p  afat06p afat06p  afat10p afat10p  afat14p  afat14p afat18p;
array    vfatp   {*}   vfat86p  vfat86p  vfat90p  vfat90p  vfat94p  vfat94p  vfat98p  vfat98p  vfat02p   vfat02p  vfat06p vfat06p  vfat10p vfat10p  vfat14p  vfat14p vfat18p;

array eendo2c    {*} eendo2c88 eendo2c88 eendo2c90 eendo2c92 eendo2c94 eendo2c96 eendo2c98 eendo2c00 eendo2c02 eendo2c04 eendo2c06 eendo2c08 eendo2c10 eendo2c12 eendo2c14 eendo2c16 eendo2c18  ;
array eendo3c    {*} eendo3c88 eendo3c88 eendo3c90    eendo3c92    eendo3c94    eendo3c96    eendo3c98    eendo3c00    eendo3c02  eendo3c04    eendo3c06    eendo3c08    eendo3c10 eendo3c12 eendo3c14 eendo3c16 eendo3c18 ;


   /* Create ltf using the original irt values */
ltf=.; do i=1 to dim(irt)-1; if irt{i}>0 then ltf=irt{i}; end;
if ltf=irt{dim(irt)-1} then ltf=.;

/* Correct irt values and calculate age */
do i=1 to dim(irt)-1;
    if (irt{i} <  (1009+24*i)) or (irt{i} >= (1033+24*i)) then irt{i}=1009+24*i;
    age{i}=int((irt{i}-dbmy09)/12);
end;

do j=2 to DIM(bmi);
      if age{j}=.      then age{j}=age{j-1}+2;
      if bmi{j}=.      then bmi{j}=bmi{j-1};
      if smok{j}=.     then smok{j}=smok{j-1};
      if cignm{j}=.    then cignm{j}=cignm{j-1};
      if alcon{j}=.    then alcon{j}=alcon{j-1};
      if act{j}=.      then act{j}=act{j-1};
      if mvit{j}=.     then mvit{j}=mvit{j-1};
      if aspi{j}=.     then aspi{j}=aspi{j-1};
      if antihp{j}=.   then antihp{j}=antihp{j-1};
      if antitc{j}=.   then antitc{j}=antitc{j-1};
      if diab{j}=.     then diab{j}=diab{j-1};
      if canc{j}=.     then canc{j}=canc{j-1};
      if nSESArr{j}=.  then nSESArr{j}=nSESArr{j-1};
end;

/* -single  s-stop updated  v-average  sv-stop updated average */
/* -s */
array calors {17} calor86ns calor86ns calor90ns calor90ns calor94ns calor94ns calor98ns calor98ns calor02ns calor02ns calor06ns calor06ns calor10ns calor10ns calor14ns calor14ns calor18ns;
array aheims {17} ahei86s ahei86s ahei90s ahei90s ahei94s ahei94s ahei98s ahei98s ahei02s ahei02s ahei06s ahei06s ahei10s ahei10s ahei14s ahei14s ahei18s;   
array lcss   {17}  lcs86s    lcs86s    lcs90s    lcs90s    lcs94s    lcs94s    lcs98s    lcs98s    lcs02s    lcs02s    lcs06s    lcs06s    lcs10s  lcs10s  lcs14s  lcs14s lcs18s;
array alcss  {17}  alcs86s   alcs86s   alcs90s   alcs90s   alcs94s   alcs94s   alcs98s   alcs98s   alcs02s   alcs02s   alcs06s   alcs06s   alcs10s alcs10s alcs14s alcs14s alcs18s;
array vlcss  {17}  vlcs86s   vlcs86s   vlcs90s   vlcs90s   vlcs94s   vlcs94s   vlcs98s   vlcs98s   vlcs02s   vlcs02s   vlcs06s   vlcs06s   vlcs10s vlcs10s vlcs14s vlcs14s vlcs18s;
array ulcss  {17}  ulcs86s   ulcs86s   ulcs90s   ulcs90s   ulcs94s   ulcs94s   ulcs98s   ulcs98s   ulcs02s   ulcs02s   ulcs06s   ulcs06s   ulcs10s ulcs10s ulcs14s ulcs14s ulcs18s;
array hlcss  {17}  hlcs86s   hlcs86s   hlcs90s   hlcs90s   hlcs94s   hlcs94s   hlcs98s   hlcs98s   hlcs02s   hlcs02s   hlcs06s   hlcs06s   hlcs10s hlcs10s hlcs14s hlcs14s hlcs18s;
array lfss   {17}  lfs86s    lfs86s    lfs90s    lfs90s    lfs94s    lfs94s    lfs98s    lfs98s    lfs02s    lfs02s    lfs06s    lfs06s    lfs10s  lfs10s  lfs14s  lfs14s  lfs18s;
array ulfss  {17}  ulfs86s   ulfs86s   ulfs90s   ulfs90s   ulfs94s   ulfs94s   ulfs98s   ulfs98s   ulfs02s   ulfs02s   ulfs06s   ulfs06s   ulfs10s ulfs10s ulfs14s ulfs14s ulfs18s;
array hlfss  {17}  hlfs86s   hlfs86s   hlfs90s   hlfs90s   hlfs94s   hlfs94s   hlfs98s   hlfs98s   hlfs02s   hlfs02s   hlfs06s   hlfs06s   hlfs10s hlfs10s hlfs14s hlfs14s hlfs18s;
array alfss {17} alfs86s alfs86s alfs90s alfs90s alfs94s alfs94s alfs98s alfs98s alfs02s alfs02s alfs06s alfs06s alfs10s alfs10s alfs14s alfs14s alfs18s;
array vlfss {17} vlfs86s vlfs86s vlfs90s vlfs90s vlfs94s vlfs94s vlfs98s vlfs98s vlfs02s vlfs02s vlfs06s vlfs06s vlfs10s vlfs10s vlfs14s vlfs14s vlfs18s;

/* -sv */
array calorsv {17}  calor86nsv calor86nsv calor90nsv calor90nsv calor94nsv calor94nsv calor98nsv calor98nsv calor02nsv calor02nsv calor06nsv calor06nsv calor10nsv calor10nsv calor14nsv calor14nsv calor18nsv;
array aheimsv {17}  ahei86sv ahei86sv ahei90sv ahei90sv ahei94sv ahei94sv ahei98sv ahei98sv ahei02sv ahei02sv ahei06sv ahei06sv ahei10sv ahei10sv ahei14sv ahei14sv ahei18sv;   
array lcssv   {17}  lcs86sv    lcs86sv    lcs90sv    lcs90sv    lcs94sv    lcs94sv    lcs98sv    lcs98sv    lcs02sv    lcs02sv    lcs06sv    lcs06sv    lcs10sv  lcs10sv  lcs14sv  lcs14sv lcs18sv;
array alcssv  {17}  alcs86sv   alcs86sv   alcs90sv   alcs90sv   alcs94sv   alcs94sv   alcs98sv   alcs98sv   alcs02sv   alcs02sv   alcs06sv   alcs06sv   alcs10sv alcs10sv alcs14sv alcs14sv alcs18sv;
array vlcssv  {17}  vlcs86sv   vlcs86sv   vlcs90sv   vlcs90sv   vlcs94sv   vlcs94sv   vlcs98sv   vlcs98sv   vlcs02sv   vlcs02sv   vlcs06sv   vlcs06sv   vlcs10sv vlcs10sv vlcs14sv vlcs14sv vlcs18sv;
array ulcssv  {17}  ulcs86sv   ulcs86sv   ulcs90sv   ulcs90sv   ulcs94sv   ulcs94sv   ulcs98sv   ulcs98sv   ulcs02sv   ulcs02sv   ulcs06sv   ulcs06sv   ulcs10sv ulcs10sv ulcs14sv ulcs14sv ulcs18sv;
array hlcssv  {17}  hlcs86sv   hlcs86sv   hlcs90sv   hlcs90sv   hlcs94sv   hlcs94sv   hlcs98sv   hlcs98sv   hlcs02sv   hlcs02sv   hlcs06sv   hlcs06sv   hlcs10sv hlcs10sv hlcs14sv hlcs14sv hlcs18sv;
array lfssv   {17}  lfs86sv    lfs86sv    lfs90sv    lfs90sv    lfs94sv    lfs94sv    lfs98sv    lfs98sv    lfs02sv    lfs02sv    lfs06sv    lfs06sv    lfs10sv  lfs10sv  lfs14sv  lfs14sv  lfs18sv;
array ulfssv  {17}  ulfs86sv   ulfs86sv   ulfs90sv   ulfs90sv   ulfs94sv   ulfs94sv   ulfs98sv   ulfs98sv   ulfs02sv   ulfs02sv   ulfs06sv   ulfs06sv   ulfs10sv ulfs10sv ulfs14sv ulfs14sv ulfs18sv;
array hlfssv  {17}  hlfs86sv   hlfs86sv   hlfs90sv   hlfs90sv   hlfs94sv   hlfs94sv   hlfs98sv   hlfs98sv   hlfs02sv   hlfs02sv   hlfs06sv   hlfs06sv   hlfs10sv hlfs10sv hlfs14sv hlfs14sv hlfs18sv;
array alfssv {17} alfs86sv alfs86sv alfs90sv alfs90sv alfs94sv alfs94sv alfs98sv alfs98sv alfs02sv alfs02sv alfs06sv alfs06sv alfs10sv alfs10sv alfs14sv alfs14sv alfs18sv;
array vlfssv {17} vlfs86sv vlfs86sv vlfs90sv vlfs90sv vlfs94sv vlfs94sv vlfs98sv vlfs98sv vlfs02sv vlfs02sv vlfs06sv vlfs06sv vlfs10sv vlfs10sv vlfs14sv vlfs14sv vlfs18sv;

array whgrn_sv {17}  whgrn_s86sv whgrn_s86sv whgrn_s90sv whgrn_s90sv whgrn_s94sv whgrn_s94sv whgrn_s98sv whgrn_s98sv whgrn_s02sv whgrn_s02sv whgrn_s06sv whgrn_s06sv whgrn_s10sv whgrn_s10sv whgrn_s14sv whgrn_s14sv whgrn_s18sv;
array rgrain_sv{17}  rgrain_s86sv rgrain_s86sv rgrain_s90sv rgrain_s90sv rgrain_s94sv rgrain_s94sv rgrain_s98sv rgrain_s98sv rgrain_s02sv rgrain_s02sv rgrain_s06sv rgrain_s06sv rgrain_s10sv rgrain_s10sv rgrain_s14sv rgrain_s14sv rgrain_s18sv ;   
array frui_sv  {17}  frui_s86sv    frui_s86sv    frui_s90sv    frui_s90sv    frui_s94sv    frui_s94sv    frui_s98sv    frui_s98sv    frui_s02sv    frui_s02sv    frui_s06sv    frui_s06sv    frui_s10sv frui_s10sv frui_s14sv frui_s14sv frui_s18sv;
array leg_sv   {17}  leg_s86sv   leg_s86sv   leg_s90sv   leg_s90sv   leg_s94sv   leg_s94sv   leg_s98sv   leg_s98sv   leg_s02sv   leg_s02sv   leg_s06sv   leg_s06sv   leg_s10sv leg_s10sv leg_s14sv leg_s14sv leg_s18sv ;
array pota_sv  {17}  pota_s86sv   pota_s86sv   pota_s90sv   pota_s90sv   pota_s94sv   pota_s94sv   pota_s98sv   pota_s98sv   pota_s02sv   pota_s02sv   pota_s06sv   pota_s06sv   pota_s10sv pota_s10sv pota_s14sv pota_s14sv pota_s18sv ;
array veg_noleg_sv  {17}  veg_noleg_s86sv   veg_noleg_s86sv   veg_noleg_s90sv   veg_noleg_s90sv   veg_noleg_s94sv   veg_noleg_s94sv   veg_noleg_s98sv   veg_noleg_s98sv   veg_noleg_s02sv   veg_noleg_s02sv   veg_noleg_s06sv   veg_noleg_s06sv   veg_noleg_s10sv veg_noleg_s10sv veg_noleg_s14sv veg_noleg_s14sv veg_noleg_s18sv ;
array fish_sv  {17}  fish_s86sv   fish_s86sv   fish_s90sv   fish_s90sv   fish_s94sv   fish_s94sv   fish_s98sv   fish_s98sv   fish_s02sv   fish_s02sv   fish_s06sv   fish_s06sv   fish_s10sv fish_s10sv fish_s14sv fish_s14sv fish_s18sv ;
array rpmeats_sv   {17}  rpmeats_s86sv    rpmeats_s86sv    rpmeats_s90sv    rpmeats_s90sv    rpmeats_s94sv    rpmeats_s94sv    rpmeats_s98sv    rpmeats_s98sv    rpmeats_s02sv    rpmeats_s02sv    rpmeats_s06sv    rpmeats_s06sv    rpmeats_s10sv rpmeats_s10sv rpmeats_s14sv rpmeats_s14sv rpmeats_s18sv ;
array ssb_sv  {17}  ssb_s86sv   ssb_s86sv   ssb_s90sv   ssb_s90sv   ssb_s94sv   ssb_s94sv   ssb_s98sv   ssb_s98sv   ssb_s02sv   ssb_s02sv   ssb_s06sv   ssb_s06sv   ssb_s10sv ssb_s10sv ssb_s14sv ssb_s14sv ssb_s18sv ;
array dairy_sv  {17}  dairy_s86sv   dairy_s86sv   dairy_s90sv   dairy_s90sv   dairy_s94sv   dairy_s94sv   dairy_s98sv   dairy_s98sv   dairy_s02sv   dairy_s02sv   dairy_s06sv   dairy_s06sv   dairy_s10sv dairy_s10sv dairy_s14sv dairy_s14sv dairy_s18sv ;
array swtsnochoc_sv {17} swtsnochoc_s86sv swtsnochoc_s86sv swtsnochoc_s90sv swtsnochoc_s90sv swtsnochoc_s94sv swtsnochoc_s94sv swtsnochoc_s98sv swtsnochoc_s98sv swtsnochoc_s02sv swtsnochoc_s02sv swtsnochoc_s06sv swtsnochoc_s06sv swtsnochoc_s10sv swtsnochoc_s10sv swtsnochoc_s14sv swtsnochoc_s14sv swtsnochoc_s18sv  ;
array poul_sv {17} poul_s86sv poul_s86sv poul_s90sv poul_s90sv poul_s94sv poul_s94sv poul_s98sv poul_s98sv poul_s02sv poul_s02sv poul_s06sv poul_s06sv poul_s10sv poul_s10sv poul_s14sv poul_s14sv poul_s18sv  ;

array carbsv {17}  carbo86sv carbo86sv carbo90sv carbo90sv carbo94sv carbo94sv carbo98sv carbo98sv carbo02sv carbo02sv carbo06sv carbo06sv carbo10sv carbo10sv carbo14sv carbo14sv carbo18nsv;
array hcarbsv {17}  hcarbo86sv hcarbo86sv hcarbo90sv hcarbo90sv hcarbo94sv hcarbo94sv hcarbo98sv hcarbo98sv hcarbo02sv hcarbo02sv hcarbo06sv hcarbo06sv hcarbo10sv hcarbo10sv hcarbo14sv hcarbo14sv hcarbo18sv;   
array sumref_carbsv   {17}  sumref_carb86sv    sumref_carb86sv    sumref_carb90sv    sumref_carb90sv    sumref_carb94sv    sumref_carb94sv    sumref_carb98sv    sumref_carb98sv    sumref_carb02sv    sumref_carb02sv    sumref_carb06sv    sumref_carb06sv    sumref_carb10sv  sumref_carb10sv  sumref_carb14sv  sumref_carb14sv sumref_carb18sv;
array protsv  {17}  prot86sv   prot86sv   prot90sv   prot90sv   prot94sv   prot94sv   prot98sv   prot98sv   prot02sv   prot02sv   prot06sv   prot06sv   prot10sv prot10sv prot14sv prot14sv prot18sv;
array aprotsv  {17}  aprot86sv   aprot86sv   aprot90sv   aprot90sv   aprot94sv   aprot94sv   aprot98sv   aprot98sv   aprot02sv   aprot02sv   aprot06sv   aprot06sv   aprot10sv aprot10sv aprot14sv aprot14sv aprot18sv;
array vprotsv  {17}  vprot86sv   vprot86sv   vprot90sv   vprot90sv   vprot94sv   vprot94sv   vprot98sv   vprot98sv   vprot02sv   vprot02sv   vprot06sv   vprot06sv   vprot10sv vprot10sv vprot14sv vprot14sv vprot18sv;
array tfatsv  {17}  tfat86sv   tfat86sv   tfat90sv   tfat90sv   tfat94sv   tfat94sv   tfat98sv   tfat98sv   tfat02sv   tfat02sv   tfat06sv   tfat06sv   tfat10sv tfat10sv tfat14sv tfat14sv tfat18sv;
array afatsv {17} afat86sv afat86sv afat90sv afat90sv afat94sv afat94sv afat98sv afat98sv afat02sv afat02sv afat06sv afat06sv afat10sv afat10sv afat14sv afat14sv afat18sv;
array vfatsv {17} vfat86sv vfat86sv vfat90sv vfat90sv vfat94sv vfat94sv vfat98sv vfat98sv vfat02sv vfat02sv vfat06sv vfat06sv vfat10sv vfat10sv vfat14sv vfat14sv vfat18sv;

/*************************************************************************
***** If an irt date is before June of that qq year or after or equal ****
***** to the next qq year it is incorrect and should be defaulted to  ****
***** June of that qq year.    Make time period indicator tvar=0.     ****
*************************************************************************/



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
    if irt{i}< dtdth <=irt{i+1} and (140<=newicda<=208) then do; dead_can=1; tdead_can=dtdth-irt{i}; dtdx_dead_can=dtdth; end;

    dead_cvd   = 0; tdead_cvd=irt{i+1}-irt{i}; dtdx_dead_cvd=.;
    if irt{i}< dtdth <=irt{i+1} and (390<=newicda<=459) then do; dead_cvd=1; tdead_cvd=dtdth-irt{i}; dtdx_dead_cvd=dtdth; end;

    dead_resp=0;  tdead_resp=irt{i+1}-irt{i};  dtdx_dead_resp=.;      
    if irt{i}<dtdth<=irt{i+1} and (460<=newicda<520) then do; dead_resp=1; tdead_resp=dtdth-irt{i}; dtdx_dead_resp=dtdth; end;

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
    lcssvcon  =lcssv{i};
    alcssvcon =alcssv{i};
    vlcssvcon =vlcssv{i};
    ulcssvcon =ulcssv{i};
    hlcssvcon =hlcssv{i};
    lfssvcon  =lfssv{i};
    ulfssvcon =ulfssv{i};
    hlfssvcon =hlfssv{i};
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
%indic3(vbl=ethnic, prefix=ethnic, reflev=1, min=2, max=4, missing=., usemiss=0,
           label1='causcasian',
           label2='african',
           label3='asian',
           label4='other');
     
%indic3(vbl=eth3g, prefix=eth3g, reflev=1, min=2, max=3, missing=., usemiss=0,
           label1='causcasian',
           label2='african',
           label3='asian & others');

/* AGE GROUP */
      agecon=age{i};
      if age{i}<60 then agegp=1;     **** Define the agegp in the i-th period;
            else if age{i}<65 then agegp=2;
            else if age{i}<70 then agegp=3;
            else if age{i}<75 then agegp=4;
            else if age{i}<80 then agegp=5;
            else agegp=6;
            if age{i}=. then agegp=.;

      %indic3(vbl=agegp, prefix=agegp, min=2, max=6, reflev=1, missing=., usemiss=0,
                   label1='age<60ys',
                   label2='age60-64ys',
                   label3='age65-69ys',
                   label4='age70-74ys',
                   label5='age75-79ys',
                   label6='age>=80ys');
     if agecon>65 then agesub=1; else agesub=0;

/* Indicator for Smoking */ 
/* smoke  $label 1.Never;  2.Past; 3.No, unknown past history; 4.Current; 5.PASSTHRU;
   cgnm $label 1.1-4/d; 2.5-14; 3.15-24; 4.25-34; 5.35-44;6.45 or more; 7.PASSTHRU; 0.not app or no qx */ 

      if smok{i}=5 then smok{i}=.;      
      if smok{i}<4 then cignm{i}=0;      
     
      if cignm{i}=7 then cignm{i}=.;      
      if cignm{i}=0 then cignm{i}=.;      
    
      if smok{i}=1 then smm(i)=1;                    /*   smkm=1, never*/
      if smok{i}=2 then smm(i)=2;                    /*   smkm=2, past*/
      if smok{i}=3 then smm(i)=2;      
      if smok{i}=. then smm(i)=.;      
      if smok{i}=4 then 
         do;     
         if 1 <= cignm{i} and cignm{i} <= 2 then smm(i)=3;            /*   smkm=3, 1-14*/
         else if cignm{i}=3 then smm(i)=4;           /*    smkm=4,15-24*/
         else if 4 <= cignm{i} and cignm{i} <= 6 then smm(i)=5;       /* smkm=5: >=25/d */
         else if cignm{i}=. then smm(i)=3;           /* smkm=3: currnt smker, amount unknown */
         end;                                        /* smm=1 never; 2 past 3 1-14 4 15-24 5 >=25/d */

     smkk=smm{i};
     %indic3(vbl=smkk, prefix=smkk, min=2,max=5, reflev=1, missing=., usemiss=0,
              label1='never smker',
              label2='past smker',
              label3='curr 1-14 smker',
              label4='curr 15-24 smker',
              label5='curr 25+ smker');
      
      if smkk=2 then smkstatus=2;                     /* past*/
      else if smkk>=3 then smkstatus=3;               /* current*/
      else smkstatus=1;                               /* never*/
      %indic3(vbl=smkstatus, prefix=smkstatus, min=2, max=3, reflev=1, missing=., usemiss=0);
      if smkstatus=1 then smknever=1; else smknever=0;

/* Indicator for BMI */
      bmicon=bmi{i};if bmi{i}=0 then bmicon=.;

      if bmi{i}>0 and bmi{i}<21 then bmicq{i}=1;
         else if bmi{i}>=21 and bmi{i}<25 then bmicq{i}=2;
         else if bmi{i}>=25 and bmi{i}<30 then bmicq{i}=3;
         else if bmi{i}>=30 and bmi{i}<32 then bmicq{i}=4;
         else if bmi{i}>=32 then bmicq{i}=5;
         else bmicq{i}=3;
      
      if 25.00<= bmicon  then bmic25=1; else  bmic25=0; /*Binary flags for 25 and 30*/
      if 30.00<= bmicon  then bmic30=1; else  bmic30=0;

       bmic=bmicq{i};                                   /*5-level categorical BMI per visit*/
      %indic3(vbl=bmic, prefix=bmic, min=2, max=5, reflev=1, missing=., usemiss=0,
                label1='bmi <21',
                label2='bmi 21-24.9',
                label3='bmi 25-29.9',
                label4='bmi 30.0-31.9',
                label5='bmi 32.0+');
           
            if   bmic <=2 then bmic3cat=1;              /*3-level categorical BMI per visit*/
     else   if   bmic <=3 then bmic3cat=2;
     else   if   4 <=bmic  then bmic3cat=3;  
     %indic3(vbl= bmic3cat, prefix= bmic3cat, min=2, max=3, reflev=1, missing=., usemiss=0);
    
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


     bmib=bmicq{1};                                     /*Baseline 5-level categorical BMI*/
      %indic3(vbl=bmib, prefix=bmib, min=2, max=5, reflev=1, missing=., usemiss=0,
                label1='bmi <21',
                label2='bmi 21-24.9',
                label3='bmi 25-29.9',
                label4='bmi 30.0-31.9',
                label5='bmi 32.0+');
     if   bmib <=2 then bmib25=0; else  bmib25=1;
     if   bmib <=3 then bmib30=0; else  bmib30=1;
           if   bmib <=2 then bmib3cat=1;               /*Baseline 3-level categorical BMI*/
     else  if   bmib <=3 then bmib3cat=2;
     else  if  4 <=bmib  then bmib3cat=3;

     %indic3(vbl= bmib3cat, prefix= bmib3cat, min=2, max=3, reflev=1, missing=., usemiss=0);

    bmibcon=bmi{1};
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


/* Indicator Physical Activity */
      actcon=act{i}; /***continue, hours/wk***/

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

/* Indicator Alcohol consumption*/
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
      
/*** multiple vitamin ***/
        if mvit{i}=2 then mv=1; else mv=0;

/*** aspirin ***/
        asp=aspi{i}; if asp ne 1 then asp=0;

/***family history ***/
       if mifh ne 1 then mifh=0;
       if dbfh ne 1 then dbfh=0;
       if cafh ne 1 then cafh=0;

/* Indicator for History of antihypertensive user */
            select(antihp{i});
                when (1)     antihbp=1;
                otherwise    antihbp=0;
            end;

            if hbp{i}=1 then highbp=1; else highbp=0;
            if hbp{1}=1 then basehighbp=1; else basehighbp=0;
      
/****** Indicator for TC lowering ********/
            select(antitc{i});
                when (1)     antihtc=1;
                otherwise    antihtc=0;
            end;

           if chol{i}=1 then hightc=1; else hightc=0;
           if chol{1}=1 then basehightc=1; else basehightc=0;
      
/****** Indicator for History of diabetes ******/
           if diab{i}=1 then diabetes=1; else diabetes=0;

/****** Indicator for History of cancer ******/
           if canc{i}=1 then cancer=1; else cancer=0;

/****** endoscopy for screening purpose ******/
    endo2c=eendo2c{i}; /* ever screening endo, 1=no screening endo, 2=screening endo */
    endo3c=eendo3c{i}; /* ever endo3c, 1=no endo, 2=endo for symptoms, 3=screening endo */
    if endo2c = . then endo2c = 1;
/****** Mammogram in last 2 yrs *******/
    mam=0;

/* nSES */
nSES=nSESArr{i};


/* Baseline Food Inatke */
basewhgrn=whgrn86d;
basergrain=rgrain86d;
basefrui=frui86d;
baseleg=leg86d;
basepota=pota86d;
baseveg_noleg=veg_noleg86d;
basefish=fish86d;
baserpmeats=rpmeats86d;
basessb=ssb86d;
basedairy=dairy86d;
basesweets=swtsnochoc86d;
basepoultry=poul86d;


/****************  BASELINE EXCLUSIONS ********************************************/
   if i=1 then do;
   %exclude(exrec eq 1);   *multiple records and not in master file;
   %exclude(dbmy09 le 0);
   
   %exclude(baseffq eq 1);
   %exclude(lcs86 eq .);
   %exclude(alcs86 eq .);
   %exclude(vlcs86 eq .);
   %exclude(ulcs86 eq .);
   %exclude(hlcs86 eq .);
   %exclude(lfs86 eq .);
   %exclude(ulfs86 eq .);
   %exclude(hlfs86 eq .);

   %exclude(ltf eq irt{1});

   %exclude(0 lt dtdth   le irt{i});  
   %exclude(0 lt dtdx_allca  le irt{i}); 
   %exclude(can86   eq 1);        
   %exclude(0 lt dtdxcvd le irt{i});   
   %exclude(mi86    eq 1);
   %exclude(str86   eq 1);
   %exclude(cabg86  eq 1);    

   %exclude(0 lt dtdxdb2 le irt{i}); 
   %exclude(db86    eq 1);  

   %exclude(age86 eq .);

   %exclude(calor86n gt 4200); 
   %exclude(calor86n lt 800);
   %exclude(calor86n eq .);  


   %output();
  end;

  else if i> 1 then do;
   %exclude(irt{i-1} lt dtdth le irt{i} );       *censor death; 
   %exclude(irt{i-1} lt dtdx_obeca le irt{i} );  *censor obesity-related cancer: stop follow-up after obesity-related cancer;
   %output();
  end;
  end;

/* END OF DO-LOOP OVER TIME PERIODs */
%endex();

proc sort data=hpfs; by period; run;

proc rank data=hpfs out=hpfs group=5;
by  period;
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

data hpfs;
set hpfs end=_end_;

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



keep 
      id  period 
      /**Time scale**/
      rtmnyr86 rtmnyr88  rtmnyr90 rtmnyr92 rtmnyr94 rtmnyr96 rtmnyr98 rtmnyr00 rtmnyr02 rtmnyr04 rtmnyr06 rtmnyr08 rtmnyr10 rtmnyr12 rtmnyr14 rtmnyr16 rtmnyr18 cutoff 
      t86 t88 t90 t92 t94 t96 t98 t00 t02 t04 t06 t08 t10 t12 t14 t16 t18
      dtdxchd dtdxmi dtdth dtdxca dtdx_allca dtdx_obeca dtdx_nonobeca fup_allca
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
      white  eth3g 
      agecon agegp  agesub  
      bmicon bmic3cat &bmic3cat_ bmib3cat &bmib3cat_  bmic25 bmic30 bmic9 &bmic9_ bmib9 &bmib9_ 
      smkstatus &smkstatus_ smknever alc  &alc_ neverdrinker  
      alcocon actcon actc &actc_ active 
      dbfh mifh cafh basehighbp basehightc highbp hightc diabetes cancer
      mv asp antihbp antihtc
      ltf
      basewhgrn basergrain basefrui  baseleg basepota baseveg_noleg  basefish baserpmeats  basessb basedairy basesweets basepoultry 
      whgrn_svcon rgrain_svcon frui_svcon leg_svcon pota_svcon veg_noleg_svcon fish_svcon rpmeats_svcon ssb_svcon 
      dairy_svcon swtsnochoc_svcon poul_svcon  
      mam endo2c endo3c 
      nSES 
; 
run;


/*
proc means data=hpfs mean std median n nmiss min max;
var 
calorscon  calorsvcon  aheimscon  aheimsvcon agecon bmicon alcocon actcon 

lcsscon    lcssvcon   alcsscon   alcssvcon  vlcsscon   vlcssvcon  ulcsscon   ulcssvcon  hlcsscon   hlcssvcon
lfsscon    lfssvcon   ulfsscon   ulfssvcon  hlfsscon   hlfssvcon;
run;

proc freq data=hpfs;
tables 
agegp white bmib3cat smkstatus alc actc mifh basehighbp basehightc mv asp
;
run;
*/