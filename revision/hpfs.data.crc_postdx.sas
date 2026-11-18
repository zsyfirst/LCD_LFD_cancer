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
               stage, duration, tumor differentiation, location.
revision:      add SES as covarites
Date:  10/2026

******************************************************/


/** Cancer **/
data allcancer;
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

        if cancerp=1 and dateca>0 and ca=1 and icdx ne 173.2 and icdx ne 173.1 then dtdxca=dateca;


    /* ---- Site groups ---- */
      /* Digestive */
      if 140 <= icdx <= 149 and ca=1 then do; oral_ca=1;      dtdx_oral      = dtdxca; end;
      if icdx=150 and ca=1 then do; esoph_ca=1;               dtdx_esoph     = dtdxca; end;
      if icdx=151 and ca=1 then do; stomach_ca=1;             dtdx_stomach   = dtdxca; end;
      if icdx=152 and ca=1 then do; smallintest_ca=1;         dtdx_smallintest=dtdxca; end;
      colorect_ca=0; if icdx in (153,154) and ca=1 then do; colorect_ca=1;   dtdx_colorect  = dtdxca; end;
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
     proc sort data=allcancer; by id dtdxca; run;

     data allcancer;
     set allcancer end=_end_;

     dtdx_obeca=.;
     if obese_ca=1 and dtdxca>0 and ca=1 and icdx in (150,151,153,154,155,156,157,174,182,183,189,193,203,225) then dtdx_obeca=dtdxca;
     dtdx_nonobeca=.;
     if nonobeca=1 and dtdxca>0 and ca=1  then dtdx_nonobeca=dtdxca;
    
     keep id conf dtdxca colorect_ca brca_ca dtdx_colorect dtdx_brca 
          /*obese_ca nonobeca smoke_ca r_smoke_ca  
          esoph_ca colorect_ca stomach_ca liver_ca gallbladder_ca panc_ca kidney_ca thyroid_ca myeloma_ca menin_ca */
          ;
  run;


/************** Colorectal Cancer **************/
data crc_hp;

infile 'CRC_8620_TNM.dat';
input 
hpfsid     1-6   /* id number for HPFS*/
colorect   7     /* all colorectal cancer-conf 11-19*/
colon      8     /* colon cancer, MOSTLY CONF=11, BUT NOT ALL*/
rectum     9     /* rectal cancer, MOSTLY CONF=11, BUT NOT ALL */
prox       10    /* proximal cancer, MOSTLY CONF=11, BUT NOT ALL */
distal     11    /* distal cancer, MOSTLY CONF=11, BUT NOT ALL */
distcr     12    /* distal and or rectum, MOSTLY CONF=11, BUT NOT ALL */
dateca     13-16 /* date of dx in NHS mo-yr format*/
coldirt    17    /*colon cancer including those without site info conf 11-19*/
dead       18    /* died of CRC either primary or secondary cause of death*/
advanced   19   /* old staging: died of CRC or lymphnode and or distant mets*/
mets       20    /* old staging: lympnode and or distant mets*/
ab         21    /* old staging: early stage, see coding sheet state 1 or 2*/
conf       22-24 /* confirmation codes medical records if dx through dx follow up*/
dconf      25-27 /* confirmation code death file if dx through death file, primary cause*/
celltype   28-29  /* Celltype 
                10    Epithelial (unspecified)
    11    Adenocarcinoma
    12    Mucinous Adenocarcinoma
    13    Signet-ring Cell
    14    Adeno-Squamous
    15    Undifferentiated*/
grade      30    /* Histological Grade
    1       Well Differentiated
    2       Grade 2 Moderately Diff..
    3       Grade 3 Poorly Diff..
    4       Unspecified*/
other_cond 31   /* OTHER CONDITIONS
    0        None
    1        Chronic Ulcerative Colitis
    2        Chron s granulomatous Colitis
    3        Unspecified Inflammatory Bowel Disease
    4        Familial Polyposis Syndromes*/
tstage     32   /* EXTENT OF DISEASE THROUGH BOWEL WALL
    1 T1 - invades submucosa
    2  T2 - invades muscularis
    3 T3 -  invades through muscularis into subserosa
    4 T4 - invades other organs structures, perforation
    5 Tis - carcinoma in situ, invasion of lamina propria
    6 Tx  - unknown*/
nstage     33   /* NUMBER OF LYMPH NODES INVOLVED
    1 N0 - no regional lymph nodes metastases
    2 N1 - metastases in 1 - 3 regional lymph nodes
    3 N2 - metastases in 4 or more regional lymph nodes
    4 Positive lymph nodes - number unknown
    5 Nx - unknown status*/   
mstage     34   /* DISTANT METASTASES
    1 M0 - no distant metastases known
    2 M1 - distant metastases*/
mwork      35   
    /* TNM Stage Grouping I-IV, 1=yes, else missing, see below*/ 
    /* TNM stage grouping:*/ 
    /*Stage IV any T, any N, M1*/
    /*Stage III any T, N1,N2, M0, also includes our NHS code tstage=4 which is number of lymphnodes unknown*/
    /*Stage II T3, N0, M0 or T4, N0, M0*/
    /*Stage I T1, N0, M0 or T2, N0, M0*/
TNM_1 36 
TNM_2 37
TNM_3 38
TNM_4 39 
overlap_conf   40     /* diagnosed in both disease and death primary cause fu 1=yes*/
finalstage     41     /* TNM_1   finalstage = 1
                   TNM_2   finalstage = 2
             TNM_3   finalstage = 3
             TNM_4   finalstage = 4*/

dconfs         42-43  /* confirmed as secondary cause of death, conf code*/
overlap_confs  44     /*diagnosed in disease and death secondary cause fu 1=yes*/
sitehp 45-46  /* subcode 5 and 6 */
crcother       47     /* confirmed colorectal cancer without information on site */
stage          48     /* DUKES stage 
      1= spread to mucosa and muscle only
      2= preirectum or pericolic and lymph
      3= lymph node and /or  distant  metastasis
      4= unspecified*/
dukes_sub      49    /* 1=yes 0=no  substituted dukes into finalstage if finalstage is missing and dukes stage s 1 or 2 */
; 

id=hpfsid;
dtdxca=dateca;
dtdx_cr = int(dateca/12 + 1900);
if conf in (11,12,13,14,15,16,17,18,19)  and colorect = 1 then cr=1; else cr=0;
conf1119=0; if 11<=conf<=19 then conf1119=1;

     if prox=1 then subsite_cr=1;
else if distal=1 then subsite_cr=2;
else if rectum=1 then subsite_cr=3;
else subsite_cr=4;

grade_cr = grade;
if grade_cr not in (1,2,3,4) then grade_cr = 4; /* unknown / unspecified */

*stage_cr = finalstage;
*if tstage = 5 and nstage = 1 and mstage = 1 then stage_cr = 1;
*if stage_cr = . then stage_cr = .;   /* unknown / unspecified */

/* ---------- AJCC stage for colorectal cancer ---------- */
stage_cr  = finalstage;                                   /* 1. official stage grouping */
stage_src = ifn(stage_cr ne ., 1, .);

if stage_cr = . then do;                                  /* 2. derive from T, N, M     */
    if      mstage = 2 or mwork = 1               then stage_cr = 4;   /* distant metastasis       */
    else if nstage in (2,3,4)                     then stage_cr = 3;   /* any node-positive, M0/Mx */
    else if nstage = 1 and tstage in (3,4)        then stage_cr = 2;   /* N0, T3-T4                */
    else if nstage = 1 and tstage in (1,2)        then stage_cr = 1;   /* N0, T1-T2                */
    else if nstage = 1 and tstage = 5             then stage_cr = 1;   /* Tis N0 -- as in your current code */
    if stage_cr ne . then stage_src = 2;
end;

if stage_cr = . then do;                                  /* 3. Dukes stage             */
    if      stage = 1                                    then stage_cr = 1;   /* Dukes A ~ stage I  */
    else if stage = 2                                    then stage_cr = 2;   /* Dukes B ~ stage II */
    else if stage = 3 and (mstage = 1 or mwork = 2)      then stage_cr = 3;   /* Dukes C, confirmed M0 */
    if stage_cr ne . then stage_src = 3;
end;

/* SEER summary stage: needs only node status, so it can be filled even when stage_cr is not */
if      stage_cr = 4 or mstage = 2 or mwork = 1     then sumstage = 3;   /* distant (excluded) */
else if nstage = 1 or stage_cr in (1,2)             then sumstage = 1;   /* localized */
else if nstage in (2,3,4) or stage_cr = 3           then sumstage = 2;   /* regional  */
else if ab = 1                                      then sumstage = 1;   /* old flag: early stage */
else if mets = 1 and mstage ne 2                    then sumstage = 2;   /* old flag: nodes/mets  */
else                                                     sumstage = .;

keep id cr colorect conf1119 stage_cr stage_src stage finalstage grade_cr subsite_cr dtdxca dtdx_cr;
run;
proc sort data=crc_hp; by id dtdxca; run;
proc sort data=crc_hp nodupkey; by id; run;

proc means data=crc_hp; var _numeric_;
proc freq; tables cr*(stage_cr grade_cr subsite_cr)/missing;run; 
proc freq data=crc_hp;tables dtdx_cr;run; /*dtdx_cr 1986 to 2020;*/

proc freq data=crc_hp;
    tables stage_src*stage_cr / missing;          /* how many recovered at each step */
    tables finalstage*stage_cr / missing;         /* step 1 must be unchanged        */
    tables stage*stage_cr / missing;              /* Dukes vs final: should broadly agree */
run;

                                                                             
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
alfs86=sum(carbo86pdd,aprot86pdd,vfat86pdd);
vlfs86=sum(carbo86pdd,vprot86pdd,afat86pdd);

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
              ); if mvt86d=0 or mvt86d=1 then mvt86=1; else if mvt86d=2 then mvt86=2; else mvt86=3; run;


%hp88(keep=hbp88 db88 chol88 trig88 mi88 ang88 str88 cabg88 pe88
              colc88  pros88  lymp88  ocan88  mel88  
              mvt88 mvtr88 asp88 mar88  livng88
              betab88 lasix88 diur88 ccblo88 ald88
              chrx88); 
              if asp88=1 then asp88=1;else asp88=0;
              if mvt88=1 then mvtr88=2; else if mvt88=2 then mvtr88=1; else mvtr88=3; 
              mvt88=mvtr88;
              run;

%hp90(keep=hbp90 db90 chol90 trig90 mi90 ang90 str90 cabg90 pe90
              colc90  pros90  lymp90  ocan90  mel90  
              asp90 mar90 livng90 
              mdb90 fdb90 sdb90 dbfh90 
              betab90 lasix90 diur90 ccblo90 ald90
              chrx90); 
              if asp90=1 then asp90=1;else asp90=0;
              if fdb90=1 or mdb90=1 or sdb90=1 then dbfh90=1; else dbfh90=0;run;

%h90_dt(keep=mvt90d mvt90   ); mvt90=mvt90d; run;


%hp92(keep=hbp92 db92 chol92 trig92 mi92 ang92 str92 cabg92 
              colc92  pros92  lymp92  ocan92  mel92  
              mvt92 asp92 mar92 livng92
              betab92 lasix92 thiaz92 ccblo92 ald92  
              chrx92
              fdb92 mdb92 sdb92 dbfh92
              mlng92 mclc92 flng92 fclc92 fpro92 slng92 sclc92 spro92 cafh92);
              if fdb92=1 or mdb92=1 or sdb92=1 then dbfh92=1; else dbfh92=0; 
              if mlng92=1 or mclc92=1 or flng92=1 or fclc92=1 or fpro92=1 or slng92=1 or sclc92=1 or spro92=1 then cafh92=1; 
              else cafh92=0;
              if asp92=1 then asp92=1;else asp92=0;run;

%hp94(keep=hbp94 db94 chol94 trig94 mi94 ang94 str94 cabg94 pe94
              colc94  pros94  lymp94  ocan94  mel94  
              asp94 mar94 livng94
              betab94 lasix94 thiaz94 ccblo94 ald94 
              chrx94);
              if asp94=1 then asp94=1;else asp94=0;run;


%h94_dt(keep=  mvt94d mvt94  ); mvt94=mvt94d; run;

%hp96(keep=hbp96 db96 chol96 trig96 mi96 ang96 str96 cabg96 pe96 
              colc96  pros96  lymp96  ocan96  mel96 
              asp96 aspd96 mvt96 mar96  livar96 
              betab96 lasix96 thiaz96 ccblo96 ald96 
              chrx96
              pclc96 sclc196 sclc296 fpro96 bpro196 bpro296 sbrc196 mbrcn96 cafh96); 
              if aspd96 in (2,3,4,5) then asp96=1; else asp96=0;
              if pclc96=1 or sclc196=1 or sclc296=1 or fpro96=1 or bpro196=1 or bpro296=1 or sbrc196=1 or mbrcn96=1 then cafh96=1; 
              else cafh96=0;
              run;

%hp98(keep=hbp98 db98 chol98 trig98 mi98 ang98 str98 cabg98 pe98 
              colc98  pros98  lymp98  ocan98  mel98 
              mar98  livng98  asp98 aspfr98
              betab98 lasix98 thiaz98 ccblo98 ald98 
              chrx98);
              if aspfr98 in (2,3,4,5,6) then asp98=1; else asp98=0;run;


%h98_dt(keep=  mvt98d mvt98 );
mvt98=mvt98d;
run;


%hp00(keep=hbp00 db00 chol00 trig00 mi00 ang00 strk00 cabg00 pe00 /* str=strk */
            colc00  pros00  lymp00  ocan00  mel00 
            asp00 mvt00 mar00 livng00
            betab00 lasix00 thiaz00 ccblo00 oanth00
            chrx00 ochrx00);

%hp02(keep=hbp02 db02 chol02 trig02 mi02 ang02 strk02 cabg02 pe02 /* str=strk */
              colc02  pros02  lymp02  ocan02  mel02 
              betab02 lasix02 thiaz02 calcb02 anthp02  
              asp02 mvt02 mar02 lalon02
              stat02 ochrx02 ); 

%hp04(keep=hbp04 db04 chol04 trig04 mi04 mi1d04 ang04 strk04 cabg04 pe04
              colc04  pros04  lymp04  ocan04  mel04 
              mar04 lalon04 asp04 mvt04 
              betab04 ace04 lasix04 thiaz04 calcb04 anthp04 
              stat04 mev04 zoc04 crest04 prav04 lip04 ostat04 ochrx04);

%hp06(keep=hbp06 hbpd06 db06 chol06 trig06 mi06 ang06 strk06 cabg06  pe06
              colc06  pros06  lymp06  ocan06  mel06 
              betab06 ace06 lasix06 thiaz06 calcb06 anthp06 
              asp06 mvt06 mar06 lalon06 
              stat06 mev06 zoc06 crest06 prav06 lip06 ostat06 ochrx06);

%hp08(keep=hbp08 db08 chol08 trig08 mi08 ang08 strk08 cabg08 pe08
              colc08  pros08  lymp08  ocan08  mel08
              sclc108 pmel08 smel08 mpnc08 fpnc08 spnc08 bpnc08 cafh08
              asp08 mvt08  mar08 lalon08 
              betab08 ace08 lasix08 thiaz08 calcb08 arb08 anthp08 
              stat08 mev08 zoc08 crest08 prav08 lip08 ostat08 ochrx08);
              if sclc108=1 or pmel08=1 or smel08=1 or mpnc08=1 or fpnc08=1 or spnc08=1 or bpnc08=1 then cafh08=1; 
              else cafh08=0;
              run;

%hp10(keep=hbp10 db10  chol10  trig10  cabg10 ang10 mi10 strk10 pe10 
              colc10  pros10  blad10  mel10   lymp10  ocan10  /* bladder cancer new */
              betab10 calcb10 ace10   arb10   thiaz10 lasix10 anthp10 
              stat10  mev10   zoc10   crest10 prav10  lip10   ostat10 ochrx10 
              mar10  lalon10   mvt10 asp10);
             
%hp12(keep= db12  hbp12 chol12 mi12  ang12 cabg12  strk12  
colc12  pros12   blad12  mel12 lymp12  ocan12    
betab12 calcb12 ace12   arb12   thiaz12 lasix12 anthp12
stat12 mev12  zoc12 crest12  prav12  lip12  ostat12 ochrx12
mvt12 asp12 mar12 lalon12);   *renal failure no longer asked asfter 2010;


%hp14(keep=  db14  hbp14  chol14 mi14  ang14  cabg14  strk14 
colc14  pros14   blad14  mel14 lymp14  ocan14  
betab14 calcb14 ace14   arb14   thiaz14 lasix14 anthp14
stat14 mev14 prav14 zoc14 crest14  lip14  ostat14 ochrx14
mvt14 asp14 mar14 lalon14);

%hp16(keep= db16    hbp16    chol16 mi16    ang16    cabg16  stk16 /* change */
colc16  pros16   blad16  mel16 lymp16  ocan16  
/* multivitamin emphysema/COPD NOT asked in 2016 */ 
bblock16 cblock16 aceinhb16  angio16  thiaz16 lasix16 antihy16
stat16 prava16 zocor16 crestor16  lipit16  othmed16 ochrx16
asp16 mar16 lalon16); 


%hp18(keep= db18    hbp18    chol18 mi18    ang18    cabg18  strk18 
colc18  pros18   blad18  mel18 lymp18  ocan18   
mvt18  asp18  mar18 lalon18
bblock18 cblock18 aceinhb18  angio18  thiaz18 lasix18 antihy18
stat18 prava18 zocor18 crestor18  lipit18  othmed18 ochrx18);


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
   merge hp_der(in=mst) hp_der_2 hp_dead crc_hp(in=crc) cvd16 diabetes hp_mets
                        hp86 hp88 hp90 hp92 hp94 hp96 hp98 hp00 hp02 hp04 hp06 hp08 hp10 hp12 hp14 hp16 hp18 
                        h86_dt h90_dt h94_dt h98_dt
                        h86_nts h90_nts h94_nts h98_nts h02_nts h06_nts h10_nts h14_nts h18_nts 
                        ahei86 ahei90 ahei94 ahei98 ahei02 ahei06 ahei10 ahei14 ahei18 
                        nses8818 end=_end_;
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

      %indic3(vbl=grade_cr,prefix=grade_cr,min=1, max=4,reflev=1,missing=.,usemiss=0,
          label1='Well Differentiated',
          label2='Grade 2 Moderately Diff',
          label3='Grade 3 Poorly Diff',
          label4='Unspecified');


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
   if cafh86=1 or cafh90=1 or cafh92=1 or cafh96=1 or cafh08=1 then cafh=1; else cafh=0;


/* cancer */
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


/* chronic disease status (flag for stopping diet updates) */
   array hbp     {*}  hbp86   hbp88   hbp90   hbp92   hbp94   hbp96   hbp98   hbp00    hbp02    hbp04    hbp06     hbp08   hbp10  hbp12  hbp14  hbp16  hbp18; 
   array chol    {*}  chol86  chol88  chol90  chol92  chol94  chol96  chol98  chol00   chol02   chol04   chol06    chol08  chol10 chol12 chol14 chol16 chol18;
   array diab    {*}  db86    db88    db90    db92    db94    db96    db98    db00     db02     db04     db06      db08    db10   db12   db14   db16   db18;
   array mi      {*}  mi86    mi88    mi90    mi92    mi94    mi96    mi98    mi00     mi02     mi04     mi06      mi08    mi10   mi12   mi14   mi16   mi18;
   array cabg    {*}  cabg86  cabg88  cabg90  cabg92  cabg94  cabg96  cabg98  cabg00   cabg02   cabg04   cabg06    cabg08  cabg10 cabg12 cabg14 cabg16 cabg18;
   array strk    {*}  str86   str88   str90   str92   str94   str96   str98   strk00   strk02   strk04   strk06    strk08  strk10 strk12 strk14 strk16 strk18;
   array canc_prev {*} can86  can88   can90   can92   can94   can96   can98   can00    can02    can04    can06     can08   can10  can12  can14 can16  can18;
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
        if mi(i) =1 or strk(i)=1 or cabg(i)=1  or diab(i)=1  /*or canc_prev(i)=1*/  then chrodis(i)=1; else chrodis(i)=0; 
    end;

    do i = 2 to dim(chrodis); if chrodis{i-1} = 1 then chrodis{i} = 1; 
    end; /* Make chrodis cumulative */

if nSES_88 = . then nSES_88m = 0.0296195; else nSES_88m = nSES_88;


/***************** Diet scores: carry forward + STOP updating after cancer/CVD *******************/
/* original, simple updated  scores */  
array irt   {*}  rtmnyr86 rtmnyr90 rtmnyr94 rtmnyr98 rtmnyr02 rtmnyr06 rtmnyr10 rtmnyr14 rtmnyr18;
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
array alcon {*}  alco86n  alco90n  alco94n  alco98n  alco02n  alco06n  alco10n alco14n alco18n;

array    carbp  {*}   carbo86p         carbo90p        carbo94p         carbo98p         carbo02p        carbo06p         carbo10p        carbo14p  carbo18p;
array    hcarbp  {*}   hcarbo86p        hcarbo90p       hcarbo94p        hcarbo98p        hcarbo02p       hcarbo06p        hcarbo10p       hcarbo14p hcarbo18p ;
array sumref_carbp {*}  sumref_carb86p   sumref_carb90p  sumref_carb94p   sumref_carb98p   sumref_carb02p  sumref_carb06p   sumref_carb10p  sumref_carb14p  sumref_carb18p;
array     protp  {*}   prot86p          prot90p         prot94p          prot98p          prot02p         prot06p          prot10p         prot14p  prot18p;
array    aprotp  {*}   aprot86p         aprot90p        aprot94p         aprot98p         aprot02p        aprot06p         aprot10p        aprot14p aprot18p ;
array    vprotp  {*}   vprot86p         vprot90p        vprot94p         vprot98p         vprot02p        vprot06p         vprot10p        vprot14p  vprot18p ;
array    tfatp   {*}   tfat86p          tfat90p         tfat94p          tfat98p          tfat02p         tfat06p          tfat10p         tfat14p  tfat18p;
array    afatp   {*}   afat86p          afat90p         afat94p          afat98p          afat02p         afat06p          afat10p         afat14p   afat18p;
array    vfatp   {*}   vfat86p          vfat90p         vfat94p          vfat98p          vfat02p         vfat06p          vfat10p         vfat14p   vfat18p;

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
array alcons {*} alco86ns alco90ns alco94ns alco98ns alco02ns alco06ns alco10ns alco14ns alco18ns;

array    carbs  {*}   carbo86s         carbo90s        carbo94s         carbo98s         carbo02s        carbo06s         carbo10s        carbo14s  carbo18s ;
array    hcarbs  {*}   hcarbo86s        hcarbo90s       hcarbo94s        hcarbo98s        hcarbo02s       hcarbo06s        hcarbo10s       hcarbo14s hcarbo18s ;
array sumref_carbs {*}  sumref_carb86s   sumref_carb90s  sumref_carb94s   sumref_carb98s   sumref_carb02s  sumref_carb06s   sumref_carb10s  sumref_carb14s  sumref_carb18s ;
array     prots  {*}   prot86s          prot90s         prot94s          prot98s          prot02s         prot06s          prot10s         prot14s  prot18s ;
array    aprots  {*}   aprot86s         aprot90s        aprot94s         aprot98s         aprot02s        aprot06s         aprot10s        aprot14s aprot18s ;
array    vprots  {*}   vprot86s         vprot90s        vprot94s         vprot98s         vprot02s        vprot06s         vprot10s        vprot14s  vprot18s ;
array    tfats   {*}   tfat86s          tfat90s         tfat94s          tfat98s          tfat02s         tfat06s          tfat10s         tfat14s  tfat18s ;
array    afats   {*}   afat86s          afat90s         afat94s          afat98s          afat02s         afat06s          afat10s         afat14s   afat18s ;
array    vfats   {*}   vfat86s          vfat90s         vfat94s          vfat98s          vfat02s         vfat06s          vfat10s         vfat14s   vfat18s ;


/* nominal FFQ dates (Jan of each FFQ year, months since 1900), used when irt is missing */
array ffqdt {9} _temporary_ (1033 1081 1129 1177 1225 1273 1321 1369 1417 );   /* 1986 1990 1994 1998 2002 2006 2010 2014 2018 */


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
array diet{9,13}
    calor86n ahei86 lcs86  alcs86 vlcs86  ulcs86  hlcs86 lfs86    ulfs86 hlfs86 alfs86 vlfs86  alco86n
    calor90n ahei90 lcs90  alcs90 vlcs90  ulcs90  hlcs90 lfs90    ulfs90 hlfs90 alfs90 vlfs90  alco90n
    calor94n ahei94 lcs94  alcs94 vlcs94  ulcs94  hlcs94 lfs94    ulfs94 hlfs94 alfs94 vlfs94  alco94n
    calor98n ahei98 lcs98  alcs98 vlcs98  ulcs98  hlcs98 lfs98    ulfs98 hlfs98 alfs98 vlfs98  alco98n
    calor02n ahei02 lcs02  alcs02 vlcs02  ulcs02  hlcs02 lfs02    ulfs02 hlfs02 alfs02 vlfs02  alco02n
    calor06n ahei06 lcs06  alcs06 vlcs06  ulcs06  hlcs06 lfs06    ulfs06 hlfs06 alfs06 vlfs06  alco06n
    calor10n ahei10 lcs10  alcs10 vlcs10  ulcs10  hlcs10 lfs10    ulfs10 hlfs10 alfs10 vlfs10  alco10n
    calor14n ahei14 lcs14  alcs14 vlcs14  ulcs14  hlcs14 lfs14    ulfs14 hlfs14 alfs14 vlfs14  alco14n
    calor18n ahei18 lcs18  alcs18 vlcs18  ulcs18  hlcs18 lfs18    ulfs18 hlfs18 alfs18 vlfs18  alco18n
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

    array dietin{9,13}
    calor86ns ahei86s lcs86s  alcs86s vlcs86s  ulcs86s  hlcs86s lfs86s    ulfs86s hlfs86s alfs86s vlfs86s  alco86ns
    calor90ns ahei90s lcs90s  alcs90s vlcs90s  ulcs90s  hlcs90s lfs90s    ulfs90s hlfs90s alfs90s vlfs90s  alco90ns
    calor94ns ahei94s lcs94s  alcs94s vlcs94s  ulcs94s  hlcs94s lfs94s    ulfs94s hlfs94s alfs94s vlfs94s  alco94ns
    calor98ns ahei98s lcs98s  alcs98s vlcs98s  ulcs98s  hlcs98s lfs98s    ulfs98s hlfs98s alfs98s vlfs98s  alco98ns
    calor02ns ahei02s lcs02s  alcs02s vlcs02s  ulcs02s  hlcs02s lfs02s    ulfs02s hlfs02s alfs02s vlfs02s  alco02ns
    calor06ns ahei06s lcs06s  alcs06s vlcs06s  ulcs06s  hlcs06s lfs06s    ulfs06s hlfs06s alfs06s vlfs06s  alco06ns
    calor10ns ahei10s lcs10s  alcs10s vlcs10s  ulcs10s  hlcs10s lfs10s    ulfs10s hlfs10s alfs10s vlfs10s  alco10ns
    calor14ns ahei14s lcs14s  alcs14s vlcs14s  ulcs14s  hlcs14s lfs14s    ulfs14s hlfs14s alfs14s vlfs14s  alco14ns
    calor18ns ahei18s lcs18s  alcs18s vlcs18s  ulcs18s  hlcs18s lfs18s    ulfs18s hlfs18s alfs18s vlfs18s  alco18ns
    ;
    array dietout{9,13}
    calor86nsv  ahei86sv  lcs86sv alcs86sv vlcs86sv  ulcs86sv  hlcs86sv  lfs86sv  ulfs86sv  hlfs86sv  alfs86sv  vlfs86sv alco86sv
    calor90nsv  ahei90sv  lcs90sv alcs90sv vlcs90sv  ulcs90sv  hlcs90sv  lfs90sv  ulfs90sv  hlfs90sv  alfs90sv  vlfs90sv alco90sv
    calor94nsv  ahei94sv  lcs94sv alcs94sv vlcs94sv  ulcs94sv  hlcs94sv  lfs94sv  ulfs94sv  hlfs94sv  alfs94sv  vlfs94sv alco94sv
    calor98nsv  ahei98sv  lcs98sv alcs98sv vlcs98sv  ulcs98sv  hlcs98sv  lfs98sv  ulfs98sv  hlfs98sv  alfs98sv  vlfs98sv alco98sv
    calor02nsv  ahei02sv  lcs02sv alcs02sv vlcs02sv  ulcs02sv  hlcs02sv  lfs02sv  ulfs02sv  hlfs02sv  alfs02sv  vlfs02sv alco02sv
    calor06nsv  ahei06sv  lcs06sv alcs06sv vlcs06sv  ulcs06sv  hlcs06sv  lfs06sv  ulfs06sv  hlfs06sv  alfs06sv  vlfs06sv alco06sv
    calor10nsv  ahei10sv  lcs10sv alcs10sv vlcs10sv  ulcs10sv  hlcs10sv  lfs10sv  ulfs10sv  hlfs10sv  alfs10sv  vlfs10sv alco10sv
    calor14nsv  ahei14sv  lcs14sv alcs14sv vlcs14sv  ulcs14sv  hlcs14sv  lfs14sv  ulfs14sv  hlfs14sv  alfs14sv  vlfs14sv alco14sv
    calor18nsv  ahei18sv  lcs18sv alcs18sv vlcs18sv  ulcs18sv  hlcs18sv  lfs18sv  ulfs18sv  hlfs18sv  alfs18sv  vlfs18sv alco18sv
    ;

/* calculate cumulative average intake */
do j=1 to 13;
  dietout{1,j}=dietin{1,j};
 do i=2 to 9;
  sumvar=0;
   n=0;
     do k=1 to i;
      if (dietin{k,j} ne .) then do;
      n=n+1; sumvar=sumvar+dietin{k,j};
      end;
     end;
    if n=0 then dietout{i,j}=.;     /* no post-dx FFQ yet -> missing */
   else dietout{i,j}=sumvar/n;
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
   if n=0 then macronv{i,j}=.;
   else macronv{i,j}=sumvar/n;
  end;
end;
run;


proc datasets nolist;
   delete hp_der hp_der_2  diabetes cvd16 hp_dead hp_mets allcancer crc_hp
   hp86 hp88 hp90 hp92 hp94 hp96 hp98 hp00 hp02 hp04 hp06 hp08 hp10 hp12 hp14 hp16 hp18
   h86_nts h90_nts h94_nts h98_nts h02_nts h06_nts h10_nts h14_nts h18_nts 
   h86_dt h90_dt h94_dt h98_dt 
   ahei86 ahei90 ahei94 ahei98 ahei02 ahei06 ahei10 ahei18 ahei14
   nses8818 ;
   run;

/*************************************************************************************/
/*                          Large section of arrays                                  */
/*************************************************************************************/
data hpfs;
     set base end=_end_;

cutoff=1441; /* 2020.01 */
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
array nSESArr{*} nSES_88m nSES_88m nSES_90 nSES_92 nSES_94 nSES_96 nSES_98 nSES_00 nSES_02 nSES_04 nSES_06 nSES_08 nSES_10 nSES_12 nSES_14 nSES_16 nSES_18 ;
array conf1119a  {17};
array act   {*} act86 act88 act90 act92 act94 act96 act98 act00 act02 act04 act06 act08 act10 act12 act14 act16 act18;
array mvit  {*} mvt86 mvt88 mvt90 mvt92 mvt94 mvt96 mvt98 mvt00 mvt02 mvt04 mvt06 mvt08 mvt10 mvt12 mvt14 mvt16 mvt18; 
array aspi  {*} asp86 asp88 asp90 asp92 asp94 asp96 asp98 asp00 asp02 asp04 asp06 asp08 asp10 asp12 asp14 asp16 asp18;
*array mar  {*} mar86 mar88  mar90  mar92  mar94  mar96  mar98  mar00  mar02  mar04  mar06  mar08 mar10 mar12 mar14 mar16 mar18;
*array alone{*} livng86 livng88 livng90 livng92 livng94 livng96 livar96 livng00  lalon02  lalon04  lalon06 lalon08 lalon10 lalon12 lalon14 lalon16 lalon18;
array hbp   {*}  hbp86   hbp88   hbp90   hbp92   hbp94   hbp96   hbp98   hbp00    hbp02    hbp04    hbp06     hbp08   hbp10  hbp12  hbp14  hbp16  hbp18; 
array chol  {*}  chol86  chol88  chol90  chol92  chol94  chol96  chol98  chol00   chol02   chol04   chol06    chol08  chol10 chol12 chol14 chol16 chol18;
array antihp {*} htnrx86 htnrx88 htnrx90 htnrx92 htnrx94 htnrx96 htnrx98 htnrx00 htnrx02 htnrx04 htnrx06 htnrx08 htnrx10 htnrx12 htnrx14 htnrx16 htnrx18;
array antitc {*} hchtx86 hchtx88 hchtx90 hchtx92 hchtx94 hchtx96 hchtx98 hchtx00 hchtx02 hchtx04 hchtx06 hchtx08 hchtx10 hchtx12 hchtx14 hchtx16 hchtx18;
array diab   {*} db86    db88    db90    db92    db94    db96    db98    db00     db02     db04     db06      db08    db10   db12   db14   db16   db18;
array mi    {*} mi86    mi88    mi90    mi92    mi94    mi96    mi98    mi00     mi02     mi04     mi06      mi08    mi10   mi12   mi14   mi16   mi18;
array cabg  {*} cabg86  cabg88  cabg90  cabg92  cabg94  cabg96  cabg98  cabg00   cabg02   cabg04   cabg06    cabg08  cabg10 cabg12 cabg14 cabg16 cabg18;
array strk  {*} str86   str88   str90   str92   str94   str96   str98   strk00   strk02   strk04   strk06    strk08  strk10 strk12 strk14 strk16 strk18;
array canc   {*} can86  can88   can90   can92   can94   can96   can98   can00    can02    can04    can06     can08   can10  can12  can14 can16  can18;

/*************************************************************************
***** If an irt date is before June of that qq year or after or equal ****
***** to the next qq year it is incorrect and should be defaulted to  ****
***** June of that qq year.    Make time period indicator tvar=0.     ****
*************************************************************************/

   /* Create ltf using the original irt values */
ltf=.; do i=1 to dim(irt)-1; if irt{i}>0 then ltf=irt{i}; end;
if ltf=irt{dim(irt)-1} then ltf=.;

/* Correct irt values and calculate age */
do i=1 to dim(irt)-1;
    if (irt{i} <  (1009+24*i)) or (irt{i} >= (1033+24*i)) then irt{i}=1009+24*i;
    age{i}=int((irt{i}-dbmy09)/12);
end;

/* There are a few cases with date of diagnosis after date of death.  Since were more sure about date of death, use that date before we assign period of diagnosis */
    if 0<dtdxca<9999 and 0<dtdth<9999 then do;
        if dtdxca>dtdth then dtdxca=dtdth;
    end;

do j=2 to DIM(bmi);
      if age{j}=.      then age{j}=age{j-1}+2;
      *if bmi{j}=.      then bmi{j}=bmi{j-1};
      *if smok{j}=.     then smok{j}=smok{j-1};
      if cignm{j}=.    then cignm{j}=cignm{j-1};
      *if alcon{j}=.    then alcon{j}=alcon{j-1};
      if act{j}=.      then act{j}=act{j-1};
      if mvit{j}=.     then mvit{j}=mvit{j-1};
      if aspi{j}=.     then aspi{j}=aspi{j-1};
      if hbp{j}=.      then hbp{j}=hbp{j-1};
      if chol{j}=.     then chol{j}=chol{j-1};
      if diab{j}=.     then diab{j}=diab{j-1};
      if nSESArr{j}=.  then nSESArr{j}=nSESArr{j-1};
end;

do i=2 to DIM(bmi);
if irt{i-1}<=dtdxca and dtdxca<irt{i} then do;
  *age_dx=age{i};
  age_dx=age{i-1} + (dtdxca - irt{i-1})/12; /* age at diagnosis*/
  bmi_dx=bmi{i-1}; /* BMI prior to diagnosis*/
  *bmi_dx=bmi{i}; 
  highbp_dx=hbp{i};
  hightc_dx=chol{i};
  db_dx=diab{i};
end;end;
if age_dx = . and dtdxca >= irt{17} then do;
  age_dx=age{17};
  bmi_dx=bmi{17};
  highbp_dx=hbp{17};
  hightc_dx=chol{17};
  db_dx=diab{17};
end;



/* -original  s-stop updated  sv-postdx cum averaged   pre-predx   post-postdx single */
/* -original*/
array calor {17}  calor86n  calor86n calor90n  calor90n calor94n  calor94n calor98n  calor98n calor02n  calor02n calor06n calor06n calor10n calor10n calor14n calor14n calor18n;
array aheim {17}  ahei86 ahei86 ahei90 ahei90 ahei94 ahei94 ahei98 ahei98 ahei02 ahei02 ahei06 ahei06 ahei10 ahei10 ahei14 ahei14 ahei18;   
array lcs   {17}  lcs86    lcs86    lcs90    lcs90   lcs94    lcs94    lcs98    lcs98    lcs02    lcs02    lcs06    lcs06    lcs10 lcs10 lcs14 lcs14 lcs18;
array alcs  {17}  alcs86   alcs86   alcs90   alcs90  alcs94   alcs94   alcs98   alcs98   alcs02   alcs02   alcs06   alcs06   alcs10 alcs10 alcs14 alcs14 alcs18;
array vlcs  {17}  vlcs86   vlcs86   vlcs90   vlcs90  vlcs94   vlcs94   vlcs98   vlcs98   vlcs02   vlcs02   vlcs06   vlcs06   vlcs10 vlcs10 vlcs14 vlcs14 vlcs18;
array ulcs  {17}  ulcs86   ulcs86   ulcs90   ulcs90  ulcs94   ulcs94   ulcs98   ulcs98   ulcs02   ulcs02   ulcs06   ulcs06   ulcs10 ulcs10 ulcs14 ulcs14 ulcs18;
array hlcs  {17}  hlcs86   hlcs86   hlcs90   hlcs90  hlcs94   hlcs94   hlcs98   hlcs98   hlcs02   hlcs02   hlcs06   hlcs06   hlcs10 hlcs10 hlcs14 hlcs14 hlcs18;
array lfs   {17}  lfs86    lfs86    lfs90    lfs90   lfs94    lfs94    lfs98    lfs98    lfs02    lfs02    lfs06    lfs06    lfs10 lfs10 lfs14 lfs14 lfs18;
array ulfs  {17}  ulfs86   ulfs86   ulfs90   ulfs90  ulfs94   ulfs94   ulfs98   ulfs98   ulfs02   ulfs02   ulfs06   ulfs06   ulfs10 ulfs10 ulfs14 ulfs14 ulfs18;
array hlfs  {17}  hlfs86   hlfs86   hlfs90   hlfs90  hlfs94   hlfs94   hlfs98   hlfs98   hlfs02   hlfs02   hlfs06   hlfs06   hlfs10 hlfs10 hlfs14 hlfs14 hlfs18;
array alfs {17} alfs86 alfs86 alfs90 alfs90 alfs94 alfs94 alfs98 alfs98 alfs02 alfs02 alfs06 alfs06 alfs10 alfs10 alfs14 alfs14 alfs18;
array vlfs {17} vlfs86 vlfs86 vlfs90 vlfs90 vlfs94 vlfs94 vlfs98 vlfs98 vlfs02 vlfs02 vlfs06 vlfs06 vlfs10 vlfs10 vlfs14 vlfs14 vlfs18;


/* -s stop updated */
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
array alcons {17} alco86ns alco86ns alco90ns alco90ns alco94ns alco94ns alco98ns alco98ns alco02ns alco02ns alco06ns alco06ns alco10ns alco10ns alco14ns alco14ns alco18ns;

array    carbs  {*}   carbo86s  carbo86s  carbo90s  carbo90s  carbo94s  carbo94s  carbo98s  carbo98s  carbo02s   carbo02s  carbo06s carbo06s  carbo10s carbo10s  carbo14s  carbo14s carbo18s ;
array    hcarbs  {*}   hcarbo86s  hcarbo86s  hcarbo90s  hcarbo90s  hcarbo94s  hcarbo94s  hcarbo98s  hcarbo98s  hcarbo02s   hcarbo02s  hcarbo06s hcarbo06s  hcarbo10s hcarbo10s  hcarbo14s  hcarbo14s hcarbo18s ;
array sumref_carbs {*}  sumref_carb86s  sumref_carb86s  sumref_carb90s  sumref_carb90s  sumref_carb94s  sumref_carb94s  sumref_carb98s  sumref_carb98s  sumref_carb02s   sumref_carb02s  sumref_carb06s sumref_carb06s  sumref_carb10s sumref_carb10s  sumref_carb14s  sumref_carb14s sumref_carb18s ;
array     prots  {*}   prot86s  prot86s  prot90s  prot90s  prot94s  prot94s  prot98s  prot98s  prot02s   prot02s  prot06s prot06s  prot10s prot10s  prot14s  prot14s prot18s ;
array    aprots  {*}   aprot86s  aprot86s  aprot90s  aprot90s  aprot94s  aprot94s  aprot98s  aprot98s  aprot02s   aprot02s  aprot06s aprot06s  aprot10s aprot10s  aprot14s  aprot14s aprot18s ;
array    vprots  {*}   vprot86s  vprot86s  vprot90s  vprot90s  vprot94s  vprot94s  vprot98s  vprot98s  vprot02s   vprot02s  vprot06s vprot06s  vprot10s vprot10s  vprot14s  vprot14s vprot18s ;
array    tfats   {*}   tfat86s  tfat86s  tfat90s  tfat90s  tfat94s  tfat94s  tfat98s  tfat98s  tfat02s   tfat02s  tfat06s tfat06s  tfat10s tfat10s  tfat14s  tfat14s tfat18s ;
array    afats   {*}   afat86s  afat86s  afat90s  afat90s  afat94s  afat94s  afat98s  afat98s  afat02s   afat02s  afat06s afat06s  afat10s afat10s  afat14s  afat14s afat18s ;
array    vfats   {*}   vfat86s  vfat86s  vfat90s  vfat90s  vfat94s  vfat94s  vfat98s  vfat98s  vfat02s   vfat02s  vfat06s vfat06s  vfat10s vfat10s  vfat14s  vfat14s vfat18s ;

/* -sv post-dx, cumulative average */
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
array alcosv {17} alco86sv alco86sv alco90sv alco90sv alco94sv alco94sv alco98sv alco98sv alco02sv alco02sv alco06sv alco06sv alco10sv alco10sv alco14sv alco14sv alco18sv;

array carbsv {17}  carbo86sv carbo86sv carbo90sv carbo90sv carbo94sv carbo94sv carbo98sv carbo98sv carbo02sv carbo02sv carbo06sv carbo06sv carbo10sv carbo10sv carbo14sv carbo14sv carbo18nsv;
array hcarbsv {17}  hcarbo86sv hcarbo86sv hcarbo90sv hcarbo90sv hcarbo94sv hcarbo94sv hcarbo98sv hcarbo98sv hcarbo02sv hcarbo02sv hcarbo06sv hcarbo06sv hcarbo10sv hcarbo10sv hcarbo14sv hcarbo14sv hcarbo18sv;   
array sumref_carbsv   {17}  sumref_carb86sv    sumref_carb86sv    sumref_carb90sv    sumref_carb90sv    sumref_carb94sv    sumref_carb94sv    sumref_carb98sv    sumref_carb98sv    sumref_carb02sv    sumref_carb02sv    sumref_carb06sv    sumref_carb06sv    sumref_carb10sv  sumref_carb10sv  sumref_carb14sv  sumref_carb14sv sumref_carb18sv;
array protsv  {17}  prot86sv   prot86sv   prot90sv   prot90sv   prot94sv   prot94sv   prot98sv   prot98sv   prot02sv   prot02sv   prot06sv   prot06sv   prot10sv prot10sv prot14sv prot14sv prot18sv;
array aprotsv  {17}  aprot86sv   aprot86sv   aprot90sv   aprot90sv   aprot94sv   aprot94sv   aprot98sv   aprot98sv   aprot02sv   aprot02sv   aprot06sv   aprot06sv   aprot10sv aprot10sv aprot14sv aprot14sv aprot18sv;
array vprotsv  {17}  vprot86sv   vprot86sv   vprot90sv   vprot90sv   vprot94sv   vprot94sv   vprot98sv   vprot98sv   vprot02sv   vprot02sv   vprot06sv   vprot06sv   vprot10sv vprot10sv vprot14sv vprot14sv vprot18sv;
array tfatsv  {17}  tfat86sv   tfat86sv   tfat90sv   tfat90sv   tfat94sv   tfat94sv   tfat98sv   tfat98sv   tfat02sv   tfat02sv   tfat06sv   tfat06sv   tfat10sv tfat10sv tfat14sv tfat14sv tfat18sv;
array afatsv {17} afat86sv afat86sv afat90sv afat90sv afat94sv afat94sv afat98sv afat98sv afat02sv afat02sv afat06sv afat06sv afat10sv afat10sv afat14sv afat14sv afat18sv;
array vfatsv {17} vfat86sv vfat86sv vfat90sv vfat90sv vfat94sv vfat94sv vfat98sv vfat98sv vfat02sv vfat02sv vfat06sv vfat06sv vfat10sv vfat10sv vfat14sv vfat14sv vfat18sv;

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
array calor_pre {17}  calor86p_pre calor86p_pre calor90p_pre calor90p_pre calor94p_pre calor94p_pre calor98p_pre calor98p_pre calor02p_pre calor02p_pre calor06p_pre calor06p_pre calor10p_pre calor10p_pre calor14p_pre calor14p_pre calor18p_pre;
array aheim_pre {17}  ahei86_pre ahei86_pre ahei90_pre ahei90_pre ahei94_pre ahei94_pre ahei98_pre ahei98_pre ahei02_pre ahei02_pre ahei06_pre ahei06_pre ahei10_pre ahei10_pre ahei14_pre ahei14_pre ahei18_pre;   
array lcs_pre   {17}  lcs86_pre    lcs86_pre    lcs90_pre    lcs90_pre    lcs94_pre    lcs94_pre    lcs98_pre    lcs98_pre    lcs02_pre    lcs02_pre    lcs06_pre    lcs06_pre    lcs10_pre   lcs10_pre  lcs14_pre lcs14_pre lcs18_pre;
array vlcs_pre  {17}  vlcs86_pre   vlcs86_pre   vlcs90_pre   vlcs90_pre   vlcs94_pre   vlcs94_pre   vlcs98_pre   vlcs98_pre   vlcs02_pre   vlcs02_pre   vlcs06_pre   vlcs06_pre   vlcs10_pre  vlcs10_pre vlcs14_pre vlcs14_pre vlcs18_pre;
array alcs_pre  {17}  alcs86_pre   alcs86_pre   alcs90_pre   alcs90_pre   alcs94_pre   alcs94_pre   alcs98_pre   alcs98_pre   alcs02_pre   alcs02_pre   alcs06_pre   alcs06_pre   alcs10_pre  alcs10_pre alcs14_pre alcs14_pre alcs18_pre;
array ulcs_pre  {17}  ulcs86_pre   ulcs86_pre   ulcs90_pre   ulcs90_pre   ulcs94_pre   ulcs94_pre   ulcs98_pre   ulcs98_pre   ulcs02_pre   ulcs02_pre   ulcs06_pre   ulcs06_pre   ulcs10_pre  ulcs10_pre ulcs14_pre ulcs14_pre ulcs18_pre; 
array hlcs_pre  {17}  hlcs86_pre   hlcs86_pre   hlcs90_pre   hlcs90_pre   hlcs94_pre   hlcs94_pre   hlcs98_pre   hlcs98_pre   hlcs02_pre   hlcs02_pre   hlcs06_pre   hlcs06_pre   hlcs10_pre  hlcs10_pre hlcs14_pre hlcs14_pre hlcs18_pre;
array lfs_pre   {17}  lfs86_pre    lfs86_pre    lfs90_pre    lfs90_pre    lfs94_pre    lfs94_pre    lfs98_pre    lfs98_pre    lfs02_pre    lfs02_pre    lfs06_pre    lfs06_pre    lfs10_pre   lfs10_pre  lfs14_pre lfs14_pre lfs18_pre;
array vlfs_pre  {17}  vlfs86_pre   vlfs86_pre   vlfs90_pre   vlfs90_pre   vlfs94_pre   vlfs94_pre   vlfs98_pre   vlfs98_pre   vlfs02_pre   vlfs02_pre   vlfs06_pre   vlfs06_pre   vlfs10_pre  vlfs10_pre vlfs14_pre vlfs14_pre vlfs18_pre;
array alfs_pre  {17}  alfs86_pre   alfs86_pre   alfs90_pre   alfs90_pre   alfs94_pre   alfs94_pre   alfs98_pre   alfs98_pre   alfs02_pre   alfs02_pre   alfs06_pre   alfs06_pre   alfs10_pre  alfs10_pre alfs14_pre alfs14_pre alfs18_pre;
array ulfs_pre  {17}  ulfs86_pre   ulfs86_pre   ulfs90_pre   ulfs90_pre   ulfs94_pre   ulfs94_pre   ulfs98_pre   ulfs98_pre   ulfs02_pre   ulfs02_pre   ulfs06_pre   ulfs06_pre   ulfs10_pre  ulfs10_pre ulfs14_pre ulfs14_pre ulfs18_pre; 
array hlfs_pre  {17}  hlfs86_pre   hlfs86_pre   hlfs90_pre   hlfs90_pre   hlfs94_pre   hlfs94_pre   hlfs98_pre   hlfs98_pre   hlfs02_pre   hlfs02_pre   hlfs06_pre   hlfs06_pre   hlfs10_pre  hlfs10_pre hlfs14_pre hlfs14_pre hlfs18_pre;
array alco_pre  {17}  alco86_pre   alco86_pre   alco90_pre   alco90_pre   alco94_pre   alco94_pre   alco98_pre   alco98_pre   alco02_pre   alco02_pre   alco06_pre   alco06_pre   alco10_pre  alco10_pre alco14_pre alco14_pre alco18_pre;

/*Post-dx, donnot use the carry forward value to ensure the cumulaitve averaged validity*/
array calor_post {17}  calor86p_post calor86p_post calor90p_post calor90p_post calor94p_post calor94p_post calor98p_post calor98p_post calor02p_post calor02p_post calor06p_post calor06p_post calor10p_post calor10p_post calor14p_post calor14p_post calor18p_post;
array aheim_post {17}  ahei86_post ahei86_post ahei90_post ahei90_post ahei94_post ahei94_post ahei98_post ahei98_post ahei02_post ahei02_post ahei06_post ahei06_post ahei10_post ahei10_post ahei14_post ahei14_post ahei18_post;   
array lcs_post   {17}  lcs86_post    lcs86_post    lcs90_post    lcs90_post    lcs94_post    lcs94_post    lcs98_post    lcs98_post    lcs02_post    lcs02_post    lcs06_post    lcs06_post    lcs10_post   lcs10_post  lcs14_post lcs14_post lcs18_post;
array vlcs_post  {17}  vlcs86_post   vlcs86_post   vlcs90_post   vlcs90_post   vlcs94_post   vlcs94_post   vlcs98_post   vlcs98_post   vlcs02_post   vlcs02_post   vlcs06_post   vlcs06_post   vlcs10_post  vlcs10_post vlcs14_post vlcs14_post vlcs18_post;
array alcs_post  {17}  alcs86_post   alcs86_post   alcs90_post   alcs90_post   alcs94_post   alcs94_post   alcs98_post   alcs98_post   alcs02_post   alcs02_post   alcs06_post   alcs06_post   alcs10_post  alcs10_post alcs14_post alcs14_post alcs18_post;
array ulcs_post  {17}  ulcs86_post   ulcs86_post   ulcs90_post   ulcs90_post   ulcs94_post   ulcs94_post   ulcs98_post   ulcs98_post   ulcs02_post   ulcs02_post   ulcs06_post   ulcs06_post   ulcs10_post  ulcs10_post ulcs14_post ulcs14_post ulcs18_post; 
array hlcs_post  {17}  hlcs86_post   hlcs86_post   hlcs90_post   hlcs90_post   hlcs94_post   hlcs94_post   hlcs98_post   hlcs98_post   hlcs02_post   hlcs02_post   hlcs06_post   hlcs06_post   hlcs10_post  hlcs10_post hlcs14_post hlcs14_post hlcs18_post;
array lfs_post   {17}  lfs86_post    lfs86_post    lfs90_post    lfs90_post    lfs94_post    lfs94_post    lfs98_post    lfs98_post    lfs02_post    lfs02_post    lfs06_post    lfs06_post    lfs10_post   lfs10_post lfs14_post lfs14_post lfs18_post;
array vlfs_post  {17}  vlfs86_post   vlfs86_post   vlfs90_post   vlfs90_post   vlfs94_post   vlfs94_post   vlfs98_post   vlfs98_post   vlfs02_post   vlfs02_post   vlfs06_post   vlfs06_post   vlfs10_post  vlfs10_post vlfs14_post vlfs14_post vlfs18_post;
array alfs_post  {17}  alfs86_post   alfs86_post   alfs90_post   alfs90_post   alfs94_post   alfs94_post   alfs98_post   alfs98_post   alfs02_post   alfs02_post   alfs06_post   alfs06_post   alfs10_post  alfs10_post alfs14_post alfs14_post alfs18_post;
array ulfs_post  {17}  ulfs86_post   ulfs86_post   ulfs90_post   ulfs90_post   ulfs94_post   ulfs94_post   ulfs98_post   ulfs98_post   ulfs02_post   ulfs02_post   ulfs06_post   ulfs06_post   ulfs10_post  ulfs10_post ulfs14_post ulfs14_post ulfs18_post; 
array hlfs_post  {17}  hlfs86_post   hlfs86_post   hlfs90_post   hlfs90_post   hlfs94_post   hlfs94_post   hlfs98_post   hlfs98_post   hlfs02_post   hlfs02_post   hlfs06_post   hlfs06_post   hlfs10_post  hlfs10_post hlfs14_post hlfs14_post hlfs18_post;
array alco_post  {17}  alco86_post   alco86_post   alco90_post   alco90_post   alco94_post   alco94_post   alco98_post   alco98_post   alco02_post   alco02_post   alco06_post   alco06_post   alco10_post  alco10_post alco14_post alco14_post alco18_post;

/* Compute pre, post, and chg */
do i=2 to DIM(calor)+1;
    if irt{i-1}<dtdxca and dtdxca <=irt{i} then do;
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


/* FFQ positions in the 17-element period arrays: 86 90 94 98 02 06 10 14 18 */
array ffqk {9} _temporary_ (1 3 5 7 9 11 13 15 17);

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

/*-crc Post-dx, within 6 month to 4 years */
array calor_crc {17}  calor86p_crc calor86p_crc calor90p_crc calor90p_crc calor94p_crc calor94p_crc calor98p_crc calor98p_crc calor02p_crc calor02p_crc calor06p_crc calor06p_crc calor10p_crc calor10p_crc calor14p_crc calor14p_crc calor18p_crc;
array aheim_crc {17}  ahei86_crc ahei86_crc ahei90_crc ahei90_crc ahei94_crc ahei94_crc ahei98_crc ahei98_crc ahei02_crc ahei02_crc ahei06_crc ahei06_crc ahei10_crc ahei10_crc ahei14_crc ahei14_crc ahei18_crc;   
array lcs_crc   {17}  lcs86_crc    lcs86_crc    lcs90_crc    lcs90_crc    lcs94_crc    lcs94_crc    lcs98_crc    lcs98_crc    lcs02_crc    lcs02_crc    lcs06_crc    lcs06_crc    lcs10_crc   lcs10_crc lcs14_crc lcs14_crc lcs18_crc;
array vlcs_crc  {17}  vlcs86_crc   vlcs86_crc   vlcs90_crc   vlcs90_crc   vlcs94_crc   vlcs94_crc   vlcs98_crc   vlcs98_crc   vlcs02_crc   vlcs02_crc   vlcs06_crc   vlcs06_crc   vlcs10_crc  vlcs10_crc vlcs14_crc vlcs14_crc vlcs18_crc;
array alcs_crc  {17}  alcs86_crc   alcs86_crc   alcs90_crc   alcs90_crc   alcs94_crc   alcs94_crc   alcs98_crc   alcs98_crc   alcs02_crc   alcs02_crc   alcs06_crc   alcs06_crc   alcs10_crc  alcs10_crc alcs14_crc alcs14_crc alcs18_crc;
array ulcs_crc  {17}  ulcs86_crc   ulcs86_crc   ulcs90_crc   ulcs90_crc   ulcs94_crc   ulcs94_crc   ulcs98_crc   ulcs98_crc   ulcs02_crc   ulcs02_crc   ulcs06_crc   ulcs06_crc   ulcs10_crc  ulcs10_crc ulcs14_crc ulcs14_crc ulcs18_crc; 
array hlcs_crc  {17}  hlcs86_crc   hlcs86_crc   hlcs90_crc   hlcs90_crc   hlcs94_crc   hlcs94_crc   hlcs98_crc   hlcs98_crc   hlcs02_crc   hlcs02_crc   hlcs06_crc   hlcs06_crc   hlcs10_crc  hlcs10_crc hlcs14_crc hlcs14_crc hlcs18_crc;
array lfs_crc   {17}  lfs86_crc    lfs86_crc    lfs90_crc    lfs90_crc    lfs94_crc    lfs94_crc    lfs98_crc    lfs98_crc    lfs02_crc    lfs02_crc    lfs06_crc    lfs06_crc    lfs10_crc   lfs10_crc lfs14_crc lfs14_crc lfs18_crc;
array vlfs_crc  {17}  vlfs86_crc   vlfs86_crc   vlfs90_crc   vlfs90_crc   vlfs94_crc   vlfs94_crc   vlfs98_crc   vlfs98_crc   vlfs02_crc   vlfs02_crc   vlfs06_crc   vlfs06_crc   vlfs10_crc  vlfs10_crc vlfs14_crc vlfs14_crc vlfs18_crc;
array alfs_crc  {17}  alfs86_crc   alfs86_crc   alfs90_crc   alfs90_crc   alfs94_crc   alfs94_crc   alfs98_crc   alfs98_crc   alfs02_crc   alfs02_crc   alfs06_crc   alfs06_crc   alfs10_crc  alfs10_crc alfs14_crc alfs14_crc alfs18_crc;
array ulfs_crc  {17}  ulfs86_crc   ulfs86_crc   ulfs90_crc   ulfs90_crc   ulfs94_crc   ulfs94_crc   ulfs98_crc   ulfs98_crc   ulfs02_crc   ulfs02_crc   ulfs06_crc   ulfs06_crc   ulfs10_crc  ulfs10_crc ulfs14_crc ulfs14_crc ulfs18_crc; 
array hlfs_crc  {17}  hlfs86_crc   hlfs86_crc   hlfs90_crc   hlfs90_crc   hlfs94_crc   hlfs94_crc   hlfs98_crc   hlfs98_crc   hlfs02_crc   hlfs02_crc   hlfs06_crc   hlfs06_crc   hlfs10_crc  hlfs10_crc hlfs14_crc hlfs14_crc hlfs18_crc;
array alco_crc  {17}  alco86_crc   alco86_crc   alco90_crc   alco90_crc   alco94_crc   alco94_crc   alco98_crc   alco98_crc   alco02_crc   alco02_crc   alco06_crc   alco06_crc   alco10_crc  alco10_crc alco14_crc alco14_crc alco18_crc;

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
        alco_crc{i}  = crc_alco_fix;
    end;
end;


        elig_any = 0;
        if dtdxca > 0 then do f = 1 to dim(ffqk);  k = ffqk{f};   /* -- FFQ positions */
           nomd = 1009+24*k ;       /* -- NHS ; NHSII 1074+24*k ; HPFS 1009+24*k */
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
    dead_can   = 0; tdead_can=irt{i+1}-irt{i}; dtdx_dead_can=.;
    if irt{i}< dtdth <=irt{i+1} and (140<=newicda<=208) then do; dead_can=1; tdead_can=dtdth-irt{i}; dtdx_dead_can=dtdth; end;

    dead_nocan=0; tdead_nocan=irt{i+1}-irt{i}; dtdx_dead_nocan=.;
    if irt{i}<dtdth<=irt{i+1} and dead_can ne 1 then do; dead_nocan=1; tdead_nocan=dtdth-irt{i};dtdx_dead_nocan=dtdth; end;

    dead_brc   = 0; tdead_brc=irt{i+1}-irt{i}; dtdx_dead_brc=.; /*breast cancer death*/
    if irt{i}< dtdth <=irt{i+1} and (newicda=174) then do; dead_brc=1; tdead_brc=dtdth-irt{i}; dtdx_dead_brc=dtdth; end;

    dead_crc   = 0; tdead_crc=irt{i+1}-irt{i}; dtdx_dead_crc=.; /*colorectal cancer death*/
    if irt{i}< dtdth <=irt{i+1} and newicda in (153,154) then do; dead_crc=1; tdead_crc=dtdth-irt{i}; dtdx_dead_crc=dtdth; end;

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
   %indic3(vbl=crcduration_c, prefix=crcduration_c, min=2, max=4, reflev=1, usemiss=0,
    label1='0-5 years', label2='5-10 years',label3='10-15 yaers',label4='>15 years');

   /* Year of diagnosis */
      if dtdxca=9999 then yeardx_crc=.;
      else yeardx_crc=round(((dtdxca/12)+1900-0.5)); /* Calendar year of diagnosis */



/***************time since diagnosis (post-dx survival time scale) ******************/
            
      if dtdxca > . and conf1119 = 1  then do;
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
    /* 1. single updated score at cycle i */
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

array   smkka   {17}    smkk86   smkk88   smkk90   smkk92   smkk94    smkk96   smkk98   smkk00   smkk02   smkk04   smkk06   smkk08  smkk10 smkk10 smkk14 smkk14 smkk18;
array   smksta  {17}    smksta86   smksta88   smksta90   smksta92   smksta94  smksta96  smksta98  smksta00   smksta02   smksta04   smksta06   smksta08  smksta10 smksta10 smksta14 smksta14 smksta18;

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
         
     length smkk 3 smkstatus 3;
     smkk=smm{i};
     /* store time-varying 5-level category */
     smkka{i} = smkk;
     %indic3(vbl=smkk, prefix=smkk, min=2,max=5, reflev=1, usemiss=0,
              label1='never smker',
              label2='past smker',
              label3='curr 1-14 smker',
              label4='curr 15-24 smker',
              label5='curr 25+ smker');
      
      /* 3-level: 1 never 2 past 3 current */
      if smkk=1 then smkstatus=1;
      else if smkk=2 then smkstatus=2;
      else if smkk>=3 then smkstatus=3;
      else smkstatus=.;
      /* store time-varying 3-level category */
      smksta{i} = smkstatus;
      %indic3(vbl=smkstatus, prefix=smkstatus, min=2, max=3, reflev=1,  usemiss=1);

      if smkstatus=1  then smknever=1; else smknever=0;  

/* Indicator for BMI */
      bmicon=bmi{i};if bmi{i}=0 then bmicon=.;

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
          %indic3(vbl=alcovc_post, prefix=alcovc_post, min=2, max=5, reflev=1, missing=., usemiss=1,
              label1='none',
              label2='0.1-4.9 g/d',
              label3='5.0-14.9 g/d',
              label4='15.0-29.9 g/d',
              label5='30+ g/d');
      
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
            if highbp_dx=1 then diaghbp=1; else diaghbp=0;
               
/****** Indicator for TC lowering ********/
            select(antitc{i});
                when (1)     antihtc=1;
                otherwise    antihtc=0;
            end;

           if chol{i}=1 then hightc=1; else hightc=0;
           if chol{1}=1 then basehightc=1; else basehightc=0;
           if hightc_dx=1 then diaghtc=1; else diaghtc=0;

/****** Indicator for History of diabetes ******/
           if diab{i}=1 then diabetes=1; else diabetes=0;

/* nSES */
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
            mvitd=mvit{b}; aspd=aspi{b};
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
        mvitp=mvit{b-1}; aspp=aspi{b-1}; 
        hbpp=hbp{b-1}; cholp=chol{b-1}; diabp=diab{b-1};
        alcp=alcq{b-1};
        nSESp=nSESArr{b-1};
        aheip=aheim{b-1}; /*lcsp=lcs{b-1}; alcsp=alcs{b-1}; vlcsp=vlcs{b-1}; ulcsp=ulcs{b-1}; hlcsp=hlcs{b-1};*/
        calorp=calor{b-1}; /*lfsp=lfs{b-1}; alfsp=alfs{b-1}; vlfsp=vlfs{b-1}; ulfsp=ulfs{b-1}; hlfsp=hlfs{b-1};*/
 /* Covariates post diagnosis */
        bmicona=bmi{b+1}; agecona=age{b+1}; actcona=act{b+1}; smkkaa=smkka{b+1}; smkstaa=smksta{b+1};
        mvita=mvit{b+1}; aspa=aspi{b+1}; 
        hbpa=hbp{b+1}; chola=chol{b+1}; diaba=diab{b+1};
        alca=alcq{b+1};
        nSESa=nSESArr{b+1};
        aheia=aheim{b+1}; /*lcsa=lcs{b+1}; alcsa=alcs{b+1}; vlcsa=vlcs{b+1}; ulcsa=ulcs{b+1}; hlcsa=hlcs{b+1};*/
        calora=calor{b+1}; /*lfsa=lfs{b+1}; alfsa=alfs{b+1}; vlfsa=vlfs{b+1}; ulfsa=ulfs{b+1}; hlfsa=hlfs{b+1};*/
 /* change  */
        bmiconch=bmi{b+1}-bmi{b-1}; ageconch=age{b+1}-age{b-1}; 
        actconch=act{b+1}-act{b-1}; 
        smkkch=smkka{b+1}-smkka{b-1}; smkstch=smksta{b+1}-smksta{b-1};
        mvitch=mvit{b+1}-mvit{b-1};   aspch=aspi{b+1}-aspi{b-1}; 
        nSESch=nSESArr{b+1}-nSESArr{b-1};
        hbpch=hbp{b+1}-hbp{b-1}; cholch=chol{b+1}-chol{b-1}; diabch=diab{b+1}-diab{b-1};
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
   %exclude(exrec eq 1);      *before 1986;
   %exclude(crc_cases eq 1);  * not in CRC file;
   %exclude(conf1119 ne 1);  

   %exclude(dbmy09 le 0);
   
   %exclude(0 lt dtdth   le irt{i});   
   %exclude(0 lt dtdxca  le irt{i});    
   %exclude(dtdxca ge cutoff); *cases diagnosed after the end of follow-up;

   %exclude(ltf eq irt{1});
   %exclude(age86 eq .);
   %exclude(stage_cr eq 4);    *stage IV CRC;
   %exclude(stage_cr eq .);    *stage missing CRC;
   %exclude(canpre eq 0, skip=t);  * skip individuals diagnosed after current time period;
   %output();
  end;

  else if i> 1 then do;
   %exclude(irt{i-1} lt dtdth le irt{i} );
   %exclude(stage_cr eq 4);    *stage IV CRC;
   %exclude(stage_cr eq .);    *stage missing CRC;  
   %exclude(canpre eq 0, skip=t);  * skip individuals diagnosed after current time period; 
   %output();
  end;
end;

/* END OF DO-LOOP OVER TIME PERIODs */
%endex();


