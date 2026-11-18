/*****************************************************
Purpose:       low carb, low fat, and post obesity-related cancer survival analysis
Study design:  prospective cohort (NHS: 1986-2022.6, NHSII: 1991-2023.6)
Endpoints:     mortality among breast cancer survivors
Exposures:     post-diagnostic cumulative lcs, alcs, vlcs, ulcs, hlcs, lfs, ulfs, hlfs (main)
               post-diagnostic simple update diet score
Inclusion:     - breast cancer / colorectal cancer (analyzed seperately)
               - ICDA codes 153 (colon cancer) and 154 (rectal cancer)
               - confirmation codes 11-19 (including cases identified from death file);
               - diagnosed between 1984 and 2020
               - stage I-III
Exclusion:     had diagnosis of cancer prior to cohort baseline, 
               had missing information on age at enrollment, had missing information on post-diagnosis diet
Covariates:    age, interval, cohort, white, energy quintile, smoking, alcohol, physical activity, family history of cancer,
               hormone use, oral contraceptive use, multivitamin use, aspirin use, 
               BMI at diagnosis, pre-diagnosis diet score, 
               cancer characteristics (clinical stage, estrogen receptor status, her2 status, radiation therapy, chemotherapy, and hormonal treatment for cancer)
revision:      add SES and surgery to covariates (1 =mastectomy, 2 =lumpectomy, 3=unknown type, 4=no surgery, . =missing)

Date:  10/2026

******************************************************/


*------------------- Read case file ---------------------*;
/************** Breast Cancer **************/
/* load breast cancer dataset: br_marker_treatment*/   
data br_marker_treatment;
    %include 'br7624.052926.cases.input';
    casefile = 1;
    dtdx = dxmonth;

    **ER PR STATUS**;
****************;   
*ER Positive*;                        if era_results in (1,3) then er_status = 1; *ER positive - if borderline, code as positive because cutoffs for positivity used to be more stringent than they are now, according to Wendy Chen*/    
*ER Negative*;                   else if era_results = 2 then er_status = 0;
*ER Unknown*;                    else                         er_status = 9;      
    

*PR Positive*;                         if pra_results in (1,3) then pr_status = 1;
*PR Negative*;                    else if pra_results = 2 then pr_status = 0;
*PR Unknown*;                     else                         pr_status = 9;

*HER2 NEU STATUS*;
*****************;
*H Positive*;                          if her2neu_results=1 then her2_status=1;
*H Negative*;                     else if her2neu_results=2 then her2_status=0;
*H Unknown*;                      else                           her2_status=9;


**TUMOR SIZE**;
**************;
*0.1 to <=1.0 cm*;       if tumor_size=7        then tmsize=1;  
*1.1 to <=2.0 cm*;  else if tumor_size in (8,2) then tmsize=2;  
*2.1 to <=4.0 cm*;  else if tumor_size=3        then tmsize=3;  
*4.1+ cm*;          else if tumor_size=4        then tmsize=4;  


**LYMPHY NODE INVOLVEMENT**;
***************************;
*No nodes  involved*;                    if lymph_nodes=2        then node=1; 
*1-3 nodes involved*;               else if lymph_nodes=3        then node=2;
*4-9 nodes involved*;               else if lymph_nodes=4        then node=3; 
*10+ nodes involved or metastatic*; else if lymph_nodes in (5,6) then node=4; 

**GRADE**;
*********;
*Grade
          1=Predominantly well-differentiated
          2=Moderately differentiated / balanced pattern
          3=Poorly differentiated
          4=Not applicable, not given / in-situ
*;

rename tumor_grade= grade_path;

**HOW DISCOVERED**;
******************;

* firstdetect              i1             20 how the cancer discovered;
* 1.self exam;
* 2.health professional exam;
* 3.husband or other non health professional exam;
* 4.routine screening (mammography);
  

if (firstdetect<1 or firstdetect>4) then firstdetect=.;


**STAGE**;
*********;

* Stage Deriviation (based on code from Rong Hu)*;

*Jae H Kang changed this code slightly so that we are not assuming missing nodes are "no nodal involvement" and missing size is no longer "0.1-1 cm")*;

*Stage 1 = tumor size<=2cm and no nodal involvement*;                                      
if   tumor_size in (2,7,8) and lymph_nodes=2            then stage=1;   
 
*Stage 2 = tumor size<=2cm & 1-3 nodes or 2-4cm & 0-3 nodes or 4+cm & 0 nodes*;       
else if ((tumor_size in (2,7,8) and lymph_nodes=3) or (tumor_size=3 and lymph_nodes in (2,3)) or 
(tumor_size=4 and lymph_nodes=2))                       then stage=2; 
  
*Stage 3 = tumor size<=2cm & 4+  nodes or 2-4cm & 4+  nodes or >4cm & 1+ nodes*;      
else if ((tumor_size in (2,7,8) and lymph_nodes in (4,5)) or (tumor_size=3 and lymph_nodes in (4,5)) or 
(tumor_size=4 and lymph_nodes in (3,4,5))) then stage=3; 
  
*Stage 4 = 10+ nodes involved or metastatic*;                                         
else if  lymph_nodes=6  then stage=4;   
else stage =5;  /* Unknown/unspecified */

*INVASIVENESS*;
**************;

*canc_in_situ;
*        1.Ducts;
*        2.Lobules;
*        3.Other;
*        4.Both intraductal and lobular;
*        6.Intraductal - no mention;
*        7.Intraductal - possible invasion;
*        8.Intraductal - early microinvasion;
*        9.Not in-situ;

*invasive:
1   no primary structure evident
2   ductal
3   lobular
4   both ductal and lobular
8   neither
5   (in-situ)
9   not specified
*;

       if canc_in_situ=1                then inv_type=1; *DCIS*;
  else if canc_in_situ=2                then inv_type=2; *LCIS*;
  else if canc_in_situ=8 and invasive=2 then inv_type=4; *invasive ductal*;
  else if invasive=2                    then inv_type=4; *invasive ductal*;
  else if invasive=3                    then inv_type=5; *invasive lobular*;
  else if invasive=4                    then inv_type=6; *invasive both*;
  else if canc_in_situ in (3,4,6,7)     then inv_type=3; *other in situ lesion*;
  else if invasive=1                    then inv_type=7; *invasive not specified*;
run;
proc sort; by id; run;
              


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

        if cohort ne 2 then delete; /*include NHS2*/
        if conf in (11,12,13,14,15,16,17,18,19) then ca=1; /* ca=1: confirmed cancer ; exclude in situ ca by taking out conf code 25*/

        dtdx_allca=.;
        if dtdxca>0 and ca=1 and (140<=icda2<=172 or 174<=icda2<=207 or icda2=225) then dtdx_allca=dtdxca;

        if icda2=174 and ca=1 then do; brca_ca=1;                dtdx_brca      = dtdxca; end;

        keep id conf dtdx_allca brca_ca dtdx_brca;
     run;
     proc sort data=allcancer nodupkey; 
          by id dtdx_allca; 
     run;




data allcancer;
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
dtdx_allca=.;
if dxmonth>0 and ca=1 and (140<=icda2<=172 or 174<=icda2<=207 or icda2 =225) then dtdx_allca=dxmonth;

if icda2=174 and ca=1 then do; brca_ca=1;                dtdx_brca      = dxmonth; end;

keep id conf dtdx_allca brca_ca dtdx_brca;
run;

proc sort data=allcancer nodupkey;
    by id dtdx_allca;
run;



                                                                            
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


/*  diabetes 2024 */
data diabetes;
infile 'db7624.110123.cases' lrecl=183 recfm=d;

input
@1    id        6.
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
    merge dead diabetes update_cabg mi tot_str;
    by id;

    /* CVD */
    dtdxchd=.;
    if 0<dtdxmi<9999 or  0<dtdxcabg<9999 then dtdxchd=min(dtdxmi,dtdxcabg, .);
    
    dtdxcvd=.;
    if 0<dtdxchd<9999 or  0<dtdxstr<9999 then dtdxcvd=min(dtdxchd,dtdxstr,.);

    keep id dtdth  newicda newicda3 dtdxdb2 dtdxcvd dtdxchd dtdxstr dtdxmi;
    run; proc sort; by id; run;

    PROC DATASETS;
    delete diabetes dead mi cabg update_cabg nur14 nur16 nur20 tot_str; run;





*----------------2.Exposures-----------------------------;
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




*----------------3.Covariables-----------------------------;
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
                 mclc88 fclc88 bclc88 sclc88 mbrcn88 sbrcn88 mcanc88 fcanc88 bcanc88 scanc88 cafh88);
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

               hxbrc76 hxbrc78 hxbrc80 hxbrc82 hxbrc84 hxbrc86 hxbrc88 hxbrc90 hxbrc92 
               hxbrc94 hxbrc96 hxbrc98 hxbrc00 hxbrc02 hxbrc04 hxbrc06 hxbrc08 hxbrc10 
               hxbrc12 hxbrc14 hxbrc16  /*Family (mother, sister) history of Breast Cancer*/
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
                db80 hbp80 chol80
                wt76 ht76  wt80);

%nur82(keep=db82 hbp82 chol82 str82 wt82 ); /* cabg82 not asked */

%nur84(keep=  db84 hbp84 chol84 cabg84 str84 wt84); /* CODING CHANGE: 0=blank, 1=no, 2=yes, range=0-2 */
       array fix db84 hbp84 chol84 cabg84 str84;
       do over fix; fix=fix-1; end;

%nur86(keep=db86 hbp86 chol86 cabg86 str86 wt86);

%n86_dt(keep=  waist86  hip86 nblnk86d  cerbr86d);

%nur88(keep=db88 hbp88 chol88 str88  cabg88 wt88
            thiaz88 ccblo88 betab88 ace88 bprx88 /* medications asked */
            clrx88
            insul88 odbrx88 /* Hypoglycemic */
            );

%nur90(keep=db90 hbp90 chol90 cabg90 str90 wt90);

%nur92(keep=db92 hbp92 chol92 cabg92 str92 alone92 wt92);

%nur94(keep=db94 hbp94 chol94 cabg94 str94 wt94
            thiaz94 lasix94 ccblo94 betab94 bprx94 
            couma94 digox94 
            clrx94 
            tag94 zantc94 /* PPIs */
            insul94 ohypo94 /* Hypoglycemic */
            );

%nur96(keep=db96 hbp96 chol96 cabg96 str96  alone96  marry96 wt96
            thiaz96  lasix96  ccblo96  betab96  ace96  bprx96  
            couma96  digox96  antia96 /* anticoagulant  heart-failure antiarrhythmic */
            chrx96
            antid96 /* antidepressant */
            tag96 h2blk96 /* PPIs */
            insul96 ohypo96 /* Hypoglycemic */
            );

%nur98(keep=db98 hbp98 chol98 cabg98 str98  wt98
            thiaz98  lasix98  ccblo98  betab98  ace98 bprx98 
            couma98   digox98   antia98 
            chrx98
            antid98 /* antidepressant */
            tag98 h2blk98 /* PPIs */
            insul98 ohypo98 /* Hypoglycemic */
            );

%nur00(keep=db00 hbp00 chol00 cabg00 str00 alone00 marry00 wt00
            thiaz00  lasix00  ccblo00  betab00  ace00 bprx00   
            couma00   digox00   antia00
            chrx00 oclrx00 /* Statin  other cholesterol-lowering */
            przc00 zol00 paxil00 celex00 antid00 /* antidepressant */
            tag00 h2blk00 prilo00 /* PPIs */
            insul00 ohypo00 /* Hypoglycemic */
            );

%nur02(keep=db02 hbp02 chol02 cabg02 str02 wt02
            thiaz02 lasix02 ace02 ccblo02 betab02 bprx02 
            przc02 zol02 paxil02 celex02 antid02 
            couma02   digox02   antia02 
            chrx02 oclrx02
            h2blk02 prilo02 /* PPIs */
            insul02 ohypo02
            cerbr02d);

%nur04(keep=db04 hbp04 chol04 cabg04 str04 alone04 marry04 wt04
            lasix04 thiaz04 ace04 ccblo04 betab04 k04 bprex04 
            mev04 zoc04 crest04 prav04 lip04 lesc04 oclrx04
            insul04 ohypo04 /* Hypoglycemic */
            ssri04 antid04  /* antidepressant */
            h2blk04 prilo04
            );

%nur06(keep=db06 hbp06 chol06 str06 cabg06 wt06 
            antid06 ssri06 
            bprx06 lasix06 thiaz06 ace06 ccblo06 betab06 k06
            mev06 zoc06 crest06 prav06 lip06 ostat06 oclrx06
            insul06 ohypo06 /* Hypoglycemic */
            ssri06 antid06  /* antidepressant */
            prilo06 h2blk06 cerbr06d
            );

%nur08(keep=db08 hbp08 chol08 cabg08 str08 alone08 marry08 wt08
            antid08 ssri08 
            bprx08 lasix08 thiaz08 ace08 ccblo08 betab08 k08 
            mev08 zoc08 crest08 prav08 lip08 ostat08 oclrx08
            insul08 metf08 ohypo08
            ssri08 antid08  /* antidepressant */
            prilo08 h2blk08
            );

%nur10(keep=db10 hbp10 chol10 str10 cabg10 wt10 
            antid10 ssri10 
            bprx10 lasix10 thiaz10 ace10 ccblo10 betab10 angio10 k10
            stat10 mev10 prav10 zoc10 lip10 crest10 ostat10 oclrx10
            insul10 metf10 avand10 ohypo10 
            prilo10 h2blk10
            cerbr10d);

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

*----------------3.Merge-----------------------------;

/**** Post-breast cancer ****/
    data base;
    merge der7620(in=mst)  br_marker_treatment(in=brc)  disease familyhis 
          n86_nts n90_nts n94_nts  n98_nts  n02_nts  n06_nts n10_nts
          n767880 nur82 nur84 nur86 nur88 nur90 nur92 nur94 nur96 nur98 nur00 nur02 nur04 nur06 nur08 nur10
          ahei2010_8410 meds8020 supp8016 act8614
          nses7616 
          end=_end_;
     by id;
     exrec=1; if first.id and mst then exrec=0; 
     brc_cases=1; if brc then brc_cases=0; 

      /* 4-year weight change */
          wtchg86c   = 0; /*baseline change = 0*/
          wtchg90c = wt90-wt86;  
          wtchg94c = wt94-wt90;
          wtchg98c = wt98-wt94;
          wtchg02c = wt02-wt98;  
          wtchg06c = wt06-wt02;  
          wtchg10c = wt10-wt06;  
          *wtchg14c = wt14-wt10;
          *wtchg18c = wt18-wt14;

      /* 4-year % weight change */
        wtchg86p = 0;
        if wt86 > 0 then wtchg90p = 100*(wt90-wt86)/wt86; else wtchg90p = .;
        if wt90 > 0 then wtchg94p = 100*(wt94-wt90)/wt90; else wtchg94p = .;
        if wt94 > 0 then wtchg98p = 100*(wt98-wt94)/wt94; else wtchg98p = .;
        if wt98 > 0 then wtchg02p = 100*(wt02-wt98)/wt98; else wtchg02p = .;
        if wt02 > 0 then wtchg06p = 100*(wt06-wt02)/wt02; else wtchg06p = .;
        if wt06 > 0 then wtchg10p = 100*(wt10-wt06)/wt06; else wtchg10p = .;

      /* 4-year BMI change */
          bmichg86c   = 0;
          bmichg90c = bmi90 - bmi86;
          bmichg94c = bmi94 - bmi90;
          bmichg98c = bmi98 - bmi94;
          bmichg02c = bmi02 - bmi98;
          bmichg06c = bmi06 - bmi02;
          bmichg10c = bmi10 - bmi06;
          *bmichg14c = bmi14 - bmi10;
          *bmichg18c = bmi18 - bmi14;
          
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


/* chronic disease status (flag for stopping diet updates) */
    *make intermediate outcomes cumulative varibles: i.e. once cancer, always cancer;
    array canc_prev{18} can76  can78  can80  can82  can84 can86  can88  can90  can92  can94 
                        can96  can98  can00  can02  can04 can06  can08  can10;
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
        if mi(i) =1 or strk(i)=1 or cabg(i)=1  or diab(i)=1  /*or canc_prev(i)=1*/  then chrodis(i)=1; else chrodis(i)=0; 
    end;

    do i = 2 to dim(chrodis); if chrodis{i-1} = 1 then chrodis{i} = 1; 
    end; /* Make chrodis cumulative */

/***************** Diet scores: carry forward + STOP updating after cancer/CVD *******************/
/* simple updated  scores */  
array irt    {*}  irt86 irt90 irt94 irt98 irt02 irt06 irt10 ;
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
array alcon {*}  alco86n  alco90n  alco94n  alco98n  alco02n  alco06n  alco10n;
array wtc4  {*}  wtchg86c wtchg90c wtchg94c wtchg98c wtchg02c wtchg06c wtchg10c;   /* 4-year change variables that do not carry forward missings */
array wtc4p {*}  wtchg86p wtchg90p wtchg94p wtchg98p wtchg02p wtchg06p wtchg10p;         
array bmichg4 {*}  bmichg86c bmichg90c bmichg94c bmichg98c bmichg02c bmichg06c bmichg10c;

array    carbp  {*}   carbo86p         carbo90p        carbo94p         carbo98p         carbo02p        carbo06p      carbo10p  ;
array    hcarbp  {*}   hcarbo86p        hcarbo90p       hcarbo94p        hcarbo98p        hcarbo02p       hcarbo06p     hcarbo10p   ;
array sumref_carbp {*}  sumref_carb86p   sumref_carb90p  sumref_carb94p   sumref_carb98p   sumref_carb02p  sumref_carb06p  sumref_carb10p ;
array     protp  {*}   prot86p          prot90p         prot94p          prot98p          prot02p         prot06p       prot10p ;
array    aprotp  {*}   aprot86p         aprot90p        aprot94p         aprot98p         aprot02p        aprot06p      aprot10p;
array    vprotp  {*}   vprot86p         vprot90p        vprot94p         vprot98p         vprot02p        vprot06p      vprot10p;
array    tfatp   {*}   tfat86p          tfat90p         tfat94p          tfat98p          tfat02p         tfat06p       tfat10p;
array    afatp   {*}   afat86p          afat90p         afat94p          afat98p          afat02p         afat06p       afat10p;
array    vfatp   {*}   vfat86p          vfat90p         vfat94p          vfat98p          vfat02p         vfat06p       vfat10p;

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
array alcons {*} alco86ns alco90ns alco94ns alco98ns alco02ns alco06ns alco10ns;
array wtc4s  {*}  wtchg86cs wtchg90cs wtchg94cs wtchg98cs wtchg02cs wtchg06cs wtchg10cs;   /* 4-year change variables that do not carry forward missings */
array wtc4ps {*}  wtchg86ps wtchg90ps wtchg94ps wtchg98ps wtchg02ps wtchg06ps wtchg10ps;         
array bmichg4s {*}  bmichg86cs bmichg90cs bmichg94cs bmichg98cs bmichg02cs bmichg06cs bmichg10cs;

array    carbs  {*}   carbo86s         carbo90s        carbo94s         carbo98s         carbo02s        carbo06s         carbo10s         ;
array    hcarbs  {*}   hcarbo86s        hcarbo90s       hcarbo94s        hcarbo98s        hcarbo02s       hcarbo06s        hcarbo10s        ;
array sumref_carbs {*}  sumref_carb86s   sumref_carb90s  sumref_carb94s   sumref_carb98s   sumref_carb02s  sumref_carb06s   sumref_carb10s   ;
array     prots  {*}   prot86s          prot90s         prot94s          prot98s          prot02s         prot06s          prot10s           ;
array    aprots  {*}   aprot86s         aprot90s        aprot94s         aprot98s         aprot02s        aprot06s         aprot10s          ;
array    vprots  {*}   vprot86s         vprot90s        vprot94s         vprot98s         vprot02s        vprot06s         vprot10s          ;
array    tfats   {*}   tfat86s          tfat90s         tfat94s          tfat98s          tfat02s         tfat06s          tfat10s           ;
array    afats   {*}   afat86s          afat90s         afat94s          afat98s          afat02s         afat06s          afat10s          ;
array    vfats   {*}   vfat86s          vfat90s         vfat94s          vfat98s          vfat02s         vfat06s          vfat10s          ;


/* nominal FFQ dates (June of each FFQ year, months since 1900), used when irt is missing */
array ffqdt {7} _temporary_ (1038 1086 1134 1182 1230 1278 1326);   /* 1986 1990 1994 1998 2002 2006 2010 */


/******** STEP 1. Post-diagnosis inputs: RAW values, post-diagnosis FFQs only ********/
/* Must run BEFORE the carry-forward in Step 2, so a missing post-dx FFQ is not      */
/* filled with a pre-diagnosis value.                                                */
do i = 1 to dim(calor);
    calors{i} = calor{i};   aheims{i} = aheim{i};
    lcss{i}   = lcs{i};     alcss{i}  = alcs{i};  vlcss{i} = vlcs{i};  ulcss{i} = ulcs{i};  hlcss{i} = hlcs{i};
    lfss{i}   = lfs{i};     alfss{i}  = alfs{i};  vlfss{i} = vlfs{i};  ulfss{i} = ulfs{i};  hlfss{i} = hlfs{i};
    alcons{i} = alcon{i};
    wtc4s{i}  = wtc4{i};    wtc4ps{i} = wtc4p{i}; bmichg4s{i} = bmichg4{i};   /* not restricted, as in the original */

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
array diet{7,13}
    calor86n ahei86 lcs86  alcs86 vlcs86  ulcs86  hlcs86 lfs86    ulfs86 hlfs86 alfs86 vlfs86  alco86n
    calor90n ahei90 lcs90  alcs90 vlcs90  ulcs90  hlcs90 lfs90    ulfs90 hlfs90 alfs90 vlfs90  alco90n
    calor94n ahei94 lcs94  alcs94 vlcs94  ulcs94  hlcs94 lfs94    ulfs94 hlfs94 alfs94 vlfs94  alco94n
    calor98n ahei98 lcs98  alcs98 vlcs98  ulcs98  hlcs98 lfs98    ulfs98 hlfs98 alfs98 vlfs98  alco98n
    calor02n ahei02 lcs02  alcs02 vlcs02  ulcs02  hlcs02 lfs02    ulfs02 hlfs02 alfs02 vlfs02  alco02n
    calor06n ahei06 lcs06  alcs06 vlcs06  ulcs06  hlcs06 lfs06    ulfs06 hlfs06 alfs06 vlfs06  alco06n
    calor10n ahei10 lcs10  alcs10 vlcs10  ulcs10  hlcs10 lfs10    ulfs10 hlfs10 alfs10 vlfs10  alco10n
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

array dietin{7,13}
    calor86ns ahei86s lcs86s  alcs86s vlcs86s  ulcs86s  hlcs86s lfs86s    ulfs86s hlfs86s alfs86s vlfs86s  alco86ns
    calor90ns ahei90s lcs90s  alcs90s vlcs90s  ulcs90s  hlcs90s lfs90s    ulfs90s hlfs90s alfs90s vlfs90s  alco90ns
    calor94ns ahei94s lcs94s  alcs94s vlcs94s  ulcs94s  hlcs94s lfs94s    ulfs94s hlfs94s alfs94s vlfs94s  alco94ns
    calor98ns ahei98s lcs98s  alcs98s vlcs98s  ulcs98s  hlcs98s lfs98s    ulfs98s hlfs98s alfs98s vlfs98s  alco98ns
    calor02ns ahei02s lcs02s  alcs02s vlcs02s  ulcs02s  hlcs02s lfs02s    ulfs02s hlfs02s alfs02s vlfs02s  alco02ns
    calor06ns ahei06s lcs06s  alcs06s vlcs06s  ulcs06s  hlcs06s lfs06s    ulfs06s hlfs06s alfs06s vlfs06s  alco06ns
    calor10ns ahei10s lcs10s  alcs10s vlcs10s  ulcs10s  hlcs10s lfs10s    ulfs10s hlfs10s alfs10s vlfs10s  alco10ns
;      
array dietout{7,13}
    calor86nsv  ahei86sv  lcs86sv  alcs86sv vlcs86sv  ulcs86sv  hlcs86sv  lfs86sv  ulfs86sv  hlfs86sv alfs86sv  vlfs86sv alco86sv
    calor90nsv  ahei90sv  lcs90sv  alcs90sv vlcs90sv  ulcs90sv  hlcs90sv  lfs90sv  ulfs90sv  hlfs90sv alfs90sv  vlfs90sv alco90sv
    calor94nsv  ahei94sv  lcs94sv  alcs94sv vlcs94sv  ulcs94sv  hlcs94sv  lfs94sv  ulfs94sv  hlfs94sv alfs94sv  vlfs94sv alco94sv
    calor98nsv  ahei98sv  lcs98sv  alcs98sv vlcs98sv  ulcs98sv  hlcs98sv  lfs98sv  ulfs98sv  hlfs98sv alfs98sv  vlfs98sv alco98sv
    calor02nsv  ahei02sv  lcs02sv  alcs02sv vlcs02sv  ulcs02sv  hlcs02sv  lfs02sv  ulfs02sv  hlfs02sv alfs02sv  vlfs02sv alco02sv
    calor06nsv  ahei06sv  lcs06sv  alcs06sv vlcs06sv  ulcs06sv  hlcs06sv  lfs06sv  ulfs06sv  hlfs06sv alfs06sv  vlfs06sv alco06sv
    calor10nsv  ahei10sv  lcs10sv  alcs10sv vlcs10sv  ulcs10sv  hlcs10sv  lfs10sv  ulfs10sv  hlfs10sv alfs10sv  vlfs10sv alco10sv
;


/* calculate cumulative average intake */
do j=1 to 13;
  dietout{1,j}=dietin{1,j};
 do i=2 to 7;
  sumvar=0;
   n=0;
     do k=1 to i;
      if (dietin{k,j} ne .) then do;
      n=n+1; sumvar=sumvar+dietin{k,j};
      end;
     end;
        if n = 0 then dietout{i,j} = .;          /* no post-dx FFQ yet -> missing */
        else          dietout{i,j} = sumvar / n;
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
        if n = 0 then macronv{i,j} = .;
        else          macronv{i,j} = sumvar / n;
    end;
end;
run;


PROC DATASETS;
delete    der7620 br_marker_treatment disease familyhis 
          n86_nts n90_nts n94_nts  n98_nts  n02_nts  n06_nts n10_nts
          n767880 nur82 nur84 nur86 nur88 nur90 nur92 nur94 nur96 nur98 nur00 nur02 nur04 nur06 nur08 nur10
          ahei2010_8410 meds8020 supp8016 act8614
          ;
RUN;


/*************************************************************************************/
/*                          Large section of arrays                                  */
/*************************************************************************************/
data nhs1;
     set base end=_end_;

cutoff=1470;    /*** Set cutoff at 2022.06 ***/
array   irt     {14}    irt86   irt88   irt90   irt92   irt94   irt96   irt98   irt00   irt02   irt04   irt06   irt08   irt10   cutoff;
array   tvar    {13}    t86     t88     t90     t92     t94     t96     t98     t00     t02     t04     t06     t08     t10 ;
array   wt      {13}    wt86    wt88    wt90    wt92    wt94    wt96    wt98    wt00    wt02    wt04    wt06    wt08    wt10    ;
array   bmi     {13}    bmi86   bmi88   bmi90   bmi92   bmi94   bmi96   bmi98   bmi00   bmi02   bmi04   bmi06   bmi08   bmi10   ;
array   age     {13}    age86   age88   age90   age92   age94   age96   age98   age00   age02   age04   age06   age08   age10   ;
array   hor     {13}    nhor86  nhor88  nhor90  nhor92  nhor94  nhor96  nhor98  nhor00  nhor02  nhor04  nhor06  nhor08  nhor10  ;
array   phms    {13}    phms86  phms88  phms90  phms92  phms94  phms96  phms98  phms00  phms02  phms04  phms06  phms08  phms10  ;
array   smok    {13}    smkdr86 smkdr88 smkdr90 smkdr92 smkdr94 smkdr96 smkdr98 smkdr00 smkdr02 smkdr04 smkdr06 smkdr08 smkdr10 ;
array   alcon   {13}    alco86n alco86n alco90n alco90n alco94n alco94n alco98n alco98n alco02n alco02n alco06n alco06n alco10n ;
array   alcq    {13}    alcc86  alcc88  alcc90  alcc92  alcc94  alcc96  alcc98  alcc00  alcc02  alcc04  alcc06  alcc08  alcc10  ;
array nSESArr   {13}    nSES_86m nSES_88m nSES_90m nSES_92m nSES_94m nSES_96m nSES_98m nSES_00m nSES_02m nSES_04m nSES_06m nSES_08m nSES_10m ;

array conf1119a {13};
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
array   strk    {13}    str86   str88   str90   str92   str94   str96   str98   str00   str02   str04   str06   str08   str10   ;
array   canc    {13}    can86   can88   can90   can92   can94   can96   can98   can00   can02   can04   can06   can08   can10;
array   wtc4    {13}    wtchg86c wtchg86c wtchg90c wtchg90c wtchg94c wtchg94c wtchg98c wtchg98c wtchg02c wtchg02c wtchg06c wtchg06c wtchg10c;
array   wtc4p   {13}    wtchg86p wtchg86p wtchg90p wtchg90p wtchg94p wtchg94p wtchg98p wtchg98p wtchg02p wtchg02p wtchg06p wtchg06p wtchg10p;         
array   bmichg4   {13}    bmichg86c bmichg86c bmichg90c bmichg90c bmichg94c bmichg94c bmichg98c bmichg98c bmichg02c bmichg02c bmichg06c bmichg06c bmichg10c;


**************************************************************************
***** If an irt date is before June of that qq year or after or equal ****
***** to the next qq year it is incorrect and should be defaulted to  ****
***** June of that qq year.    Make time period indicator tvar=0.       ****
**************************************************************************;

/* From 1986: create ltf using the original irt values */
ltf=.;
do i=1 to dim(irt)-1; if irt{i}>0 then ltf=irt{i}; end;
if ltf=irt{dim(irt)-1} then ltf=.;

/* Correct irt values */
do i=1 to dim(irt)-1;
    if (irt{i} <  (1014+24*i)) or (irt{i} >= (1038+24*i)) then irt{i}=1014+24*i;
end;

/* There are a few cases with date of diagnosis after date of death.  Since were more sure about date of death, use that date before we assign period of diagnosis */
    if 0<dtdxca<9999 and 0<dtdth<9999 then do;
        if dtdxca>dtdth then dtdxca=dtdth;
    end;




do j=2 to DIM(bmi); /*carry forward*/
      if age{j}=.      then age{j}=age{j-1}+2;
      *if bmi{j}=.      then bmi{j}=bmi{j-1};
      if hor{j}=.      then hor{j}=hor{j-1}; /* New PMH use, $label 1.Pre/Missing Mnp; 2.Never Used; 3.Current User; 4.Past User; 5.Status Unknown; 6.Missing */
      *if smok{j}=.     then smok{j}=smok{j-1}; /* smoking status range 0-15 0.missing 1.never 2-8.past 9-15.curr*/
      *if alcon{j}=.    then alcon{j}=alcon{j-1};
      if act{j}=.      then act{j}=act{j-1};
      if mvit{j}=.     then mvit{j}=mvit{j-1};
      if aspi{j}=.     then aspi{j}=aspi{j-1};
      *if marriage{j}=. then marriage{j}=marriage{j-1};
      *if alone{j}=.    then alone{j}=alone{j-1};
      if antihp{j}=.   then antihp{j}=antihp{j-1};
      if antitc{j}=.   then antitc{j}=antitc{j-1};
      if canc{j}=.     then canc{j}=canc{j-1};
      if hbp{j}=.       then hbp{j}=hbp{j-1};
      if chol{j}=.      then chol{j}=chol{j-1};
      if nSESArr{j}=.  then nSESArr{j}=nSESArr{j-1}; 
end;

do i=2 to DIM(bmi);
if irt{i-1}<=dtdxca<irt{i} then do;
  age_dx=age{i-1} + (dtdxca - irt{i-1})/12; /* age at diagnosis*/
  bmi_dx=bmi{i-1}; /* BMI prior to diagnosis*/
  *bmi_dx=bmi{i}; 
  highbp_dx=hbp{i};
  hightc_dx=chol{i};
  db_dx=diab{i};
end;end;
if age_dx = . and dtdxca >= irt{13} then do;
  age_dx=age{13};
  bmi_dx=bmi{13};
  highbp_dx=hbp{13};
  hightc_dx=chol{13};
  db_dx=diab{13};
end;

/* -single  s-stop updated  v-average  sv-stop updated average*/

array calor {13}  calor86n  calor86n calor90n  calor90n calor94n  calor94n calor98n  calor98n calor02n  calor02n calor06n calor06n calor10n;
array aheim {13}  ahei86 ahei86 ahei90 ahei90 ahei94 ahei94 ahei98 ahei98 ahei02 ahei02 ahei06 ahei06 ahei10 ;   
array lcs   {13}  lcs86    lcs86    lcs90    lcs90   lcs94    lcs94    lcs98    lcs98    lcs02    lcs02    lcs06    lcs06    lcs10;
array alcs  {13}  alcs86   alcs86   alcs90   alcs90  alcs94   alcs94   alcs98   alcs98   alcs02   alcs02   alcs06   alcs06   alcs10;
array vlcs  {13}  vlcs86   vlcs86   vlcs90   vlcs90  vlcs94   vlcs94   vlcs98   vlcs98   vlcs02   vlcs02   vlcs06   vlcs06   vlcs10;
array ulcs  {13}  ulcs86   ulcs86   ulcs90   ulcs90  ulcs94   ulcs94   ulcs98   ulcs98   ulcs02   ulcs02   ulcs06   ulcs06   ulcs10;
array hlcs  {13}  hlcs86   hlcs86   hlcs90   hlcs90  hlcs94   hlcs94   hlcs98   hlcs98   hlcs02   hlcs02   hlcs06   hlcs06   hlcs10;
array lfs   {13}  lfs86    lfs86    lfs90    lfs90   lfs94    lfs94    lfs98    lfs98    lfs02    lfs02    lfs06    lfs06    lfs10;
array ulfs  {13}  ulfs86   ulfs86   ulfs90   ulfs90  ulfs94   ulfs94   ulfs98   ulfs98   ulfs02   ulfs02   ulfs06   ulfs06   ulfs10;
array hlfs  {13}  hlfs86   hlfs86   hlfs90   hlfs90  hlfs94   hlfs94   hlfs98   hlfs98   hlfs02   hlfs02   hlfs06   hlfs06   hlfs10;
array alfs {13} alfs86 alfs86 alfs90 alfs90 alfs94 alfs94 alfs98 alfs98 alfs02 alfs02 alfs06 alfs06 alfs10 ;
array vlfs {13} vlfs86 vlfs86 vlfs90 vlfs90 vlfs94 vlfs94 vlfs98 vlfs98 vlfs02 vlfs02 vlfs06 vlfs06 vlfs10 ;


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
array alcons {13} alco86ns alco86ns alco90ns alco90ns alco94ns alco94ns alco98ns alco98ns alco02ns alco02ns alco06ns alco06ns alco10ns ;

array    carbs  {*}   carbo86s    carbo86s     carbo90s    carbo90s    carbo94s    carbo94s     carbo98s   carbo98s      carbo02s   carbo02s     carbo06s   carbo06s   carbo10s  ;
array    hcarbs  {*}   hcarbo86s   hcarbo86s    hcarbo90s   hcarbo90s   hcarbo94s   hcarbo94s    hcarbo98s  hcarbo98s     hcarbo02s  hcarbo02s    hcarbo06s  hcarbo06s  hcarbo10s   ;
array sumref_carbs {*}  sumref_carb86s   sumref_carb86s sumref_carb90s  sumref_carb90s sumref_carb94s  sumref_carb94s sumref_carb98s  sumref_carb98s sumref_carb02s  sumref_carb02s sumref_carb06s  sumref_carb06s sumref_carb10s ;
array     prots  {*}   prot86s     prot86s      prot90s     prot90s     prot94s     prot94s      prot98s     prot98s     prot02s     prot02s    prot06s    prot06s   prot10s ;
array    aprots  {*}   aprot86s     aprot86s      aprot90s     aprot90s     aprot94s     aprot94s      aprot98s     aprot98s     aprot02s     aprot02s    aprot06s    aprot06s   aprot10s ;
array    vprots  {*}   vprot86s     vprot86s      vprot90s     vprot90s     vprot94s     vprot94s      vprot98s     vprot98s     vprot02s     vprot02s    vprot06s    vprot06s   vprot10s ;
array    tfats   {*}   tfat86s     tfat86s      tfat90s     tfat90s     tfat94s     tfat94s      tfat98s     tfat98s     tfat02s     tfat02s    tfat06s    tfat06s   tfat10s ;
array    afats   {*}   afat86s     afat86s      afat90s     afat90s     afat94s     afat94s      afat98s     afat98s     afat02s     afat02s    afat06s    afat06s   afat10s ;
array    vfats   {*}   vfat86s     vfat86s      vfat90s     vfat90s     vfat94s     vfat94s      vfat98s     vfat98s     vfat02s     vfat02s    vfat06s    vfat06s   vfat10s ;

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
array alcosv {13} alco86sv alco86sv alco90sv alco90sv alco94sv alco94sv alco98sv alco98sv alco02sv alco02sv alco06sv alco06sv alco10sv ;

array carbsv {13}  carbo86sv carbo86sv carbo90sv carbo90sv carbo94sv carbo94sv carbo98sv carbo98sv carbo02sv carbo02sv carbo06sv carbo06sv carbo10sv  ;
array hcarbsv {13}  hcarbo86sv hcarbo86sv hcarbo90sv hcarbo90sv hcarbo94sv hcarbo94sv hcarbo98sv hcarbo98sv hcarbo02sv hcarbo02sv hcarbo06sv hcarbo06sv hcarbo10sv  ;   
array sumref_carbsv   {13}  sumref_carb86sv    sumref_carb86sv    sumref_carb90sv    sumref_carb90sv    sumref_carb94sv    sumref_carb94sv    sumref_carb98sv    sumref_carb98sv    sumref_carb02sv    sumref_carb02sv    sumref_carb06sv    sumref_carb06sv    sumref_carb10sv  ;
array protsv  {13}  prot86sv   prot86sv   prot90sv   prot90sv   prot94sv   prot94sv   prot98sv   prot98sv   prot02sv   prot02sv   prot06sv   prot06sv   prot10sv  ;
array aprotsv  {13}  aprot86sv   aprot86sv   aprot90sv   aprot90sv   aprot94sv   aprot94sv   aprot98sv   aprot98sv   aprot02sv   aprot02sv   aprot06sv   aprot06sv   aprot10sv  ;
array vprotsv  {13}  vprot86sv   vprot86sv   vprot90sv   vprot90sv   vprot94sv   vprot94sv   vprot98sv   vprot98sv   vprot02sv   vprot02sv   vprot06sv   vprot06sv   vprot10sv  ;
array tfatsv  {13}  tfat86sv   tfat86sv   tfat90sv   tfat90sv   tfat94sv   tfat94sv   tfat98sv   tfat98sv   tfat02sv   tfat02sv   tfat06sv   tfat06sv   tfat10sv  ;
array afatsv {13} afat86sv afat86sv afat90sv afat90sv afat94sv afat94sv afat98sv afat98sv afat02sv afat02sv afat06sv afat06sv afat10sv  ;
array vfatsv {13} vfat86sv vfat86sv vfat90sv vfat90sv vfat94sv vfat94sv vfat98sv vfat98sv vfat02sv vfat02sv vfat06sv vfat06sv vfat10sv  ;

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
array calor_pre {13}  calor86p_pre calor86p_pre calor90p_pre calor90p_pre calor94p_pre calor94p_pre calor98p_pre calor98p_pre calor02p_pre calor02p_pre calor06p_pre calor06p_pre calor10p_pre ;
array aheim_pre {13}  ahei86_pre ahei86_pre ahei90_pre ahei90_pre ahei94_pre ahei94_pre ahei98_pre ahei98_pre ahei02_pre ahei02_pre ahei06_pre ahei06_pre ahei10_pre ;   
array lcs_pre   {13}  lcs86_pre    lcs86_pre    lcs90_pre    lcs90_pre    lcs94_pre    lcs94_pre    lcs98_pre    lcs98_pre    lcs02_pre    lcs02_pre    lcs06_pre    lcs06_pre    lcs10_pre    ;
array vlcs_pre  {13}  vlcs86_pre   vlcs86_pre   vlcs90_pre   vlcs90_pre   vlcs94_pre   vlcs94_pre   vlcs98_pre   vlcs98_pre   vlcs02_pre   vlcs02_pre   vlcs06_pre   vlcs06_pre   vlcs10_pre   ;
array alcs_pre  {13}  alcs86_pre   alcs86_pre   alcs90_pre   alcs90_pre   alcs94_pre   alcs94_pre   alcs98_pre   alcs98_pre   alcs02_pre   alcs02_pre   alcs06_pre   alcs06_pre   alcs10_pre   ;
array ulcs_pre  {13}  ulcs86_pre   ulcs86_pre   ulcs90_pre   ulcs90_pre   ulcs94_pre   ulcs94_pre   ulcs98_pre   ulcs98_pre   ulcs02_pre   ulcs02_pre   ulcs06_pre   ulcs06_pre   ulcs10_pre   ; 
array hlcs_pre  {13}  hlcs86_pre   hlcs86_pre   hlcs90_pre   hlcs90_pre   hlcs94_pre   hlcs94_pre   hlcs98_pre   hlcs98_pre   hlcs02_pre   hlcs02_pre   hlcs06_pre   hlcs06_pre   hlcs10_pre   ;
array lfs_pre   {13}  lfs86_pre    lfs86_pre    lfs90_pre    lfs90_pre    lfs94_pre    lfs94_pre    lfs98_pre    lfs98_pre    lfs02_pre    lfs02_pre    lfs06_pre    lfs06_pre    lfs10_pre    ;
array vlfs_pre  {13}  vlfs86_pre   vlfs86_pre   vlfs90_pre   vlfs90_pre   vlfs94_pre   vlfs94_pre   vlfs98_pre   vlfs98_pre   vlfs02_pre   vlfs02_pre   vlfs06_pre   vlfs06_pre   vlfs10_pre   ;
array alfs_pre  {13}  alfs86_pre   alfs86_pre   alfs90_pre   alfs90_pre   alfs94_pre   alfs94_pre   alfs98_pre   alfs98_pre   alfs02_pre   alfs02_pre   alfs06_pre   alfs06_pre   alfs10_pre   ;
array ulfs_pre  {13}  ulfs86_pre   ulfs86_pre   ulfs90_pre   ulfs90_pre   ulfs94_pre   ulfs94_pre   ulfs98_pre   ulfs98_pre   ulfs02_pre   ulfs02_pre   ulfs06_pre   ulfs06_pre   ulfs10_pre   ; 
array hlfs_pre  {13}  hlfs86_pre   hlfs86_pre   hlfs90_pre   hlfs90_pre   hlfs94_pre   hlfs94_pre   hlfs98_pre   hlfs98_pre   hlfs02_pre   hlfs02_pre   hlfs06_pre   hlfs06_pre   hlfs10_pre   ;
array alco_pre  {13}  alco86_pre   alco86_pre   alco90_pre   alco90_pre   alco94_pre   alco94_pre   alco98_pre   alco98_pre   alco02_pre   alco02_pre   alco06_pre   alco06_pre   alco10_pre   ;

/*Post-dx, donnot use the carry forward value to ensure the cumulaitve averaged validity*/
array calor_post {13}  calor86p_post calor86p_post calor90p_post calor90p_post calor94p_post calor94p_post calor98p_post calor98p_post calor02p_post calor02p_post calor06p_post calor06p_post calor10p_post ;
array aheim_post {13}  ahei86_post ahei86_post ahei90_post ahei90_post ahei94_post ahei94_post ahei98_post ahei98_post ahei02_post ahei02_post ahei06_post ahei06_post ahei10_post ;   
array lcs_post   {13}  lcs86_post    lcs86_post    lcs90_post    lcs90_post    lcs94_post    lcs94_post    lcs98_post    lcs98_post    lcs02_post    lcs02_post    lcs06_post    lcs06_post    lcs10_post    ;
array vlcs_post  {13}  vlcs86_post   vlcs86_post   vlcs90_post   vlcs90_post   vlcs94_post   vlcs94_post   vlcs98_post   vlcs98_post   vlcs02_post   vlcs02_post   vlcs06_post   vlcs06_post   vlcs10_post   ;
array alcs_post  {13}  alcs86_post   alcs86_post   alcs90_post   alcs90_post   alcs94_post   alcs94_post   alcs98_post   alcs98_post   alcs02_post   alcs02_post   alcs06_post   alcs06_post   alcs10_post   ;
array ulcs_post  {13}  ulcs86_post   ulcs86_post   ulcs90_post   ulcs90_post   ulcs94_post   ulcs94_post   ulcs98_post   ulcs98_post   ulcs02_post   ulcs02_post   ulcs06_post   ulcs06_post   ulcs10_post   ; 
array hlcs_post  {13}  hlcs86_post   hlcs86_post   hlcs90_post   hlcs90_post   hlcs94_post   hlcs94_post   hlcs98_post   hlcs98_post   hlcs02_post   hlcs02_post   hlcs06_post   hlcs06_post   hlcs10_post   ;
array lfs_post   {13}  lfs86_post    lfs86_post    lfs90_post    lfs90_post    lfs94_post    lfs94_post    lfs98_post    lfs98_post    lfs02_post    lfs02_post    lfs06_post    lfs06_post    lfs10_post    ;
array vlfs_post  {13}  vlfs86_post   vlfs86_post   vlfs90_post   vlfs90_post   vlfs94_post   vlfs94_post   vlfs98_post   vlfs98_post   vlfs02_post   vlfs02_post   vlfs06_post   vlfs06_post   vlfs10_post   ;
array alfs_post  {13}  alfs86_post   alfs86_post   alfs90_post   alfs90_post   alfs94_post   alfs94_post   alfs98_post   alfs98_post   alfs02_post   alfs02_post   alfs06_post   alfs06_post   alfs10_post   ;
array ulfs_post  {13}  ulfs86_post   ulfs86_post   ulfs90_post   ulfs90_post   ulfs94_post   ulfs94_post   ulfs98_post   ulfs98_post   ulfs02_post   ulfs02_post   ulfs06_post   ulfs06_post   ulfs10_post   ; 
array hlfs_post  {13}  hlfs86_post   hlfs86_post   hlfs90_post   hlfs90_post   hlfs94_post   hlfs94_post   hlfs98_post   hlfs98_post   hlfs02_post   hlfs02_post   hlfs06_post   hlfs06_post   hlfs10_post   ;
array alco_post  {13}  alco86_post   alco86_post   alco90_post   alco90_post   alco94_post   alco94_post   alco98_post   alco98_post   alco02_post   alco02_post   alco06_post   alco06_post   alco10_post   ;

/* Compute pre, post, and chg */
do i=2 to DIM(calor)+1;
    if irt{i-1}<dtdxca<=irt{i} then do;
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

/* FFQ positions in the 17-element period arrays: 86 90 94 98 02 06 10   */
array ffqk {7} _temporary_ (1 3 5 7 9 11 13  );

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
array calor_crc {13}  calor86p_crc calor86p_crc calor90p_crc calor90p_crc calor94p_crc calor94p_crc calor98p_crc calor98p_crc calor02p_crc calor02p_crc calor06p_crc calor06p_crc calor10p_crc ;
array aheim_crc {13}  ahei86_crc ahei86_crc ahei90_crc ahei90_crc ahei94_crc ahei94_crc ahei98_crc ahei98_crc ahei02_crc ahei02_crc ahei06_crc ahei06_crc ahei10_crc ;   
array lcs_crc   {13}  lcs86_crc    lcs86_crc    lcs90_crc    lcs90_crc    lcs94_crc    lcs94_crc    lcs98_crc    lcs98_crc    lcs02_crc    lcs02_crc    lcs06_crc    lcs06_crc    lcs10_crc    ;
array vlcs_crc  {13}  vlcs86_crc   vlcs86_crc   vlcs90_crc   vlcs90_crc   vlcs94_crc   vlcs94_crc   vlcs98_crc   vlcs98_crc   vlcs02_crc   vlcs02_crc   vlcs06_crc   vlcs06_crc   vlcs10_crc   ;
array alcs_crc  {13}  alcs86_crc   alcs86_crc   alcs90_crc   alcs90_crc   alcs94_crc   alcs94_crc   alcs98_crc   alcs98_crc   alcs02_crc   alcs02_crc   alcs06_crc   alcs06_crc   alcs10_crc   ;
array ulcs_crc  {13}  ulcs86_crc   ulcs86_crc   ulcs90_crc   ulcs90_crc   ulcs94_crc   ulcs94_crc   ulcs98_crc   ulcs98_crc   ulcs02_crc   ulcs02_crc   ulcs06_crc   ulcs06_crc   ulcs10_crc   ; 
array hlcs_crc  {13}  hlcs86_crc   hlcs86_crc   hlcs90_crc   hlcs90_crc   hlcs94_crc   hlcs94_crc   hlcs98_crc   hlcs98_crc   hlcs02_crc   hlcs02_crc   hlcs06_crc   hlcs06_crc   hlcs10_crc   ;
array lfs_crc   {13}  lfs86_crc    lfs86_crc    lfs90_crc    lfs90_crc    lfs94_crc    lfs94_crc    lfs98_crc    lfs98_crc    lfs02_crc    lfs02_crc    lfs06_crc    lfs06_crc    lfs10_crc    ;
array vlfs_crc  {13}  vlfs86_crc   vlfs86_crc   vlfs90_crc   vlfs90_crc   vlfs94_crc   vlfs94_crc   vlfs98_crc   vlfs98_crc   vlfs02_crc   vlfs02_crc   vlfs06_crc   vlfs06_crc   vlfs10_crc   ;
array alfs_crc  {13}  alfs86_crc   alfs86_crc   alfs90_crc   alfs90_crc   alfs94_crc   alfs94_crc   alfs98_crc   alfs98_crc   alfs02_crc   alfs02_crc   alfs06_crc   alfs06_crc   alfs10_crc   ;
array ulfs_crc  {13}  ulfs86_crc   ulfs86_crc   ulfs90_crc   ulfs90_crc   ulfs94_crc   ulfs94_crc   ulfs98_crc   ulfs98_crc   ulfs02_crc   ulfs02_crc   ulfs06_crc   ulfs06_crc   ulfs10_crc   ; 
array hlfs_crc  {13}  hlfs86_crc   hlfs86_crc   hlfs90_crc   hlfs90_crc   hlfs94_crc   hlfs94_crc   hlfs98_crc   hlfs98_crc   hlfs02_crc   hlfs02_crc   hlfs06_crc   hlfs06_crc   hlfs10_crc   ;
array alco_crc  {13}  alco86_crc   alco86_crc   alco90_crc   alco90_crc   alco94_crc   alco94_crc   alco98_crc   alco98_crc   alco02_crc   alco02_crc   alco06_crc   alco06_crc   alco10_crc   ;

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
        if dtdxca > 0 then do f = 1 to dim(ffqk);  k = ffqk{f};  /* -- FFQ positions */
           nomd = 1014 + 24*k;        /* -- NHS ; NHSII 1074+24*k ; HPFS 1009+24*k */
           if nomd > dtdxca and nomd < cutoff and (dtdth in (.,0) or dtdth > nomd)
              then elig_any = 1;     /* -- alive when an FFQ was mailed after dx */
        end;



*%beginex(idfile=exclude.nhs1.out,timevar=i);
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
    if irt{i} < dtdth <= irt{i+1} and conf1119=1 then do; dead = 1; tdead = dtdth - irt{i}; end;
    /* Cause-specific death */
    dead_can   = 0; tdead_can=irt{i+1}-irt{i}; dtdx_dead_can=.;
    if irt{i}< dtdth <=irt{i+1} and (140<=newicda3<=208) and conf1119=1 then do; dead_can=1; tdead_can=dtdth-irt{i}; dtdx_dead_can=dtdth; end;

    dead_brc   = 0; tdead_brc=irt{i+1}-irt{i}; dtdx_dead_brc=.; /*breast cancer death*/
    if irt{i}< dtdth <=irt{i+1} and (newicda3=174) and conf1119=1 then do; dead_brc=1; tdead_brc=dtdth-irt{i}; dtdx_dead_brc=dtdth; end;

    dead_crc   = 0; tdead_crc=irt{i+1}-irt{i}; dtdx_dead_crc=.; /*colorectal cancer death*/
    if irt{i}< dtdth <=irt{i+1} and newicda3 in (153,154) and conf1119=1 then do; dead_crc=1; tdead_crc=dtdth-irt{i}; dtdx_dead_crc=dtdth; end;

    dead_nocan=0; tdead_nocan=irt{i+1}-irt{i}; dtdx_dead_nocan=.;
    if irt{i}<dtdth<=irt{i+1} and dead_can ne 1 then do; dead_nocan=1; dtdx_dead_nocan=dtdth; end;

    dead_cvd   = 0; tdead_cvd=irt{i+1}-irt{i}; dtdx_dead_cvd=.;
    if irt{i}< dtdth <=irt{i+1} and (390<=newicda3<=459) and conf1119=1 then do; dead_cvd=1; tdead_cvd=dtdth-irt{i}; dtdx_dead_cvd=dtdth; end;

    dead_resp=0;  tdead_resp=irt{i+1}-irt{i}; dtdx_dead_resp=.;
    if irt{i}<dtdth<=irt{i+1} and (460<=newicda3<520) and conf1119=1  then do; dead_resp=1; tdead_resp=dtdth-irt{i}; dtdx_dead_resp=dtdth; end;

    dead_oth=0;tdead_oth=irt{i+1}-irt{i}; dtdx_dead_oth=.;
    if irt{i}<dtdth<=irt{i+1} and (dead_can=0 and dead_cvd=0 and dead_resp=0 ) and conf1119=1 then do; dead_oth=1; tdead_oth=dtdth-irt{i}; dtdx_dead_oth=dtdth; end;

    /* Calendar-time person-time for reference */
    ptime=irt{i+1}-irt{i}; 
    if irt{i}<dtdth<=irt{i+1}  and conf1119=1  then  ptime=dtdth-irt{i};

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
   
/****** breast cancer duration ******/
        bcduration = (irt{i} - dtdxca)/12;
        if bcduration < 1 then bcduration = 1;
      *Define bc duration category; 
      if 0<bcduration<=5 then bcduration_c=1;
      else if  5<bcduration<=10 then bcduration_c=2;
      else if 10<bcduration<=15 then bcduration_c=3;
      else if 15<bcduration then bcduration_c=4;
 %indic3(vbl=bcduration_c, prefix=bcduration_c, min=2, max=4, reflev=1, usemiss=0,
    label1='0-5 years', label2='5-10 years',label3='10-15 yaers',label4='>15 years');

 /* Year of diagnosis */
    if dtdxca=9999 then yeardx_breast=.;
    else yeardx_breast=round(((dtdxca/12)+1900-0.5)); /* Calendar year of diagnosis */

/***************time since diagnosis (post-dx survival time scale) ******************/
            
      if dtdxca > . and conf1119 = 1 then do;
            if irt{i+1} <= dtdxca then do; /* pre-dx time: no post-dx contribution */
                tsd = .; time_begin = .; time_end = .; person_time = 0;
            end;
            else do;
                if irt{i} le dtdxca lt irt{i+1} then do; /* Interval contains the diagnosis date */
                    time_begin=0;
                end;          
                else do;                                 /* Intervals after diagnosis */
                time_begin = irt{i} - dtdxca; end;
                tsd = time_begin;

                end_cal = min(irt{i+1}, cutoff);         /*calendar end of the interval*/
                if 0 < dtdth <= end_cal then end_cal = dtdth;
                time_end = end_cal - dtdxca;

                person_time = max(0, time_end - time_begin);
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

    /* 3. -svchg: change from pre-dx: cumulative post - single pre*/
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
%indic3(vbl=eth3g, prefix=eth3g, reflev=1, min=2, max=3, missing=., usemiss=0,
           label1='causcasian',
           label2='african',
           label3='asian & others');

/* AGE GROUP */
     agecon=age{i};
     age_bl = age{1};             /* baseline age */
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
array   smkka   {13}    smkk86   smkk88   smkk90   smkk92   smkk94    smkk96   smkk98   smkk00   smkk02   smkk04   smkk06   smkk08  smkk10 ;
array   smksta  {13}    smksta86   smksta88   smksta90   smksta92   smksta94  smksta96  smksta98  smksta00   smksta02   smksta04   smksta06   smksta08  smksta10 ;

length smkk 3 smkstatus 3;
        if      smok{i} = 1                      then smkk = 1;  * never;
        else if smok{i} in (2,3,4,5,6,7,8)       then smkk = 2;  * past smoker;
        else if smok{i} in (9,10)                then smkk = 3;  * current, 1-14 cigs;
        else if smok{i} = 11                     then smkk = 4;  * current, 15-24 cigs;
        else if smok{i} in (12,13,14)            then smkk = 5;  * current, 25+ cigs;
        else                                          smkk = .;  * missing;


        /* store time-varying 5-level category */
     smkka{i} = smkk;
     %indic3(vbl=smkk, prefix=smkk, min=2,max=5, reflev=1,  usemiss=1,
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
      %indic3(vbl=smkstatus, prefix=smkstatus, min=2, max=3, reflev=1,  usemiss=1);

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

/*** Indicator Alcohol consumption ***/
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

        if alc in (2,3,4,5,6) then neverdrinker=0; else neverdrinker=1;

        /* post dx alcohol intake, cumulative from i+1*/
        alcosvcon = alcosv{i};
              if      alcosv{i}=0       then alcovc_post=1;
            else if 0 <  alcosv{i}<5    then alcovc_post=2;
            else if 5<=  alcosv{i}<15   then alcovc_post=3;
            else if 15<=  alcosv{i}<30  then alcovc_post=4;
            else if 30<= alcosv{i}      then alcovc_post=5;
            else                             alcovc_post=9;                    
          
          label alcovc_post='Post-dx Cumulative average daily alcohol intake, g/d, categorical';
          %indic3(vbl=alcovc_post, prefix=alcovc_post, min=2, max=5, reflev=1, missing=9, usemiss=1,
              label1='none',
              label2='0.1-4.9 g/d',
              label3='5.0-14.9 g/d',
              label4='15.0-29.9 g/d',
              label5='30+ g/d');


/*** Indicator for Weight Change ***/
    wtcont     =  wt{i};        *baseline weight in lbs (continuous);
    wtcontkg   =  wt{i}/2.2;    *baseline weight in kg (continuous);
    wtc4cont   =  wtc4{i};      *4-yr change in weight in lbs (continuous);
    wtc4contkg =  wtc4{i}/2.2;  *4-yr change in weight in kg (continuous);
    wtc4contp   =  wtc4p{i};    *4-yr change in % weight(continuous);


/*** Indicator for BMI ***/
         bmicon=bmi{i};         *baseline bmi (continuous);
         bmichg4con =  bmichg4{i};  *4-yr change in bmi (continuous);

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

/***family history of MI***/
       if mifh ne 1 then mifh=0; 
       
/****** Indicator for History of antihypertensive user *******/
            select(antihp{i});
                when (1)     antihbp=1;
                otherwise    antihbp=0;
            end;

            if hbp{i}=1  then highbp=1; else highbp=0;
            if hbp{1}=1 then basehighbp=1; else basehighbp=0;

            if highbp_dx=1 then diaghbp=1; else diaghbp=0;
         

/****** Indicator for TC lowering *******/
            select(antitc{i});
                when (1)     antihtc=1;
                otherwise    antihtc=0;
            end;

           if chol{i}=1 then hightc=1; else hightc=0;
           if chol{1}=1 then basehightc=1; else basehightc=0;

           if hightc_dx=1 then diaghtc=1; else diaghtc=0;

/****** Indicator for History of diabetes *******/
           if diab{i}=1 then diabetes=1; else diabetes=0;

/** Indicator for nSES **/
nSES = nSESArr{i};

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

if 0 <= actconp and actconp < 7  then preact=1;
      else if actconp>=7 and actconp<15 then preact=2;
      else if actconp>=15 and actconp<25 then preact=3;
      else if actconp>=25  then preact=4;
      else  preact=2;

%indic3(vbl=preact, prefix=preact, min=2, max=4, reflev=1, missing=., usemiss=0,
        label1='<7 total mets/wk',
        label2='7 to < 15 mets/wk',
        label3='15 to < 25 mets/wk',
        label4='25+ mets/wk');


%indic3(vbl=smkstp, prefix=smkstp, min=2, max=3, reflev=1, missing=., usemiss=1);
if aspp = . then aspp = 0;
if mvitp = . then mvitp = 0;

est_status2 = est_status + 1; if est_status2 = . then est_status2 = 3;
%indic3(vbl=est_status2,       prefix=est_status2,    min=1, max=3, reflev=1, usemiss=0, label1='ER-', label2='ER+', label3='Missing'); 

/****************  BASELINE EXCLUSIONS ********************************************/

if i=1 then do;
   %exclude(exrec eq 1);                 *multiple records and not in master file;
   %exclude(brc_cases eq 1);             *not in brc file;
   %exclude(conf1119 ne 1);

   %exclude(0 lt dtdth   le irt{i});
   %exclude(0 lt dtdxca  le irt{i}); 
   %exclude(dtdxca ge cutoff); *cases diagnosed after the end of follow-up;

   %exclude(ltf eq irt{1});
   %exclude(yobf le 20);
   %exclude(age86 eq .);
   %exclude(stage eq 4);
   %exclude(stage eq 5);
   %exclude(canpre eq 0, skip=t);  * skip individuals diagnosed after current time period;
   %output();
  end;

/* EXCLUSIONS DURING FOLLOW-UP */
   if i>1 then do;
   %exclude(irt{i-1} lt dtdth le irt{i} );
   %exclude(stage eq 4);
   %exclude(stage eq 5);
   %exclude(canpre eq 0, skip=t);  * skip individuals diagnosed after current time period;
   *%exclude(lcssvcon eq ., skip=t);
   %output();
  end;
end;        

/* END OF DO-LOOP OVER TIME PERIODs */
%endex();
