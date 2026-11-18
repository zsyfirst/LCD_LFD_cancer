/*
# Updated: 10/2026
# Purpose of this program: preparing data for sensitivity analysis: 4-year low carb change and weight change
# content: Table S18
----------------------------------------------------------

Exposures:  4-year changes in low carb score (5 types) per quintile/SD increase:

1. lcs - Low carb diet (per quintile: LCSCQUIN; per SD: LCSCSD)
2. alcs - animal based low carb diet (per quintile: ALCSCQUIN; per SD: ALCSCSD)
3. vlcs - plant based low carb diet (per quintile: VLCSCQUIN; per SD: VLCSCSD)
4. ulcs - low-quality low carb diet (per quintile: ULCSCQUIN; per SD: ULCSCSD)
5. hlcs - high-quality low carb diet (per quintile: HLCSCQUIN; per SD: HLCSCSD)
----------------------------------------------------------

Outcome: 4-year weight change in kilogram (WTC4CONTKG)
--------------------
Follow-up period: NHS: 1984-2010, NHSII: 1991-2019, HPFS: 1986-2018
HPFS FFQ cycle: 1986, 1990, 1994, 1998, 2002, 2006, 2010, 2014, 2018
NHS FFQ cycle:  1986, 1990, 1994, 1998, 2002, 2006, 2010
NHS2 FFQ cycle: 1991, 1995, 1999, 2003, 2007, 2011, 2015, 2019
----------------------------------------------------------

Covariates: (all time-varying unless described as baseline covariates)
  AGECONT (age, continuous), 
  RACE (race: white, African American, Asian, others), 
  FHDB (family history of diabetes, yes/no), 
  HBPBASE (baseline hypertension, yes/no), 
  CHOLBASE (baseline hypercholesterolemia, yes/no), 
  CALORB (baseline total caloric intake, continuous, kcal/day), 
  BMIB (baseline body mass index, <25.0, 25.0-29.9, or ≥30.0 kg/m2), 
  SMKC4CAT (change in smoking status: stayed never smoker, stayed former smoker, stayed current smoker, change from former to current smoker, change from never to current smoker, and change from current to former smoker), 
  ACTB (baseline in physical activity, METs-hour/week), 
  ACTSC4CONT (change in physical activity, METs-hour/week),
  ALCOC4 (change in alcohol consumption, continuous, g/day), 
  PPMH (postmenopausal hormone use, women only, premenopausal, never, former, current, or missing), 
  POCU (oral contraceptive use, NHSII only, yes/no),
  STUDY (data origin, NHS, NHSII or HPFS)

*/

%include 'lowcarb_cancer/DietChg_WtChg/lcs_data_pooled.sas';

data all;
  set all end=_end_;

%indic3(vbl=   lcscquin ,prefix=  lcscquin,  min=0, max=4, reflev=2, missing=., usemiss=0);  
%indic3(vbl=  alcscquin ,prefix=  alcscquin, min=0, max=4, reflev=2, missing=., usemiss=0);  
%indic3(vbl=  vlcscquin ,prefix=  vlcscquin, min=0, max=4, reflev=2, missing=., usemiss=0);  
%indic3(vbl=  hlcscquin ,prefix=  hlcscquin, min=0, max=4, reflev=2, missing=., usemiss=0);  
%indic3(vbl=  ulcscquin ,prefix=  ulcscquin, min=0, max=4, reflev=2, missing=., usemiss=0);  


%indic3(vbl=race,prefix=race,min=2,max=4,reflev=1, missing=., usemiss=0,
                label1='caucasian',
                label2='asian',
                label3='hispanic',
                label4='african american');
%indic3(vbl=smkc4cat, prefix=smkc4cat, min=2, max=6, reflev=1, missing=., usemiss=0,
          label1='never-never',
          label2='current-past',
          label3='past-current',
          label4='never-current',
          label5='past-past',
          label6='current-current');
%indic3(vbl=bmib, prefix=bmib, min=2, max=3, reflev=1, missing=., usemiss=0,
            label1='bmi  <25.0',
            label2='bmi =>25.0 & <30.0',
            label3='bmi =>30.0 ');
%indic3(vbl=age55, prefix=age55, min=1, max=1, reflev=0, missing=., usemiss=0,
          label0='<55 yrs',
          label1='>=55 yrs'); 
%indic3(vbl=study, prefix=study, min=1, max=2, reflev=0, missing=., usemiss=0,
            label0='NHS',
            label1='NHS2',
            label2='HPFS');


 %let ppmh_=ppmh1 ppmh2 ppmh3 ppmh4;
 %let pocu_=pocu1 pocu2 pocu3;

run;

proc means data=all n nmiss mean std median min max;
run;  



/* descriptive stats of score quintiles */

proc means data=all n nmiss mean std median min max;
class lcscquin ;
  var lcscont
      lcscontc 
      lcscsd 
      baselcscont 
      wtc4contkg;
run;  

proc means data=all n nmiss mean std median min max;
class alcscquin ;
  var alcscont
      alcscontc 
      alcscsd 
      basealcscont 
      wtc4contkg;
run;  

proc means data=all n nmiss mean std median min max;
class vlcscquin ;
  var vlcscont
      vlcscontc 
      vlcscsd 
      basevlcscont 
      wtc4contkg;
run;  
proc means data=all n nmiss mean std median min max;
class hlcscquin ;
  var hlcscont
      hlcscontc 
      hlcscsd 
      basehlcscont 
      wtc4contkg;
run;  
proc means data=all n nmiss mean std median min max;
class ulcscquin ;
  var ulcscont
      ulcscontc 
      ulcscsd 
      baseulcscont 
      wtc4contkg;
run;


proc means data=all n nmiss mean std median min max;
class lfscquin ;
  var lfscont
      lfscontc 
      lfscsd 
      baselfscont 
      wtc4contkg;
run;  

proc means data=all n nmiss mean std median min max;
class alfscquin ;
  var alfscont
      alfscontc 
      alfscsd 
      basealfscont 
      wtc4contkg;
run;  

proc means data=all n nmiss mean std median min max;
class vlfscquin ;
  var vlfscont
      vlfscontc 
      vlfscsd 
      basevlfscont 
      wtc4contkg;
run;  
proc means data=all n nmiss mean std median min max;
class hlfscquin ;
  var hlfscont
      hlfscontc 
      hlfscsd 
      basehlfscont 
      wtc4contkg;
run;  
proc means data=all n nmiss mean std median min max;
class ulfscquin ;
  var ulfscont
      ulfscontc 
      ulfscsd 
      baseulfscont 
      wtc4contkg;
run;    
proc means data=all mean median p5 p10 q1 q3 p90 min max n nmiss; class lcscquin  ;var lcscontc ;run;
proc means data=all mean median p5 p10 q1 q3 p90 min max n nmiss; class alcscquin  ;var alcscontc ;run;
proc means data=all mean median p5 p10 q1 q3 p90 min max n nmiss; class vlcscquin  ;var vlcscontc ;run;
proc means data=all mean median p5 p10 q1 q3 p90 min max n nmiss; class hlcscquin  ;var hlcscontc ;run;
proc means data=all mean median p5 p10 q1 q3 p90 min max n nmiss; class ulcscquin  ;var ulcscontc ;run;

proc means data=all mean median p5 p10 q1 q3 p90 min max n nmiss; class lfscquin  ;var lfscontc ;run;
proc means data=all mean median p5 p10 q1 q3 p90 min max n nmiss; class alfscquin  ;var alfscontc ;run;
proc means data=all mean median p5 p10 q1 q3 p90 min max n nmiss; class vlfscquin  ;var vlfscontc ;run;
proc means data=all mean median p5 p10 q1 q3 p90 min max n nmiss; class hlfscquin  ;var hlfscontc ;run;
proc means data=all mean median p5 p10 q1 q3 p90 min max n nmiss; class ulfscquin  ;var ulfscontc ;run;

/**************/
/* Table 2  */

/* Age-adjusted  */

/* per 1-SD increase */
proc genmod data = all;
class id;
model wtc4contkg =  lcscsd agecont;
repeated subject=id/corr= un modelse;
run; 
proc genmod data = all;
class id;
model wtc4contkg =  alcscsd agecont;
repeated subject=id/corr= un modelse;
run; 
proc genmod data = all;
class id;
model wtc4contkg =  vlcscsd agecont;
repeated subject=id/corr= un modelse;
run; 
proc genmod data = all;
class id;
model wtc4contkg =  hlcscsd agecont;
repeated subject=id/corr= un modelse;
run; 
proc genmod data = all;
class id;
model wtc4contkg =  ulcscsd agecont;
repeated subject=id/corr= un modelse;
run; 

proc genmod data = all;
class id;
model wtc4contkg =  lfscsd agecont;
repeated subject=id/corr= un modelse;
run; 
proc genmod data = all;
class id;
model wtc4contkg =  alfscsd agecont;
repeated subject=id/corr= un modelse;
run; 
proc genmod data = all;
class id;
model wtc4contkg =  vlfscsd agecont;
repeated subject=id/corr= un modelse;
run; 
proc genmod data = all;
class id;
model wtc4contkg =  hlfscsd agecont;
repeated subject=id/corr= un modelse;
run; 
proc genmod data = all;
class id;
model wtc4contkg =  ulfscsd agecont;
repeated subject=id/corr= un modelse;
run; 

/* per quintile increase */
proc genmod data = all;
class id;
model wtc4contkg =  &lcscquin_ agecont;
repeated subject=id/corr= un modelse;
run; 
proc genmod data = all;
class id;
model wtc4contkg =  &alcscquin_  agecont;
repeated subject=id/corr= un modelse;
run; 
proc genmod data = all;
class id;
model wtc4contkg =  &vlcscquin_  agecont;
repeated subject=id/corr= un modelse;
run; 
proc genmod data = all;
class id;
model wtc4contkg =  &hlcscquin_  agecont;
repeated subject=id/corr= un modelse;
run; 
proc genmod data = all;
class id;
model wtc4contkg =  &ulcscquin_  agecont;
repeated subject=id/corr= un modelse;
run; 

proc genmod data = all;
class id;
model wtc4contkg =  &lfscquin_ agecont;
repeated subject=id/corr= un modelse;
run; 
proc genmod data = all;
class id;
model wtc4contkg =  &alfscquin_  agecont;
repeated subject=id/corr= un modelse;
run; 
proc genmod data = all;
class id;
model wtc4contkg =  &vlfscquin_  agecont;
repeated subject=id/corr= un modelse;
run; 
proc genmod data = all;
class id;
model wtc4contkg =  &hlfscquin_  agecont;
repeated subject=id/corr= un modelse;
run; 
proc genmod data = all;
class id;
model wtc4contkg =  &ulfscquin_  agecont;
repeated subject=id/corr= un modelse;
run; 

/* p for trend */
proc genmod data = all;
class id;
model wtc4contkg =  lcscquintt agecont;
repeated subject=id/corr= un modelse;
run; 
proc genmod data = all;
class id;
model wtc4contkg =  alcscquintt  agecont;
repeated subject=id/corr= un modelse;
run; 
proc genmod data = all;
class id;
model wtc4contkg =  vlcscquintt  agecont;
repeated subject=id/corr= un modelse;
run; 
proc genmod data = all;
class id;
model wtc4contkg =  hlcscquintt agecont;
repeated subject=id/corr= un modelse;
run; 
proc genmod data = all;
class id;
model wtc4contkg =  ulcscquintt  agecont;
repeated subject=id/corr= un modelse;
run; 
proc genmod data = all;
class id;
model wtc4contkg =  lfscquintt agecont;
repeated subject=id/corr= un modelse;
run; 
proc genmod data = all;
class id;
model wtc4contkg =  alfscquintt  agecont;
repeated subject=id/corr= un modelse;
run; 
proc genmod data = all;
class id;
model wtc4contkg =  vlfscquintt  agecont;
repeated subject=id/corr= un modelse;
run; 
proc genmod data = all;
class id;
model wtc4contkg =  hlfscquintt agecont;
repeated subject=id/corr= un modelse;
run; 
proc genmod data = all;
class id;
model wtc4contkg =  ulfscquintt  agecont;
repeated subject=id/corr= un modelse;
run; 


/** Multivariable model **/
* adjusted for age, race, family history of T2D, baseline HBP, hypercholesterolemia, 
baseline BMI, change in smoking status, baseline  total caloric intake,
baseline and change in physical activity, change in alcohol consumption, postmenopausal hormone use, and oral contraceptive use
data origin ;

/* per 1-SD increase */
proc genmod data = all;
class id;
model wtc4contkg =  lcscsd  
          hbpbase cholbase fhdb calorb    
          agecont &bmib_ actsc4cont actb  alcoc4  &race_ 
          &smkc4cat_  &ppmh_ &pocu_ study;
repeated subject=id/corr= un modelse;
run; 

proc genmod data = all;
class id;
model wtc4contkg =  alcscsd  
          hbpbase cholbase fhdb calorb    
          agecont &bmib_ actsc4cont actb  alcoc4  &race_ 
          &smkc4cat_  &ppmh_ &pocu_ study;
repeated subject=id/corr= un modelse;
run; 

proc genmod data = all;
class id;
model wtc4contkg =  vlcscsd  
          hbpbase cholbase fhdb calorb    
          agecont &bmib_ actsc4cont actb  alcoc4  &race_ 
          &smkc4cat_  &ppmh_ &pocu_ study;
repeated subject=id/corr= un modelse;
run; 

proc genmod data = all;
class id;
model wtc4contkg =  hlcscsd  
          hbpbase cholbase fhdb calorb    
          agecont &bmib_ actsc4cont actb  alcoc4  &race_ 
          &smkc4cat_  &ppmh_ &pocu_ study;
repeated subject=id/corr= un modelse;
run; 

proc genmod data = all;
class id;
model wtc4contkg =  ulcscsd  
          hbpbase cholbase fhdb calorb    
          agecont &bmib_ actsc4cont actb  alcoc4  &race_ 
          &smkc4cat_  &ppmh_ &pocu_ study;
repeated subject=id/corr= un modelse;
run; 

proc genmod data = all;
class id;
model wtc4contkg =  lfscsd  
          hbpbase cholbase fhdb calorb    
          agecont &bmib_ actsc4cont actb  alcoc4  &race_ 
          &smkc4cat_  &ppmh_ &pocu_ study;
repeated subject=id/corr= un modelse;
run; 

proc genmod data = all;
class id;
model wtc4contkg =  alfscsd  
          hbpbase cholbase fhdb calorb    
          agecont &bmib_ actsc4cont actb  alcoc4  &race_ 
          &smkc4cat_  &ppmh_ &pocu_ study;
repeated subject=id/corr= un modelse;
run; 

proc genmod data = all;
class id;
model wtc4contkg =  vlfscsd  
          hbpbase cholbase fhdb calorb    
          agecont &bmib_ actsc4cont actb  alcoc4  &race_ 
          &smkc4cat_  &ppmh_ &pocu_ study;
repeated subject=id/corr= un modelse;
run; 

proc genmod data = all;
class id;
model wtc4contkg =  hlfscsd  
          hbpbase cholbase fhdb calorb    
          agecont &bmib_ actsc4cont actb  alcoc4  &race_ 
          &smkc4cat_  &ppmh_ &pocu_ study;
repeated subject=id/corr= un modelse;
run; 

proc genmod data = all;
class id;
model wtc4contkg =  ulfscsd  
          hbpbase cholbase fhdb calorb    
          agecont &bmib_ actsc4cont actb  alcoc4  &race_ 
          &smkc4cat_  &ppmh_ &pocu_ study;
repeated subject=id/corr= un modelse;
run; 


/* per quintile increase */
proc genmod data = all;
class id;
model wtc4contkg =  &lcscquin_  
          hbpbase cholbase fhdb calorb    
          agecont &bmib_ actsc4cont actb  alcoc4  &race_ 
          &smkc4cat_  &ppmh_ &pocu_ study;
repeated subject=id/corr= un modelse;
run; 

proc genmod data = all;
class id;
model wtc4contkg =  &alcscquin_  
          hbpbase cholbase fhdb calorb    
          agecont &bmib_ actsc4cont actb  alcoc4  &race_ 
          &smkc4cat_  &ppmh_ &pocu_ study;
repeated subject=id/corr= un modelse;
run; 

proc genmod data = all;
class id;
model wtc4contkg =  &vlcscquin_  
          hbpbase cholbase fhdb calorb    
          agecont &bmib_ actsc4cont actb  alcoc4  &race_ 
          &smkc4cat_  &ppmh_ &pocu_ study;
repeated subject=id/corr= un modelse;
run; 

proc genmod data = all;
class id;
model wtc4contkg =  &hlcscquin_  
          hbpbase cholbase fhdb calorb    
          agecont &bmib_ actsc4cont actb  alcoc4  &race_ 
          &smkc4cat_  &ppmh_ &pocu_ study;
repeated subject=id/corr= un modelse;
run; 

proc genmod data = all;
class id;
model wtc4contkg =  &ulcscquin_  
          hbpbase cholbase fhdb calorb    
          agecont &bmib_ actsc4cont actb  alcoc4  &race_ 
          &smkc4cat_  &ppmh_ &pocu_ study;
repeated subject=id/corr= un modelse;
run; 

proc genmod data = all;
class id;
model wtc4contkg =  &lfscquin_  
          hbpbase cholbase fhdb calorb    
          agecont &bmib_ actsc4cont actb  alcoc4  &race_ 
          &smkc4cat_  &ppmh_ &pocu_ study;
repeated subject=id/corr= un modelse;
run; 

proc genmod data = all;
class id;
model wtc4contkg =  &alfscquin_  
          hbpbase cholbase fhdb calorb    
          agecont &bmib_ actsc4cont actb  alcoc4  &race_ 
          &smkc4cat_  &ppmh_ &pocu_ study;
repeated subject=id/corr= un modelse;
run; 

proc genmod data = all;
class id;
model wtc4contkg =  &vlfscquin_  
          hbpbase cholbase fhdb calorb    
          agecont &bmib_ actsc4cont actb  alcoc4  &race_ 
          &smkc4cat_  &ppmh_ &pocu_ study;
repeated subject=id/corr= un modelse;
run; 

proc genmod data = all;
class id;
model wtc4contkg =  &hlfscquin_  
          hbpbase cholbase fhdb calorb    
          agecont &bmib_ actsc4cont actb  alcoc4  &race_ 
          &smkc4cat_  &ppmh_ &pocu_ study;
repeated subject=id/corr= un modelse;
run; 

proc genmod data = all;
class id;
model wtc4contkg =  &ulfscquin_  
          hbpbase cholbase fhdb calorb    
          agecont &bmib_ actsc4cont actb  alcoc4  &race_ 
          &smkc4cat_  &ppmh_ &pocu_ study;
repeated subject=id/corr= un modelse;
run; 


/* p for trend */
proc genmod data = all;
class id;
model wtc4contkg =  lcscquintt
          hbpbase cholbase fhdb calorb    
          agecont &bmib_ actsc4cont actb  alcoc4  &race_ 
          &smkc4cat_  &ppmh_ &pocu_ study;
repeated subject=id/corr= un modelse;
run; 

proc genmod data = all;
class id;
model wtc4contkg =  alcscquintt
          hbpbase cholbase fhdb calorb    
          agecont &bmib_ actsc4cont actb  alcoc4  &race_ 
          &smkc4cat_  &ppmh_ &pocu_ study;
repeated subject=id/corr= un modelse;
run; 

proc genmod data = all;
class id;
model wtc4contkg =  vlcscquintt
          hbpbase cholbase fhdb calorb    
          agecont &bmib_ actsc4cont actb  alcoc4  &race_ 
          &smkc4cat_  &ppmh_ &pocu_ study;
repeated subject=id/corr= un modelse;
run; 

proc genmod data = all;
class id;
model wtc4contkg = hlcscquintt
          hbpbase cholbase fhdb calorb    
          agecont &bmib_ actsc4cont actb  alcoc4  &race_ 
          &smkc4cat_  &ppmh_ &pocu_ study;
repeated subject=id/corr= un modelse;
run; 

proc genmod data = all;
class id;
model wtc4contkg = ulcscquintt
          hbpbase cholbase fhdb calorb    
          agecont &bmib_ actsc4cont actb  alcoc4  &race_ 
          &smkc4cat_  &ppmh_ &pocu_ study;
repeated subject=id/corr= un modelse;
run; 

proc genmod data = all;
class id;
model wtc4contkg =  lfscquintt
          hbpbase cholbase fhdb calorb    
          agecont &bmib_ actsc4cont actb  alcoc4  &race_ 
          &smkc4cat_  &ppmh_ &pocu_ study;
repeated subject=id/corr= un modelse;
run; 

proc genmod data = all;
class id;
model wtc4contkg =  alfscquintt
          hbpbase cholbase fhdb calorb    
          agecont &bmib_ actsc4cont actb  alcoc4  &race_ 
          &smkc4cat_  &ppmh_ &pocu_ study;
repeated subject=id/corr= un modelse;
run; 

proc genmod data = all;
class id;
model wtc4contkg =  vlfscquintt
          hbpbase cholbase fhdb calorb    
          agecont &bmib_ actsc4cont actb  alcoc4  &race_ 
          &smkc4cat_  &ppmh_ &pocu_ study;
repeated subject=id/corr= un modelse;
run; 

proc genmod data = all;
class id;
model wtc4contkg = hlfscquintt
          hbpbase cholbase fhdb calorb    
          agecont &bmib_ actsc4cont actb  alcoc4  &race_ 
          &smkc4cat_  &ppmh_ &pocu_ study;
repeated subject=id/corr= un modelse;
run; 

proc genmod data = all;
class id;
model wtc4contkg = ulfscquintt
          hbpbase cholbase fhdb calorb    
          agecont &bmib_ actsc4cont actb  alcoc4  &race_ 
          &smkc4cat_  &ppmh_ &pocu_ study;
repeated subject=id/corr= un modelse;
run; 
