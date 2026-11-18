/*****************************************************
Purpose:       low carb, low fat, and risk of Obesity-related cancer
Study design:  prospective cohort (NHS: 1986-2022.6, NHSII: 1991-2023.6, HPFS: 1986-2024.1)

Contents:      Table S20 --- PY and HR (Obesity-related cancer as outcome)
Exposure:      Using Continuously Updated Diet Score (No Stopping of Dietary Updates After Diagnosis of CVD, Diabetes, or Cancer)
Date:          10/2026
******************************************************/


/************************ Call in Cleaned Datasets**********************************/
%include 'lowcarb_cancer/incidence/without_stop/data.sas';


/************************ CASE-PERSON YEARS*****************************************/
%macro casepy(cohort,exp,case);
proc sql;
create table temp1 as
select *, 
sum(_ppm) as sumppm,sum(&case) as &cohort._sum_&case
from two
group by &exp;
;
quit;
data &cohort._&exp._&case;
set temp1;
where &exp ne .;
by &exp;
if first.&exp;
sumpy=round(sumppm/12,1);
cohort="&cohort";
casepy_&case._&cohort=(put(&cohort._sum_&case,comma20.) || '/' || put(sumpy,comma20. -L));
casepy_&case=&cohort._sum_&case ;
keep cohort &exp casepy_&case._&cohort casepy_&case sumpy;
run;
proc transpose data=&cohort._&exp._&case out=&cohort._&exp._&case.print prefix=&exp;
by cohort;
var casepy_&case._&cohort casepy_&case sumpy;
run;
title"**************************&cohort. for &case. Case/PY**************************";
proc print data=&cohort._&exp._&case.print;
run;
title"";
data &cohort._&exp._&case;
length index $6.;
set &cohort._&exp._&case;
index="&exp";
run;
%mend;

/*==== NHS1 ====*/
%pre_pm(data=nhs1,out=TWO,case=can_obese,timevar=period,
irt=irt86 irt88 irt90 irt92 irt94 irt96 irt98 irt00 irt02 irt04 irt06 irt08 irt10,cutoff=cutoff,
dtdx=dtdx_obeca, dtdth=dtdth, var=agegp lcssvq); 
%pm(data=TWO, case=can_obese, exposure=lcssvq, strata=agegp);
%casepy(cohort=nhs1, exp=lcssvq, case=can_obese);

%pre_pm(data=nhs1,out=TWO,case=can_obese,timevar=period,
irt=irt86 irt88 irt90 irt92 irt94 irt96 irt98 irt00 irt02 irt04 irt06 irt08 irt10,cutoff=cutoff,
dtdx=dtdx_obeca, dtdth=dtdth, var=agegp alcssvq); 
%pm(data=TWO, case=can_obese, exposure=alcssvq, strata=agegp);
%casepy(cohort=nhs1, exp=alcssvq, case=can_obese);

%pre_pm(data=nhs1,out=TWO,case=can_obese,timevar=period,
irt=irt86 irt88 irt90 irt92 irt94 irt96 irt98 irt00 irt02 irt04 irt06 irt08 irt10,cutoff=cutoff,
dtdx=dtdx_obeca, dtdth=dtdth, var=agegp vlcssvq); 
%pm(data=TWO, case=can_obese, exposure=vlcssvq, strata=agegp);
%casepy(cohort=nhs1, exp=vlcssvq, case=can_obese);

%pre_pm(data=nhs1,out=TWO,case=can_obese,timevar=period,
irt=irt86 irt88 irt90 irt92 irt94 irt96 irt98 irt00 irt02 irt04 irt06 irt08 irt10,cutoff=cutoff,
dtdx=dtdx_obeca, dtdth=dtdth, var=agegp ulcssvq); 
%pm(data=TWO, case=can_obese, exposure=ulcssvq, strata=agegp);
%casepy(cohort=nhs1, exp=ulcssvq, case=can_obese);

%pre_pm(data=nhs1,out=TWO,case=can_obese,timevar=period,
irt=irt86 irt88 irt90 irt92 irt94 irt96 irt98 irt00 irt02 irt04 irt06 irt08 irt10,cutoff=cutoff,
dtdx=dtdx_obeca, dtdth=dtdth, var=agegp hlcssvq); 
%pm(data=TWO, case=can_obese, exposure=hlcssvq, strata=agegp);
%casepy(cohort=nhs1, exp=hlcssvq, case=can_obese);

%pre_pm(data=nhs1,out=TWO,case=can_obese,timevar=period,
irt=irt86 irt88 irt90 irt92 irt94 irt96 irt98 irt00 irt02 irt04 irt06 irt08 irt10,cutoff=cutoff,
dtdx=dtdx_obeca, dtdth=dtdth, var=agegp lfssvq); 
%pm(data=TWO, case=can_obese, exposure=lfssvq, strata=agegp);
%casepy(cohort=nhs1, exp=lfssvq, case=can_obese);

%pre_pm(data=nhs1,out=TWO,case=can_obese,timevar=period,
irt=irt86 irt88 irt90 irt92 irt94 irt96 irt98 irt00 irt02 irt04 irt06 irt08 irt10,cutoff=cutoff,
dtdx=dtdx_obeca, dtdth=dtdth, var=agegp alfssvq); 
%pm(data=TWO, case=can_obese, exposure=alfssvq, strata=agegp);
%casepy(cohort=nhs1, exp=alfssvq, case=can_obese);

%pre_pm(data=nhs1,out=TWO,case=can_obese,timevar=period,
irt=irt86 irt88 irt90 irt92 irt94 irt96 irt98 irt00 irt02 irt04 irt06 irt08 irt10,cutoff=cutoff,
dtdx=dtdx_obeca, dtdth=dtdth, var=agegp vlfssvq); 
%pm(data=TWO, case=can_obese, exposure=vlfssvq, strata=agegp);
%casepy(cohort=nhs1, exp=vlfssvq, case=can_obese);

%pre_pm(data=nhs1,out=TWO,case=can_obese,timevar=period,
irt=irt86 irt88 irt90 irt92 irt94 irt96 irt98 irt00 irt02 irt04 irt06 irt08 irt10,cutoff=cutoff,
dtdx=dtdx_obeca, dtdth=dtdth, var=agegp ulfssvq); 
%pm(data=TWO, case=can_obese, exposure=ulfssvq, strata=agegp);
%casepy(cohort=nhs1, exp=ulfssvq, case=can_obese);

%pre_pm(data=nhs1,out=TWO,case=can_obese,timevar=period,
irt=irt86 irt88 irt90 irt92 irt94 irt96 irt98 irt00 irt02 irt04 irt06 irt08 irt10,cutoff=cutoff,
dtdx=dtdx_obeca, dtdth=dtdth, var=agegp hlfssvq); 
%pm(data=TWO, case=can_obese, exposure=hlfssvq, strata=agegp);
%casepy(cohort=nhs1, exp=hlfssvq, case=can_obese);



/*==== NHS2 ====*/
%pre_pm(data=nhs2,out=TWO,case=can_obese,timevar=period,
irt=retmo91 retmo93 retmo95 retmo97 retmo99 retmo01 retmo03 retmo05 retmo07 retmo09 retmo11 retmo13 retmo15 retmo17 retmo19,
cutoff=cutoff, dtdx=dtdx_obeca, dtdth=dtdth, var=agegp lcssvq); 
%pm(data=TWO, case=can_obese, exposure=lcssvq, strata=agegp);
%casepy(cohort=nhs2, exp=lcssvq, case=can_obese);

%pre_pm(data=nhs2,out=TWO,case=can_obese,timevar=period,
irt=retmo91 retmo93 retmo95 retmo97 retmo99 retmo01 retmo03 retmo05 retmo07 retmo09 retmo11 retmo13 retmo15 retmo17 retmo19,
cutoff=cutoff, dtdx=dtdx_obeca, dtdth=dtdth, var=agegp alcssvq); 
%pm(data=TWO, case=can_obese, exposure=alcssvq, strata=agegp);
%casepy(cohort=nhs2, exp=alcssvq, case=can_obese);

%pre_pm(data=nhs2,out=TWO,case=can_obese,timevar=period,
irt=retmo91 retmo93 retmo95 retmo97 retmo99 retmo01 retmo03 retmo05 retmo07 retmo09 retmo11 retmo13 retmo15 retmo17 retmo19,
cutoff=cutoff, dtdx=dtdx_obeca, dtdth=dtdth, var=agegp vlcssvq); 
%pm(data=TWO, case=can_obese, exposure=vlcssvq, strata=agegp);
%casepy(cohort=nhs2, exp=vlcssvq, case=can_obese);

%pre_pm(data=nhs2,out=TWO,case=can_obese,timevar=period,
irt=retmo91 retmo93 retmo95 retmo97 retmo99 retmo01 retmo03 retmo05 retmo07 retmo09 retmo11 retmo13 retmo15 retmo17 retmo19,
cutoff=cutoff, dtdx=dtdx_obeca, dtdth=dtdth, var=agegp ulcssvq); 
%pm(data=TWO, case=can_obese, exposure=ulcssvq, strata=agegp);
%casepy(cohort=nhs2, exp=ulcssvq, case=can_obese);

%pre_pm(data=nhs2,out=TWO,case=can_obese,timevar=period,
irt=retmo91 retmo93 retmo95 retmo97 retmo99 retmo01 retmo03 retmo05 retmo07 retmo09 retmo11 retmo13 retmo15 retmo17 retmo19,
cutoff=cutoff, dtdx=dtdx_obeca, dtdth=dtdth, var=agegp hlcssvq); 
%pm(data=TWO, case=can_obese, exposure=hlcssvq, strata=agegp);
%casepy(cohort=nhs2, exp=hlcssvq, case=can_obese);

%pre_pm(data=nhs2,out=TWO,case=can_obese,timevar=period,
irt=retmo91 retmo93 retmo95 retmo97 retmo99 retmo01 retmo03 retmo05 retmo07 retmo09 retmo11 retmo13 retmo15 retmo17 retmo19,
cutoff=cutoff, dtdx=dtdx_obeca, dtdth=dtdth, var=agegp lfssvq); 
%pm(data=TWO, case=can_obese, exposure=lfssvq, strata=agegp);
%casepy(cohort=nhs2, exp=lfssvq, case=can_obese);

%pre_pm(data=nhs2,out=TWO,case=can_obese,timevar=period,
irt=retmo91 retmo93 retmo95 retmo97 retmo99 retmo01 retmo03 retmo05 retmo07 retmo09 retmo11 retmo13 retmo15 retmo17 retmo19,
cutoff=cutoff, dtdx=dtdx_obeca, dtdth=dtdth, var=agegp alfssvq); 
%pm(data=TWO, case=can_obese, exposure=alfssvq, strata=agegp);
%casepy(cohort=nhs2, exp=alfssvq, case=can_obese);

%pre_pm(data=nhs2,out=TWO,case=can_obese,timevar=period,
irt=retmo91 retmo93 retmo95 retmo97 retmo99 retmo01 retmo03 retmo05 retmo07 retmo09 retmo11 retmo13 retmo15 retmo17 retmo19,
cutoff=cutoff, dtdx=dtdx_obeca, dtdth=dtdth, var=agegp vlfssvq); 
%pm(data=TWO, case=can_obese, exposure=vlfssvq, strata=agegp);
%casepy(cohort=nhs2, exp=vlfssvq, case=can_obese);

%pre_pm(data=nhs2,out=TWO,case=can_obese,timevar=period,
irt=retmo91 retmo93 retmo95 retmo97 retmo99 retmo01 retmo03 retmo05 retmo07 retmo09 retmo11 retmo13 retmo15 retmo17 retmo19,
cutoff=cutoff, dtdx=dtdx_obeca, dtdth=dtdth, var=agegp ulfssvq); 
%pm(data=TWO, case=can_obese, exposure=ulfssvq, strata=agegp);
%casepy(cohort=nhs2, exp=ulfssvq, case=can_obese);

%pre_pm(data=nhs2,out=TWO,case=can_obese,timevar=period,
irt=retmo91 retmo93 retmo95 retmo97 retmo99 retmo01 retmo03 retmo05 retmo07 retmo09 retmo11 retmo13 retmo15 retmo17 retmo19,
cutoff=cutoff, dtdx=dtdx_obeca, dtdth=dtdth, var=agegp hlfssvq); 
%pm(data=TWO, case=can_obese, exposure=hlfssvq, strata=agegp);
%casepy(cohort=nhs2, exp=hlfssvq, case=can_obese);

/*==== HPFS ====*/
%pre_pm(data=hpfs,out=TWO,case=can_obese,timevar=period,
irt=rtmnyr86 rtmnyr88 rtmnyr90 rtmnyr92 rtmnyr94 rtmnyr96 rtmnyr98 rtmnyr00 rtmnyr02 rtmnyr04 rtmnyr06 rtmnyr08 rtmnyr10 rtmnyr12 rtmnyr14 rtmnyr16 rtmnyr18,cutoff=cutoff,
dtdx=dtdx_obeca, dtdth=dtdth, var=agegp lcssvq); 
%pm(data=TWO, case=can_obese, exposure=lcssvq, strata=agegp);
%casepy(cohort=hpfs, exp=lcssvq, case=can_obese);

%pre_pm(data=hpfs,out=TWO,case=can_obese,timevar=period,
irt=rtmnyr86 rtmnyr88 rtmnyr90 rtmnyr92 rtmnyr94 rtmnyr96 rtmnyr98 rtmnyr00 rtmnyr02 rtmnyr04 rtmnyr06 rtmnyr08 rtmnyr10 rtmnyr12 rtmnyr14 rtmnyr16 rtmnyr18,cutoff=cutoff,
dtdx=dtdx_obeca, dtdth=dtdth, var=agegp alcssvq); 
%pm(data=TWO, case=can_obese, exposure=alcssvq, strata=agegp);
%casepy(cohort=hpfs, exp=alcssvq, case=can_obese);

%pre_pm(data=hpfs,out=TWO,case=can_obese,timevar=period,
irt=rtmnyr86 rtmnyr88 rtmnyr90 rtmnyr92 rtmnyr94 rtmnyr96 rtmnyr98 rtmnyr00 rtmnyr02 rtmnyr04 rtmnyr06 rtmnyr08 rtmnyr10 rtmnyr12 rtmnyr14 rtmnyr16 rtmnyr18,cutoff=cutoff,
dtdx=dtdx_obeca, dtdth=dtdth, var=agegp vlcssvq); 
%pm(data=TWO, case=can_obese, exposure=vlcssvq, strata=agegp);
%casepy(cohort=hpfs, exp=vlcssvq, case=can_obese);

%pre_pm(data=hpfs,out=TWO,case=can_obese,timevar=period,
irt=rtmnyr86 rtmnyr88 rtmnyr90 rtmnyr92 rtmnyr94 rtmnyr96 rtmnyr98 rtmnyr00 rtmnyr02 rtmnyr04 rtmnyr06 rtmnyr08 rtmnyr10 rtmnyr12 rtmnyr14 rtmnyr16 rtmnyr18,cutoff=cutoff,
dtdx=dtdx_obeca, dtdth=dtdth, var=agegp ulcssvq); 
%pm(data=TWO, case=can_obese, exposure=ulcssvq, strata=agegp);
%casepy(cohort=hpfs, exp=ulcssvq, case=can_obese);

%pre_pm(data=hpfs,out=TWO,case=can_obese,timevar=period,
irt=rtmnyr86 rtmnyr88 rtmnyr90 rtmnyr92 rtmnyr94 rtmnyr96 rtmnyr98 rtmnyr00 rtmnyr02 rtmnyr04 rtmnyr06 rtmnyr08 rtmnyr10 rtmnyr12 rtmnyr14 rtmnyr16 rtmnyr18,cutoff=cutoff,
dtdx=dtdx_obeca, dtdth=dtdth, var=agegp hlcssvq); 
%pm(data=TWO, case=can_obese, exposure=hlcssvq, strata=agegp);
%casepy(cohort=hpfs, exp=hlcssvq, case=can_obese);

%pre_pm(data=hpfs,out=TWO,case=can_obese,timevar=period,
irt=rtmnyr86 rtmnyr88 rtmnyr90 rtmnyr92 rtmnyr94 rtmnyr96 rtmnyr98 rtmnyr00 rtmnyr02 rtmnyr04 rtmnyr06 rtmnyr08 rtmnyr10 rtmnyr12 rtmnyr14 rtmnyr16 rtmnyr18,cutoff=cutoff,
dtdx=dtdx_obeca, dtdth=dtdth, var=agegp lfssvq); 
%pm(data=TWO, case=can_obese, exposure=lfssvq, strata=agegp);
%casepy(cohort=hpfs, exp=lfssvq, case=can_obese);

%pre_pm(data=hpfs,out=TWO,case=can_obese,timevar=period,
irt=rtmnyr86 rtmnyr88 rtmnyr90 rtmnyr92 rtmnyr94 rtmnyr96 rtmnyr98 rtmnyr00 rtmnyr02 rtmnyr04 rtmnyr06 rtmnyr08 rtmnyr10 rtmnyr12 rtmnyr14 rtmnyr16 rtmnyr18,cutoff=cutoff,
dtdx=dtdx_obeca, dtdth=dtdth, var=agegp alfssvq); 
%pm(data=TWO, case=can_obese, exposure=alfssvq, strata=agegp);
%casepy(cohort=hpfs, exp=alfssvq, case=can_obese);

%pre_pm(data=hpfs,out=TWO,case=can_obese,timevar=period,
irt=rtmnyr86 rtmnyr88 rtmnyr90 rtmnyr92 rtmnyr94 rtmnyr96 rtmnyr98 rtmnyr00 rtmnyr02 rtmnyr04 rtmnyr06 rtmnyr08 rtmnyr10 rtmnyr12 rtmnyr14 rtmnyr16 rtmnyr18,cutoff=cutoff,
dtdx=dtdx_obeca, dtdth=dtdth, var=agegp vlfssvq); 
%pm(data=TWO, case=can_obese, exposure=vlfssvq, strata=agegp);
%casepy(cohort=hpfs, exp=vlfssvq, case=can_obese);

%pre_pm(data=hpfs,out=TWO,case=can_obese,timevar=period,
irt=rtmnyr86 rtmnyr88 rtmnyr90 rtmnyr92 rtmnyr94 rtmnyr96 rtmnyr98 rtmnyr00 rtmnyr02 rtmnyr04 rtmnyr06 rtmnyr08 rtmnyr10 rtmnyr12 rtmnyr14 rtmnyr16 rtmnyr18,cutoff=cutoff,
dtdx=dtdx_obeca, dtdth=dtdth, var=agegp ulfssvq); 
%pm(data=TWO, case=can_obese, exposure=ulfssvq, strata=agegp);
%casepy(cohort=hpfs, exp=ulfssvq, case=can_obese);

%pre_pm(data=hpfs,out=TWO,case=can_obese,timevar=period,
irt=rtmnyr86 rtmnyr88 rtmnyr90 rtmnyr92 rtmnyr94 rtmnyr96 rtmnyr98 rtmnyr00 rtmnyr02 rtmnyr04 rtmnyr06 rtmnyr08 rtmnyr10 rtmnyr12 rtmnyr14 rtmnyr16 rtmnyr18,cutoff=cutoff,
dtdx=dtdx_obeca, dtdth=dtdth, var=agegp hlfssvq); 
%pm(data=TWO, case=can_obese, exposure=hlfssvq, strata=agegp);
%casepy(cohort=hpfs, exp=hlfssvq, case=can_obese);

/* Print */
data printcasepy;
    set 
    NHS1_LCSSVQ_can_obesePRINT
    NHS1_ALCSSVQ_can_obesePRINT
    NHS1_VLCSSVQ_can_obesePRINT
    NHS1_ULCSSVQ_can_obesePRINT
    NHS1_HLCSSVQ_can_obesePRINT
    NHS1_LFSSVQ_can_obesePRINT
    NHS1_ALFSSVQ_can_obesePRINT
    NHS1_VLFSSVQ_can_obesePRINT
    NHS1_ULFSSVQ_can_obesePRINT
    NHS1_HLFSSVQ_can_obesePRINT

    NHS2_LCSSVQ_can_obesePRINT
    NHS2_ALCSSVQ_can_obesePRINT
    NHS2_VLCSSVQ_can_obesePRINT
    NHS2_ULCSSVQ_can_obesePRINT
    NHS2_HLCSSVQ_can_obesePRINT
    NHS2_LFSSVQ_can_obesePRINT
    NHS2_ALFSSVQ_can_obesePRINT
    NHS2_VLFSSVQ_can_obesePRINT
    NHS2_ULFSSVQ_can_obesePRINT
    NHS2_HLFSSVQ_can_obesePRINT

    HPFS_LCSSVQ_can_obesePRINT
    HPFS_ALCSSVQ_can_obesePRINT
    HPFS_VLCSSVQ_can_obesePRINT
    HPFS_ULCSSVQ_can_obesePRINT
    HPFS_HLCSSVQ_can_obesePRINT
    HPFS_LFSSVQ_can_obesePRINT
    HPFS_ALFSSVQ_can_obesePRINT
    HPFS_VLFSSVQ_can_obesePRINT
    HPFS_ULFSSVQ_can_obesePRINT
    HPFS_HLFSSVQ_can_obesePRINT
    ;
run;

proc export data=printcasepy
    dbms=csv
    outfile='lowcarb_cancer/incidence/without_stop/table/casePY_bycohort_can_obese.csv'
    replace;
    run;


/*==== pool ====*/
%macro cpy(exp,case);
data casepy;
    set 
    hpfs_&exp._&case.
    nhs1_&exp._&case.
    nhs2_&exp._&case. 
    ;
run;

proc sql;
create table temp1 as
select *, 
sum(casepy_&case.) as sumcase_all,sum(sumpy) as sumpy_all 
from casepy
group by &exp;
;
quit;

proc sort data=temp1;
by &exp cohort ;
run;
data cp_&exp._&case.;
set temp1 ;
by &exp;
if first.&exp;
casepy_&case._all=(put(sumcase_all,comma20.) || '/' || put(sumpy_all,comma20. -L));
run;

proc transpose data=cp_&exp._&case. out=cp_&exp._&case.;
id &exp;
var casepy_&case._all;
run;

data cp_&exp._&case.;
length index $6. cohort $6.;
retain cohort index ;
set  cp_&exp._&case.;
index="&exp";
cohort="Pool";
run;
proc print data=cp_&exp._&case.;
run;
%mend;

%cpy(lcssvq,    can_obese);
%cpy(alcssvq,   can_obese);
%cpy(vlcssvq,   can_obese);
%cpy(ulcssvq,   can_obese);
%cpy(hlcssvq,   can_obese);
%cpy(lfssvq,    can_obese);
%cpy(alfssvq,   can_obese);
%cpy(vlfssvq,   can_obese);
%cpy(ulfssvq,   can_obese);
%cpy(hlfssvq,   can_obese);

data casepy_can_obese;
    length index $6.;
    set  
    cp_lcssvq_can_obese
    cp_alcssvq_can_obese
    cp_vlcssvq_can_obese
    cp_ulcssvq_can_obese
    cp_hlcssvq_can_obese
	cp_lfssvq_can_obese
    cp_alfssvq_can_obese
    cp_vlfssvq_can_obese
    cp_ulfssvq_can_obese
    cp_hlfssvq_can_obese
    ;
run;

proc export data=casepy_can_obese
    dbms=csv
    outfile='lowcarb_cancer/incidence/without_stop/table/casePY_pool_can_obese.csv'
    replace;
    run;

/************************ HR *****************************************/
%include 'macro/res.sas';

%res(analysis=avestop.adj0, catlin=T,   
         data=pool, event=can_obese, time=tobeca, 
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         agevar=agecon,
         strata=agecon cohort intervalp, 
         labels=F, 
         exp1=&lcssvq_   ,    adjust1=   ,
         exp2=&alcssvq_  ,    adjust2=   ,
         exp3=&vlcssvq_  ,    adjust3=   ,
         exp4=&ulcssvq_  ,    adjust4=   ,
         exp5=&hlcssvq_  ,    adjust5=   , 
         exp6=&lfssvq_   ,    adjust6=   ,
         exp7=&ulfssvq_  ,    adjust7=   ,
         exp8=&hlfssvq_  ,    adjust8=   ,
         exp9=&alfssvq_  ,    adjust9=   ,
         exp10=&vlfssvq_  ,    adjust10=   );
         Run;

%res(analysis=avestop.adj1,catlin=T,   
         data=pool, event=can_obese, time=tobeca, 
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         agevar=agecon,
         strata=agecon cohort intervalp, 
         labels=F, 
         exp1=&lcssvq_   ,   adjust1= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp2=&alcssvq_  ,   adjust2= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp3=&vlcssvq_  ,   adjust3= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp4=&ulcssvq_  ,   adjust4= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp5=&hlcssvq_  ,   adjust5= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  , 
         exp6=&lfssvq_   ,   adjust6= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp7=&ulfssvq_  ,   adjust7= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp8=&hlfssvq_  ,   adjust8= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp9=&alfssvq_  ,   adjust9= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp10=&vlfssvq_  ,  adjust10=&calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  );
         Run;

%res(analysis=unit.adj0, catlin=T,   
         data=pool, event=can_obese, time=tobeca, 
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         agevar=agecon,
         strata=agecon cohort intervalp, 
         labels=F, 
         exp1=lcs_sd   ,    adjust1=  ,
         exp2=alcs_sd  ,    adjust2=  ,
         exp3=vlcs_sd  ,    adjust3=  ,
         exp4=ulcs_sd  ,    adjust4=  ,
         exp5=hlcs_sd  ,    adjust5=  , 
         exp6=lfs_sd   ,    adjust6=  ,
         exp7=ulfs_sd  ,    adjust7=  ,
         exp8=hlfs_sd  ,    adjust8=  ,
         exp9=alfs_sd  ,    adjust9=  ,
         exp10=vlfs_sd  ,   adjust10=  );      
         Run;

%res(analysis=unit.adj1,catlin=T,   
         data=pool, event=can_obese, time=tobeca, 
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         agevar=agecon,
         strata=agecon cohort intervalp, 
         labels=F, 
         exp1=lcs_sd   ,   adjust1= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp2=alcs_sd  ,   adjust2= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp3=vlcs_sd  ,   adjust3= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp4=ulcs_sd  ,   adjust4= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp5=hlcs_sd  ,   adjust5= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  , 
         exp6=lfs_sd   ,   adjust6= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp7=ulfs_sd  ,   adjust7= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp8=hlfs_sd  ,   adjust8= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp9=alfs_sd  ,   adjust9= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp10=vlfs_sd  ,  adjust10=&calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ); 
         Run;

%res(analysis=trend.adj0, catlin=T,   
         data=pool, event=can_obese, time=tobeca, 
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         agevar=agecon,
         strata=agecon cohort intervalp, 
         labels=F, 
         exp1=lcssvtrend  ,     adjust1=  ,
         exp2=alcssvtrend ,     adjust2=  ,
         exp3=vlcssvtrend ,     adjust3=  ,
         exp4=ulcssvtrend ,     adjust4=  ,
         exp5=hlcssvtrend ,     adjust5=  , 
         exp6=lfssvtrend  ,     adjust6=  ,
         exp7=ulfssvtrend ,     adjust7=  ,
         exp8=hlfssvtrend ,     adjust8=  ,
         exp9=alfssvtrend ,     adjust9=  ,
         exp10=vlfssvtrend ,    adjust10=  );
         Run;

%res(analysis=trend.adj1,catlin=T,   
         data=pool, event=can_obese, time=tobeca, 
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         agevar=agecon,
         strata=agecon cohort intervalp, 
         labels=F, 
         exp1=lcssvtrend  ,    adjust1= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp2=alcssvtrend ,    adjust2= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp3=vlcssvtrend ,    adjust3= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp4=ulcssvtrend ,    adjust4= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp5=hlcssvtrend ,    adjust5= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  , 
         exp6=lfssvtrend  ,    adjust6= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp7=ulfssvtrend ,    adjust7= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp8=hlfssvtrend ,    adjust8= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp9=alfssvtrend ,    adjust9= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  ,
         exp10=vlfssvtrend ,   adjust10=&calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam  );
         Run;

%include 'macro/resout.sas';
/************************ HR *****************************************/

%resout(outfile=table/tableS20_withoutstop, dbms=tab, opt=replace);                                                                                                                             
endsas;