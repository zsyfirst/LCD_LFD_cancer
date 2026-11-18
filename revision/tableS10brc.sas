/*****************************************************
Purpose:       low carb, low fat, and risk of Obesity-related cancer
Study design:  prospective cohort (NHS: 1986-2022.6, NHSII: 1991-2023.6, HPFS: 1986-2024.1)

Contents:      Table S10 --- PY and HR (individual cancer as outcome: Breast cancer)
Date:          10/2026
******************************************************/


/************************ Call in Cleaned Datasets**********************************/
%include '/lowcarb_cancer/revision/data.inc.sas';

/* Breast cancer */ 
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
set &cohort._&exp._&case;
index="&exp";
run;
%mend;

/*==== NHS1 ====*/
%pre_pm(data=nhs1,out=TWO,case=can_brca,timevar=period,
irt=irt86 irt88 irt90 irt92 irt94 irt96 irt98 irt00 irt02 irt04 irt06 irt08 irt10,cutoff=cutoff,
dtdx=dtdx_brca, dtdth=dtdth, var=agegp lcssvq); 
%pm(data=TWO, case=can_brca, exposure=lcssvq, strata=agegp);
%casepy(cohort=nhs1, exp=lcssvq, case=can_brca);

%pre_pm(data=nhs1,out=TWO,case=can_brca,timevar=period,
irt=irt86 irt88 irt90 irt92 irt94 irt96 irt98 irt00 irt02 irt04 irt06 irt08 irt10,cutoff=cutoff,
dtdx=dtdx_brca, dtdth=dtdth, var=agegp alcssvq); 
%pm(data=TWO, case=can_brca, exposure=alcssvq, strata=agegp);
%casepy(cohort=nhs1, exp=alcssvq, case=can_brca);

%pre_pm(data=nhs1,out=TWO,case=can_brca,timevar=period,
irt=irt86 irt88 irt90 irt92 irt94 irt96 irt98 irt00 irt02 irt04 irt06 irt08 irt10,cutoff=cutoff,
dtdx=dtdx_brca, dtdth=dtdth, var=agegp vlcssvq); 
%pm(data=TWO, case=can_brca, exposure=vlcssvq, strata=agegp);
%casepy(cohort=nhs1, exp=vlcssvq, case=can_brca);

%pre_pm(data=nhs1,out=TWO,case=can_brca,timevar=period,
irt=irt86 irt88 irt90 irt92 irt94 irt96 irt98 irt00 irt02 irt04 irt06 irt08 irt10,cutoff=cutoff,
dtdx=dtdx_brca, dtdth=dtdth, var=agegp ulcssvq); 
%pm(data=TWO, case=can_brca, exposure=ulcssvq, strata=agegp);
%casepy(cohort=nhs1, exp=ulcssvq, case=can_brca);

%pre_pm(data=nhs1,out=TWO,case=can_brca,timevar=period,
irt=irt86 irt88 irt90 irt92 irt94 irt96 irt98 irt00 irt02 irt04 irt06 irt08 irt10,cutoff=cutoff,
dtdx=dtdx_brca, dtdth=dtdth, var=agegp hlcssvq); 
%pm(data=TWO, case=can_brca, exposure=hlcssvq, strata=agegp);
%casepy(cohort=nhs1, exp=hlcssvq, case=can_brca);

%pre_pm(data=nhs1,out=TWO,case=can_brca,timevar=period,
irt=irt86 irt88 irt90 irt92 irt94 irt96 irt98 irt00 irt02 irt04 irt06 irt08 irt10,cutoff=cutoff,
dtdx=dtdx_brca, dtdth=dtdth, var=agegp lfssvq); 
%pm(data=TWO, case=can_brca, exposure=lfssvq, strata=agegp);
%casepy(cohort=nhs1, exp=lfssvq, case=can_brca);

%pre_pm(data=nhs1,out=TWO,case=can_brca,timevar=period,
irt=irt86 irt88 irt90 irt92 irt94 irt96 irt98 irt00 irt02 irt04 irt06 irt08 irt10,cutoff=cutoff,
dtdx=dtdx_brca, dtdth=dtdth, var=agegp alfssvq); 
%pm(data=TWO, case=can_brca, exposure=alfssvq, strata=agegp);
%casepy(cohort=nhs1, exp=alfssvq, case=can_brca);

%pre_pm(data=nhs1,out=TWO,case=can_brca,timevar=period,
irt=irt86 irt88 irt90 irt92 irt94 irt96 irt98 irt00 irt02 irt04 irt06 irt08 irt10,cutoff=cutoff,
dtdx=dtdx_brca, dtdth=dtdth, var=agegp vlfssvq); 
%pm(data=TWO, case=can_brca, exposure=vlfssvq, strata=agegp);
%casepy(cohort=nhs1, exp=vlfssvq, case=can_brca);

%pre_pm(data=nhs1,out=TWO,case=can_brca,timevar=period,
irt=irt86 irt88 irt90 irt92 irt94 irt96 irt98 irt00 irt02 irt04 irt06 irt08 irt10,cutoff=cutoff,
dtdx=dtdx_brca, dtdth=dtdth, var=agegp ulfssvq); 
%pm(data=TWO, case=can_brca, exposure=ulfssvq, strata=agegp);
%casepy(cohort=nhs1, exp=ulfssvq, case=can_brca);

%pre_pm(data=nhs1,out=TWO,case=can_brca,timevar=period,
irt=irt86 irt88 irt90 irt92 irt94 irt96 irt98 irt00 irt02 irt04 irt06 irt08 irt10,cutoff=cutoff,
dtdx=dtdx_brca, dtdth=dtdth, var=agegp hlfssvq); 
%pm(data=TWO, case=can_brca, exposure=hlfssvq, strata=agegp);
%casepy(cohort=nhs1, exp=hlfssvq, case=can_brca);



/*==== NHS2 ====*/
%pre_pm(data=nhs2,out=TWO,case=can_brca,timevar=period,
irt=retmo91 retmo93 retmo95 retmo97 retmo99 retmo01 retmo03 retmo05 retmo07 retmo09 retmo11 retmo13 retmo15 retmo17 retmo19,
cutoff=cutoff, dtdx=dtdx_brca, dtdth=dtdth, var=agegp lcssvq); 
%pm(data=TWO, case=can_brca, exposure=lcssvq, strata=agegp);
%casepy(cohort=nhs2, exp=lcssvq, case=can_brca);

%pre_pm(data=nhs2,out=TWO,case=can_brca,timevar=period,
irt=retmo91 retmo93 retmo95 retmo97 retmo99 retmo01 retmo03 retmo05 retmo07 retmo09 retmo11 retmo13 retmo15 retmo17 retmo19,
cutoff=cutoff, dtdx=dtdx_brca, dtdth=dtdth, var=agegp alcssvq); 
%pm(data=TWO, case=can_brca, exposure=alcssvq, strata=agegp);
%casepy(cohort=nhs2, exp=alcssvq, case=can_brca);

%pre_pm(data=nhs2,out=TWO,case=can_brca,timevar=period,
irt=retmo91 retmo93 retmo95 retmo97 retmo99 retmo01 retmo03 retmo05 retmo07 retmo09 retmo11 retmo13 retmo15 retmo17 retmo19,
cutoff=cutoff, dtdx=dtdx_brca, dtdth=dtdth, var=agegp vlcssvq); 
%pm(data=TWO, case=can_brca, exposure=vlcssvq, strata=agegp);
%casepy(cohort=nhs2, exp=vlcssvq, case=can_brca);

%pre_pm(data=nhs2,out=TWO,case=can_brca,timevar=period,
irt=retmo91 retmo93 retmo95 retmo97 retmo99 retmo01 retmo03 retmo05 retmo07 retmo09 retmo11 retmo13 retmo15 retmo17 retmo19,
cutoff=cutoff, dtdx=dtdx_brca, dtdth=dtdth, var=agegp ulcssvq); 
%pm(data=TWO, case=can_brca, exposure=ulcssvq, strata=agegp);
%casepy(cohort=nhs2, exp=ulcssvq, case=can_brca);

%pre_pm(data=nhs2,out=TWO,case=can_brca,timevar=period,
irt=retmo91 retmo93 retmo95 retmo97 retmo99 retmo01 retmo03 retmo05 retmo07 retmo09 retmo11 retmo13 retmo15 retmo17 retmo19,
cutoff=cutoff, dtdx=dtdx_brca, dtdth=dtdth, var=agegp hlcssvq); 
%pm(data=TWO, case=can_brca, exposure=hlcssvq, strata=agegp);
%casepy(cohort=nhs2, exp=hlcssvq, case=can_brca);

%pre_pm(data=nhs2,out=TWO,case=can_brca,timevar=period,
irt=retmo91 retmo93 retmo95 retmo97 retmo99 retmo01 retmo03 retmo05 retmo07 retmo09 retmo11 retmo13 retmo15 retmo17 retmo19,
cutoff=cutoff, dtdx=dtdx_brca, dtdth=dtdth, var=agegp lfssvq); 
%pm(data=TWO, case=can_brca, exposure=lfssvq, strata=agegp);
%casepy(cohort=nhs2, exp=lfssvq, case=can_brca);

%pre_pm(data=nhs2,out=TWO,case=can_brca,timevar=period,
irt=retmo91 retmo93 retmo95 retmo97 retmo99 retmo01 retmo03 retmo05 retmo07 retmo09 retmo11 retmo13 retmo15 retmo17 retmo19,
cutoff=cutoff, dtdx=dtdx_brca, dtdth=dtdth, var=agegp alfssvq); 
%pm(data=TWO, case=can_brca, exposure=alfssvq, strata=agegp);
%casepy(cohort=nhs2, exp=alfssvq, case=can_brca);

%pre_pm(data=nhs2,out=TWO,case=can_brca,timevar=period,
irt=retmo91 retmo93 retmo95 retmo97 retmo99 retmo01 retmo03 retmo05 retmo07 retmo09 retmo11 retmo13 retmo15 retmo17 retmo19,
cutoff=cutoff, dtdx=dtdx_brca, dtdth=dtdth, var=agegp vlfssvq); 
%pm(data=TWO, case=can_brca, exposure=vlfssvq, strata=agegp);
%casepy(cohort=nhs2, exp=vlfssvq, case=can_brca);

%pre_pm(data=nhs2,out=TWO,case=can_brca,timevar=period,
irt=retmo91 retmo93 retmo95 retmo97 retmo99 retmo01 retmo03 retmo05 retmo07 retmo09 retmo11 retmo13 retmo15 retmo17 retmo19,
cutoff=cutoff, dtdx=dtdx_brca, dtdth=dtdth, var=agegp ulfssvq); 
%pm(data=TWO, case=can_brca, exposure=ulfssvq, strata=agegp);
%casepy(cohort=nhs2, exp=ulfssvq, case=can_brca);

%pre_pm(data=nhs2,out=TWO,case=can_brca,timevar=period,
irt=retmo91 retmo93 retmo95 retmo97 retmo99 retmo01 retmo03 retmo05 retmo07 retmo09 retmo11 retmo13 retmo15 retmo17 retmo19,
cutoff=cutoff, dtdx=dtdx_brca, dtdth=dtdth, var=agegp hlfssvq); 
%pm(data=TWO, case=can_brca, exposure=hlfssvq, strata=agegp);
%casepy(cohort=nhs2, exp=hlfssvq, case=can_brca);

/*==== HPFS ====*/
%pre_pm(data=hpfs,out=TWO,case=can_brca,timevar=period,
irt=rtmnyr86 rtmnyr88 rtmnyr90 rtmnyr92 rtmnyr94 rtmnyr96 rtmnyr98 rtmnyr00 rtmnyr02 rtmnyr04 rtmnyr06 rtmnyr08 rtmnyr10 rtmnyr12 rtmnyr14 rtmnyr16 rtmnyr18,cutoff=cutoff,
dtdx=dtdx_brca, dtdth=dtdth, var=agegp lcssvq); 
%pm(data=TWO, case=can_brca, exposure=lcssvq, strata=agegp);
%casepy(cohort=hpfs, exp=lcssvq, case=can_brca);

%pre_pm(data=hpfs,out=TWO,case=can_brca,timevar=period,
irt=rtmnyr86 rtmnyr88 rtmnyr90 rtmnyr92 rtmnyr94 rtmnyr96 rtmnyr98 rtmnyr00 rtmnyr02 rtmnyr04 rtmnyr06 rtmnyr08 rtmnyr10 rtmnyr12 rtmnyr14 rtmnyr16 rtmnyr18,cutoff=cutoff,
dtdx=dtdx_brca, dtdth=dtdth, var=agegp alcssvq); 
%pm(data=TWO, case=can_brca, exposure=alcssvq, strata=agegp);
%casepy(cohort=hpfs, exp=alcssvq, case=can_brca);

%pre_pm(data=hpfs,out=TWO,case=can_brca,timevar=period,
irt=rtmnyr86 rtmnyr88 rtmnyr90 rtmnyr92 rtmnyr94 rtmnyr96 rtmnyr98 rtmnyr00 rtmnyr02 rtmnyr04 rtmnyr06 rtmnyr08 rtmnyr10 rtmnyr12 rtmnyr14 rtmnyr16 rtmnyr18,cutoff=cutoff,
dtdx=dtdx_brca, dtdth=dtdth, var=agegp vlcssvq); 
%pm(data=TWO, case=can_brca, exposure=vlcssvq, strata=agegp);
%casepy(cohort=hpfs, exp=vlcssvq, case=can_brca);

%pre_pm(data=hpfs,out=TWO,case=can_brca,timevar=period,
irt=rtmnyr86 rtmnyr88 rtmnyr90 rtmnyr92 rtmnyr94 rtmnyr96 rtmnyr98 rtmnyr00 rtmnyr02 rtmnyr04 rtmnyr06 rtmnyr08 rtmnyr10 rtmnyr12 rtmnyr14 rtmnyr16 rtmnyr18,cutoff=cutoff,
dtdx=dtdx_brca, dtdth=dtdth, var=agegp ulcssvq); 
%pm(data=TWO, case=can_brca, exposure=ulcssvq, strata=agegp);
%casepy(cohort=hpfs, exp=ulcssvq, case=can_brca);

%pre_pm(data=hpfs,out=TWO,case=can_brca,timevar=period,
irt=rtmnyr86 rtmnyr88 rtmnyr90 rtmnyr92 rtmnyr94 rtmnyr96 rtmnyr98 rtmnyr00 rtmnyr02 rtmnyr04 rtmnyr06 rtmnyr08 rtmnyr10 rtmnyr12 rtmnyr14 rtmnyr16 rtmnyr18,cutoff=cutoff,
dtdx=dtdx_brca, dtdth=dtdth, var=agegp hlcssvq); 
%pm(data=TWO, case=can_brca, exposure=hlcssvq, strata=agegp);
%casepy(cohort=hpfs, exp=hlcssvq, case=can_brca);

%pre_pm(data=hpfs,out=TWO,case=can_brca,timevar=period,
irt=rtmnyr86 rtmnyr88 rtmnyr90 rtmnyr92 rtmnyr94 rtmnyr96 rtmnyr98 rtmnyr00 rtmnyr02 rtmnyr04 rtmnyr06 rtmnyr08 rtmnyr10 rtmnyr12 rtmnyr14 rtmnyr16 rtmnyr18,cutoff=cutoff,
dtdx=dtdx_brca, dtdth=dtdth, var=agegp lfssvq); 
%pm(data=TWO, case=can_brca, exposure=lfssvq, strata=agegp);
%casepy(cohort=hpfs, exp=lfssvq, case=can_brca);

%pre_pm(data=hpfs,out=TWO,case=can_brca,timevar=period,
irt=rtmnyr86 rtmnyr88 rtmnyr90 rtmnyr92 rtmnyr94 rtmnyr96 rtmnyr98 rtmnyr00 rtmnyr02 rtmnyr04 rtmnyr06 rtmnyr08 rtmnyr10 rtmnyr12 rtmnyr14 rtmnyr16 rtmnyr18,cutoff=cutoff,
dtdx=dtdx_brca, dtdth=dtdth, var=agegp alfssvq); 
%pm(data=TWO, case=can_brca, exposure=alfssvq, strata=agegp);
%casepy(cohort=hpfs, exp=alfssvq, case=can_brca);

%pre_pm(data=hpfs,out=TWO,case=can_brca,timevar=period,
irt=rtmnyr86 rtmnyr88 rtmnyr90 rtmnyr92 rtmnyr94 rtmnyr96 rtmnyr98 rtmnyr00 rtmnyr02 rtmnyr04 rtmnyr06 rtmnyr08 rtmnyr10 rtmnyr12 rtmnyr14 rtmnyr16 rtmnyr18,cutoff=cutoff,
dtdx=dtdx_brca, dtdth=dtdth, var=agegp vlfssvq); 
%pm(data=TWO, case=can_brca, exposure=vlfssvq, strata=agegp);
%casepy(cohort=hpfs, exp=vlfssvq, case=can_brca);

%pre_pm(data=hpfs,out=TWO,case=can_brca,timevar=period,
irt=rtmnyr86 rtmnyr88 rtmnyr90 rtmnyr92 rtmnyr94 rtmnyr96 rtmnyr98 rtmnyr00 rtmnyr02 rtmnyr04 rtmnyr06 rtmnyr08 rtmnyr10 rtmnyr12 rtmnyr14 rtmnyr16 rtmnyr18,cutoff=cutoff,
dtdx=dtdx_brca, dtdth=dtdth, var=agegp ulfssvq); 
%pm(data=TWO, case=can_brca, exposure=ulfssvq, strata=agegp);
%casepy(cohort=hpfs, exp=ulfssvq, case=can_brca);

%pre_pm(data=hpfs,out=TWO,case=can_brca,timevar=period,
irt=rtmnyr86 rtmnyr88 rtmnyr90 rtmnyr92 rtmnyr94 rtmnyr96 rtmnyr98 rtmnyr00 rtmnyr02 rtmnyr04 rtmnyr06 rtmnyr08 rtmnyr10 rtmnyr12 rtmnyr14 rtmnyr16 rtmnyr18,cutoff=cutoff,
dtdx=dtdx_brca, dtdth=dtdth, var=agegp hlfssvq); 
%pm(data=TWO, case=can_brca, exposure=hlfssvq, strata=agegp);
%casepy(cohort=hpfs, exp=hlfssvq, case=can_brca);

/* Print */
data printcasepy;
    set 
    NHS1_LCSSVQ_can_brcaPRINT
    NHS1_ALCSSVQ_can_brcaPRINT
    NHS1_VLCSSVQ_can_brcaPRINT
    NHS1_ULCSSVQ_can_brcaPRINT
    NHS1_HLCSSVQ_can_brcaPRINT
    NHS1_LFSSVQ_can_brcaPRINT
    NHS1_ALFSSVQ_can_brcaPRINT
    NHS1_VLFSSVQ_can_brcaPRINT
    NHS1_ULFSSVQ_can_brcaPRINT
    NHS1_HLFSSVQ_can_brcaPRINT

    NHS2_LCSSVQ_can_brcaPRINT
    NHS2_ALCSSVQ_can_brcaPRINT
    NHS2_VLCSSVQ_can_brcaPRINT
    NHS2_ULCSSVQ_can_brcaPRINT
    NHS2_HLCSSVQ_can_brcaPRINT
    NHS2_LFSSVQ_can_brcaPRINT
    NHS2_ALFSSVQ_can_brcaPRINT
    NHS2_VLFSSVQ_can_brcaPRINT
    NHS2_ULFSSVQ_can_brcaPRINT
    NHS2_HLFSSVQ_can_brcaPRINT

    HPFS_LCSSVQ_can_brcaPRINT
    HPFS_ALCSSVQ_can_brcaPRINT
    HPFS_VLCSSVQ_can_brcaPRINT
    HPFS_ULCSSVQ_can_brcaPRINT
    HPFS_HLCSSVQ_can_brcaPRINT
    HPFS_LFSSVQ_can_brcaPRINT
    HPFS_ALFSSVQ_can_brcaPRINT
    HPFS_VLFSSVQ_can_brcaPRINT
    HPFS_ULFSSVQ_can_brcaPRINT
    HPFS_HLFSSVQ_can_brcaPRINT
    ;
run;

proc export data=printcasepy
    dbms=csv
    outfile='lowcarb_cancer/revision/site_specific/table/casePY_bycohort_can_brca.csv'
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
retain cohort index ;
length index $32 cohort $8;
set  cp_&exp._&case.;
index="&exp";
cohort="Pool";
run;
proc print data=cp_&exp._&case.;
run;
%mend;

%cpy(lcssvq,    can_brca);
%cpy(alcssvq,   can_brca);
%cpy(vlcssvq,   can_brca);
%cpy(ulcssvq,   can_brca);
%cpy(hlcssvq,   can_brca);
%cpy(lfssvq,    can_brca);
%cpy(alfssvq,   can_brca);
%cpy(vlfssvq,   can_brca);
%cpy(ulfssvq,   can_brca);
%cpy(hlfssvq,   can_brca);

data casepy_can_brca;
 length index $32 cohort $8;
    set  
    cp_lcssvq_can_brca
    cp_alcssvq_can_brca
    cp_vlcssvq_can_brca
    cp_ulcssvq_can_brca
    cp_hlcssvq_can_brca
    cp_lfssvq_can_brca
    cp_alfssvq_can_brca
    cp_vlfssvq_can_brca
    cp_ulfssvq_can_brca
    cp_hlfssvq_can_brca
    ;
run;

proc export data=casepy_can_brca
    dbms=csv
    outfile='lowcarb_cancer/revision/site_specific/table/casePY_pool_can_brca.csv'
    replace;
    run;

/************************ HR *****************************************/
%include 'macro/res.sas';

%res(analysis=avestop.adj0, catlin=T,   
         data=pool, event=can_brca, time=tbrca, 
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
         data=pool, event=can_brca, time=tbrca, 
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         agevar=agecon,
         strata=agecon cohort intervalp, 
         labels=F, 
         exp1=&lcssvq_   ,   adjust1= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam nSES  ,
         exp2=&alcssvq_  ,   adjust2= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam nSES  ,
         exp3=&vlcssvq_  ,   adjust3= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam nSES  ,
         exp4=&ulcssvq_  ,   adjust4= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam nSES  ,
         exp5=&hlcssvq_  ,   adjust5= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam nSES  , 
         exp6=&lfssvq_   ,   adjust6= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam nSES  ,
         exp7=&ulfssvq_  ,   adjust7= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam nSES  ,
         exp8=&hlfssvq_  ,   adjust8= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam nSES  ,
         exp9=&alfssvq_  ,   adjust9= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam nSES  ,
         exp10=&vlfssvq_  ,  adjust10=&calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam nSES  );
         Run;

%res(analysis=unit.adj0, catlin=T,   
         data=pool, event=can_brca, time=tbrca, 
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
         data=pool, event=can_brca, time=tbrca, 
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         agevar=agecon,
         strata=agecon cohort intervalp, 
         labels=F, 
         exp1=lcs_sd   ,   adjust1= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam nSES  ,
         exp2=alcs_sd  ,   adjust2= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam nSES  ,
         exp3=vlcs_sd  ,   adjust3= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam nSES  ,
         exp4=ulcs_sd  ,   adjust4= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam nSES  ,
         exp5=hlcs_sd  ,   adjust5= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam nSES  , 
         exp6=lfs_sd   ,   adjust6= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam nSES  ,
         exp7=ulfs_sd  ,   adjust7= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam nSES  ,
         exp8=hlfs_sd  ,   adjust8= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam nSES  ,
         exp9=alfs_sd  ,   adjust9= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam nSES  ,
         exp10=vlfs_sd  ,  adjust10=&calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam nSES  ); 
         Run;

%res(analysis=trend.adj0, catlin=T,   
         data=pool, event=can_brca, time=tbrca, 
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
         data=pool, event=can_brca, time=tbrca, 
         timevar=tp01 tp02 tp03 tp04 tp05 tp06 tp07 tp08 tp09 tp10 tp11 tp12 tp13 tp14 tp15 tp16 tp17 tp18,
         id=newid, 
         tvar=intervalp,
         agevar=agecon,
         strata=agecon cohort intervalp, 
         labels=F, 
         exp1=lcssvtrend  ,    adjust1= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam nSES  ,
         exp2=alcssvtrend ,    adjust2= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam nSES  ,
         exp3=vlcssvtrend ,    adjust3= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam nSES  ,
         exp4=ulcssvtrend ,    adjust4= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam nSES  ,
         exp5=hlcssvtrend ,    adjust5= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam nSES  , 
         exp6=lfssvtrend  ,    adjust6= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam nSES  ,
         exp7=ulfssvtrend ,    adjust7= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam nSES  ,
         exp8=hlfssvtrend ,    adjust8= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam nSES  ,
         exp9=alfssvtrend ,    adjust9= &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam nSES  ,
         exp10=vlfssvtrend ,   adjust10=&calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam nSES  );
         Run;

%include 'macro/resout.sas';
/************************ HR *****************************************/

%resout(outfile=table/tableS10_can_breast, dbms=tab, opt=replace);                                                                                                                             
endsas;