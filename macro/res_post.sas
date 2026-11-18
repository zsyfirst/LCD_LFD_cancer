*options mprint symbolgen;
%macro res(analysis=,catlin=T,
           data=,event=,time=,timevar=,id=,agevar=,tvar=,contvar=,
           labels=F,outest=,covout=F,byvar=,where=,strata=,modopt=,
           tbegtend=F,
           tbeg=,
           tend=,
           exp1=,adjust1=,
           exp2=,adjust2=,
           exp3=,adjust3=,
           exp4=,adjust4=,
           exp5=,adjust5=,
           exp6=,adjust6=,
           exp7=,adjust7=,
           exp8=,adjust8=,
           exp9=,adjust9=,
           exp10=,adjust10=,
           exp11=,adjust11=,
           exp12=,adjust12=,
           exp13=,adjust13=,
           exp14=,adjust14=,
           exp15=,adjust15=,
           exp16=,adjust16=,
           exp17=,adjust17=,
           exp18=,adjust18=,
           exp19=,adjust19=,
           exp20=,adjust20=);

%local blank num adjust covlist;
%global resno resxx tresxx;
%let blank=%quote( );
%if %length(&resno)=0 %then %let resno=1;
%else %let resno=%eval(&resno.+1);
%if %length(&resxx)=0 %then %let resxx=res&resno;
%else %let resxx=&resxx.&blank.res&resno.;
%do i=1 %to 20;
    %let num&i=%eval(%sysfunc(count(trim(left(&&adjust&i)),&blank.))+1);
%end;

%mphreg9(data=&data, event=&event, 
         timevar=&timevar,
         id=&id, tvar=&tvar,
         contvar=&contvar,byvar=&byvar,where=&where,strata=&strata,modopt=&modopt,
         outest=&outest.,covout=&covout.,
         outdat=res&resno.,
         labels=&labels,
         tbegtend=&tbegtend,
         tbeg=&tbeg,
         tend=&tend,
         model1= &exp1. &adjust1.,
         model2= &exp2. &adjust2.,
         model3= &exp3. &adjust3.,
         model4= &exp4. &adjust4.,
         model5= &exp5. &adjust5.,
         model6= &exp6. &adjust6.,
         model7= &exp7. &adjust7.,
         model8= &exp8. &adjust8.,
         model9= &exp9. &adjust9.,
         model10= &exp10. &adjust10.,
         model11= &exp11. &adjust11.,
         model12= &exp12. &adjust12.,
         model13= &exp13. &adjust13.,
         model14= &exp14. &adjust14.,
         model15= &exp15. &adjust15.,
         model16= &exp16. &adjust16.,
         model17= &exp17. &adjust17.,
         model18= &exp18. &adjust18.,
         model19= &exp19. &adjust19.,
         model20= &exp20. &adjust20.);
run;
data res&resno.;
    set res&resno. end=_end_;
    drop hazardratio;
run;
data res&resno.;
    set res&resno. end=_end_;
    %do j=1 %to 20;
        %if %length(&&adjust&j)>0 %then %do;
            %do k=1 %to &&num&j;
                if modelno=&j and variable="%scan(&&adjust&j.,&k,&blank.)" then delete;
            %end;
        %end;
    %end;
    format no $ 8. event $ 8. analysis $ 20.;
    no="&resno.";
    event="&event";
    analysis="&analysis.";
    hazardratio=exp(estimate);
/*    a=int(sum(hazardratio,0.005)/10**2)*10**2;
    b=int(sum(lcl,0.005)/10**2)*10**2;
    c=int(sum(ucl,0.005)/10**2)*10**2;*/
    hr=put(round(hazardratio,0.01),4.2);
    ll=put(round(lcl,0.01),4.2);
    ul=put(round(ucl,0.01),4.2);
    result=trim(hr)||" ("||trim(ll)||", "||trim(ul)||")";

    %if &catlin=T %then %do;
        if mod(modelno,2)=0 then modelno=modelno-1;
    %end;
    prev=substr(variable,1,length(variable)-1);
    postv=substr(variable,length(variable),1);
keep no modelno analysis &byvar variable estimate stderr probchisq hazardratio lcl ucl result event prev postv;
run;
data res&resno.;
    format no modelno analysis &byvar variable estimate stderr probchisq hazardratio lcl ucl result event prev postv;
    set res&resno. end=_end_;
run;

%if &covout ne T %then %do;
%if %length(&tresxx)=0 %then %let tresxx=tres&resno;
%else %let tresxx=&tresxx.&blank.tres&resno.;
proc sort data=res&resno.;
    by analysis &byvar modelno prev;
run;
proc transpose data=res&resno. out=tres&resno. prefix=level;
    by analysis &byvar modelno prev;
    id postv;
    var result;
run;
%end;

%if &covout=T %then %do;
    %do i=1 %to 20;
        %if %length(&&exp&i)=0 %then %goto skip;
        %let covnum&i=%eval(%sysfunc(count(trim(left(&&exp&i)),&blank.))+2);
        %if %length(&covlist)=0 %then %let covlist=cov&i;
        %else %let covlist=&covlist.&blank.cov&i.;
        data cov&i;
            set cov&i(rename=(_name_=variable) firstobs=2 obs=&&covnum&i.) end=_end_;
            modelno=&i.;
            keep variable modelno &&exp&i;
        run;
        %skip:
    %end;
    data cov;
        set &covlist end=_end_;
    run;
    proc sort data=res&resno.;by modelno variable;run;
    proc sort data=cov;by modelno variable;run;
    data res&resno.;
        merge res&resno. cov;
        by modelno variable;
    run;
    proc datasets noprint;
        delete cov &covlist;
    run;
%end;
%mend;
