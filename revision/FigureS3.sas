/*****************************************************
Purpose:       low carb, low fat, and risk of Obesity-related cancer
Study design:  prospective cohort (NHS: 1986-2022.6, NHSII: 1991-2023.6, HPFS: 1986-2024.1)

Contents:      Figure S3 --- Cubic Spline (Obesity-related cancer as outcome)
Revision:      Add nSES for cavariates
Date:          10/2026
******************************************************/


/************************ Call in Cleaned Datasets**********************************/
%include 'lowcarb_cancer/revision/data.inc.sas';


/************ 4 knots ************/

%LGTPHCURV9(    data    = pool, 
        exposure= lcssvcon, 
        case    = can_obese, 
        model   = cox, 
        time    = tobeca,
        strata  = agecon cohort intervalp ,
        adj = &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam nSES    , 
        refval  = MIN,  
        lpct    = 1, 
        hpct    = 99,  
        nk  = 4, 
        SELECT  = 1,     /* force all spline variables */
        /*SELECT   = 3,*/  /*use automatic stepwise selection*/    
        testrep  = long,     
        outplot = JPEG,
        pwhich  = spline,  
        footer  = NONE,
        GRAPHTIT= NONE,
        plot    = 4,
        e   = T,     
        ci  = 2,     
        displayx= T,  
        klines  = T,   
        plotdec = T,
        PLOTDATA = pool.lcs.force_min.txt);  /* Output filename */

%LGTPHCURV9(    data    = pool, 
        exposure= alcssvcon, 
        case    = can_obese, 
        model   = cox, 
        time    = tobeca,
        strata  = agecon cohort intervalp ,
        adj = &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam nSES    , 
        refval  = MIN,  
        lpct    = 1, 
        hpct    = 99,  
        nk  = 4, 
        SELECT  = 1,     /* force all spline variables */      
        testrep  = long,
        outplot = JPEG,
        pwhich  = spline,  
        footer  = NONE,
        GRAPHTIT= NONE,
        plot    = 4,
        e   = T,     
        ci  = 2,     
        displayx= T,  
        klines  = T,   
        plotdec = T,
        PLOTDATA = pool.alcs.force_min.txt);  /* Output filename */

%LGTPHCURV9(    data    = pool, 
        exposure= vlcssvcon, 
        case    = can_obese, 
        model   = cox, 
        time    = tobeca,
        strata  = agecon cohort intervalp ,
        adj = &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam nSES    , 
        refval  = MIN,  
        lpct    = 1, 
        hpct    = 99,  
        nk  = 4, 
        SELECT  = 1,     /* force all spline variables */      
        testrep  = long,
        outplot = JPEG,
        pwhich  = spline,  
        footer  = NONE,
        GRAPHTIT= NONE,
        plot    = 4,
        e   = T,     
        ci  = 2,     
        displayx= T,  
        klines  = T,   
        plotdec = T,
        PLOTDATA = pool.vlcs.force_min.txt);  /* Output filename */



%LGTPHCURV9(    data    = pool, 
        exposure= ulcssvcon, 
        case    = can_obese, 
        model   = cox, 
        time    = tobeca,
        strata  = agecon cohort intervalp ,
        adj = &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam nSES    , 
        refval  = MIN,  
        lpct    = 1, 
        hpct    = 99,  
        nk  = 4, 
        SELECT  = 1,     /* force all spline variables */      
        testrep  = long,
        outplot = JPEG,
        pwhich  = spline,  
        footer  = NONE,
        GRAPHTIT= NONE,
        plot    = 4,
        e   = T,     
        ci  = 2,     
        displayx= T,  
        klines  = T,   
        plotdec = T,
        PLOTDATA = pool.ulcs.force_min.txt);  /* Output filename */

%LGTPHCURV9(    data    = pool, 
        exposure= hlcssvcon, 
        case    = can_obese, 
        model   = cox, 
        time    = tobeca,
        strata  = agecon cohort intervalp ,
        adj = &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam nSES    , 
        refval  = MIN,  
        lpct    = 1, 
        hpct    = 99,  
        nk  = 4, 
        SELECT  = 1,     /* force all spline variables */      
        testrep  = long,
        outplot = JPEG,
        pwhich  = spline,  
        footer  = NONE,
        GRAPHTIT= NONE,
        plot    = 4,
        e   = T,     
        ci  = 2,     
        displayx= T,  
        klines  = T,   
        plotdec = T,
        PLOTDATA = pool.hlcs.force_min.txt);  /* Output filename */


%LGTPHCURV9(    data    = pool, 
        exposure= lfssvcon, 
        case    = can_obese, 
        model   = cox, 
        time    = tobeca,
        strata  = agecon cohort intervalp ,
        adj = &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam nSES    , 
        refval  = MIN,  
        lpct    = 1, 
        hpct    = 99,  
        nk  = 4, 
        SELECT  = 1,     /* force all spline variables */      
        testrep  = long,
        outplot = JPEG,
        pwhich  = spline,  
        footer  = NONE,
        GRAPHTIT= NONE,
        plot    = 4,
        e   = T,     
        ci  = 2,     
        displayx= T,  
        klines  = T,   
        plotdec = T,
        PLOTDATA = pool.lfs.force_min.txt);  /* Output filename */


%LGTPHCURV9(    data    = pool, 
        exposure= alfssvcon, 
        case    = can_obese, 
        model   = cox, 
        time    = tobeca,
        strata  = agecon cohort intervalp ,
        adj = &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam nSES    , 
        refval  = MIN,  
        lpct    = 1, 
        hpct    = 99,  
        nk  = 4, 
        SELECT  = 1,     /* force all spline variables */      
        testrep  = long,
        outplot = JPEG,
        pwhich  = spline,  
        footer  = NONE,
        GRAPHTIT= NONE,
        plot    = 4,
        e   = T,     
        ci  = 2,     
        displayx= T,  
        klines  = T,   
        plotdec = T,
        PLOTDATA = pool.alfs.force_min.txt);  /* Output filename */


%LGTPHCURV9(    data    = pool, 
        exposure= vlfssvcon, 
        case    = can_obese, 
        model   = cox, 
        time    = tobeca,
        strata  = agecon cohort intervalp ,
        adj = &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam nSES    , 
        refval  = MIN,  
        lpct    = 1, 
        hpct    = 99,  
        nk  = 4, 
        SELECT  = 1,     /* force all spline variables */      
        testrep  = long,
        outplot = JPEG,
        pwhich  = spline,  
        footer  = NONE,
        GRAPHTIT= NONE,
        plot    = 4,
        e   = T,     
        ci  = 2,     
        displayx= T,  
        klines  = T,   
        plotdec = T,
        PLOTDATA = pool.vlfs.force_min.txt);  /* Output filename */


%LGTPHCURV9(    data    = pool, 
        exposure= ulfssvcon, 
        case    = can_obese, 
        model   = cox, 
        time    = tobeca,
        strata  = agecon cohort intervalp ,
        adj = &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam nSES    , 
        refval  = MIN,  
        lpct    = 1, 
        hpct    = 99,  
        nk  = 4, 
        SELECT  = 1,     /* force all spline variables */      
        testrep  = long,
        outplot = JPEG,
        pwhich  = spline,  
        footer  = NONE,
        GRAPHTIT= NONE,
        plot    = 4,
        e   = T,     
        ci  = 2,     
        displayx= T,  
        klines  = T,   
        plotdec = T,
        PLOTDATA = pool.ulfs.force_min.txt);  /* Output filename */


%LGTPHCURV9(    data    = pool, 
        exposure= hlfssvcon, 
        case    = can_obese, 
        model   = cox, 
        time    = tobeca,
        strata  = agecon cohort intervalp ,
        adj = &calorsvq_ &aheimsvq_ &smkstatus_ &alc_ &actc_ &bmib9_ white mv asp cafh pmhever ocever endo2c mam nSES    , 
        refval  = MIN,  
        lpct    = 1, 
        hpct    = 99,  
        nk  = 4, 
        SELECT  = 1,     /* force all spline variables */      
        testrep  = long,
        outplot = JPEG,
        pwhich  = spline,  
        footer  = NONE,
        GRAPHTIT= NONE,
        plot    = 4,
        e   = T,     
        ci  = 2,     
        displayx= T,  
        klines  = T,   
        plotdec = T,
        PLOTDATA = pool.hlfs.force_min.txt);  /* Output filename */





