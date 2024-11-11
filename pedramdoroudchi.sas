/* Create dataset in SAS */
data dat;
    /* Initialize counters */
    /* Total number of observations for each variable */
    do i = 1 to 300; /* Total of 300 observations */
        
        /* Populate FCD */
        if i <= 145 then FCD = 'Positive'; 
        else if i <= 148 then FCD = 'Negative'; 
        else if i = 149 then FCD = 'Positive'; 
        else if i = 150 then FCD = 'Negative'; 
        else if i <= 155 then FCD = 'Positive'; 
        else if i = 156 then FCD = 'Negative'; 
        else if i <= 158 then FCD = 'Positive'; 
        else FCD = 'Negative'; /* Remaining 142 Negative */
        
        /* Populate CCD1 */
        if i <= 150 then CCD1 = 'Positive'; 
        else CCD1 = 'Negative'; 
        
        /* Populate CCD2 */
        if i <= 148 then CCD2 = 'Positive'; 
        else if i <= 150 then CCD2 = 'Negative'; 
        else if i <= 156 then CCD2 = 'Positive'; 
        else CCD2 = 'Negative'; /* Remaining 144 Negative */
        
        output; /* Output the current observation */
    end;

    /* Keep the variables in the dataset */
    keep FCD CCD1 CCD2;
run;


/* Create a confusion matrix for CCD1 Positive */
proc sql;
    create table confusion_matrix_ccd_pos as
    select 
        FCD,
        sum(CCD2 = 'Positive') as CCD1_plus_CCD2_plus,
        sum(CCD2 = 'Negative') as CCD1_plus_CCD2_minus,
        count(*) as Total_CCD1_plus
    from dat
    where CCD1 = 'Positive'
    group by FCD
    order by FCD desc;
quit;

/* Create total for CCD1 Positive */
proc sql;
    create table confusion_matrix_ccd_pos_tot as
    select 
        'Total' as FCD,
        sum(CCD2 = 'Positive') as CCD1_plus_CCD2_plus,
        sum(CCD2 = 'Negative') as CCD1_plus_CCD2_minus,
        count(*) as Total_CCD1_plus
    from dat
    where CCD1 = 'Positive';
quit;

/* Combine CCD1 Positive and its total */
data confusion_matrix_ccd_pos_w_tot;
    set confusion_matrix_ccd_pos confusion_matrix_ccd_pos_tot;
run;

/* Create a confusion matrix for CCD1 Negative */
proc sql;
    create table confusion_matrix_ccd_neg as
    select 
        FCD,
        sum(CCD2 = 'Positive') as CCD1_minus_CCD2_plus,
        sum(CCD2 = 'Negative') as CCD1_minus_CCD2_minus,
        count(*) as Total_CCD1_minus
    from dat
    where CCD1 = 'Negative'
    group by FCD
    order by FCD desc;
quit;

/* Create total for CCD1 Negative */
proc sql;
    create table confusion_matrix_ccd_neg_tot as
    select 
        'Total' as FCD,
        sum(CCD2 = 'Positive') as CCD1_minus_CCD2_plus,
        sum(CCD2 = 'Negative') as CCD1_minus_CCD2_minus,
        count(*) as Total_CCD1_minus
    from dat
    where CCD1 = 'Negative';
quit;

/* Combine CCD1 Negative and its total */
data confusion_matrix_ccd_neg_w_tot;
    set confusion_matrix_ccd_neg confusion_matrix_ccd_neg_tot;
run;

/* Merge the confusion matrices for CCD1 Positive and Negative */
proc sql;
    create table confusion_matrix as
    select 
        coalesce(a.FCD, b.FCD) as FCD,
        a.CCD1_plus_CCD2_plus,
        a.CCD1_plus_CCD2_minus,
        a.Total_CCD1_plus,
        b.CCD1_minus_CCD2_plus,
        b.CCD1_minus_CCD2_minus,
        b.Total_CCD1_minus
    from confusion_matrix_ccd_pos_w_tot as a
    full join confusion_matrix_ccd_neg_w_tot as b
    on a.FCD = b.FCD;
quit;

/* Final formatting of FCD */
data confusion_matrix;
    set confusion_matrix;
    length FCD $ 10;
    if FCD = 'Positive' then FCD = 'FCD+';
    else if FCD = 'Negative' then FCD = 'FCD-';
    else if FCD = 'Total' then FCD = 'Total';
run;

/* Final confusion matrix name formatting */
proc sql;
    create table confusion_matrix_final as
    select 
        FCD,
        CCD1_plus_CCD2_plus as 'CCD1+, CCD2+'n,
        CCD1_plus_CCD2_minus as 'CCD1+, CCD2-'n,
        Total_CCD1_plus as 'CCD1+ Total'n,
        CCD1_minus_CCD2_plus as 'CCD1-, CCD2+'n,
        CCD1_minus_CCD2_minus as 'CCD1-, CCD2-'n,
        Total_CCD1_minus as 'CCD1- Total'n
    from confusion_matrix
    order by case 
        when FCD = 'FCD+' then 1
        when FCD = 'FCD-' then 2
        when FCD = 'Total' then 3
    end;
quit;


/* Calculate zPPA1_hat */
proc sql;
    create table zPPA1_hat as
    select 
        ( 
            (select CCD1_plus_CCD2_plus from confusion_matrix where FCD = 'Total') -
            (select Total_CCD1_plus from confusion_matrix where FCD = 'FCD+')
        ) / 
        (select Total_CCD1_plus from confusion_matrix where FCD = 'Total') as zPPA1_hat
    from 
        confusion_matrix
    where 
        FCD = 'Total';
quit;


/* Set up bootstrap parameters */
%let n_samples = 1000; /* Number of bootstrap samples */
%let sample_size = 150; /* Sample size for each bootstrap */

/* Create a filtered dataset for positive CCD1 observations */
data positive_CCD1;
    set dat;
    if CCD1 = 'Positive'; /* Keep only observations where CCD1 is Positive */
run;

/* Create a dataset for bootstrap samples */
proc surveyselect data=positive_CCD1 out=bootstrap_samples
    method=urs /* Use URS method for bootstrap sampling */
    n=&sample_size /* Sample size for each bootstrap sample */
    reps=&n_samples; /* Number of bootstrap samples */
run;

/* Calculate zPPA1_hats using PROC SQL, incorporating NumberHits */
proc sql;
    create table zPPA1_hats as
    select 
        sum(case when CCD2 = 'Positive' then NumberHits else 0 end) as count_CCD2_Pos,
        sum(case when FCD = 'Positive' then NumberHits else 0 end) as count_FCD_Pos,
        (&sample_size) as sample_size,
        (calculated count_CCD2_Pos - calculated count_FCD_Pos) / &sample_size as zPPA1_hat
    from bootstrap_samples
    group by Replicate; /* Group by each bootstrap sample */
quit;

/* Calculate the 95% confidence intervals */
proc univariate data=zPPA1_hats noprint;
    var zPPA1_hat; /* Use the zPPA1_hat values */
    output out=zPPA1_ci pctlpre=CI_ pctlpts=2.5 97.5; /* Get the 2.5 and 97.5 percentiles */
run;


/* Set the prevalence value */
%let prev = 0.3;

/* Calculate zPPA2_hat */
proc sql;
    create table zPPA2_hat as
    select 
        (
            (&prev * (select CCD1_plus_CCD2_plus from confusion_matrix where FCD = 'Total') * 
                (select Total_CCD1_minus from confusion_matrix where FCD = 'Total')) - 
            (&prev * (select CCD1_plus_CCD2_plus from confusion_matrix where FCD = 'FCD+') * 
                (select Total_CCD1_minus from confusion_matrix where FCD = 'Total')) - 
            ((1 - &prev) * (select CCD1_minus_CCD2_plus from confusion_matrix where FCD = 'FCD+') * 
                (select Total_CCD1_plus from confusion_matrix where FCD = 'Total'))
        ) / 
        (
            (&prev * (select CCD1_plus_CCD2_plus from confusion_matrix where FCD = 'Total') * 
                (select Total_CCD1_minus from confusion_matrix where FCD = 'Total')) + 
            ((1 - &prev) * (select CCD1_minus_CCD2_plus from confusion_matrix where FCD = 'Total') * 
                (select Total_CCD1_plus from confusion_matrix where FCD = 'Total'))
        ) as zPPA2_hat
    from 
        confusion_matrix
    where 
        FCD = 'Total';
quit;

/* Set up bootstrap parameters */
%let n_samples = 1000; /* Number of bootstrap samples */
%let sample_size = 300; /* Sample size for each bootstrap */

/* Create a dataset for bootstrap samples */
proc surveyselect data=dat out=bootstrap_samples
    method=urs /* Use URS method for bootstrap sampling */
    n=&sample_size /* Sample size for each bootstrap sample */
    reps=&n_samples; /* Number of bootstrap samples */
run;

/* Calculate zPPA2_hats using PROC SQL, incorporating NumberHits */
proc sql;
    create table zPPA2_hats as
    select 
        sum(case when CCD1 = 'Positive' and CCD2 = 'Positive' then NumberHits else 0 end) as count_CCD1_Pos_CCD2_Pos,
        sum(case when CCD1 = 'Negative' and CCD2 = 'Positive' then NumberHits else 0 end) as count_CCD1_Neg_CCD2_Pos,
        sum(case when FCD = 'Positive' and CCD1 = 'Positive' and CCD2 = 'Positive' then NumberHits else 0 end) as count_FCD_Pos_CCD1_Pos_CCD2_Pos,
        sum(case when FCD = 'Positive' and CCD1 = 'Negative' and CCD2 = 'Positive' then NumberHits else 0 end) as count_FCD_Pos_CCD1_Neg_CCD2_Pos,
        sum(case when CCD1 = 'Negative' then NumberHits else 0 end) as count_CCD1_Neg,
        sum(case when CCD1 = 'Positive' then NumberHits else 0 end) as count_CCD1_Pos,
        (&prev * calculated count_CCD1_Pos_CCD2_Pos * calculated count_CCD1_Neg - 
         &prev * calculated count_FCD_Pos_CCD1_Pos_CCD2_Pos * calculated count_CCD1_Neg - 
         (1 - &prev) * calculated count_FCD_Pos_CCD1_Neg_CCD2_Pos * calculated count_CCD1_Pos) /
        (&prev * calculated count_CCD1_Pos_CCD2_Pos * calculated count_CCD1_Neg + 
         (1 - &prev) * calculated count_CCD1_Neg_CCD2_Pos * calculated count_CCD1_Pos) as zPPA2_hat
    from bootstrap_samples
    group by Replicate; /* Group by each bootstrap sample */
quit;

/* Calculate the 95% confidence intervals */
proc univariate data=zPPA2_hats noprint;
    var zPPA2_hat; /* Use the zPPA2_hat values */
    output out=zPPA2_ci pctlpre=CI_ pctlpts=2.5 97.5; /* Get the 2.5 and 97.5 percentiles */
run;


/* Calculate zNPA1_hat */
proc sql;
    create table zNPA1_hat as
    select 
        ( 
            (select CCD1_minus_CCD2_minus from confusion_matrix where FCD = 'Total') -
            (select Total_CCD1_minus from confusion_matrix where FCD = 'FCD-')
        ) / 
        (select Total_CCD1_minus from confusion_matrix where FCD = 'Total') as zNPA1_hat
    from 
        confusion_matrix
    where 
        FCD = 'Total';
quit;


/* Set up bootstrap parameters */
%let n_samples = 1000; /* Number of bootstrap samples */
%let sample_size = 150; /* Sample size for each bootstrap */

/* Create a filtered dataset for negative CCD1 observations */
data negative_CCD1;
    set dat;
    if CCD1 = 'Negative'; /* Keep only observations where CCD1 is Negative */
run;

/* Create a dataset for bootstrap samples */
proc surveyselect data=negative_CCD1 out=bootstrap_samples
    method=urs /* Use URS method for bootstrap sampling */
    n=&sample_size /* Sample size for each bootstrap sample */
    reps=&n_samples; /* Number of bootstrap samples */
run;

/* Calculate zNPA1_hats using PROC SQL, incorporating NumberHits */
proc sql;
    create table zNPA1_hats as
    select 
        sum(case when CCD2 = 'Negative' then NumberHits else 0 end) as count_CCD2_Neg,
        sum(case when FCD = 'Negative' then NumberHits else 0 end) as count_FCD_Neg,
        (&sample_size) as sample_size,
        (calculated count_CCD2_Neg - calculated count_FCD_Neg) / &sample_size as zNPA1_hat
    from bootstrap_samples
    group by Replicate; /* Group by each bootstrap sample */
quit;

/* Calculate the 95% confidence intervals */
proc univariate data=zNPA1_hats noprint;
    var zNPA1_hat; /* Use the zNPA1_hat values */
    output out=zNPA1_ci pctlpre=CI_ pctlpts=2.5 97.5; /* Get the 2.5 and 97.5 percentiles */
run;


/* Calculate zNPA2_hat */
proc sql;
    create table zNPA2_hat as
    select 
        (
            ((1 - &prev) * (select CCD1_minus_CCD2_minus from confusion_matrix where FCD = 'Total') * 
                (select Total_CCD1_plus from confusion_matrix where FCD = 'Total')) - 
            ((1 - &prev) * (select CCD1_minus_CCD2_minus from confusion_matrix where FCD = 'FCD-') * 
                (select Total_CCD1_plus from confusion_matrix where FCD = 'Total')) - 
            (&prev * (select CCD1_plus_CCD2_minus from confusion_matrix where FCD = 'FCD-') * 
                (select Total_CCD1_minus from confusion_matrix where FCD = 'Total'))
        ) / 
        (
            ((1 - &prev) * (select CCD1_minus_CCD2_minus from confusion_matrix where FCD = 'Total') * 
                (select Total_CCD1_plus from confusion_matrix where FCD = 'Total')) + 
            (&prev * (select CCD1_plus_CCD2_minus from confusion_matrix where FCD = 'Total') * 
                (select Total_CCD1_minus from confusion_matrix where FCD = 'Total'))
        ) as zNPA2_hat
    from 
        confusion_matrix
    where 
        FCD = 'Total';
quit;

/* Set up bootstrap parameters */
%let n_samples = 1000; /* Number of bootstrap samples */
%let sample_size = 300; /* Sample size for each bootstrap */

/* Create a dataset for bootstrap samples */
proc surveyselect data=dat out=bootstrap_samples
    method=urs /* Use URS method for bootstrap sampling */
    n=&sample_size /* Sample size for each bootstrap sample */
    reps=&n_samples; /* Number of bootstrap samples */
run;

/* Calculate zNPA2_hats using PROC SQL, incorporating NumberHits */
proc sql;
    create table zNPA2_hats as
    select 
        sum(case when CCD1 = 'Negative' and CCD2 = 'Negative' then NumberHits else 0 end) as count_CCD1_Neg_CCD2_Neg,
        sum(case when CCD1 = 'Positive' and CCD2 = 'Negative' then NumberHits else 0 end) as count_CCD1_Pos_CCD2_Neg,
        sum(case when FCD = 'Negative' and CCD1 = 'Negative' and CCD2 = 'Negative' then NumberHits else 0 end) as count_FCD_Neg_CCD1_Neg_CCD2_Neg,
        sum(case when FCD = 'Negative' and CCD1 = 'Positive' and CCD2 = 'Negative' then NumberHits else 0 end) as count_FCD_Neg_CCD1_Pos_CCD2_Neg,
        sum(case when CCD1 = 'Negative' then NumberHits else 0 end) as count_CCD1_Neg,
        sum(case when CCD1 = 'Positive' then NumberHits else 0 end) as count_CCD1_Pos,
        ((1 - &prev) * calculated count_CCD1_Neg_CCD2_Neg * calculated count_CCD1_Pos - 
         (1 - &prev) * calculated count_FCD_Neg_CCD1_Neg_CCD2_Neg * calculated count_CCD1_Pos - 
         &prev * calculated count_FCD_Neg_CCD1_Pos_CCD2_Neg * calculated count_CCD1_Neg) /
        ((1 - &prev) * calculated count_CCD1_Neg_CCD2_Neg * calculated count_CCD1_Pos + 
         &prev * calculated count_CCD1_Pos_CCD2_Neg * calculated count_CCD1_Neg) as zNPA2_hat
    from bootstrap_samples
    group by Replicate; /* Group by each bootstrap sample */
quit;

/* Calculate the 95% confidence intervals */
proc univariate data=zNPA2_hats noprint;
    var zNPA2_hat; /* Use the zNPA2_hat values */
    output out=zNPA2_ci pctlpre=CI_ pctlpts=2.5 97.5; /* Get the 2.5 and 97.5 percentiles */
run;


/* Create the analysis table */
data analysis_tbl;
    /* Set the prevalence value */
    Prevalence = '30.0%';

    /* ZPPA1 Estimate and CI */
    set zPPA1_hat; /* Assuming zPPA1_hat contains the estimated value */
    ZPPA1_Estimate = put(zPPA1_hat, 8.3); /* Format to three decimal places */
    set zPPA1_ci; /* Assuming zPPA1_ci contains the CI */
    ZPPA1_95_CI = catx(' ', '(', put(CI_2_5, 8.3), ',', put(CI_97_5, 8.3), ')');

    /* ZPPA2 Estimate and CI */
    set zPPA2_hat; /* Assuming zPPA2_hat contains the estimated value */
    ZPPA2_Estimate = put(zPPA2_hat, 8.3); /* Format to three decimal places */
    set zPPA2_ci; /* Assuming zPPA2_ci contains the CI */
    ZPPA2_95_CI = catx(' ', '(', put(CI_2_5, 8.3), ',', put(CI_97_5, 8.3), ')');

    /* ZNPA1 Estimate and CI */
    set zNPA1_hat; /* Assuming zNPA1_hat contains the estimated value */
    ZNPA1_Estimate = put(zNPA1_hat, 8.3); /* Format to three decimal places */
    set zNPA1_ci; /* Assuming zNPA1_ci contains the CI */
    ZNPA1_95_CI = catx(' ', '(', put(CI_2_5, 8.3), ',', put(CI_97_5, 8.3), ')');

    /* ZNPA2 Estimate and CI */
    set zNPA2_hat; /* Assuming zNPA2_hat contains the estimated value */
    ZNPA2_Estimate = put(zNPA2_hat, 8.3); /* Format to three decimal places */
    set zNPA2_ci; /* Assuming zNPA2_ci contains the CI */
    ZNPA2_95_CI = catx(' ', '(', put(CI_2_5, 8.3), ',', put(CI_97_5, 8.3), ')');

    /* Output the observation */
    output; 
run;

proc sql;
    create table analysis_tbl_final as
    select 
        Prevalence as Prev,
        ZPPA1_Estimate as zPPA1,
        ZPPA1_95_CI as 'ZPPA1 CI'n,
        ZPPA2_Estimate as zPPA2,
        ZPPA2_95_CI as 'ZPPA2 CI'n,
        ZNPA1_Estimate as zNPA1,
        ZNPA1_95_CI as 'ZNPA1 CI'n,
        ZNPA2_Estimate as zNPA2,
        ZNPA2_95_CI as 'ZNPA2 CI'n
    from analysis_tbl
quit;


/* Define output path */
ods word file="/home/u63705201/analysis_tables_SAS_PD.docx";

/* Output confusion matrix */
proc print data=confusion_matrix_final noobs;
    title "Confusion Matrix";
run;

/* Output analysis table */
proc print data=analysis_tbl_final noobs;
    title "Analysis Table";
run;

/* Close ODS destination */
ods word close;

