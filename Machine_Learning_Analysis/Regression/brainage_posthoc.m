cd('/volume/projects/JF_ECNP/Analysis/Analysis/R_models_150224/Regression/')
aal3_HCnorm=readtable('aal3_HCnorm_formanuscript_new.csv')
allatlas_HCnorm=readtable('allatlas_HCnorm_formanuscript_new.csv')
schaefer_HCnorm=readtable('schaefer_HCnorm_formanuscript_new.csv')
hammers_HCnorm=readtable('hammers_HCnorm_formanuscript_new.csv')

aal3_HCrem=readtable('aal3_HCremaining_formanuscript_new.csv');
allatlas_HCrem=readtable('allatlas_HCremaining_formanuscript_new.csv');
schaefer_HCrem=readtable('schaefer_HCremaining_formanuscript_new.csv');
hammers_HCrem=readtable('hammers_HCremaining_formanuscript_new.csv');

aal3_SCZ=readtable('aal3_SCZ_formanuscript_new.csv');
allatlas_SCZ=readtable('allatlas_SCZ_formanuscript_new.csv');
schaefer_SCZ=readtable('schaefer_SCZ_formanuscript_new.csv');
hammers_SCZ=readtable('hammers_SCZ_formanuscript_new.csv');


%%HC NORM
absolute_differences = abs(aal3_HCnorm.predicted_age_corr - aal3_HCnorm.original_age);
MAE=mean(absolute_differences)
[r,p]=corr(aal3_HCnorm.original_age,aal3_HCnorm.predicted_age_corr)
model = fitlm(aal3_HCnorm.predicted_age_corr, aal3_HCnorm.original_age);
model.Rsquared
brainage=mean(aal3_HCnorm.predicted_age_corr-aal3_HCnorm.original_age)
stddev=std(aal3_HCnorm.predicted_age_corr-aal3_HCnorm.original_age)

absolute_differences = abs(allatlas_HCnorm.predicted_age_corr - allatlas_HCnorm.original_age);
MAE=mean(absolute_differences)
[r,p]=corr(allatlas_HCnorm.original_age,allatlas_HCnorm.predicted_age_corr)
model = fitlm(allatlas_HCnorm.predicted_age_corr, allatlas_HCnorm.original_age);
model.Rsquared
brainage=mean(allatlas_HCnorm.predicted_age_corr-allatlas_HCnorm.original_age)
stddev=std(allatlas_HCnorm.predicted_age_corr-allatlas_HCnorm.original_age)


absolute_differences = abs(schaefer_HCnorm.predicted_age_corr - schaefer_HCnorm.original_age);
MAE=mean(absolute_differences)
[r,p]=corr(schaefer_HCnorm.original_age,schaefer_HCnorm.predicted_age_corr)
model = fitlm(schaefer_HCnorm.predicted_age_corr, schaefer_HCnorm.original_age);
model.Rsquared
brainage=mean(schaefer_HCnorm.predicted_age_corr-schaefer_HCnorm.original_age)
stddev=std(schaefer_HCnorm.predicted_age_corr-schaefer_HCnorm.original_age)


absolute_differences = abs(hammers_HCnorm.predicted_age_corr - hammers_HCnorm.original_age);
MAE=mean(absolute_differences)
[r,p]=corr(hammers_HCnorm.original_age,hammers_HCnorm.predicted_age_corr)
model = fitlm(hammers_HCnorm.predicted_age_corr, hammers_HCnorm.original_age);
model.Rsquared
brainage=mean(hammers_HCnorm.predicted_age_corr-hammers_HCnorm.original_age)
stddev=std(hammers_HCnorm.predicted_age_corr-hammers_HCnorm.original_age)


%%HC REM
absolute_differences = abs(aal3_HCrem.predicted_age_corr - aal3_HCrem.original_age);
MAE=mean(absolute_differences)
[r,p]=corr(aal3_HCrem.original_age,aal3_HCrem.predicted_age_corr)
model = fitlm(aal3_HCrem.predicted_age_corr, aal3_HCrem.original_age);
model.Rsquared
brainage=mean(aal3_HCrem.predicted_age_corr-aal3_HCrem.original_age)
stddev=std(aal3_HCrem.predicted_age_corr-aal3_HCrem.original_age)


absolute_differences = abs(allatlas_HCrem.predicted_age_corr - allatlas_HCrem.original_age);
MAE=mean(absolute_differences)
[r,p]=corr(allatlas_HCrem.original_age,allatlas_HCrem.predicted_age_corr)
model = fitlm(allatlas_HCrem.predicted_age_corr, allatlas_HCrem.original_age);
model.Rsquared
brainage=mean(allatlas_HCrem.predicted_age_corr-allatlas_HCrem.original_age)
stddev=std(allatlas_HCrem.predicted_age_corr-allatlas_HCrem.original_age)


absolute_differences = abs(schaefer_HCrem.predicted_age_corr - schaefer_HCrem.original_age);
MAE=mean(absolute_differences)
[r,p]=corr(schaefer_HCrem.original_age,schaefer_HCrem.predicted_age_corr)
model = fitlm(schaefer_HCrem.predicted_age_corr, schaefer_HCrem.original_age);
model.Rsquared
brainage=mean(schaefer_HCrem.predicted_age_corr-schaefer_HCrem.original_age)
stddev=std(schaefer_HCrem.predicted_age_corr-schaefer_HCrem.original_age)


absolute_differences = abs(hammers_HCrem.predicted_age_corr - hammers_HCrem.original_age);
MAE=mean(absolute_differences)
[r,p]=corr(hammers_HCrem.original_age,hammers_HCrem.predicted_age_corr)
model = fitlm(hammers_HCrem.predicted_age_corr, hammers_HCrem.original_age);
model.Rsquared
brainage=mean(hammers_HCrem.predicted_age_corr-hammers_HCrem.original_age)
stddev=std(hammers_HCrem.predicted_age_corr-hammers_HCrem.original_age)


%%SCZ
absolute_differences = abs(aal3_SCZ.predicted_age_corr - aal3_SCZ.original_age);
MAE=mean(absolute_differences)
[r,p]=corr(aal3_SCZ.original_age,aal3_SCZ.predicted_age_corr)
model = fitlm(aal3_SCZ.predicted_age_corr, aal3_SCZ.original_age);
model.Rsquared
brainage=mean(aal3_SCZ.predicted_age_corr-aal3_SCZ.original_age)
stddev=std(aal3_SCZ.predicted_age_corr-aal3_SCZ.original_age)


absolute_differences = abs(allatlas_SCZ.predicted_age_corr - allatlas_SCZ.original_age);
MAE=mean(absolute_differences)
[r,p]=corr(allatlas_SCZ.original_age,allatlas_SCZ.predicted_age_corr)
model = fitlm(allatlas_SCZ.predicted_age_corr, allatlas_SCZ.original_age);
model.Rsquared
brainage=mean(allatlas_SCZ.predicted_age_corr-allatlas_SCZ.original_age)
stddev=std(allatlas_SCZ.predicted_age_corr-allatlas_SCZ.original_age)


absolute_differences = abs(schaefer_SCZ.predicted_age_corr - schaefer_SCZ.original_age);
MAE=mean(absolute_differences)
[r,p]=corr(schaefer_SCZ.original_age,schaefer_SCZ.predicted_age_corr)
model = fitlm(schaefer_SCZ.predicted_age_corr, schaefer_SCZ.original_age);
model.Rsquared
brainage=mean(schaefer_SCZ.predicted_age_corr-schaefer_SCZ.original_age)
stddev=std(schaefer_SCZ.predicted_age_corr-schaefer_SCZ.original_age)


absolute_differences = abs(hammers_SCZ.predicted_age_corr - hammers_SCZ.original_age);
MAE=mean(absolute_differences)
[r,p]=corr(hammers_SCZ.original_age,hammers_SCZ.predicted_age_corr)
model = fitlm(hammers_SCZ.predicted_age_corr, hammers_SCZ.original_age);
model.Rsquared
brainage=mean(hammers_SCZ.predicted_age_corr-hammers_SCZ.original_age)
stddev=std(hammers_SCZ.predicted_age_corr-hammers_SCZ.original_age)
BrainAGE_HCremaining_hammers_corrected=hammers_HCrem.BrainAGE_corr


[h,p,ci,stats]=ttest2(BrainAGE_SCZ_hammers_corrected,BrainAGE_HCremaining_hammers_corrected)



data_neuroimg_regr_unif=readtable('data_neuroimg_regr_unif.csv');
data_neuroimg_regr_rem=readtable('data_neuroimg_regr_rem.csv');


female_unif_count = sum(data_neuroimg_regr_unif.SEX == 1);
female_rem_count = sum(data_neuroimg_regr_rem.SEX == 1);

% Concatenate the SEX variables from both datasets
sex_concatenated = [data_neuroimg_regr_unif.SEX; data_neuroimg_regr_rem.SEX];

% Create a grouping variable indicating the dataset
grouping_variable = [3*ones(size(data_neuroimg_regr_unif.SEX)); 2*ones(size(data_neuroimg_regr_rem.SEX))];

% Create the contingency table
[tbl, chi2,p] = crosstab(sex_concatenated, grouping_variable);

% Calculate degrees of freedom
[n_rows, n_cols] = size(tbl);
degrees_of_freedom = (n_rows - 1) * (n_cols - 1);

unique(data_neuroimg_regr_unif_muc.SITE_ID)
data_neuroimg_regr_unif_muc=data_neuroimg_regr_unif.AGE(data_neuroimg_regr_unif.SITE_ID==1)
data_neuroimg_regr_rem_muc=data_neuroimg_regr_unif.AGE(data_neuroimg_regr_rem.SITE_ID==1)
data_neuroimg_regr_unif_osl=data_neuroimg_regr_unif.AGE(data_neuroimg_regr_unif.SITE_ID==9)
data_neuroimg_regr_rem_osl=data_neuroimg_regr_unif.AGE(data_neuroimg_regr_rem.SITE_ID==9)


