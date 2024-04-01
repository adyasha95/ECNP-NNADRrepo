HC_norm_allatlas=readtable('/volume/projects/JF_ECNP/Analysis/Analysis/R_models_150224/Regression/allatlas_HCnorm_new.csv');
HC_rem_allatlas=readtable('/volume/projects/JF_ECNP/Analysis/Analysis/R_models_150224/Regression/allatlas_HCrem_new.csv');
SCZ_allatlas=readtable('/volume/projects/JF_ECNP/Analysis/Analysis/R_models_150224/Regression/allatlas_scz_new.csv');

HC_norm_schaefer=readtable('/volume/projects/JF_ECNP/Analysis/Analysis/R_models_150224/Regression/schaefer200_gmv_HCnorm_new.csv');
HC_rem_schaefer=readtable('/volume/projects/JF_ECNP/Analysis/Analysis/R_models_150224/Regression/schaefer200_gmv_HCrem_new.csv');
SCZ_schaefer=readtable('/volume/projects/JF_ECNP/Analysis/Analysis/R_models_150224/Regression/schaefer200_gmv_scz_new.csv');

HC_norm_hammers=readtable('/volume/projects/JF_ECNP/Analysis/Analysis/R_models_150224/Regression/hammers_gmv_HCnorm_new.csv');
HC_rem_hammers=readtable('/volume/projects/JF_ECNP/Analysis/Analysis/R_models_150224/Regression/hammers_gmv_HCrem_new.csv');
SCZ_hammers=readtable('/volume/projects/JF_ECNP/Analysis/Analysis/R_models_150224/Regression/hammers_gmv_scz_new.csv');

HC_norm_aal3=readtable('/volume/projects/JF_ECNP/Analysis/Analysis/R_models_150224/Regression/aal3_HCnorm_new.csv');
HC_rem_aal3=readtable('/volume/projects/JF_ECNP/Analysis/Analysis/R_models_150224/Regression/aal3_HCrem_new.csv');
SCZ_aal3=readtable('/volume/projects/JF_ECNP/Analysis/Analysis/R_models_150224/Regression/aal3_scz_new.csv');

out_dir='/volume/projects/JF_ECNP/Analysis/Analysis/R_models_150224/Regression/'

%% all atlas together

AGE_HC_bal_allatlas = HC_norm_allatlas.original_ages;
AGE_SCZ_allatlas=SCZ_allatlas.original_age_scz;
AGE_HC_remaining_alltlas=HC_rem_allatlas.original_age_HC_rem;
predicted_AGE_HC_allatlas = HC_norm_allatlas.average_predictions;
predicted_AGE_Scz_allatlas = SCZ_allatlas.average_predictions_patients_scz;
predicted_AGE_HCremaining_allatlas = HC_rem_allatlas.average_predictions_HC_rem;

BrainAGE_HC_allatlas=predicted_AGE_HC_allatlas-AGE_HC_bal_allatlas;
meanBrainAGE_allatlas=mean(BrainAGE_HC_allatlas)
BrainAGE_HCrem_allatlas = predicted_AGE_HCremaining_allatlas-AGE_HC_remaining_alltlas;
meanBrainAGE_HC_apply=mean(BrainAGE_HCrem_allatlas)
BrainAGE_Scz_allatlas = predicted_AGE_Scz_allatlas-AGE_SCZ_allatlas;
meanBrainAGE_Scz=mean(BrainAGE_Scz_allatlas)

IN=struct;
Y= predicted_AGE_HC_allatlas;
X= AGE_HC_bal_allatlas;
IN.TrPred=predicted_AGE_HC_allatlas;
IN.TrObs=AGE_HC_bal_allatlas;
 [ sY, Y, IN ] = nk_CorrectPredTails(Y, X, IN)
 predicted_age_HC_allatlas_corrected=sY;
 BrainAGE_HC_allatlas_corrected=predicted_age_HC_allatlas_corrected - X;
 meanBrainAGE_HC_allatlas_corrected=mean(BrainAGE_HC_allatlas_corrected)
 stdBrainAGE_HC_allatlas_corrected=std(BrainAGE_HC_allatlas_corrected)
 scatter(predicted_age_HC_allatlas_corrected,AGE_HC_bal_allatlas)

table_HC_norm_allatlas=table(AGE_HC_bal_allatlas,predicted_AGE_HC_allatlas,BrainAGE_HC_allatlas,predicted_age_HC_allatlas_corrected,BrainAGE_HC_allatlas_corrected,'VariableNames',{'original_age','predicted_age','BrainAGE','predicted_age_corr','BrainAGE_corr'})
writetable(table_HC_norm_allatlas,fullfile(out_dir,'allatlas_HCnorm_formanuscript_new.csv'));

Y= predicted_AGE_HCremaining_allatlas;
X= AGE_HC_remaining_alltlas
IN.TrPred=predicted_AGE_HC_allatlas;
IN.TrObs=AGE_HC_bal_allatlas;

 [ sY, Y, IN ] = nk_CorrectPredTails(Y, X, IN)
 predicted_age_HCremaining_allatlas_corrected=sY;
 BrainAGE_HCremaining_allatlas_corrected=predicted_age_HCremaining_allatlas_corrected - X;
 meanBrainAGE_HCremaining_allatlas_corrected=mean(BrainAGE_HCremaining_allatlas_corrected)
  stdBrainAGE_HCremaining_allatlas_corrected=std(BrainAGE_HCremaining_allatlas_corrected)

table_HC_remaining_allatlas=table(AGE_HC_remaining_alltlas,predicted_AGE_HCremaining_allatlas,BrainAGE_HCrem_allatlas,predicted_age_HCremaining_allatlas_corrected,BrainAGE_HCremaining_allatlas_corrected,'VariableNames',{'original_age','predicted_age','BrainAGE','predicted_age_corr','BrainAGE_corr'})
scatter(predicted_age_HCremaining_allatlas_corrected,AGE_HC_remaining_alltlas)
writetable(table_HC_remaining_allatlas,fullfile(out_dir,'allatlas_HCremaining_formanuscript_new.csv'));

 
Y= predicted_AGE_Scz_allatlas;
X= AGE_SCZ_allatlas;
IN.TrPred=predicted_AGE_HC_allatlas;
IN.TrObs=AGE_HC_bal_allatlas;

 [ sY, Y, IN ] = nk_CorrectPredTails(Y, X, IN)
 predicted_age_SCZ_allatlas_corrected=sY;
 BrainAGE_SCZ_allatlas_corrected=predicted_age_SCZ_allatlas_corrected - X;
 meanBrainAGE_Scz_all_atlas_corrected=mean(BrainAGE_SCZ_allatlas_corrected)
  stdBrainAGE_Scz_all_atlas_corrected=std(BrainAGE_SCZ_allatlas_corrected)

table_SCZ_allatlas=table(AGE_SCZ_allatlas,predicted_AGE_Scz_allatlas,BrainAGE_Scz_allatlas,predicted_age_SCZ_allatlas_corrected,BrainAGE_SCZ_allatlas_corrected,'VariableNames',{'original_age','predicted_age','BrainAGE','predicted_age_corr','BrainAGE_corr'})
scatter(predicted_AGE_Scz_allatlas,AGE_SCZ_allatlas)
writetable(table_SCZ_allatlas,fullfile(out_dir,'allatlas_SCZ_formanuscript_new.csv'));

%% Schaefer 200 gmv
AGE_HC_bal_schaefer = HC_norm_schaefer.original_ages;
AGE_SCZ_schaefer=SCZ_schaefer.original_age_scz;
AGE_HC_remaining_alltlas=HC_rem_schaefer.original_age_HC_rem;
predicted_AGE_HC_schaefer = HC_norm_schaefer.average_predictions;
predicted_AGE_Scz_schaefer = SCZ_schaefer.average_predictions_patients_scz;
predicted_AGE_HCremaining_schaefer = HC_rem_schaefer.average_predictions_HC_rem;

BrainAGE_HC_schaefer=predicted_AGE_HC_schaefer-AGE_HC_bal_schaefer;
meanBrainAGE_schaefer=mean(BrainAGE_HC_schaefer)
BrainAGE_HCrem_schaefer = predicted_AGE_HCremaining_schaefer-AGE_HC_remaining_alltlas;
meanBrainAGE_HC_apply=mean(BrainAGE_HCrem_schaefer)
BrainAGE_Scz_schaefer = predicted_AGE_Scz_schaefer-AGE_SCZ_schaefer;
meanBrainAGE_Scz=mean(BrainAGE_Scz_schaefer)

IN=struct;
Y= predicted_AGE_HC_schaefer;
X= AGE_HC_bal_schaefer;
IN.TrPred=predicted_AGE_HC_schaefer;
IN.TrObs=AGE_HC_bal_schaefer;
 [ sY, Y, IN ] = nk_CorrectPredTails(Y, X, IN)
 predicted_age_HC_schaefer_corrected=sY;
 BrainAGE_HC_schaefer_corrected=predicted_age_HC_schaefer_corrected - X;
 meanBrainAGE_HC_schaefer_corrected=mean(BrainAGE_HC_schaefer_corrected)
  stdBrainAGE_HC_schaefer_corrected=std(BrainAGE_HC_schaefer_corrected)

scatter(predicted_age_HC_schaefer_corrected,AGE_HC_bal_schaefer)

table_HC_norm_schaefer=table(AGE_HC_bal_schaefer,predicted_AGE_HC_schaefer,BrainAGE_HC_schaefer,predicted_age_HC_schaefer_corrected,BrainAGE_HC_schaefer_corrected,'VariableNames',{'original_age','predicted_age','BrainAGE','predicted_age_corr','BrainAGE_corr'})
writetable(table_HC_norm_schaefer,fullfile(out_dir,'schaefer_HCnorm_formanuscript_new.csv'));

Y= predicted_AGE_HCremaining_schaefer;
X= AGE_HC_remaining_alltlas
IN.TrPred=predicted_AGE_HC_schaefer;
IN.TrObs=AGE_HC_bal_schaefer;

 [ sY, Y, IN ] = nk_CorrectPredTails(Y, X, IN)
 predicted_age_HCremaining_schaefer_corrected=sY;
 BrainAGE_HCremaining_schaefer_corrected=predicted_age_HCremaining_schaefer_corrected - X;
 meanBrainAGE_HCremaining_schaefer_corrected=mean(BrainAGE_HCremaining_schaefer_corrected)
 stdBrainAGE_HCremaining_schaefer_corrected=std(BrainAGE_HCremaining_schaefer_corrected)

 table_HC_remaining_schaefer=table(AGE_HC_remaining_alltlas,predicted_AGE_HCremaining_schaefer,BrainAGE_HCrem_schaefer,predicted_age_HCremaining_schaefer_corrected,BrainAGE_HCremaining_schaefer_corrected,'VariableNames',{'original_age','predicted_age','BrainAGE','predicted_age_corr','BrainAGE_corr'})
scatter(predicted_age_HCremaining_schaefer_corrected,AGE_HC_remaining_alltlas)
writetable(table_HC_remaining_schaefer,fullfile(out_dir,'schaefer_HCremaining_formanuscript_new.csv'));

 
Y= predicted_AGE_Scz_schaefer;
X= AGE_SCZ_schaefer;
IN.TrPred=predicted_AGE_HC_schaefer;
IN.TrObs=AGE_HC_bal_schaefer;

 [ sY, Y, IN ] = nk_CorrectPredTails(Y, X, IN)
 predicted_age_SCZ_schaefer_corrected=sY;
 BrainAGE_SCZ_schaefer_corrected=predicted_age_SCZ_schaefer_corrected - X;
 meanBrainAGE_Scz_schaefer_corrected=mean(BrainAGE_SCZ_schaefer_corrected)
 stdBrainAGE_Scz_schaefer_corrected=std(BrainAGE_SCZ_schaefer_corrected)

 table_SCZ_schaefer=table(AGE_SCZ_schaefer,predicted_AGE_Scz_schaefer,BrainAGE_Scz_schaefer,predicted_age_SCZ_schaefer_corrected,BrainAGE_SCZ_schaefer_corrected,'VariableNames',{'original_age','predicted_age','BrainAGE','predicted_age_corr','BrainAGE_corr'})
scatter(predicted_AGE_Scz_schaefer,AGE_SCZ_schaefer)
writetable(table_SCZ_schaefer,fullfile(out_dir,'schaefer_SCZ_formanuscript_new.csv'));

%% aal3

AGE_HC_bal_aal3 = HC_norm_aal3.original_ages;
AGE_SCZ_aal3=SCZ_aal3.original_age_scz;
AGE_HC_remaining_alltlas=HC_rem_aal3.original_age_HC_rem;
predicted_AGE_HC_aal3 = HC_norm_aal3.average_predictions;
predicted_AGE_Scz_aal3 = SCZ_aal3.average_predictions_patients_scz;
predicted_AGE_HCremaining_aal3 = HC_rem_aal3.average_predictions_HC_rem;

BrainAGE_HC_aal3=predicted_AGE_HC_aal3-AGE_HC_bal_aal3;
meanBrainAGE_aal3=mean(BrainAGE_HC_aal3)
BrainAGE_HCrem_aal3 = predicted_AGE_HCremaining_aal3-AGE_HC_remaining_alltlas;
meanBrainAGE_HC_apply=mean(BrainAGE_HCrem_aal3)
BrainAGE_Scz_aal3 = predicted_AGE_Scz_aal3-AGE_SCZ_aal3;
meanBrainAGE_Scz=mean(BrainAGE_Scz_aal3)

IN=struct;
Y= predicted_AGE_HC_aal3;
X= AGE_HC_bal_aal3;
IN.TrPred=predicted_AGE_HC_aal3;
IN.TrObs=AGE_HC_bal_aal3;
 [ sY, Y, IN ] = nk_CorrectPredTails(Y, X, IN)
 predicted_age_HC_aal3_corrected=sY;
 BrainAGE_HC_aal3_corrected=predicted_age_HC_aal3_corrected - X;
 meanBrainAGE_HC_aal3_corrected=mean(BrainAGE_HC_aal3_corrected)
  stdBrainAGE_HC_aal3_corrected=std(BrainAGE_HC_aal3_corrected)

scatter(predicted_age_HC_aal3_corrected,AGE_HC_bal_aal3)

table_HC_norm_aal3=table(AGE_HC_bal_aal3,predicted_AGE_HC_aal3,BrainAGE_HC_aal3,predicted_age_HC_aal3_corrected,BrainAGE_HC_aal3_corrected,'VariableNames',{'original_age','predicted_age','BrainAGE','predicted_age_corr','BrainAGE_corr'})
writetable(table_HC_norm_aal3,fullfile(out_dir,'aal3_HCnorm_formanuscript_new.csv'));

Y= predicted_AGE_HCremaining_aal3;
X= AGE_HC_remaining_alltlas
IN.TrPred=predicted_AGE_HC_aal3;
IN.TrObs=AGE_HC_bal_aal3;

 [ sY, Y, IN ] = nk_CorrectPredTails(Y, X, IN)
 predicted_age_HCremaining_aal3_corrected=sY;
 BrainAGE_HCremaining_aal3_corrected=predicted_age_HCremaining_aal3_corrected - X;
 meanBrainAGE_HCremaining_aal3_corrected=mean(BrainAGE_HCremaining_aal3_corrected)
stdBrainAGE_HCremaining_aal3_corrected=std(BrainAGE_HCremaining_aal3_corrected)

table_HC_remaining_aal3=table(AGE_HC_remaining_alltlas,predicted_AGE_HCremaining_aal3,BrainAGE_HCrem_aal3,predicted_age_HCremaining_aal3_corrected,BrainAGE_HCremaining_aal3_corrected,'VariableNames',{'original_age','predicted_age','BrainAGE','predicted_age_corr','BrainAGE_corr'})
scatter(predicted_age_HCremaining_aal3_corrected,AGE_HC_remaining_alltlas)
writetable(table_HC_remaining_aal3,fullfile(out_dir,'aal3_HCremaining_formanuscript_new.csv'));

 
Y= predicted_AGE_Scz_aal3;
X= AGE_SCZ_aal3;
IN.TrPred=predicted_AGE_HC_aal3;
IN.TrObs=AGE_HC_bal_aal3;

 [ sY, Y, IN ] = nk_CorrectPredTails(Y, X, IN)
 predicted_age_SCZ_aal3_corrected=sY;
 BrainAGE_SCZ_aal3_corrected=predicted_age_SCZ_aal3_corrected - X;
 meanBrainAGE_Scz_aal3_corrected=mean(BrainAGE_SCZ_aal3_corrected)
  stdBrainAGE_Scz_aal3_corrected=std(BrainAGE_SCZ_aal3_corrected)

table_SCZ_aal3=table(AGE_SCZ_aal3,predicted_AGE_Scz_aal3,BrainAGE_Scz_aal3,predicted_age_SCZ_aal3_corrected,BrainAGE_SCZ_aal3_corrected,'VariableNames',{'original_age','predicted_age','BrainAGE','predicted_age_corr','BrainAGE_corr'})
scatter(predicted_AGE_Scz_aal3,AGE_SCZ_aal3)
writetable(table_SCZ_aal3,fullfile(out_dir,'aal3_SCZ_formanuscript_new.csv'));


%% Hammers


AGE_HC_bal_hammers = HC_norm_hammers.original_ages;
AGE_SCZ_hammers=SCZ_hammers.original_age_scz;
AGE_HC_remaining_alltlas=HC_rem_hammers.original_age_HC_rem;
predicted_AGE_HC_hammers = HC_norm_hammers.average_predictions;
predicted_AGE_Scz_hammers = SCZ_hammers.average_predictions_patients_scz;
predicted_AGE_HCremaining_hammers = HC_rem_hammers.average_predictions_HC_rem;

BrainAGE_HC_hammers=predicted_AGE_HC_hammers-AGE_HC_bal_hammers;
meanBrainAGE_hammers=mean(BrainAGE_HC_hammers)
BrainAGE_HCrem_hammers = predicted_AGE_HCremaining_hammers-AGE_HC_remaining_alltlas;
meanBrainAGE_HC_apply=mean(BrainAGE_HCrem_hammers)
BrainAGE_Scz_hammers = predicted_AGE_Scz_hammers-AGE_SCZ_hammers;
meanBrainAGE_Scz=mean(BrainAGE_Scz_hammers)

IN=struct;
Y= predicted_AGE_HC_hammers;
X= AGE_HC_bal_hammers;
IN.TrPred=predicted_AGE_HC_hammers;
IN.TrObs=AGE_HC_bal_hammers;
 [ sY, Y, IN ] = nk_CorrectPredTails(Y, X, IN)
 predicted_age_HC_hammers_corrected=sY;
 BrainAGE_HC_hammers_corrected=predicted_age_HC_hammers_corrected - X;
 meanBrainAGE_HC_hammers_corrected=mean(BrainAGE_HC_hammers_corrected)
stdBrainAGE_HC_hammers_corrected=std(BrainAGE_HC_hammers_corrected)
scatter(predicted_age_HC_hammers_corrected,AGE_HC_bal_hammers)

table_HC_norm_hammers=table(AGE_HC_bal_hammers,predicted_AGE_HC_hammers,BrainAGE_HC_hammers,predicted_age_HC_hammers_corrected,BrainAGE_HC_hammers_corrected,'VariableNames',{'original_age','predicted_age','BrainAGE','predicted_age_corr','BrainAGE_corr'})
writetable(table_HC_norm_hammers,fullfile(out_dir,'hammers_HCnorm_formanuscript_new.csv'));

Y= predicted_AGE_HCremaining_hammers;
X= AGE_HC_remaining_alltlas
IN.TrPred=predicted_AGE_HC_hammers;
IN.TrObs=AGE_HC_bal_hammers;

 [ sY, Y, IN ] = nk_CorrectPredTails(Y, X, IN)
 predicted_age_HCremaining_hammers_corrected=sY;
 BrainAGE_HCremaining_hammers_corrected=predicted_age_HCremaining_hammers_corrected - X;
 meanBrainAGE_HCremaining_hammers_corrected=mean(BrainAGE_HCremaining_hammers_corrected)
  stdBrainAGE_HCremaining_hammers_corrected=std(BrainAGE_HCremaining_hammers_corrected)

table_HC_remaining_hammers=table(AGE_HC_remaining_alltlas,predicted_AGE_HCremaining_hammers,BrainAGE_HCrem_hammers,predicted_age_HCremaining_hammers_corrected,BrainAGE_HCremaining_hammers_corrected,'VariableNames',{'original_age','predicted_age','BrainAGE','predicted_age_corr','BrainAGE_corr'})
scatter(predicted_age_HCremaining_hammers_corrected,AGE_HC_remaining_alltlas)
writetable(table_HC_remaining_hammers,fullfile(out_dir,'hammers_HCremaining_formanuscript_new.csv'));

 
Y= predicted_AGE_Scz_hammers;
X= AGE_SCZ_hammers;
IN.TrPred=predicted_AGE_HC_hammers;
IN.TrObs=AGE_HC_bal_hammers;

 [ sY, Y, IN ] = nk_CorrectPredTails(Y, X, IN)
 predicted_age_SCZ_hammers_corrected=sY;
 BrainAGE_SCZ_hammers_corrected=predicted_age_SCZ_hammers_corrected - X;
 meanBrainAGE_Scz_hammers_corrected=mean(BrainAGE_SCZ_hammers_corrected)
 stdmeanBrainAGE_Scz_hammers_corrected=std(BrainAGE_SCZ_hammers_corrected)
table_SCZ_hammers=table(AGE_SCZ_hammers,predicted_AGE_Scz_hammers,BrainAGE_Scz_hammers,predicted_age_SCZ_hammers_corrected,BrainAGE_SCZ_hammers_corrected,'VariableNames',{'original_age','predicted_age','BrainAGE','predicted_age_corr','BrainAGE_corr'})
scatter(predicted_AGE_Scz_hammers,AGE_SCZ_hammers)
writetable(table_SCZ_hammers,fullfile(out_dir,'hammers_SCZ_formanuscript_new.csv'));

%pairwise comparision of the brainage differences between the HC norm and SCZ groups
[h,p,ci,stats]=ttest2(BrainAGE_HC_allatlas_corrected, BrainAGE_SCZ_allatlas_corrected)
[h,p,ci,stats]=ttest2(BrainAGE_HC_schaefer_corrected, BrainAGE_SCZ_schaefer_corrected)
[h,p,ci,stats]=ttest2(BrainAGE_HC_aal3_corrected, BrainAGE_SCZ_aal3_corrected)
[h,p,ci,stats]=ttest2(BrainAGE_HC_hammers_corrected, BrainAGE_SCZ_hammers_corrected)


[h,p,ci,stats]=ttest2(BrainAGE_HC_allatlas_corrected, BrainAGE_HCremaining_allatlas_corrected)
[h,p,ci,stats]=ttest2(BrainAGE_HC_schaefer_corrected, BrainAGE_HCremaining_schaefer_corrected)
[h,p,ci,stats]=ttest2(BrainAGE_HC_aal3_corrected, BrainAGE_HCremaining_aal3_corrected)
[h,p,ci,stats]=ttest2(BrainAGE_HC_hammers_corrected, BrainAGE_HCremaining_hammers_corrected)

[h,p,ci,stats]=ttest2(BrainAGE_HCremaining_allatlas_corrected,BrainAGE_SCZ_allatlas_corrected)
[h,p,ci,stats]=ttest2(BrainAGE_HCremaining_schaefer_corrected,BrainAGE_SCZ_schaefer_corrected)
[h,p,ci,stats]=ttest2(BrainAGE_HCremaining_aal3_corrected,BrainAGE_SCZ_aal3_corrected)
[h,p,ci,stats]=ttest2(BrainAGE_HCremaining_hammers_corrected,BrainAGE_SCZ_hammers_corrected)

[h, p, ci, stats] = vartest2(AGE_HC_bal_allatlas, AGE_HC_remaining_alltlas)
