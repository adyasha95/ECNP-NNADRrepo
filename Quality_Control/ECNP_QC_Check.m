%% Script for checking all the QC measures for ECNP
clear all
% Add paths to the functions that will run the QC check
addpath('');

% Path to MRI and Clinical dictionaries
% Ex. dictionary_MRI = readtable('/volume/projects/JF_ECNP/MUC/ECNP/Dictionary/MRIdictionary_ECNP_multiAtlas_Nov22.xlsx');
dictionary_MRI = readtable('');
% Ex. dictionary_clinical = readtable('/volume/projects/JF_ECNP/MUC/ECNP/Dictionary/version8_Dictionary_ECNP_Clinical_October_7_2022.xlsx');
dictionary_clinical = readtable('');
% Merge dictionaries
dictionary_all = [dictionary_clinical;dictionary_MRI];

% Path to the data table
% Ex. filename = '/volume/projects/JF_ECNP/MUC/ECNP/allAtlases_clinical_MUC_14_Dec_2022_14-Dec-2022_all_atlases.csv';
filename = '';
opts = detectImportOptions(filename);
opts = setvartype(opts, opts.SelectedVariableNames, 'char');
Table = readtable(filename, opts);                                   
%% Round all ROI values to 2 decimals according to the dictionary and save the corrected table to be used for the next steps and the ViPAR upload
outputpath = ''; % please give the path to the folder where to save the corrected table (with the ROIs rounded to 2 decimals)
Table_corr = round(cellfun(@str2num,(strrep(table2array(Table(:, ismember(Table.Properties.VariableNames, {'TIV','NCR','ICR'}) | strfind(Table.Properties.VariableNames, {'schaefer', 'aal3', 'hammers'}))),'\N','NaN'))),2);
Table(:,ismember(Table.Properties.VariableNames, {'TIV','NCR','ICR'}) | strfind(Table.Properties.VariableNames, {'schaefer', 'aal3', 'hammers'})) = strrep(arrayfun(@num2str, Table_corr, 'UniformOutput',false),'NaN','\N');

% Add extra variables and adjust clinical variables
% Renaming PID column and removing additional PID column that was a result of table merging
Table = renamevars(Table, 'PID_ClinData', 'PID');
Table(:, 'PID_Table') = [];

% Skip running this function if your data table is in accordance with current clinical dictionary
Table = variable_check_AK(dictionary_clinical,dictionary_MRI,Table);

% Save table with updated variables and indicate file name
writetable(Table,[outputpath,'']);
%% QC Check
%Table = readtable('/volume/AR_ECNP/Analysis/ECNPRun2/ECNP/Data/table_ROIs_QC_clinical_nanCorr.csv','string');

%values and range check
[Error_values_range,Error_details_val_range] = value_range_check(Table,dictionary_all);
%MRI QC
[Error_MRI_QC,Error_details_MRI_QC] = MRI_QC_check(Table,dictionary_MRI);
%AGE, SEX check
%[Error_age_sex] = age_sex_check(Table,dictionary_clinical);
%CTQ check
[Error_clinical] = clinical_check(Table); 
%% Write the errors to a log file
Errors = [Error_values_range;Error_MRI_QC;Error_age_sex;Error_clinical];

filePh = fopen(fullfile(outputpath,['Error_log_ECNP',datestr(now),'.txt']),'w');
fprintf(filePh,'%s\n',Errors{:});
fclose(filePh);