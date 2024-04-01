%% This script is used to create a struct extracting ROI values from different atlases and quality measures
%%TASK: TO CHECK WITH THE MRI DICTIONARY ITEM NAMES
close all
clear all
clc
tic
%% User input

% Please indicate if you are using the standalone version (1) or the MATLAB-based version (0).

standalone = 1;

% Please enter the full path to the directory, where the CAT12.8 folder wascreated after processing.
% It is in you ECNP folder [/Your/Path]/ECNP/Data/raw;
% Example : OutputFolder = '/volume/AR_ECNP/Analysis/ECNPRun2/ECNP/Data/raw/';

in_dir = '';

% Please indicate the full path to the directory, where you would like to save your output file.
% can be [/Your/Path]/ECNP/Data/;

out_dir = '';

% PLease enter an output filename.

out_filename = '';

% Please enter the path incl. the filename, in which the image paths and participant
% IDs (PID) are stored:  "[/Your/Path]/ECNP/NameOfCSV.csv"
% ID_file = '/volume/AR_ECNP/Analysis/ECNPRun2/ECNP/MUC_ImagingPaths_ECNP.xlsx';
% If the file does not contain a column called 'PID', please change the
% M.PID in line 66 to the accurate variable name

ID_file = '';


% Please indicate how your files are named after the CAT processing.
% pattern_ROI: '.mat'-files in the 'label' folder (inside the 'CAT12.8' folder).
% After a prefix or pattern, the PID should follow (e.g. 'catROI_10001').
% Please enter the pattern/prefix before the PID.
% pattern_QM: '.mat'-files in the 'report' folder (inside the 'CAT12.8'folder).
% After a prefix or pattern, the PID should follow (e.g. 'cat_10001').
% Please enter the pattern before the PID.
% Note: If you are using the standalone version, the 'report' folder will
% be created outside the 'CAT12.8' folder in the input directory (in_dir).

pattern_ROI = 'catROI_';

pattern_QM = 'cat_';

% Please indicate, which atlas data (ROIs) you want to save as output
% Enter 'all' if you would like to save the output for all atlases or enter a
% specific atlas (e.g. 'Schaefer2018_200Parcels_17Networks_order').
% Please make sure you enter the full atlas name as given in your ROI files.
% Ex. atlas_to_save = {'Schaefer2018_200Parcels_17Networks_order';'aal3';'hammers'};
atlas_to_save = {};

% Please define the path to the csv file containing the clinical data.
ClinDataFile = '';

%% read IDs from file

M = readtable(ID_file);

% Please change M.PID to correct variable name in case the PID has a different name in ID_file.

PID = M.PID_MRI;

%% automatically read out the atlases from ROI struct (using first PID)

atlases = atlas_to_save;
% Please provide a short name for the atlas you have included in atlas_to_save.
% Ex. atlases_short = {'schaefer200';'aal3';'hammers'};
atlases_short = {};

%% struct creation for ROIs
for i = 1:size(PID,1)
    %disp(['ECNP Extracting ROIS for PID ',num2str(PID(i,1))]);
    
    ROI_struct = [];
    
    image_filename = [pattern_ROI,num2str(PID(i,1)),'.mat'];
    
    try
        ROI_struct = importdata(fullfile(in_dir,'CAT12.8','label',image_filename));
    catch
        warning(['The file "', image_filename, '" does not exist, please check. Assigning "/N".']);
    end
    
    for j = 1:size(atlases,1)
        
        if isfile(fullfile(in_dir,'CAT12.8','label',image_filename))
            ROI_struct.(atlases{j,1}).data.Vgm = num2cell(ROI_struct.(atlases{j,1}).data.Vgm);
            if isfield(ROI_struct.(atlases{j,1}).data,'Vwm')
                ROI_struct.(atlases{j,1}).data.Vwm = num2cell(ROI_struct.(atlases{j,1}).data.Vwm);
            end
        else
            ROI_struct.(atlases{j,1}).data.Vgm = cellstr(repmat("\N",size( S.atlases.(atlases{j,1}).Vgm,2)-1,1));
            ROI_struct.(atlases{j,1}).data.Vwm = cellstr(repmat("\N", size(S.atlases.(atlases{j,1}).Vgm,2)-1,1));
%             if isfield(ROI_struct.(atlases{j,1}).data,'Vwm')
%                 if i ~=1 && ~isfield(ROI_struct.(atlases{j,1}).data,'Vwm')
%                     ROI_struct.(atlases{j,1}).data.Vwm = cellstr(repmat("\N", size(S.atlases.(atlases{j,1}).Vgm,2)-1,1));
%                 end
%             end
        end
        if isfield(ROI_struct.(atlases{j,1}).data,'Vgm')
            S.atlases.(atlases{j,1}).Vgm(i,:) = splitvars(table([num2cell(PID(i,1)),ROI_struct.(atlases{j,1}).data.Vgm']));
        end
        if isfield(ROI_struct.(atlases{j,1}).data,'Vwm')
            S.atlases.(atlases{j,1}).Vwm(i,:) = splitvars(table([num2cell(PID(i,1)),ROI_struct.(atlases{j,1}).data.Vwm']));
        end
        if i == 1
            if isfield(ROI_struct.(atlases{j,1}).data,'Vgm')
                variable_name_prefix_gmv = strcat(repmat(atlases_short(j,1),size(ROI_struct.(atlases{j,1}).names,1),1),repmat(cellstr('_gmv_'),size(ROI_struct.(atlases{j,1}).names,1),1));
                
                S.atlases.(atlases{j,1}).Vgm.Properties.VariableNames = [cellstr('PID'); matlab.lang.makeValidName(strcat(variable_name_prefix_gmv,ROI_struct.(atlases{j,1}).names))]';
            end
            if isfield(ROI_struct.(atlases{j,1}).data,'Vwm')
                variable_name_prefix_wmv = strcat(repmat(atlases_short(j,1),size(ROI_struct.(atlases{j,1}).names,1),1),repmat(cellstr('_wmv_'),size(ROI_struct.(atlases{j,1}).names,1),1));
                S.atlases.(atlases{j,1}).Vwm.Properties.VariableNames = [cellstr('PID'); matlab.lang.makeValidName(strcat(variable_name_prefix_wmv,ROI_struct.(atlases{j,1}).names))]';
            end
        end
        
    end
    
end


%% struct creation for quality measures
for i = 1:size(PID,1)
    %disp(['ECNP Extracting QC measures for PID ',num2str(PID(i,1))]);
    
    QM_struct = [];
    
    QM_filename = [pattern_QM,num2str(PID(i,1)),'.mat'];
    
    try
        switch standalone
            case 0
                QM_struct = importdata(fullfile(in_dir,'CAT12.8','report',QM_filename));
            case 1
                QM_struct = importdata(fullfile(in_dir,'report',QM_filename));
        end
        
        QM_struct.subjectmeasures.vol_TIV = num2cell(QM_struct.subjectmeasures.vol_TIV);
        QM_struct.qualityratings.NCR = num2cell(QM_struct.qualityratings.NCR);
        QM_struct.qualityratings.ICR = num2cell(QM_struct.qualityratings.ICR);
        if QM_struct.qualityratings.IQR > 5.5
            QM_struct.qualityratings.QC_grade = num2cell(6);
        elseif QM_struct.qualityratings.IQR <= 5.5 && QM_struct.qualityratings.IQR > 4.5
            QM_struct.qualityratings.QC_grade = num2cell(5);
        elseif QM_struct.qualityratings.IQR <= 4.5 && QM_struct.qualityratings.IQR > 3.5
            QM_struct.qualityratings.QC_grade = num2cell(4);
        elseif QM_struct.qualityratings.IQR <= 3.5 && QM_struct.qualityratings.IQR > 2.5
            QM_struct.qualityratings.QC_grade = num2cell(3);
        elseif QM_struct.qualityratings.IQR <= 2.5 && QM_struct.qualityratings.IQR > 1.5
            QM_struct.qualityratings.QC_grade = num2cell(2);
        elseif QM_struct.qualityratings.IQR <= 1.5 && QM_struct.qualityratings.IQR > 0.5
            QM_struct.qualityratings.QC_grade = num2cell(1);
        end
    catch
        warning(['The file "', QM_filename, '" does not exist, please check. Assigning "\N".']);
        QM_struct.subjectmeasures.vol_TIV = cellstr('\N');
        QM_struct.qualityratings.NCR = cellstr('\N');
        QM_struct.qualityratings.ICR = cellstr('\N');
        QM_struct.qualityratings.QC_grade = cellstr('\N');
    end
    
    S.qualitymeasures(i,:) = splitvars(table([num2cell(PID(i,1)),QM_struct.subjectmeasures.vol_TIV, QM_struct.qualityratings.NCR, QM_struct.qualityratings.ICR, QM_struct.qualityratings.QC_grade]));
    if i == 1
        S.qualitymeasures.Properties.VariableNames = {'PID', 'TIV', 'NCR', 'ICR', 'QC_grade'};
    end
end

%% save file

% Uncomment for saving as struct

% if strcmp(atlas_to_save,'all')
%      save([out_dir,out_filename,'.mat'],'S');
% else
%     S_to_save.atlases.(atlas_to_save) = S.atlases.(atlas_to_save);
%     S_to_save.qualitymeasures = S.qualitymeasures;
%     S = S_to_save;
%     save([out_dir,'/',out_filename,'.mat'],'S');
% end
S_to_save=struct;
if ~strcmp(atlas_to_save,'all')
    for jj=1:size(atlas_to_save,1)
        S_to_save.atlases.atlas_to_save{jj} = S.atlases.(atlas_to_save{jj});
        S_to_save.qualitymeasures = S.qualitymeasures;
    end
end


% merge the different measures
QC = S_to_save.qualitymeasures;

for ii =1:size(atlas_to_save,1)
    if isfield(S_to_save.atlases.atlas_to_save{1,ii},'Vgm') && ~isempty(S_to_save.atlases.atlas_to_save{1,ii}.Vgm)
        GM{ii,:} = S_to_save.atlases.atlas_to_save{ii}.Vgm;
    end
    %GM_new=join(GM{ii})
    if isfield(S_to_save.atlases.atlas_to_save{ii},'Vwm') && ~isempty(S_to_save.atlases.atlas_to_save{1,ii}.Vwm)
        WM{ii,:} = S_to_save.atlases.atlas_to_save{ii}.Vwm;
    end
end
GM=GM(~cellfun('isempty',GM));
WM = WM(~cellfun('isempty',WM));
%needs to be soft coded -- the indices depend on the ordering of the atlases
GM=[GM{1} GM{2}(:,2:end) GM{3}(:,2:end)];
WM=[WM{1} WM{3}(:,2:end)];

Table = [QC,GM(:,~ismember(GM.Properties.VariableNames,{'PID'})),...
    WM(:,~ismember(WM.Properties.VariableNames,{'PID'}))];

% Table.Properties.RowNames = cellfun(@num2str,Table.PID,'UniformOutput',false);
Table.PID =  cell2mat(Table.PID);

% load the clinical data
opts = detectImportOptions(ClinDataFile);
opts = setvartype(opts, opts.SelectedVariableNames, 'char');
ClinData = readtable(ClinDataFile, opts);

% merge the clinical and mri data
Table.PID = cellstr(num2str(Table.PID));
TableMerged = outerjoin(ClinData,Table,'Keys','PID');

% Round Schaefer variables to 2 digits
%idxMRI = contains(TableMerged.Properties.VariableNames, atlases_short);
%TableMerged(:,idxMRI) = round(cell2mat(TableMerged{:,idxMRI}),2);

% replace the missing values by '\N'
var_names = TableMerged.Properties.VariableNames;
TableMerged_empty = cellfun(@isempty,table2cell(TableMerged),'UniformOutput',false);
TableMerged_cell = table2cell(TableMerged);
TableMerged_cell(cell2mat(TableMerged_empty)) = cellstr('\N');
TableMerged = splitvars(table(TableMerged_cell));
TableMerged.Properties.VariableNames = var_names;

writetable(TableMerged,fullfile(out_dir,[out_filename,'_',date,'_','all_atlases.csv']));
toc