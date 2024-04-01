function [Errors] = clinical_check(Table) 
%to check whether the data table has more than one site ID and whether any
%site ID is missing (i.e. \N)
%to check whether CTQ items were corretly reversed (where needed), and
%binarized (where needed)
%to check that the total and subscale scores (CTQ, PANSS, MADRS, HDRS)
%match the individual items
%
% Inputs
%   Table   : table containing all data that are prepared for ViPAR
%
% Output
%   Errors  : error messages, if some

%initialize array for potential errors and its counter c
Errors = {};
c = 1;

%vcheck if table contains only one site ID 
sites = unique(Table.SITE_ID);  
if (length(sites) > 1)
    error_str = strcat('Your data table contains data from more than one site ID.');
    %warning(error_str)
    Errors(c) = {error_str};
    c = c + 1;
end
%check whether site ID is missing for any subject
if (sum(isnan(str2double(sites))) > 0)
    error_str = 'Your data table contains \N as site ID.';
    %warning(error_str)
    Errors(c) = {error_str};
    c = c + 1;
end

%check reversed scores of CTQ 
disp('Checking reversed CTQ items...')
rev_items = {'CTQ_02','CTQ_05','CTQ_07','CTQ_13','CTQ_19','CTQ_26','CTQ_28'};
for i=1:length(rev_items)
    org = Table.(rev_items{i});
    org(cellfun(@isempty,org)) = {"NaN"};
    rev = Table.([rev_items{i} '_reversed']);
    rev(cellfun(@isempty,rev)) = {"NaN"};
    for j=1:length(org)
        if org{j} ~= reverse_score(rev{j})
            error_str = [rev_items{i} ' is incorrectly reversed for PID ' Table.PID{j}];
            %warning(error_str)
            Errors(c) = {error_str};
            c = c + 1;
        end
    end
end

%check binarized scores of CTQ
disp('Checking binarized CTQ items...')
binar_items = {'CTQ_10','CTQ_16','CTQ_22'};
for i=1:length(binar_items)
    org = Table.(binar_items{i});
    org(cellfun(@isempty,org)) = {"NaN"};
    binar = Table.([binar_items{i} '_binarized']);
    binar(cellfun(@isempty,binar)) = {"NaN"};
    for j=1:length(org)
        if binar{j} ~= binarize_score(org{j})
            error_str = [binar_items{i} ' is incorrectly binarized for PID ' Table.PID{j}];
            %warning(error_str)
            Errors(c) = {error_str};
            c = c + 1;
        end
    end
end

%check CTQ sum scores
disp('Checking CTQ scales...')
sum_scores = {'CTQ_eabu','CTQ_pabu','CTQ_sabu','CTQ_eneg','CTQ_pneg','CTQ_vscale','CTQ_total_score'};
eabu_items = {'CTQ_03','CTQ_08','CTQ_14','CTQ_18','CTQ_25'};
pabu_items = {'CTQ_09','CTQ_11','CTQ_12','CTQ_15','CTQ_17'};
sabu_items = {'CTQ_20','CTQ_21','CTQ_23','CTQ_24','CTQ_27'};
eneg_items = {'CTQ_05_reversed','CTQ_07_reversed','CTQ_13_reversed','CTQ_19_reversed','CTQ_28_reversed'};
pneg_items = {'CTQ_01','CTQ_02_reversed','CTQ_04','CTQ_06','CTQ_26_reversed'};
vscale_items = {'CTQ_10_binarized','CTQ_16_binarized','CTQ_22_binarized'};
score_items = {eabu_items, pabu_items, sabu_items, eneg_items, pneg_items, vscale_items, sum_scores(1:5)};
Errors = check_scores(Table, sum_scores, score_items, Errors);

%check PANSS sum scores
disp('Checking PANSS scales...')
sum_scores = {'PANSS_general','PANSS_negative','PANSS_positive','PANSS_total_score'};
panss_gen = Table.Properties.VariableNames(contains(Table.Properties.VariableNames, 'PANSS_G'));
panss_gen = setdiff(panss_gen, sum_scores(1));
panss_neg = Table.Properties.VariableNames(contains(Table.Properties.VariableNames, 'PANSS_N'));
panss_neg = setdiff(panss_neg, sum_scores(2));
panss_pos = Table.Properties.VariableNames(contains(Table.Properties.VariableNames, 'PANSS_P'));
panss_pos = setdiff(panss_pos, sum_scores(3));
score_items = {panss_gen, panss_neg, panss_pos, sum_scores(1:3)};
Errors = check_scores(Table, sum_scores, score_items, Errors);

%check MADRS sum score
disp('Checking MADRS scale...')
sum_scores = {'MADRS_total_score'};
madrs = Table.Properties.VariableNames(contains(Table.Properties.VariableNames, 'MADRS'));
madrs = setdiff(madrs, sum_scores(1));
score_items = {madrs};
Errors = check_scores(Table, sum_scores, score_items, Errors);

%check HDRS17 sum score
disp('Checking HDRS17 scale...')
sum_scores = {'HDRS17_total_score'};
hdrs17 = Table.Properties.VariableNames(contains(Table.Properties.VariableNames, 'HDRS17'));
hdrs17 = setdiff(hdrs17, sum_scores);
score_items = {hdrs17};
Errors = check_scores(Table, sum_scores, score_items, Errors);

%check HDRS21 sum score
disp('Checking HDRS21 scale...')
sum_scores = {'HDRS21_total_score'};
hdrs21 = Table.Properties.VariableNames(contains(Table.Properties.VariableNames, 'HDRS21'));
hdrs21 = setdiff(hdrs21, sum_scores);
score_items = {hdrs21};
Errors = check_scores(Table, sum_scores, score_items, Errors);

%transpose error messages
Errors = Errors';

end


% ==================== %
% Additional functions %
% ==================== %

function [rev_score] = reverse_score(score)
%reverses CTQ scores
if score == '1'
    rev_score = '5';
elseif score == '2'
    rev_score = '4';
elseif score == '3'
    rev_score = '3';
elseif score == '4'
    rev_score = '2';
elseif score == '5'
    rev_score = '1';
else
    rev_score = '\N';
end
end


function [binar_score] = binarize_score(score)
%binarizes CTQ scores
if score == '1'
    binar_score = '0';
elseif score == '2'
    binar_score = '0';
elseif score == '3'
    binar_score = '0';
elseif score == '4'
    binar_score = '0';
elseif score == '5'
    binar_score = '1';
else
    binar_score = '\N';
end
end


function Errors = check_scores(Table,sum_scores, score_items, Errors)
%check whether individual items add up to total and subscale scores
%
% Inputs
%   Table       : table containing all data that are prepared for ViPAR
%   sum_scores  : cell array containing variable names of total and subscale
%                 scores
%   score_items : nested cell array containing individual item's variable
%                 names used for the different sum_scores
%   Errors      : cell array of error messages created by preceding script
%
% Output
%   Errors      : cell array of extended (if applicable) error messages

%set counter for Error messages
c = size(Errors, 2) + 1;

for i=1:length(sum_scores)
    %counter whether individual items were completely missing for given scale
    ind_items_missing = 0;
    %extract data for each subscale
    sub_table = Table(:,score_items{i});
    %convert subscale data from strings to double, loop over variables, i.e. columns
    for j=1:size(sub_table,2)
        sub_table.(score_items{i}{j}) = str2double(sub_table.(score_items{i}{j}));
    end
    %loop over rows, i.e. subjects
    for k=1:size(sub_table, 1)
        %sum up items
        ind_sum = sum(table2array(sub_table(k,:)));
        %get subscale score
        score = str2double(table2array(Table(k, sum_scores{i})));
        %compare whether scores match
        if ~isequaln(ind_sum, score)
            %check whether individual items were available
            if (sum(~isnan(table2array(sub_table(k,:)))) == 0)
                ind_items_missing = ind_items_missing + 1;
            else
                error_str = [sum_scores{i} ' does not match the sum of the individual items for PID ' Table.PID{k}];
                %warning(error_str)
                Errors(c) = {error_str};
                c = c + 1;
            end
        end
    end
    %if individual items for particular scale were missing for one or more
    %subjects, add error message
    if ind_items_missing > 0
        error_str = ['For ' sum_scores{i} ' individual items were missing for ' ...
            num2str(ind_items_missing) ' subjects (which is ok if you do not have them available).'];
        %warning(error_str)
        Errors(c) = {error_str};
        c = c + 1;
    end
end
end