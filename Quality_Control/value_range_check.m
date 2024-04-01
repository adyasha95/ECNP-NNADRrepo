function [Error_values_range,Error_details_val_range] = value_range_check(Table,Dictionary)
%% DESCRIPTION
% This function check if data follow the specifications defined it the
% dictionary in terms of item existence, range and values
%% INPUT
% Table = table; table containing the data
% Dictionary = table; dictionary, expected headers =
% 'ItemName',VariableScale, VariableType and Value
%% OUTPUT
% Error_values_range = cell array of string containing Error_values_range messages for
% the whole item
% Error_values_range_details_val_range = cell array of string containing Error_values_range messages for
% the individual values
%% CREDENTIALS
% Anne Ruef & Madalina Buciuman @Section Precision Psychiatry
%%
% initialise the inputs
Error_values_range = [];
Error_details_val_range = [];
% check the inputs
if ~istable(Table)
    Error_values_range = [Error_values_range;{'Input Table: not a table'}];
    return
end
if ~istable(Dictionary)
    Error_values_range = [Error_values_range;{'Input Dictionary: not a table'}];
    return
end

%% read the dictionnary
% check expected headers
ExpectedHeaders = {'ItemName','VariableType','VariableScale','Value'};
if any(~ismember(ExpectedHeaders,Dictionary.Properties.VariableNames))
    Error_values_range = [Error_values_range;{'Dictionnary: incorrect headers ItemName,VariableType,VariableScale,Value are expected'}];
    disp('Dictionnary: incorrect headers ItemName,VariableType,VariableScale,Value are expected');
    return
end
% check that the headers in Table are in dictionnary
if any(~ismember(Table.Properties.VariableNames,Dictionary.ItemName));
    % define the list of the items that are not in the dictionary
    ItemNotInDict = Table.Properties.VariableNames(~ismember(Table.Properties.VariableNames,Dictionary.ItemName));
    % export the items that is not in the dictionary 
    for i = 1:numel(ItemNotInDict)
        Error_values_range = [Error_values_range;{[ItemNotInDict{i},' is not in the dictionnary please remove the item from the table']}];
        disp([ItemNotInDict{i},' is not in the dictionnary please remove the item from the table'])
    end
end
% create a variable PID that is supposed to be numeric and check its format
PID = Table.PID;
if isnumeric(PID)
elseif iscell(PID)
    try
        PID = cellfun(@str2double,Table.PID);
    catch err
        Error_values_range = [Error_values_range;{['PID is in the wrong format, numeric values are expected']}];
        disp('PID is in the wrong format, numeric values are expected')
        for i = 1:numel(PID)
            try
                PID = cellfun(@str2double,Table.PID(i));
            catch err
                Error_details_val_range = [Error_details_val_range;...
                    {['PID values for row ',num2str(i),' has wrong format']}];
            end
        end
        disp('PID done')
        return
    end
end

% check all items in dictionnary
for i = 1:size(Dictionary,1)
    try
        % export item name variable type, variable scale and values range
        ItemName        = Dictionary.ItemName{i};
        disp(['Working on item ',ItemName]);
        VariableType    = Dictionary.VariableType{i};
        VariableScale   = Dictionary.VariableScale{i};
        Value           = Dictionary.Value{i};
        
        % check if the item is in the table    
        if ismember({ItemName},Table.Properties.VariableNames)
        else
              Error_values_range = [Error_values_range;[ItemName,' is not in Table please add a column with the values \N']];
              disp([ItemName,' is not in Table please add a column with the values \N']);
              disp([ItemName,' done'])
              continue
        end
        if ismember({VariableType},{'numeric'})
            if ismember({VariableScale},{'categorical'})
                % expected format = number="CategoryName"
                TempVal = strsplit(Value,{'",','="'},'CollapseDelimiters',true);              
                NumValues = [];               
                for t = 1:2:numel(TempVal)
                    % extract the values and check the format of the
                    % dictionnary
                    if isnan(str2double(TempVal{t}))
                        Error_values_range = [Error_values_range;['dictionnary values for item ',ItemName,' are not in the right format ValueNumber="ValueName" is expected but was ',Value]];
                        NumValues = [];
                        disp(['dictionnary values for item ',ItemName,' are not in the right format ValueNumber="ValueName" is expected but was ',Value])
                        break
                    else
                        NumValues = [NumValues,str2double(TempVal{t})];
                    end
                end
                
                if isempty(NumValues) % cannot check the data 
                    Error_values_range = [Error_values_range;['values for item ',ItemName,' no values definition found for this item: please check the dictionary values']];
                    disp(['values for item ',ItemName,' no values definition found for this item: please check the dictionary values'])
                else
                    % check if the value is a cell
                    if iscell(Table.(ItemName))
                        % it means we have a nan \N
                        TempValues = unique(cellfun(@str2double,Table.(ItemName)(~ismember(Table.(ItemName),{'\N'}))));
                        % check if the values in the item correspond to the
                        % values defined in dictionary
                        if any(~ismember(TempValues,NumValues))
                            % export the error
                            Error_values_range = [Error_values_range;['values for item ',ItemName,' are not coded correctly values for this item: ',num2str(TempValues'), ' these values are expected: ',num2str(NumValues)]];
                            disp(['values for item ',ItemName,' are not coded correctly values for this item: ',num2str(TempValues'), ' these values are expected: ',num2str(NumValues)]);
                            % find the individual values that are not correct
                            ValToFind = TempValues(~ismember(TempValues,NumValues));
                            for ii = 1:numel(ValToFind)
                                for iii = 1:numel(Table.(ItemName))
                                    if isequal(Table.(ItemName){iii},num2str(ValToFind(ii)))
                                        Error_details_val_range = [Error_details_val_range;...
                                            {['values for item ',ItemName,' for PID: ',num2str(PID(iii)),' Row number ',num2str(iii),' Value ',Table.(ItemName){iii},' is incorrect, these values are expected: ',num2str(NumValues)]}];
                                    end
                                end
                            end
                        else
                            
                            if ismember({ItemName},Table.Properties.VariableNames)
                            else
                                % this is unecessary
                                Error_values_range = [Error_values_range;['item ',ItemName,' not found in table please add the item to the table']];
                                disp(['item ',ItemName,' not found in table please add the item to the table'])
                            end
                        end
                    elseif isnumeric(Table.(ItemName))
                        % check if there is any nans
                        if any(isnan(Table.(ItemName)))
                            % it means that it has been wrongly coded
                            % export the error code 
                            Error_values_range = [Error_values_range;['values for item ',ItemName,'contains nans, the nans should be replaced by \N ']];
                            disp(['values for item ',ItemName,'contains nans, the nans should be replaced by \N '])
                            % find the individual values that are incorrect
                            IDX = find(isnan(Table.(ItemName)));
                            for ii = 1:numel(IDX)
                                Error_details_val_range = [Error_details_val_range;...
                                    {[' for PID: ',num2str(PID(IDX(ii))),' Row number ',num2str(IDX(ii)),' Value ',num2str(Table.(ItemName)(IDX)),' is incorrect, these values are expected: \N']}];
                            end
                            disp([ItemName,' done'])
                            continue
                        else
                            % list all values
                            TempValues = unique(Table.(ItemName));
                            % check if any values are not correct
                            if any(~ismember(TempValues,NumValues))
                                % export the error message
                                Error_values_range = [Error_values_range;['values for item ',ItemName,' are not coded correctly, ',num2str(TempValues), ' these values are expected: ',num2str(NumValues)]];
                                disp(['values for item ',ItemName,' are not coded correctly, ',num2str(TempValues), ' these values are expected: ',num2str(NumValues)])
                                IDX = find(~ismember(Table.(ItemName),TempValues));
                                for ii = 1:numel(IDX)
                                    Error_details_val_range = [Error_details_val_range;...
                                        {[' for PID: ',num2str(PID(IDX(ii))),' Row number ',num2str(IDX(ii)),' Value ',num2str(Table.(ItemName)(IDX)),' is incorrect, these values are expected: ',num2str(NumValues)]}];
                                end
                            else
                                % maybe unecessary
                                if ismember({ItemName},Table.Properties.VariableNames)
                                else
                                    Error_values_range = [Error_values_range;['item ',ItemName,' not found in table']];
                                    disp(['item ',ItemName,' not found in table'])
                                end
                            end
                        end
                    else
                        % export the error message
                        Error_values_range = [Error_values_range;['for item ',ItemName,' values are not numeric']];
                        disp(['for item ',ItemName,' values are not numeric'])
                    end
                end
            elseif ismember({VariableScale},{'continuous'})
                % export the values from the dictionnary
                NumValues = strsplit(Value,',');
                if numel(NumValues) ~= 3
                    Error_values_range = [Error_values_range;['dictionnary values for item ',ItemName,' are not in the right format Min,Max,Decimal is expected but was ',Value,'']];
                    disp(['dictionnary values for item ',ItemName,' are not in the right format Min,Max,Decimal is expected but was ',Value,''])
                    disp([ItemName,' done'])
                    continue
                else
                    % nothing to do
                end
                % check if the item is in table
                if ismember({ItemName},Table.Properties.VariableNames)
                else
                    Error_values_range = [Error_values_range;['item ',ItemName,' not found in table']];
                    disp(['item ',ItemName,' not found in table'])
                    disp([ItemName,' done'])
                    continue
                end
                
                
                if iscell(Table.(ItemName))
                    TempValues = unique(cellfun(@str2double,Table.(ItemName)(~ismember(Table.(ItemName),{'\N'}))));
                elseif isnumeric(Table.(ItemName))
                    % check if it is nan 
                    if any(isnan(Table.(ItemName)))
                        Error_values_range = [Error_values_range;['values for item ',ItemName,'contains nans, the nans should be replaced by \N ']];
                        disp(['values for item ',ItemName,'contains nans, the nans should be replaced by \N '])
                        IDX = find(isnan(Table.(ItemName)));
                        for ii = 1:numel(IDX)
                            Error_details_val_range = [Error_details_val_range;...
                                {[' for PID: ',num2str(PID(IDX(ii))),' Row number ',num2str(IDX(ii)),' Value ',num2str(Table.(ItemName)(IDX)),' is incorrect, these values are expected: \N']}];
                        end   
                        disp([ItemName,' done'])
                        continue
                    else
                        TempValues = unique(Table.(ItemName));
                    end
                end
                % check the min
                if min(TempValues) < str2double(NumValues{1})
                    Error_values_range = [Error_values_range;['values for item ',ItemName,' min value ',num2str(min(TempValues)),' is smaller than the expected value ',NumValues{1},' please check your data']];
                    disp(['values for item ',ItemName,' min value ',num2str(min(TempValues)),' is smaller than the expected value ',NumValues{1},' please check your data'])
                    if iscell(Table.(ItemName))
                        ValToFind = TempValues(TempValues < str2double(NumValues{1}));
                        for ii = 1:numel(ValToFind)
                            for iii = 1:numel(Table.(ItemName))
                                if isequal(Table.(ItemName){iii},num2str(ValToFind(ii)))
                                    Error_details_val_range = [Error_details_val_range;...
                                        {['values for item ',ItemName, ' for PID: ',num2str(PID(iii)),' Row number ',num2str(iii),' Value ',Table.(ItemName){iii},' is incorrect, these values should not be smaller than: ',NumValues{1}]}];
                                end
                            end
                        end
                    elseif isnumeric(Table.(ItemName))
                        IDX = find(Table.(ItemName) < str2double(NumValues{1}));
                        for ii = 1:numel(IDX)
                            Error_details_val_range = [Error_details_val_range;...
                                {['values for item ',ItemName, ' for PID: ',num2str(PID(iii)),' Row number ',num2str(iii),' Value ',Table.(ItemName)(iii),' is incorrect, these values should not be smaller than: ',NumValues{1}]}];
                        end
                    else
                        % this is an error
                    end
                end
                % check the max
                if max(TempValues) > str2double(NumValues{2})
                    Error_values_range = [Error_values_range;['values for item ',ItemName,' max value ',num2str(max(TempValues)),' is higher than the expected value ',NumValues{2},' please check your data']];
                    disp(['values for item ',ItemName,' max value ',num2str(max(TempValues)),' is higher than the expected value ',NumValues{2},' please check your data'])
                    if iscell(Table.(ItemName))
                        ValToFind = TempValues(TempValues > str2double(NumValues{2}));
                        for ii = 1:numel(ValToFind)
                            for iii = 1:numel(Table.(ItemName))
                                if isequal(Table.(ItemName){iii},num2str(ValToFind(ii)))
                                    Error_details_val_range = [Error_details_val_range;...
                                        {['values for item ',ItemName, ' for PID: ',num2str(PID(iii)),' Row number ',num2str(iii),' Value ',Table.(ItemName){iii},' is incorrect, these values should not be higher than: ',NumValues{2}]}];
                                end
                            end
                        end
                    elseif isnumeric(Table.(ItemName))
                        IDX = find(Table.(ItemName) > str2double(NumValues{2}));
                        for ii = 1:numel(IDX)
                            Error_details_val_range = [Error_details_val_range;...
                                {['values for item ',ItemName, ' for PID: ',num2str(PID(iii)),' Row number ',num2str(iii),' Value ',Table.(ItemName)(IDX(ii)),' is incorrect, these values should not be higher than: ',NumValues{1}]}];
                        end
                    else
                        % this is an error
                    end
                    
                end
                % check the decimal
                if sum(double(int64(TempValues*10^(str2double(NumValues{3})))) - TempValues*10^(str2double(NumValues{3}))) > 10^(-str2double(NumValues{3})-1)
                    Error_values_range = [Error_values_range;['values for item ',ItemName,' the number of decimals ',NumValues{3},' is expected please check your data']];
                    disp(['values for item ',ItemName,' the number of decimals ',NumValues{3},' is expected please check your data'])
                    if iscell(Table.(ItemName))
                        ValToFind = TempValues(double(int64(TempValues*10^(str2double(NumValues{3})))) - TempValues*10^(str2double(NumValues{3})) ~=0);
                        for iii = 1:numel(Table.(ItemName))
                            for ii = 1:numel(ValToFind)                               
                                if isequal(Table.(ItemName){iii},num2str(ValToFind(ii)))
                                    Error_details_val_range = [Error_details_val_range;...
                                        {['values for item ',ItemName, ' for PID: ',num2str(PID(iii)),' Row number ',num2str(iii),' Value ',Table.(ItemName){iii},' is incorrect, the number of decimal expected is: ',NumValues{3}]}];
                                    break
                                end
                            end
                        end
                    elseif isnumeric(Table.(ItemName))
                        IDX = find(Table.(ItemName)((double(int64(Table.(ItemName)*10^(str2double(NumValues{3})))) - Table.(ItemName)*10^(str2double(NumValues{3}))) ~=0));
                        for ii = 1:numel(IDX)
                            Error_details_val_range = [Error_details_val_range;...
                                {['values for item ',ItemName, ' for PID: ',num2str(PID(iii)),' Row number ',num2str(iii),' Value ',Table.(ItemName)(IDX(ii)),' is incorrect, the number of decimal expected is: ',NumValues{3}]}];
                        end
                    else
                        % this is an error
                    end
                    
                end
                
            end
        else
            Error_values_range = [Error_values_range;{[ItemName,' is expected to be numeric']}];
            disp([ItemName,' is expected to be numeric'])
            disp([ItemName,' done'])
        end
        % display that it is done
        disp([ItemName,' done'])
    catch err
        % this is for any other potential errors
        Error_values_range = [Error_values_range;{[ItemName,' unknown error in value_range_check :',err.message,'']}];
        disp([ItemName,' unknown error in value_range_check :',err.message,''])
        disp([ItemName,' done'])
    end
end
end