function [expData_PL_relative_strength, expData_switching_frq] = readExperimentalData(file_experimental_data)
% Function to read and process experimental data of PL promoter strength and switching frequency.
% INPUT: Filename to experimental data. The data format should be "long format".
% OUTPUT: 
%     Two Nx2 matrices (2 column matrices, N is the number of unique promoter IDs)
%           one matrix for (normalized?) promoter strength
%           one matrix for switching frequency
%     The first column is the mean value of the data.

%file_experimental_data = 'experimental_data_v1.xlsx'; % only for testing

if exist(file_experimental_data) == 0
    error(sprintf('Cannot load experimental data. File does not exists: %s', file_experimental_data))
end

%% Read data

% xlsread('experimental_data_v1.xlsx', 'Full_data_MK_curated_2016-03')
tbl_raw_experimental_data = readtable(file_experimental_data, ...
    'FileType','spreadsheet', ...
    'Sheet', 'Full_data_MK_curated_2016-03',...
    'ReadRowNames', false,...
    'ReadVariableNames', true);


%% Summarize data

tbl_stats = grpstats(tbl_raw_experimental_data, {'Seq_name', 'Variable'}, {'mean','std'}, 'DataVars',{'Value'});
    % 'VarNames' controls the output names, but the default constructed names are good.

% Alternative approach
%res = varfun(@mean,tbl_raw_experimental_data,'InputVariables','Value', 'GroupingVariables', {'Seq_ID', 'Variable'})


%% Reshape data: unstack data from single variable into multiple variables

tbl_stats_wide = unstack(tbl_stats, {'mean_Value', 'std_Value'}, 'Variable', 'GroupingVariables', 'Seq_name');
    % vars (first argument): analogous to measure.vars
    % ivar (second argument): analogous to id.vars
    % 'GroupingVariables': The default is all the variables in T not listed in vars or ivar.
    % 'ConstantVariables': The values for these variables in W are taken from the first row in each group in T.

% Set the RowNames for convenience of indexing later on...
tbl_stats_wide.Properties.RowNames = tbl_stats_wide.Seq_name;

%% Normalize the PL promoter strength data
% We are not interested in relative PL activity values.
% Hence we use the wild-type "C1" as a normalization reference for 100% PL activity.
% This gives us RELATIVE PL activity MEAN values.
% To obtain RELATIVE standard error values, we normalize the standard error
% by the mean wild-type activity. Hence the relative std says how much the measurements varies
% compared to the wild-type activity. 
% THIS IS ALSO HOW JESPER TVENGE CALCULATED THE % STD.

wt_mean_pl_activity = tbl_stats_wide{'C1', 'mean_Value_Specific_PL_activity'};

% Copy data
tbl_stats_wide_norm = tbl_stats_wide;
% Normalize
tbl_stats_wide_norm{:,{'mean_Value_Specific_PL_activity', 'std_Value_Specific_PL_activity'}} = tbl_stats_wide_norm{:,{'mean_Value_Specific_PL_activity', 'std_Value_Specific_PL_activity'}}/wt_mean_pl_activity;


% ALTERNATIVE: here is how to do it if you do not have indexing on rowname
% tbl_stats_wide.mean_Value_Specific_PL_activity(strcmp(tbl_stats_wide.Seq_name, 'C1'));

%% Creating output variables

% TWO column matrices
expData_PL_relative_strength = tbl_stats_wide_norm{:,{'mean_Value_Specific_PL_activity', 'std_Value_Specific_PL_activity'}};
expData_switching_frq = tbl_stats_wide_norm{:,{'mean_Value_SwitchFrequency', 'std_Value_SwitchFrequency'}};

% ....
% expData_PL_relative_strength_mean = tbl_stats_wide_norm{:,{'mean_Value_Specific_PL_activity'}};
% expData_switching_frq_mean = tbl_stats_wide_norm{:,{'mean_Value_SwitchFrequency'}};
% expData_PL_relative_strength_std = tbl_stats_wide_norm{:,{'std_Value_Specific_PL_activity'}};
% expData_switching_frq_std = tbl_stats_wide_norm{:,{'std_Value_SwitchFrequency'}};

%% ----- END -----
end % end function



%% =============================== OLD STUF ==========================================

%% Smoothing

%yy = smooth(y,method)
% yy = smooth(expData_switching_frq_mean, expData_PL_relative_strength_mean, 'lowess')
    % The ouput is a vector with the same dimensions as the input vector
% http://www.mathworks.com/help/curvefit/smooth.html
% 'lowess': Local regression using weighted linear least squares and a 1st degree polynomial model
% 'loess': Local regression using weighted linear least squares and a 2nd degree polynomial model



%% UITABLE
%f = figure;
%t = uitable(f,'Data',table2cell(tbl_stats_wide_norm));


%% READING 'data summary'
% [num,txt,raw] = xlsread('experimental_data_v1.xlsx', 'Data summary - manual calc')
%expData_PL_relative_strength = xlsread(file_experimental_data, 'Data summary - manual calc', 'C:C');
%expData_switching_frq = xlsread(file_experimental_data, 'Data summary - manual calc', 'D:D');

