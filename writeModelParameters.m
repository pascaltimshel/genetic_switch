
switch MODE
    case 'mode_PSR'
        model_parameters_MODE_SPECIFIC_cell = {
            'range_PSR_PL2PR', range_PSR_PL2PR;
            };
    case 'mode_startConcentration'
        model_parameters_MODE_SPECIFIC_cell = {}; % empty
        % range_CI_initial_concentration = [XXX];
        % range_MOR_initial_concentration = [XXX];
    case 'mode_modelInspection'
        model_parameters_MODE_SPECIFIC_cell = {}; % empty
    otherwise
        error('Received unexpected MODE variable. Change the MODE variable to a valid value.')
end


model_parameters_GENERAL_cell = {
'PSR_PL2PR', PSR_PL2PR;
'CI_initial_concentration', CI_initial_concentration;
'MOR_initial_concentration', MOR_initial_concentration;
'n_simulations', n_simulations;
'p_include_O_M_site', p_include_O_M_site;
'p_protein_production_noise_level', p_protein_production_noise_level;
'n_time_steps', n_time_steps;
'threshold_time_switch_decision', threshold_time_switch_decision;
'threshold_winning_ratio_MORtoCI', threshold_winning_ratio_MORtoCI;
'kbT', kbT;
'K_or', K_or;
'K_ol', K_ol;
'K_od', K_od;
'K_om', K_om;
'G_or', G_or;
'G_ol', G_ol;
'G_od', G_od;
'G_om', G_om;
'G_pair', G_pair;
'G_trip', G_trip;
'K_d', K_d;
'K_CIMOR', K_CIMOR;
'T_generation', T_generation;
'N_CI_total', N_CI_total;
'n_1', n_1;
'r_CI', r_CI;
'tau_CI', tau_CI;
'n_2', n_2;
'N_MOR_total', N_MOR_total;
'n_3', n_3;
'r_MOR', r_MOR;
'tau_MOR', tau_MOR;
'n_4', n_4
};

%% Concatenating cells

model_parameters_cell = vertcat({'MODE', MODE}, model_parameters_MODE_SPECIFIC_cell, model_parameters_GENERAL_cell);

%% Writing out model parameters

model_parameters_table = cell2table(model_parameters_cell,'VariableNames',{'Parameter','Value'});


% writetable determines the file format based on the specified extension. The extension must be one of the following:
% 	.txt, .dat, or .csv for delimited text files
% 	.xls, .xlsm, or .xlsx for Excel spreadsheet files
% 	.xlsb for Excel spreadsheet files supported on systems with Excel for Windows

file_model_parameters = fullfile(dir_output_path, 'model_parameters.txt'); % filename | defaultdelimiter is comma IF the FILE EXTENSION is .txt
writetable(model_parameters_table, file_model_parameters, 'Delimiter', 'tab') % tab seperated


display(sprintf('Wrote model parameters to file: %s', file_model_parameters))

%% ----------------------------- OLD STUFF --------------------------------------------
 

% model_parameters_name = {'MODE';
% 'range_promoter_strength_ratio_PL2PR';
% 'CI_initial_concentration';
% 'MOR_initial_concentration';
% 'p_display_real_time_plots';
% 'p_quantify_variance';
% 'n_simulations';
% 'p_include_O_M_site';
% 'p_protein_production_noise_level';
% 'n_time_steps';
% 'threshold_time_switch_decision';
% 'threshold_winning_ratio_MORtoCI';
% 'kbT';
% 'K_or';
% 'K_ol';
% 'K_od';
% 'K_om';
% 'G_or';
% 'G_ol';
% 'G_od';
% 'G_om';
% 'G_pair';
% 'G_trip';
% 'K_d';
% 'K_CIMOR';
% 'T_generation';
% 'N_CI_total';
% 'n_1';
% 'r_CI';
% 'tau_CI';
% 'n_2';
% 'N_MOR_total';
% 'n_3';
% 'r_MOR';
% 'tau_MOR';
% 'n_4'};

%T = table(model_parameters_value,'RowNames',model_parameters_name)




% fid = fopen(myfilepath, 'wt');
% fprintf(fid, '% this is the beginning of my file');
% fprintf(fid, 'var1 = %d ',evalin('base',var1));  % if var1 is double
% fclose(fid);

 
 
% save(file_model_parameters, ...
%     'MODE',...
%     'range_promoter_strength_ratio_PL2PR',...
%     'CI_initial_concentration',...
%     'MOR_initial_concentration',...
%     'p_display_real_time_plots',...
%     'p_quantify_variance',...
%     'n_simulations',...
%     'p_include_O_M_site',...
%     'p_protein_production_noise_level',...
%     'n_time_steps',...
%     'threshold_time_switch_decision',...
%     'threshold_winning_ratio_MORtoCI',...
%     'kbT',...
%     'K_or',...
%     'K_ol',...
%     'K_od',...
%     'K_om',...
%     'G_or',...
%     'G_ol',...
%     'G_od',...
%     'G_om',...
%     'G_pair',...
%     'G_trip',...
%     'K_d',...
%     'K_CIMOR',...
%     'T_generation',...
%     'N_CI_total',...
%     'n_1',...
%     'r_CI',...
%     'tau_CI',...
%     'n_2',...
%     'N_MOR_total',...
%     'n_3',...
%     'r_MOR',...
%     'tau_MOR',...
%     'n_4', '-ascii')



