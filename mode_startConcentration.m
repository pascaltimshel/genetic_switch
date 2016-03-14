% ------------------------------------------------------------------------------------
%% ----------------------------- MODE SPECIFIC PARAMETERS ----------------------------
% ------------------------------------------------------------------------------------

% range_CI_initial_concentration = [XXX];
% range_MOR_initial_concentration = [XXX];

%filename_startConcentration_export = fullfile(dir_output_path, 'data_export_startConcentration.txt');


%% ----------------------------- INITIALIZE ----------------------------

% ...
e_msg = sprintf('MODE *%s* not implemented yet', mfilename);
error(e_msg)




%% MAIN LOOP | Initial concentrations (begyndelsesbetingelser)

% %Here is an example of how to loop over CI/MOR initial concentrations
% %for MOR_initial_concentration=range_MOR_initial_concentration
% for CI_initial_concentration=range_CI_initial_concentration
%     tic
%     % RUN SIMULATION(S)
%     % simulatePhageLifeCycle_GillespieAlgo(..., CI_initial_concentration,...)
%     % ...
%     toc
% end
