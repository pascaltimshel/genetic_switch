
% ------------------------------------------------------------------------------------
%% ----------------------------- MODE SPECIFIC PARAMETERS ----------------------------
% ------------------------------------------------------------------------------------
% The range of values for the RELATIVE promoter strenght between PL and PR
% E.g. a value of 0.5 means that P_L is 50% of the strength of P_R.
% E.g. a value of 1.2 means that P_L is 120% of the strength of P_R.
range_PSR_PL2PR = [0:0.1:1.4]; % default is [0:0.1:1.4], i.e. 10%, 20%, ..., 140%

filename_PSR_export = fullfile(dir_output_path, 'data_export_PSR.txt');


%% ----------------------------- INITIALIZE ----------------------------

PL_win_array=zeros(1,length(range_PSR_PL2PR));
PR_win_array=zeros(1,length(range_PSR_PL2PR));

%% ----------------------------- START SIMULATION LOOP ----------------------------
tic_value1=tic; % Measure the runtime
for i=1:length(range_PSR_PL2PR)
%parfor i=1:length(range_PSR_PL2PR) % parallel loop | TODO: add PARALLEL loop
    PSR_PL2PR = range_PSR_PL2PR(i);
    display(sprintf(['Simulating P_L/P_R Promoter Strength Ratio (PSR): PSR=', num2str(PSR_PL2PR)]))
    ticTime_PSRsimulation = tic;
    for simulation=1:n_simulations
        display(sprintf(['\tSimulation #', num2str(simulation)]))
        
        [time_steps, time_real, ...
        CI_total, MOR_total, ...
        CI_m, CI_d, MOR_m, CIMOR, ...
        RATIO_MORtoCI, ...
        PR_activity, PL_activity, ...
        decision_time_step, decision_time_real, ...
        PL_win, PR_win] = simulatePhageLifeCycle_GillespieAlgo(n_time_steps, threshold_time_switch_decision, threshold_winning_ratio_MORtoCI, CI_initial_concentration, MOR_initial_concentration, PSR_PL2PR, p_include_O_M_site);
        
        PL_win_array(i) = PL_win_array(i)+PL_win; % incrementing by PL_win (0 or 1, depending on what winning state)
        PR_win_array(i) = PR_win_array(i)+PR_win; % incrementing by PR_win (0 or 1, depending on what winning state)
    
    end% End for-loop (simulations).
    
    % Displaying lytic win ratio
    lytic_win_ratio=PL_win_array(i)/(PL_win_array(i)+PR_win_array(i));
    display(sprintf('Lytic (anti-immune) win ratio: %.2f', lytic_win_ratio))
    
    timeElapsed_PSRsimlation = toc(ticTime_PSRsimulation);
    display(sprintf('PSR simulation elapsed time: %.2f seconds.', timeElapsed_PSRsimlation))
end
disp(['Total simulation time:' num2str(toc(tic_value1))]);
%-----------------------------END FOR SIMULATION ----------------------------

%% POST processing of simulation
% PL wins --> Lytic life cycle (blue colony) | Anti-immune
% PR wins --> Lysogen life cycle (white colony) | Immune state
lytic_switch_frq_array=PL_win_array(:)./(PL_win_array(:)+PR_win_array(:)); %% Calculating frequencies of PL_win

%% EXPORT results to excel

% OBS: transposing some of the elements
tmp_tbl_export = table(range_PSR_PL2PR', round(lytic_switch_frq_array*100,4), PL_win_array', PR_win_array',...
    'VariableNames', {'PSR' 'Switch_frequency_pct' 'P_L_wins' 'P_R_wins'});

if save_output
    writetable(tmp_tbl_export, filename_PSR_export, 'FileType', 'text', 'Delimiter', 'tab')
    %writetable(tmp_tbl_export, filename_PSR_export, 'FileType', 'spreadsheet', 'Sheet', 'PSR_export')
        % --> Pascal OSX El Capitain got the error: "Could not start Excel server for export.  Export to a text file instead."
end

display(sprintf('Wrote PSR mode data export to file: %s', filename_PSR_export))

%% ===================== Lytic/lysogen (anti-immune/immune) main PLOT ============================

% Load experimental data file
[expData_PL_relative_strength, expData_switching_frq] = readExperimentalData(file_experimental_data); % load experimental data
    % column 1: mean values
    % column 2: std error

% ======= Initial plot setup =======
fig_lytic_switch_frequency = figure('Name','Anti-immune switch frequency','NumberTitle','off');

% REF: Cubic Smoothing Splines | http://www.mathworks.com/help/curvefit/examples/cubic-smoothing-splines.html
    % spaps(x,y,tol)
    % csaps(x,y,p) |  for p = 0, f is the least-squares straight line fit to the data,
                    % for p = 1, f is the variational, or `natural' cubic spline interpolant.
p_smoothing = 0.999; % smoothing parameter
XX_points_to_evaluate_spline = min(range_PSR_PL2PR):0.01:max(range_PSR_PL2PR); % X range to evalute the smoothing spline.

% Using new MATALB 2014b color codes: http://hoeckerson.de/notes/2016/01/default-matlab-plot-colors/
pplot_color_experimental = [0, 0.4470, 0.7410]; % Color code for 'blue'
pplot_color_simulated = [0.8500, 0.3250, 0.0980]; % Color code for 'orange'

%% ======= EXPERIMENTATAL DATA =======                   
% Y errorbar
h_expData_y_errorbar = errorbar(expData_PL_relative_strength(:,1), ... 
     expData_switching_frq(:,1), ...
     expData_switching_frq(:,2));
hold on;
% X errorbar - function downloaded from MATLAB File Exchange
h_expData_x_errorbar = herrorbar(expData_PL_relative_strength(:,1), ... 
     expData_switching_frq(:,1), ...
     expData_PL_relative_strength(:,2));
hold on;
% Spline
tmp_inverse_std_error_switching_frq = 1./expData_switching_frq(:,2); % Inverse of standard error
weights_error_measure = tmp_inverse_std_error_switching_frq/min(tmp_inverse_std_error_switching_frq); % Scale the weights.
    % POINTS WITH LARGE WEIGHTS: CLOSE fit to points
    % POINTS WITH SMALL WEIGHTS: 'Loose' fit to points
    % By dividing with min(..) we obtain error weights in the interval [1; inf]
    % ^^ *OBS*: there is no correct way of choosing the error weights
cs_expData_values = csaps(expData_PL_relative_strength(:,1), expData_switching_frq(:,1), p_smoothing, XX_points_to_evaluate_spline, weights_error_measure); % fit spline
h_expData_spline = plot(XX_points_to_evaluate_spline, cs_expData_values, 'LineWidth', 2);
hold on;

% Property Values common to all experimental plots
set([h_expData_y_errorbar, h_expData_x_errorbar, h_expData_spline], ...
    'Color', pplot_color_experimental)

set([h_expData_y_errorbar], ...
    'LineStyle', 'none',...
    'Marker', 'o')
 

%% ======= SIMULATED DATA =======
h_simData = plot(range_PSR_PL2PR, lytic_switch_frq_array);
hold on;
cs_simData_values = csaps(range_PSR_PL2PR, lytic_switch_frq_array, p_smoothing, XX_points_to_evaluate_spline); % fit spline
    % *OBS*: currently no error weights is used for the simulated data

h_simData_spline = plot(XX_points_to_evaluate_spline, cs_simData_values, 'LineWidth', 2);
hold on;

% Property Values common to all simulated plots
set([h_simData, h_simData_spline], ...
    'Color', pplot_color_simulated)

set([h_simData], ...
    'LineStyle', 'none',...
    'Marker', 'o')
 


%% ======= Plot adjustments =======
axis([0 max(range_PSR_PL2PR) 0 inf]) % set axis limits
% Axis labels:
xlabel('P_L/P_R relative promoter strength');
ylabel('Anti-immune switch frequency'); % ylabel('Lytic switch frequency');
% Legend
h_legend = legend([h_expData_spline, h_simData_spline], 'Experimental data', 'Simulated data');


%% ------------------------------ SAVE FIGURE ----------------------------
if save_output
    % Save figure
    filename_fig = fullfile(dir_output_path, 'fig_lytic_switch_frequency'); % If the file name does not include an extension, then 'print' appends the appropriate one.
    print(fig_lytic_switch_frequency, filename_fig, '-dpng') % png
    print(fig_lytic_switch_frequency, filename_fig, '-dpdf') %
    %print(fig_lytic_switch_frequency, filename_fig, '-dsvg') %
end



