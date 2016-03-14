%% Master script for modelling the TP901-1 genetic switch

%% DEPENDENCIES
% This script was tested and developed for the following MATLAB versions:
% MATLAB 2014b (area plot function will not work if running this version)
% MATLAB 2015b (everything should work)

%% Remarks
% The development of this script started September 2015 (September 4th)
% The script is a RESTRUCTERED and modified version of "TP901_work d. 3. april - CLEAN god version.m".

%% TODO - improvements
%1) Collect the normalization constants into variables, to make the code easier to change
%2) Speed-up code: parallel loops

%x) Clean up: clear variable when they are no longer needed
%x) Translate EVERYTHING to English
%x) Export more plots...

%% Initialization
clear all; % clean workspace
clc;

%% Paths and Data

addpath('Library') % path for additional functions (e.g. MATLAB File Exchange)

% Experimental data file to load
file_experimental_data = 'experimental_data_v1.xlsx';

%% ================================= SCRIPT MODE =================================

MODE = 'mode_modelInspection'; % provides an interface to inspecting the model (e.g. promoter activity plots). 
    % Requires full specification of the model.

%MODE = 'mode_PSR'; % sweeps the P_L/P_R promoter strength ratio

%MODE = 'mode_startConcentration'; % sweeps initial start concentrations

% MODE = 'mode-PLACEHOLDER'; % placeholder for sweeping other parameters (e.g. KOR).

dir_main_output = 'model_output'; % top level directory for output | default is 'model_output'
save_output = true; % if true, then the output (plots and parameters will be saved)
%verbose_script = true

%% ============================ PARAMETERS TO ADJUST ===============================
% NB: "true" is same as "logical(1)". "false" is the same as logical(0).

% Promoter strength ratio (PL/PR ratio)
PSR_PL2PR = 1;

% Initial values for CI and MOR
CI_initial_concentration = 0; % default is 0 initial concentration
MOR_initial_concentration = 0; % default is 0 initial concentration


% Number of simulations for each "setting" (e.g. number of simulations for each P_L promoter activity)
% The more simulations, the better the "parameter estimate", e.g. the switching frequency for a given P_L activity
%n_simulations = 25;
n_simulations = 10;


% Inclusion of the O_M site in the model?
% If the O_M is included in the model, the CI:MOR complex can repress PR (*CHECK THIS*)
p_include_O_M_site = true; % default is true
% p_include_O_M_site = false; % default is true

% Set the noise level in the model
% The value MUST be greater than zero. 
% Recommended range is [0.1-5.0]. 
% A high value results in more stochasticity (noise) in the model
p_protein_production_noise_level = 0.5; % default is 0.5

% Number of time steps in simulation. 
% You can etimate the "real time" simulation from the number of time steps.
% NB: n_time_steps should be high enough to reach the time for the decision ("threshold_time_switch_decision").
n_time_steps=1000; % default is 1000 steps

% Threshold for when the decision is made and the lytic/lysogen cyclus is determined
threshold_time_switch_decision = 15*60; %[seconds] default is 15min*60sec/min=15 min. 
threshold_winning_ratio_MORtoCI = 0.8; % PL "wins" if RATIO_MORtoCI(t) >= threshold_winning_ratio_MORtoCI



%% ================================= Physical constants =================================
T = 298.15; % [K]
R = 8.3144621/(4.184*10^3); % [kcal/(K*mol)]
global kbT
kbT = R*T; % [J/mol] DETTE ER FAKTISK R*T (så dårligt variabelnavn)

% BINDINGSKONSTANTER (bestemt eksperimentelt) MOLÆR
K_or = 2000*10^-9; %[M] NORMALT: 2000
K_ol = 28*10^-9; %[M]
K_od = 28*10^-9; %[M]
K_om = 10*10^-9; %[M] DENNE ER IKKE VIST EKSP.

% BINDINGSENERGIER (i JOULE eller [kcal/mol]?)
% DO NOT CHANGE ME - CHANGE K_or, K_ol, K_od, K_om instead!!
global G_or
G_or = -kbT*log(1/K_or);
global G_ol
G_ol = -kbT*log(1/K_ol);
global G_od
G_od = -kbT*log(1/K_od);
global G_om
G_om = -kbT*log(1/K_om);

% COOPERATIVE BINDING ENERGIES
% G_pair is the "additional" binding energy associated with binding of *two* CI dimers
% G_trip is the "additional" binding energy associated with binding of *three* CI dimers
global G_pair
G_pair = -20000/(4.184*10^3); % [kcal/mol] | G_pair numerator calculated to 20000 (*unsure about the calculation*)
global G_trip
G_trip = -30000/(4.184*10^3); % [kcal/mol] | G_pair numerator calculated to 30000 (*unsure about the calculation*)
% Try also -10000 and -15000 in the numerator for G_pair and G_trip respectively
% ??En faktor 10 ændring i G_kooperativtet giver en faktor ca. 1200 i vægten??

% Protein complex dissociation constants 
global K_d
K_d = (20*10^-9); % Bindingskonstant mellem CI monomerer [M]
% K_d = (20*10^-2); % Decent fit; Sigmoidal?
% K_d = (20*10^-50); % ---> All immune
global K_CIMOR
K_CIMOR = (10^-1*10^-9);  % Bindingskonstant mellem CI og MOR [M]


%% ====================== Production constants - Gillespie fra A.A ======================
T_generation=3600; % [seconds] Cell generation time. 
% Lactococcus lactis is ~ 52 min
% E.coli is ~20 min

%PARAMETRE FOR CI-PRODUKTION
%r_CItrsc=0.22;% Transcriptionsrate, Egan(2004). Bruges til at udregne raten, a_1=a_CIprod
N_CI_total=200; % Steady state "concentrations"
global n_1
n_1=round(60*p_protein_production_noise_level); % the value 30 seems good. choose values between [10-300]
global r_CI
r_CI=N_CI_total/(T_generation*n_1); % Skalering r_CItrsc/F for at opfylde lignignen N_total/T_generation ca lig med n_CIprod*r_CI
global tau_CI
tau_CI=1/3000; % [s^-1] rate for nedbrydning
global n_2
n_2=1; % Ændring i CI_total pr. nedbrydningsevent

%PARAMETRE FOR MOR-PRODUKTION
%r_MORtrsc=0.22;% Transcriptionsrate, {0.1/s;1/s} Bruges til at udregne raten, a_3=a_MORprod
N_MOR_total=200; % steady state "concentrations" [old value was 800]
global n_3
n_3=round(60*p_protein_production_noise_level); % the value 30 seems good. choose values between [10-300]
global r_MOR
r_MOR=N_MOR_total/(T_generation*n_3); % Skalering r_CItrsc/F for at opfylde lignignen N_total/T_generation ca lig med n_CIprod*r_CI
global tau_MOR
tau_MOR=1/3000;
global n_4
n_4=1; % Ændring i MOR_total pr. nedbrydningsevent


%% ====================== Create directory to save modelling results ========================
% We create an output directory that each MODE can write its output to.

if save_output
    % Create dir_main_output if it does not exist
    if exist(dir_main_output, 'dir') == 0
        mkdir(dir_main_output)
        display(sprintf('Created dir_main_output in current working directory: %s', dir_main_output))
    end

    % Create NAME for directory for current model output
    timestamp = datestr(now,'yyyy-mm-dd_HH.MM.SS'); % try also: datetime('now','TimeZone','local','Format','yyyy-MM-dd_HH:mm:ss')
    dir_output_name = sprintf('%s_%s', timestamp, MODE);
    dir_output_path = fullfile(dir_main_output, dir_output_name); % create path specification
    % Create directory
    mkdir(dir_output_path);
    display(sprintf('Created output directory: %s', dir_output_path))
end

%% ============================ Master switch ============================

try
   switch MODE
    case 'mode_PSR'
        run('mode_PSR.m')
    case 'mode_startConcentration'
        run('mode_startConcentration.m')
    case 'mode_modelInspection'
        run('mode_modelInspection.m')
    otherwise
        error('Received unexpected MODE variable. Change the MODE variable to a valid value.')
    end
catch ME % Something bad happend - we will clean up the directory we created...
   pause(0.5) % quick sleep to make sure things have "settled"
   if isEmptyDirectory(dir_output_path)
       rmdir(dir_output_path)
       display(sprintf('Removed the dir_output_path (%s) because it was empty.', dir_output_path))
   end
   display('Uuups, something unexpected happend and an error occured. Will rethrow expection.')
   rethrow(ME)
end



%% 

if save_output
    % Write model parameters
    run('writeModelParameters.m')
end



