function TP901_running %[win_ratio_pl]=

%% DEPENDENCIES
% This script was tested and developed for the following MATLAB versions:
% MATLAB 2014b (area plot function will not work if running this version)
% MATLAB 2015b (everything should work)


%% Remarks
% OBS: This script was developed in September 2015 (4. september), and is a modified version of "TP901_work d. 3. april - CLEAN god version.m".

%% TODO - improvements
%1) Collect the normalization constants into variables, to make the code
%easier to change
%2) 

%% PARAMETERS TO ADJUST
% NB: "true" is same as "logical(1)". "false" is the same as logical(0).

% Inclusion of the O_M site in the model?
% If the O_M is included in the model, the CI:MOR complex can repress PR (*CHECK THIS*)
p_include_O_M_site = true; % default is true

% Set the noise level in the model
% The value MUST be greater than zero. 
% Recommended range is [0.1-5.0]. 
% A high value results in more stochasticity (noise) in the model
p_protein_production_noise_level = 0.5; % default is 0.5

% The range of values for the RELATIVE promoter strenght between PL and PR
% E.g. a value of 0.5 means that P_L is 50% of the strength of P_R.
% E.g. a value of 1.2 means that P_L is 120% of the strength of P_R.
range_promoter_strength_ratio_PL2PR = [0:0.1:1.2]; % default is [0:0.1:1.2], i.e. 10%, 20%, ..., 120%

% number of simulations for each "setting" (e.g. number of simulations for each P_L promoter activity)
n_simulations = 10;

% Number of time steps in simulation. 
% You can etimate the "real time" simulation from the number of time steps.
% NB: n_time_steps should be high enough to reach the time for the decision ("threshold_time_switch_decision").
n_time_steps=1000; % default is 1000 steps

% Threshold for when the decision is made and the lytic/lysogen cyclus is determined
threshold_time_switch_decision = 15*60; %[seconds] default is 15min*60sec/min=15 min. 

% Initial values for CI and MOR
CI_initial_concentration = 0; % default is 0 initial concentration
MOR_initial_concentration = 0; % default is 0 initial concentration

% Additional options
p_display_real_time_plots = true; % if true, then the script will display plots with the time in estimated "real time" instead of time steps.
p_quantify_variance = true; % if true, then report the variance ("noise") of the CI and MOR protein levels. recommended (this is easily calculated)


%% Physical constants
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
global K_CIMOR
K_CIMOR = (10^-1*10^-9);  % Bindingskonstant mellem CI og MOR [M]


%% Produktionskonstanter - Gillespie fra A.A
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

%% MAIN LOOP | Initial concentrations (begyndelsesbetingelser)

% %Here is an example of how to loop over CI/MOR initial concentrations
% %for MOR_initial_concentration=range_MOR_initial_concentration
% for CI_initial_concentration=range_CI_initial_concentration
%     tic
%     % RUN SIMULATION(S)
%     % run_simulation(..., CI_initial_concentration,...)
%     % ...
%     toc
% end

%% MAIN LOOP | PROMOTER STRENGTH

PL_win=zeros(1,length(range_promoter_strength_ratio_PL2PR));
PR_win=zeros(1,length(range_promoter_strength_ratio_PL2PR));
%-----------------------------START SIMULERING----------------------------
tic_value1=tic; % Measure the runtime
for promoter_strength_ratio_PL2PR=range_promoter_strength_ratio_PL2PR
    tic
    for simulation=1:n_simulations
        display(['Simulation:', num2str(simulation)])
        
        [time_steps, time_real, ...
        CI_total, MOR_total, ...
        CI_m, CI_d, MOR_m, CIMOR, ...
        FORHOLD_MOR_CI, ...
        PR_activity, PL_activity, ...
        kritisk_tid_trin, kritisk_tid_real, ...
        PL_win, PR_win] = run_simulation(n_time_steps, threshold_time_switch_decision, CI_initial_concentration, MOR_initial_concentration, promoter_strength_ratio_PL2PR, p_include_O_M_site);


    % DISPLAYER KRITISKE TIDER:
    %display(['Kritisk_tid_trin (step):', num2str(kritisk_tid_trin)])
    %display(['Kritisk_tid_real (s):', num2str(kritisk_tid_real)])
    
    end% End for-loop (simulations).
    toc
end
%-----------------------------END FOR SIMULERING----------------------------

%% Converting to units of nano molar (nM)
CI_total_array(:,:)=CI_total_array(:,:)*10^9;
MOR_total_array(:,:)=MOR_total_array(:,:)*10^9;
CI_m(:)=CI_m(:)*10^9;
CI_d(:)=CI_d(:)*10^9;
MOR_m(:)=MOR_m(:)*10^9;
CIMOR(:)=CIMOR(:)*10^9;
CI_total(:)=CI_total(:)*10^9;
MOR_total(:)=MOR_total(:)*10^9;

%% Variance quantification
% OBS: this calculates the standard deviation for the *LAST simulation*.
if p_quantify_variance == 1
    boundery = 0.5*n_time_steps; %2*10^4;
    std_CI = std(CI_total_array([boundery:end],end),0,1);
    std_MOR = std(MOR_total_array([boundery:end],end),0,1);
    mean_CI = mean(CI_total_array([boundery:end],end));
    mean_MOR = mean(MOR_total_array([boundery:end],end));
    display(['Standard deviation for the last "simulation" for CI_total was: ', num2str(std_CI), ...
            ' (', num2str(std_CI/mean_CI*100), '% of the mean_CI)'])
    display(['Standard deviation for the last "simulation" for MOR_total was: ', num2str(std_MOR), ...
            ' (', num2str(std_MOR/mean_MOR*100), '% of the mean_MOR)'])
end % end if for bestemmelse af varians



%% Lytisk/lysogen main plot
win_ratio_pl=PL_win(:)./(PL_win(:)+PR_win(:));
%{
figure('Name','MAINPLOT: %lytiske cykler','NumberTitle','off')
plot(CI_concentations(:),win_ratio_pl(:),'--bs','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','b',...
                'MarkerSize',5);
xlabel('CI startkoncentration (nM)');
ylabel('% lytiske cykler');
axis([0 length(CI_concentations) 0 1]);
%}

disp(['The total time was approximately:' num2str(toc(tic_value1))]);
disp(['PR er lig: ' num2str(PR_win) '. PL er lig: ' num2str(PL_win) '']);
disp(['Fordelling er lig: ' num2str(win_ratio_pl) '']);


%% PLOT
% PLOT - NYT d. 15/03/2012 - Pascal
figure('Name','CI[total] vs MOR[total] plot','NumberTitle','off')

% Area plot
p_area = area([1:1:1.5*N_MOR_total], [1:1:1.5*N_MOR_total]);
set(p_area,'FaceColor',[1,0.6,0.6],'LineStyle', 'none'); % 'FaceAlpha' does not work for MATLAB2014b
hold on;

cmap = hsv(n_simulations); % loading colormap
for l=1:n_simulations
    h(l)=plot(CI_total_array(:,l),MOR_total_array(:,l), 'Color', cmap(l,:));
    hold on
    s{l}=sprintf('%s%g','simulation=',l); % defining title for legend
end
legend(h, s) % setting legend

%Plotter punkt for kritisk tid (15 min time_real)
for l=1:n_simulations
    plot(CI_total_array(kritisk_tid_trin(l),l),MOR_total_array(kritisk_tid_trin(l),l), '*', 'MarkerSize', 10, 'MarkerFaceColor', cmap(l,:))
    hold on
end
% Plot af ret linje til indikation af forhold mellem CI og MOR
plot([1:1:1.5*N_MOR_total],[1:1:1.5*N_MOR_total], 'Linewidth', 2); % Prøv med scaling constrained, og plot en ret linje
hold off; 
% Labels
xlabel('CItotal (nM)');
ylabel('MORtotal (nM)');

%% PLOT
figure('Name','Protein concentrations','NumberTitle','off')
plot(time_steps,CI_m(1,:),'g','LineWidth',1)
hold on
plot(time_steps,CI_d(1,:),'c','LineWidth',1)
hold on
plot(time_steps,MOR_m(1,:),'k','LineWidth',1)
hold on
plot(time_steps,CIMOR(1,:),'m','LineWidth',1)
hold on
plot(time_steps,CI_total(1,:),'b','LineWidth',1)
hold on
plot(time_steps,MOR_total(1,:),'r','LineWidth',1)
hold on
hold off
xlabel('Time step');
ylabel('Concentration (nM)');
legend('CI_m','CI_d','MOR_m','CI:MOR', 'CI total', 'MOR total')
title('Protein concentrations')

%% PLOT
figure('Name','Promotor activity probability','NumberTitle','off')
plot(time_steps,PL_activity(1,:),'r','LineWidth',0.01)
hold on
plot(time_steps,PR_activity(1,:),'b','LineWidth',0.01)
xlabel('Time step');
ylabel('Promotor activity probability');
legend('PL', 'PR')
axis([1 (length(time_steps)) 0 1]);
title('Promotor activity probability');


%% PLOTS | Real time 
if p_display_real_time_plots == true
    %% PLOT PROMOTER ACTIVITY
    figure('Name','Promotor activity probability (real time)','NumberTitle','off')
    plot(time_real(1,:),PL_activity(1,:),'r','LineWidth',0.01)
    hold on
    plot(time_real(1,:),PR_activity(1,:),'b','LineWidth',0.01)
    xlabel('Tid (s)');
    ylabel('Promotor activity probability');
    legend('PL', 'PR')
    axis([1 time_real(1,(length(time_steps))) 0 1]);
    title('Promotor activity probability');

    % Forhold mellem CI og MOR plot
    figure('Name','Forhold mellem MOR og CI','NumberTitle','off')
    plot(time_real(:),FORHOLD_MOR_CI(1,:),'m','LineWidth',2);
    xlabel('Tid (s)');
    ylabel('Forhold mellem MOR_t og CI_t');

    %% PLOT PROTEINS
    figure('Name','Protein concentrations i realtid','NumberTitle','off')
    plot(time_real(1,:),CI_m(1,:),'g','LineWidth',1)
    hold on
    plot(time_real(1,:),CI_d(1,:),'c','LineWidth',1)
    hold on
    plot(time_real(1,:),MOR_m(1,:),'k','LineWidth',1)
    hold on
    plot(time_real(1,:),CIMOR(1,:),'m','LineWidth',1)
    hold on
    plot(time_real(1,:),CI_total(1,:),'b','LineWidth',1)
    hold on
    plot(time_real(1,:),MOR_total(1,:),'r','LineWidth',1)
    hold on
    % PLOTTING VERTICAL LINE d. 2. april
    plot([T_generation, T_generation], [100, 500], 'r', 'LineWidth',2);
    hold on
    plot([time_real(boundery), time_real(boundery)], [100, 500], 'r', 'LineWidth',2);
    hold on
    plot([time_real(boundery), time_real(end)], [mean_CI, mean_CI], 'r', 'LineWidth',2);
    %time_steps(end)
    hold off
    xlabel('Tid (s)');
    ylabel('Concentration (nM)');
    legend('CI_m','CI_d','MOR_m','CI:MOR', 'CI total', 'MOR total')
    %set(gca, 'XTick', [1:delta_t:n_time_steps*60*delta_t]);
    %axis([1 n_time_steps 0 10^3])
    title('Protein concentrations')
end

end %End Function


