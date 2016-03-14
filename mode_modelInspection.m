
% ------------------------------------------------------------------------------------
%% ----------------------------- MODE SPECIFIC PARAMETERS ----------------------------
% ------------------------------------------------------------------------------------


% Selected simulation for INSPECTION AND PLOTTING:
idxSim = n_simulations; %using last simulation 
%idxSim = randi([1 n_simulations]); % random simulation


p_display_real_time_plots = true; % if true, then the script will display plots with the time in estimated "real time" instead of time steps.


% When reporting the variance ("noise") of the CI and MOR protein levels, you need to specify the TIME STEP THAT defines the "steady-state" used for calculation.
variance_calc_boundary = 0.5*n_time_steps; 
    % ^ *IMPORTANT* the concentrations in the last half of the simulation is used to calculate the variance


%filename_modelInspection_export = fullfile(dir_output_path, 'data_export_modelInspection.txt');



%% ----------------------------- START SIMULATION LOOP ----------------------------

for simulation=1:n_simulations
    display(sprintf(['\tSimulation #', num2str(simulation)]))

    [time_steps(simulation,:), time_real(simulation,:), ...
    CI_total(simulation,:), MOR_total(simulation,:), ...
    CI_m(simulation,:), CI_d(simulation,:), MOR_m(simulation,:), CIMOR(simulation,:), ...
    RATIO_MORtoCI(simulation,:), ...
    PR_activity(simulation,:), PL_activity(simulation,:), ...
    decision_time_step(simulation), ... % "step number" (integer)
    decision_time_real(simulation), ... % in units of seconds
    PL_win(simulation), PR_win(simulation)] = simulatePhageLifeCycle_GillespieAlgo(n_time_steps, threshold_time_switch_decision, threshold_winning_ratio_MORtoCI, CI_initial_concentration, MOR_initial_concentration, PSR_PL2PR, p_include_O_M_site);

end% End for-loop (simulations).


%% Converting to units of nano molar (nM)

CI_total=CI_total*10^9;
MOR_total=MOR_total*10^9;
CI_m=CI_m*10^9;
CI_d=CI_d*10^9;
MOR_m=MOR_m*10^9;
CIMOR=CIMOR*10^9;

%% "Steady-state" variance quantification - used to quantify the "noise level" in the simulation
% OBS: this calculates the standard deviation for the *LAST simulation*.
std_CI = std(CI_total(idxSim, [variance_calc_boundary:end]),0,2); % row-wise std
std_MOR = std(MOR_total(idxSim, [variance_calc_boundary:end]),0,2); % row-wise std
mean_CI = mean(CI_total(idxSim, [variance_calc_boundary:end]), 2); % row-wise mean
mean_MOR = mean(MOR_total(idxSim, [variance_calc_boundary:end]), 2);  % row-wise mean
display(sprintf('Selected simulation #%d/#%d | Standard deviation for CI_total was: %s (%.2f %% of the mean_CI)', idxSim, n_simulations, std_CI, std_CI/mean_CI*100));
display(sprintf('Selected simulation #%d/#%d | Standard deviation for MOR_total was: %s (%.2f %% of the mean_MOR)', idxSim, n_simulations, std_MOR, std_MOR/mean_MOR*100));



%% ======================== PLOT - CI vs MOR trajectories ========================
% Written d. 15/03/2012 - Pascal
% This plot uses all the simulations

figure('Name','CI[total] vs MOR[total] plot','NumberTitle','off')

% Area plot
p_area = area([1:1:1.5*N_MOR_total], [1:1:1.5*N_MOR_total]);
set(p_area,'FaceColor',[1,0.6,0.6],'LineStyle', 'none'); % 'FaceAlpha' does not work for MATLAB2014b
hold on;

cmap = hsv(n_simulations); % loading colormap
for l=1:n_simulations
    h(l)=plot(CI_total(l,:),MOR_total(l,:), 'Color', cmap(l,:));
    hold on
    s{l}=sprintf('%s%g','simulation=',l); % defining title for legend
end
legend(h, s) % setting legend

%Plotter punkt for kritisk tid (15 min time_real)
for l=1:n_simulations
    plot(CI_total(l,decision_time_step(l)),MOR_total(l,decision_time_step(l)), '*', 'MarkerSize', 10, 'MarkerFaceColor', cmap(l,:))
    hold on
end
% Plot af ret linje til indikation af forhold mellem CI og MOR
plot([1:1:1.5*N_MOR_total],[1:1:1.5*N_MOR_total], 'Linewidth', 2); % Prøv med scaling constrained, og plot en ret linje
hold off; 
% Labels
xlabel('CItotal (nM)');
ylabel('MORtotal (nM)');


%% =========================== Time Step plots ================================
figure('Name',sprintf('Sim #%d (time steps)',idxSim),'NumberTitle','off')
subplot(2,1,1)
%figure('Name','Protein concentrations','NumberTitle','off')
h_prot(1) = plot(time_steps(idxSim,:),CI_m(idxSim,:),'g','LineWidth',1);
hold on;
h_prot(2) = plot(time_steps(idxSim,:),CI_d(idxSim,:),'c','LineWidth',1);
hold on;
h_prot(3) = plot(time_steps(idxSim,:),MOR_m(idxSim,:),'k','LineWidth',1);
hold on;
h_prot(4) = plot(time_steps(idxSim,:),CIMOR(idxSim,:),'m','LineWidth',1);
hold on;
h_prot(5) = plot(time_steps(idxSim,:),CI_total(idxSim,:),'b','LineWidth',1);
hold on;
h_prot(6) = plot(time_steps(idxSim,:),MOR_total(idxSim,:),'r','LineWidth',1);
hold on;
% Annotation 
h_anno(1) = plot(decision_time_step(idxSim),CI_total(idxSim,decision_time_step(idxSim)),'b*');
hold on;
h_anno(2) = plot(decision_time_step(idxSim),MOR_total(idxSim,decision_time_step(idxSim)),'r*');
hold off;
xlabel('Time step');
ylabel('Concentration (nM)');
legend(h_prot, 'CI_m','CI_d','MOR_m','CI:MOR', 'CI total', 'MOR total')
title('Protein concentrations')

subplot(2,1,2)
% figure('Name','Promotor activity probability','NumberTitle','off')
h_promoterAct(1)=plot(time_steps(idxSim,:),PL_activity(idxSim,:),'r','LineWidth',1);
hold on;
h_promoterAct(2)=plot(time_steps(idxSim,:),PR_activity(idxSim,:),'b','LineWidth',1);
hold off
xlabel('Time step');
ylabel('Promoter activity probability');
legend(h_promoterAct, 'PL', 'PR')
axis([0 inf 0 inf]);
title('Promoter activity probability');


%% ======================== PLOTS | Real time ============================
if p_display_real_time_plots
    figure('Name',sprintf('Sim #%d (real time)',idxSim),'NumberTitle','off')
    
    % ---------------- PLOT PROMOTER ACTIVITY ---------------
    subplot(3,1,1)
    %figure('Name','Promoter activity probability (real time)','NumberTitle','off')
    plot(time_real(idxSim,:),PL_activity(idxSim,:),'r','LineWidth',1)
    hold on
    plot(time_real(idxSim,:),PR_activity(idxSim,:),'b','LineWidth',1)
    xlabel('Time (s)');
    ylabel('Promoter activity probability');
    legend('PL', 'PR')
    axis([0 inf 0 inf]);
    title('Promoter activity probability');


    % ---------------- PROTEINS LEVELS ---------------
    subplot(3,1,2)
    %figure('Name','Protein concentrations (real time)','NumberTitle','off')
    h_prot(1)=plot(time_real(idxSim,:),CI_m(idxSim,:),'g','LineWidth',1);
    hold on
    h_prot(2)=plot(time_real(idxSim,:),CI_d(idxSim,:),'c','LineWidth',1);
    hold on
    h_prot(3)=plot(time_real(idxSim,:),MOR_m(idxSim,:),'k','LineWidth',1);
    hold on
    h_prot(4)=plot(time_real(idxSim,:),CIMOR(idxSim,:),'m','LineWidth',1);
    hold on
    h_prot(5)=plot(time_real(idxSim,:),CI_total(idxSim,:),'b','LineWidth',1);
    hold on
    h_prot(6)=plot(time_real(idxSim,:),MOR_total(idxSim,:),'r','LineWidth',1);
    hold on
    
    % PLOTTING "ANNOTATION" LINEs
    % Generation time
    h_anno(1)=plot([T_generation, T_generation], [0, 500], 'g', 'LineWidth',1.5, 'LineStyle', '--');
    hold on
    % Boundary time for quantify variance
    h_anno(2)=plot([time_real(variance_calc_boundary), time_real(variance_calc_boundary)], [0, 500], 'c', 'LineWidth',1.5, 'LineStyle', '--');
    hold on
    h_anno(3)=plot([decision_time_real(idxSim), decision_time_real(idxSim)], [0, 500], 'r', 'LineWidth',1.5, 'LineStyle', '--');
    hold on
    % Mean CI levels
    h_anno(4)=plot([time_real(variance_calc_boundary), time_real(end)], [mean_CI, mean_CI], 'b', 'LineWidth', 1.5, 'LineStyle', '--');
    hold off
    xlabel('Time (s)');
    ylabel('Concentration (nM)');
    legend([h_prot, h_anno], 'CI_m','CI_d','MOR_m','CI:MOR', 'CI total', 'MOR total', 'Generation time', 'Time boundary quantify variance', 'Decision time', 'Mean CI levels (after boundary)')
    %set(gca, 'XTick', [1:delta_t:n_time_steps*60*delta_t]);
    axis([0 inf 0 inf])
    title('Protein concentrations')
    
    % ---------------- CI-MOR RATIO ---------------
    subplot(3,1,3)
    %figure('Name','RATIO-MORtoCI (real time)','NumberTitle','off')
    plot(time_real(idxSim,:),RATIO_MORtoCI(idxSim,:),'m','LineWidth',2);
    xlabel('Time (s)');
    ylabel('RATIO-MORtoCI');
end
