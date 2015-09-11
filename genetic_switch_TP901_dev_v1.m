function TP901_running %[win_ratio_pl]=
% Denne funktion udfører beregninger af promotoraktiviteten for PR og PL i TP901-1 bakteriofagen. 
% Modellens differentialligninger er løst i Mapple. Herefter opdateres Mapple løsningerne iterativt i Matlap. 
% Konstanter vurderet ud fra A. K. Alsings arbejde. se side 43 og 48.

% NEW: Denne function er kopieret fra "TP901_work d. 3. april - CLEAN god version.m" d. 4. september 2015. 

%% TODO - improvements
%1) Collect the normalization constants into variables, to make the code
%easier to change
%2) 

%% PARAMETERS TO ADJUST
% NB: "true" is same as "logical(1)". "false" is the same as logical(0).

p_include_O_M_site = true; % if true, then include the O_M site in the model
p_protein_production_noise_level = 0.5; % a value greater than zero. Recommended range is [0.1-5.0]. A high value results in more stochasticity (noise) in the model

% range P_L activity

p_display_real_time_plots = true; % if true, then display plots with the time in estimated "real time" instead of time steps.
p_quantify_variance = true; % if true, then report the variance ("noise") of the CI and MOR protein levels. recommended (this is easily calculated)

runs = 10; % Antal trails: antal kørsler for hver CI_startkoncentration

% Valg af tidsintervallet for kørsel (default = 20 min)
Tid=1000; % tidstrin (1000 tidstrin ~ 20 min| 3600 tidstrin ~ 3 timer)

%% Physical constants
T = 298.15; % [K]
R = 8.3144621/(4.184*10^3); % [kcal/(K*mol)]
kbT = R*T; % [J/mol] DETTE ER FAKTISK R*T (så dårligt variabelnavn)

% BINDINGSKONSTANTER (bestemt eksperimentelt) MOLÆR
K_or = 2000*10^-9; %[M] NORMALT: 2000
K_ol = 28*10^-9; %[M]
K_od = 28*10^-9; %[M]
K_om = 10*10^-9; %[M] DENNE ER IKKE VIST EKSP.

% BINDINGSENERGIER (i JOULE eller [kcal/mol]?)
% DO NOT CHANGE ME - CHANGE K_or, K_ol, K_od, K_om instead!!
G_or = -kbT*log(1/K_or);
G_ol = -kbT*log(1/K_ol);
G_od = -kbT*log(1/K_od);
G_om = -kbT*log(1/K_om);


% COOPERATIVE BINDING ENERGIES
% G_pair is the "additional" binding energy associated with binding of *two* CI dimers
% G_trip is the "additional" binding energy associated with binding of *three* CI dimers
G_pair = -20000/(4.184*10^3); % [kcal/mol] | G_pair numerator calculated to 20000 (*unsure about the calculation*)
G_trip = -30000/(4.184*10^3); % [kcal/mol] | G_pair numerator calculated to 30000 (*unsure about the calculation*)
% Try also -10000 and -15000 in the numerator for G_pair and G_trip respectively
% ??En faktor 10 ændring i G_kooperativtet giver en faktor ca. 1200 i vægten??

% Protein complex dissociation constants 
K_d = (20*10^-9); % Bindingskonstant mellem CI monomerer [M]
K_CIMOR=(10^-1*10^-9);  % Bindingskonstant mellem CI og MOR [M]


%Estimeret antal tidstrin
tids_array=[1:1:Tid];% lidt mærkeligt system, sorry...
Tid_real=zeros(1,length(tids_array));

% Begyndelsesbetingelser
CI_total(1)=0;
MOR_total(1)=0;

% VALG AF start CI-koncentrationer (overskriver vist nok "begyndelsesbetingelser" for CI_total)
CI_concentations=0; %[0:1:10];
%CI_concentations=[0:1:10];

%% Initiering af arrays
CI_total=zeros(1,length(tids_array));
MOR_total=zeros(1,length(tids_array));
CI_m=zeros(1,length(tids_array));
CI_d=zeros(1,length(tids_array));
MOR_m=zeros(1,length(tids_array));
CIMOR=zeros(1,length(tids_array));
PR=zeros(1,length(tids_array));
PL=zeros(1,length(tids_array));
FORHOLD_MOR_CI=zeros(1,length(tids_array)); % Til plotting af forhold
PL_win=zeros(1,length(CI_concentations));
PR_win=zeros(1,length(CI_concentations));

%% Produktionskonstanter - Gillespie fra A.A
T_generation=3600; % [seconds] Cell generation time. 
% Lactococcus lactis is ~ 52 min
% E.coli is ~20 min

%PARAMETRE FOR CI-PRODUKTION
r_CItrsc=0.22;% Transcriptionsrate, Egan(2004). Bruges til at udregne raten, a_1=a_CIprod
n_1=round(60*p_protein_production_noise_level); % the value 30 seems good. choose values between [10-300]
N_CI_total=200; % Steady state "concentrations"
r_CI=N_CI_total/(T_generation*n_1); % Skalering r_CItrsc/F for at opfylde lignignen N_total/T_generation ca lig med n_CIprod*r_CI
tau_CI=1/3000; % [s^-1] rate for nedbrydning DENNE VÆRDI ANGIVER 
n_2=1; % Ændring i CI_total pr. nedbrydningsevent

%PARAMETRE FOR MOR-PRODUKTION
r_MORtrsc=0.22;% Transcriptionsrate, {0.1/s;1/s} Bruges til at udregne raten, a_3=a_MORprod
n_3=round(60*p_protein_production_noise_level); % the value 30 seems good. choose values between [10-300]
N_MOR_total=200; % steady state "concentrations" [old value was 800]
r_MOR=N_MOR_total/(T_generation*n_3); % Skalering r_CItrsc/F for at opfylde lignignen N_total/T_generation ca lig med n_CIprod*r_CI
tau_MOR=1/3000;
n_4=1 ;% Ændring i MOR_total pr. nedbrydningsevent


%% MAIN LOOP 
% Måling af kørselstid
tic_value1=tic;

tic
for j=[1:1:length(CI_concentations)]
    CI_total(1)=CI_concentations(j);
    tic
    % Antal trials (antal kørsler for hver CI_startkoncentration)
    for trial=1:runs % For hver 10. trial tager det 3-4 sekunder for en CI-koncentration.
        % INITIALYZING FLAGS:
        % Idea is that initial value is 0 and when found then setting value to 1
        been_here_before = 0; % this is for udregning af "VINDER"
        flag_kritisk_tid(trial) = 0; % this if for saving critic times
        
        for t=tids_array % Hver trial tager ca. 0.26sek
            
            % UDREGNING AF CI, MOR OG CIMOR KONCENTRATIONER
            %-------------------------------------------------------%
            %MOR isoleret af lign2 og indsat i lign1. Derefter solving for CI_m.
            % OK LØSNING, men giver komplekse tal (ikke noget problem)
            CI_m(t)=(1/6)*(-9*CI_total(t)*K_d^2+36*CI_total(t)*K_d*K_CIMOR+3*K_d^2*K_CIMOR+6*K_d*K_CIMOR^2+9*MOR_total(t)*K_d^2+18*MOR_total(t)*K_d*K_CIMOR-K_d^3-8*K_CIMOR^3+3*sqrt(-12*CI_total(t)*K_d^3*K_CIMOR*MOR_total(t)+240*CI_total(t)*K_d^2*MOR_total(t)*K_CIMOR^2-24*CI_total(t)^3*K_d^3-3*CI_total(t)^2*K_d^4+12*K_d^3*K_CIMOR^3-3*K_d^4*K_CIMOR^2-12*K_d^2*K_CIMOR^4+24*MOR_total(t)^3*K_d^3-3*MOR_total(t)^2*K_d^4+96*CI_total(t)^2*K_d^2*K_CIMOR^2-12*CI_total(t)*K_d^3*K_CIMOR^2-96*CI_total(t)*K_d*K_CIMOR^4-6*CI_total(t)*K_d^4*K_CIMOR+96*CI_total(t)*K_d^2*K_CIMOR^3-72*CI_total(t)*K_d^3*MOR_total(t)^2+6*CI_total(t)*K_d^4*MOR_total(t)-48*CI_total(t)^2*K_d^3*K_CIMOR+72*CI_total(t)^2*K_d^3*MOR_total(t)+60*K_d^3*K_CIMOR*MOR_total(t)^2-6*K_d^4*K_CIMOR*MOR_total(t)+48*K_d^3*K_CIMOR^2*MOR_total(t)-24*K_d^2*K_CIMOR^3*MOR_total(t)-12*MOR_total(t)^2*K_d^2*K_CIMOR^2))^(1/3)-(6*(-(1/6)*CI_total(t)*K_d+(1/18)*K_d*K_CIMOR+(1/6)*MOR_total(t)*K_d-(1/36)*K_d^2-(1/9)*K_CIMOR^2))/(-9*CI_total(t)*K_d^2+36*CI_total(t)*K_d*K_CIMOR+3*K_d^2*K_CIMOR+6*K_d*K_CIMOR^2+9*MOR_total(t)*K_d^2+18*MOR_total(t)*K_d*K_CIMOR-K_d^3-8*K_CIMOR^3+3*sqrt(-12*CI_total(t)*K_d^3*K_CIMOR*MOR_total(t)+240*CI_total(t)*K_d^2*MOR_total(t)*K_CIMOR^2-24*CI_total(t)^3*K_d^3-3*CI_total(t)^2*K_d^4+12*K_d^3*K_CIMOR^3-3*K_d^4*K_CIMOR^2-12*K_d^2*K_CIMOR^4+24*MOR_total(t)^3*K_d^3-3*MOR_total(t)^2*K_d^4+96*CI_total(t)^2*K_d^2*K_CIMOR^2-12*CI_total(t)*K_d^3*K_CIMOR^2-96*CI_total(t)*K_d*K_CIMOR^4-6*CI_total(t)*K_d^4*K_CIMOR+96*CI_total(t)*K_d^2*K_CIMOR^3-72*CI_total(t)*K_d^3*MOR_total(t)^2+6*CI_total(t)*K_d^4*MOR_total(t)-48*CI_total(t)^2*K_d^3*K_CIMOR+72*CI_total(t)^2*K_d^3*MOR_total(t)+60*K_d^3*K_CIMOR*MOR_total(t)^2-6*K_d^4*K_CIMOR*MOR_total(t)+48*K_d^3*K_CIMOR^2*MOR_total(t)-24*K_d^2*K_CIMOR^3*MOR_total(t)-12*MOR_total(t)^2*K_d^2*K_CIMOR^2))^(1/3)-(1/6)*K_d-(1/3)*K_CIMOR;
            CI_m(t)=real(CI_m(t)); % Dont worry about this :)
            
            % De følgende to måde at regne MOR_m ud på, giver tilsyneladende sammme plot
            %MOR_m=-(CI_m(1,j)*K_d+2*CI_m(1,j)^2-CI_total*K_d)*K_CIMOR/(K_d*CI_m(1,j)); % isolering af MOR_m fra lign1
            CI_d(t)=(CI_m(t)^2)/K_d;
            MOR_m(t)=MOR_total(t)*K_CIMOR/(K_CIMOR+CI_m(t)); % BRUG DENNE, da dette udtryk er indsat i monosolveren. Isolering af MOR_m fra lign2
            CIMOR(t)=CI_m(t)*MOR_m(t)/K_CIMOR;
            %-------------------------------------------------------%
            
            %VÆGTE
            w(1)  = 1;
            w(2)  = CI_d(t)*exp(-G_or/kbT);
            w(3)  = CI_d(t)*exp(-G_ol/kbT);
            w(4)  = CI_d(t)*exp(-G_od/kbT);
            w(5)  = w(2)*w(3);
            w(6)  = w(2)*w(4);
            w(7)  = w(3)*w(4);
            w(8)  = w(2)*w(3)*w(4);
            w(9)  = w(2)*w(3)*exp(-G_pair/kbT);
            w(10) = w(2)*w(4)*exp(-G_pair/kbT);
            w(11) = w(3)*w(4)*exp(-G_pair/kbT);
            w(12) = w(8)*exp(-G_pair/kbT);
            w(13) = w(8)*exp(-G_pair/kbT);
            w(14) = w(8)*exp(-G_pair/kbT);
            w(15) = w(8)*exp(-G_trip/kbT);
            w(16) = CIMOR(t)*exp(-G_om/kbT); % all weights that includes w(16) have O_M binding.
            w(17) = w(2)*w(16); 
            w(18) = w(3)*w(16);
            w(19) = w(4)*w(16);
            w(20) = w(5)*w(16);
            w(21) = w(6)*w(16);
            w(22) = w(7)*w(16);
            w(23) = w(8)*w(16);
            w(24) = w(9)*w(16);
            w(25) = w(10)*w(16);
            w(26) = w(11)*w(16);
            w(27) = w(12)*w(16);
            w(28) = w(13)*w(16);
            w(29) = w(14)*w(16);
            w(30) = w(15)*w(16);
            
            if p_include_O_M_site;
                w_length=30; % sæt til 30 ==> MED O_M
            else
                w_length=15; % sæt til 15 ==> UDEN O_M
            end
            
            % Making a weight table?
            %pr=zeros(1,15); % NO O_M
            pr=zeros(1,30);
            pr_on=[1 3 4 7 11];
            pr(1,pr_on)=1;
            
            %pl=zeros(1,15); % NO O_M
            %pl_on=[1 4]; % NO O_M
            pl=zeros(1,30);
            pl_on=[1 4 16 19];
            pl(1,pl_on)=1;
            
            % Finding the weights
            w_total=[0;0]; % First row is for PL and second is for PR
            for k1=[1:w_length];
                w_total(1)=pl(k1)*w(k1)+w_total(1);
                w_total(2)=pr(k1)*w(k1)+w_total(2);
            end
            % Udregning af promoteraktiviteten
            Z=sum(w(1:w_length));
            PL(t)=(w_total(1)/Z); %*beta_PL
            PR(t)=(w_total(2)/Z); %*beta_PR
            
            if t+1<=length(tids_array) % avoiding out of bounds errors? check if we can update once more
                % Gillespie algorimte
                a_1=PR(t)*r_CI;
                a_2=CI_total(t)/(0.831*10^-9)*tau_CI;  % For E. coli: N_molekyler = koncentration(nM)/0.831
                a_3=PL(t)*r_MOR;
                a_4=MOR_total(t)/(0.831*10^-9)*tau_MOR; % For E. coli: N_molekyler = koncentration(nM)/0.831
                
                a_0=a_1+a_2+a_3+a_4;
                r1=random('unif', 0, 1); % random number from uniform distrib.
                r2=random('unif', 0, 1); % random number from uniform distrib.
                time_step=(1/a_0)*log(1/r1); % determining next event
                Tid_real(t+1)=Tid_real(t)+time_step; % setting time for next event
                % Hilberts way of determining which j-value to use
                if r2*a_0<=a_1
                    event=1;
                elseif r2*a_0<=a_1+a_2
                    event=2;
                elseif r2*a_0<=a_1+a_2+a_3
                    event=3;
                elseif r2*a_0<=a_1+a_2+a_3+a_4
                    event=4;
                else
                    error 'EXCEEDED LIMIT - line 313 approx'
                end
                % Setting consequence of j={1..4}
                switch event
                    case 1
                        MOR_update=0;
                        CI_update=n_1;
                    case 2
                        MOR_update=0;
                        CI_update=-1*n_2;
                    case 3
                        MOR_update=n_3;
                        CI_update=0;
                    case 4
                        MOR_update=-1*n_4;
                        CI_update=0;
                end
                
                % Opdatering af koncentrationer:
                % ***Husk at antal proteiner produceret konverteres til koncentration!***
                % For E. coli: molekyler * 0.831 = koncentration (nM)
                CI_total(t+1)=CI_total(t)+CI_update*0.831*10^-9;
                MOR_total(t+1)=MOR_total(t)+MOR_update*0.831*10^-9;
                
            end % end if tjek af tidsarray
            
            
            % GEMMER CI_total og MOR_total til plotting - NYT d. 15/03/2012 - Pascal
            % Tiderne er rækker, trials er columns
            CI_total_array(t,trial)=CI_total(t);
            MOR_total_array(t,trial)=MOR_total(t);
            
            
            % GEMMER tiden for tid_array for real tider omkring 15 min
            % Another way of finding out the "critical times" would be to use a threshold for the protein levels.
            % E.g. 295 <= CI_total <= 305
            if Tid_real(t) >= 60*15-30 && flag_kritisk_tid(trial) == 0
                kritisk_tid_trin(trial) = t;
                kritisk_tid_real(trial) = Tid_real(t); % Not really needed, but fun to save.
                flag_kritisk_tid(trial) = 1;
            end
            
            % Udregning af FORHOLD_MOR_CI forhold (TIL PLOTTING)
            FORHOLD_MOR_CI(t)=MOR_total(t)/CI_total(t);
            
            %UDREGNING AF "VINDER" Lytisk eller lysogen
            % PL_win(j): hvor j er for en given start CI koncentration.
            % Der overskrives IKKE selvom der er flere hits i intervallet for Tid_real
            if Tid_real(t) >= 60*15-15 && been_here_before == 0 %&& Tid_real(t) <= 60*15+15
                been_here_before = 1;
                FORHOLD_MOR_CI(t)=MOR_total(t)/CI_total(t);
                if FORHOLD_MOR_CI(t) >= 0.8
                    PL_win(j)=PL_win(j)+1;
                    
                else %FORHOLD_MOR_CI(t) < 0.8
                    PR_win(j)=PR_win(j)+1;
                end
                
            end % end for bestemmelse af VINDER
            
        end % End for-loop t=tidsarray.
        
    % DISPLAYER KRITISKE TIDER:
    display (['Kritisk_tid_trin (step):', num2str(kritisk_tid_trin(trial))])
    display (['Kritisk_tid_real (s):', num2str(kritisk_tid_real(trial))])
    
    end% End for-loop (trials).
    toc
end % For CI koncentrationer (j)
toc

%-----------------------------END FOR SIMULERING----------------------------
%Her følger plots og analyse af simulering
% - Varians bestemmelse
% - Plots
% - Fordeling

%% Converting to units of nano molar (nM)
CI_total_array(:,:)=CI_total_array(:,:)*10^9;
MOR_total_array(:,:)=MOR_total_array(:,:)*10^9;
CI_m(:)=CI_m(:)*10^9;
CI_d(:)=CI_d(:)*10^9;
MOR_m(:)=MOR_m(:)*10^9;
CIMOR(:)=CIMOR(:)*10^9;
CI_total(:)=CI_total(:)*10^9;
MOR_total(:)=MOR_total(:)*10^9;


% Bestemmelse af standard afvigelse d. 2. april - Pascal:
% OBS: this calculates the standard deviation for the *LAST RUN*.
if p_quantify_variance == 1
    boundery = 0.5*Tid; %2*10^4;
    std_CI = std(CI_total_array([boundery:end],end),0,1);
    std_MOR = std(MOR_total_array([boundery:end],end),0,1);
    mean_CI = mean(CI_total_array([boundery:end],end));
    mean_MOR = mean(MOR_total_array([boundery:end],end));
    
    % Pascal's implementation of calculating the standard deviation
    % NB: the results are identical.
%     SS = 0;
%     for n=[boundery:tids_array(end)]
%         SS = SS+sum((CI_total_array(n,end)-mean_CI)^2);
%     end
%     std_own = sqrt(SS/(length(CI_total_array([boundery:end],end))-1));
%     display(['Standard deviation for CI_total_OWN was: ', num2str(std_own), ...
%             ' (', num2str(std_own/mean_CI*100), '% of the mean_CI)'])

    display(['Standard deviation for the last "run" for CI_total was: ', num2str(std_CI), ...
            ' (', num2str(std_CI/mean_CI*100), '% of the mean_CI)'])
    display(['Standard deviation for the last "run" for MOR_total was: ', num2str(std_MOR), ...
            ' (', num2str(std_MOR/mean_MOR*100), '% of the mean_MOR)'])
end % end if for bestemmelse af varians



%% PLOT
% PLOT - NYT d. 15/03/2012 - Pascal
figure('Name','CI[total] vs MOR[total] plot','NumberTitle','off')

% Area plot
p_area = area([1:1:1.5*N_MOR_total], [1:1:1.5*N_MOR_total]);
set(p_area,'FaceColor',[1,0.6,0.6],'LineStyle', 'none'); % 'FaceAlpha' does not work for MATLAB2014b
hold on;

cmap = hsv(runs); % loading colormap
for l=1:runs
    h(l)=plot(CI_total_array(:,l),MOR_total_array(:,l), 'Color', cmap(l,:));
    hold on
    s{l}=sprintf('%s%g','Run=',l); % defining title for legend
end
legend(h, s) % setting legend

%Plotter punkt for kritisk tid (15 min Tid_real)
for l=1:runs
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
figure('Name','Protein koncentrationer','NumberTitle','off')
plot(tids_array,CI_m(1,:),'g','LineWidth',1)
hold on
plot(tids_array,CI_d(1,:),'c','LineWidth',1)
hold on
plot(tids_array,MOR_m(1,:),'k','LineWidth',1)
hold on
plot(tids_array,CIMOR(1,:),'m','LineWidth',1)
hold on
plot(tids_array,CI_total(1,:),'b','LineWidth',1)
hold on
plot(tids_array,MOR_total(1,:),'r','LineWidth',1)
hold on
hold off
xlabel('TIDSTRIN');
ylabel('Koncentration (nM)');
legend('CI_m','CI_d','MOR_m','CI:MOR', 'CI total', 'MOR total')
% set (gca, 'XTick', [1:delta_t:Tid*60*delta_t]);
%axis([1 Tid 0 10^3])
title('Protein koncentrationer')

%% PLOT
figure('Name','Promoteraktivitetssandsynlighed','NumberTitle','off')
plot(tids_array,PL(1,:),'r','LineWidth',0.01)
hold on
plot(tids_array,PR(1,:),'b','LineWidth',0.01)
xlabel('TIDSTRIN');
ylabel('Promoteraktivitetssandsynlighed');
legend('PL', 'PR')
axis([1 (length(tids_array)) 0 1]);
title('Promoteraktivitetssandsynlighed');


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


%% PLOTS | Real time 
if p_display_real_time_plots == true
    %% PLOT PROMOTER ACTIVITY
    figure('Name','Promoteraktivitetssandsynlighed med realtid','NumberTitle','off')
    plot(Tid_real(1,:),PL(1,:),'r','LineWidth',0.01)
    hold on
    plot(Tid_real(1,:),PR(1,:),'b','LineWidth',0.01)
    xlabel('Tid (s)');
    ylabel('Promoteraktivitetssandsynlighed');
    legend('PL', 'PR')
    axis([1 Tid_real(1,(length(tids_array))) 0 1]);
    title('Promoteraktivitetssandsynlighed');

    % Forhold mellem CI og MOR plot
    figure('Name','Forhold mellem MOR og CI','NumberTitle','off')
    plot(Tid_real(:),FORHOLD_MOR_CI(1,:),'m','LineWidth',2);
    xlabel('Tid (s)');
    ylabel('Forhold mellem MOR_t og CI_t');

    %% PLOT PROTEINS
    figure('Name','Protein koncentrationer i realtid','NumberTitle','off')
    plot(Tid_real(1,:),CI_m(1,:),'g','LineWidth',1)
    hold on
    plot(Tid_real(1,:),CI_d(1,:),'c','LineWidth',1)
    hold on
    plot(Tid_real(1,:),MOR_m(1,:),'k','LineWidth',1)
    hold on
    plot(Tid_real(1,:),CIMOR(1,:),'m','LineWidth',1)
    hold on
    plot(Tid_real(1,:),CI_total(1,:),'b','LineWidth',1)
    hold on
    plot(Tid_real(1,:),MOR_total(1,:),'r','LineWidth',1)
    hold on
    % PLOTTING VERTICAL LINE d. 2. april
    plot([T_generation, T_generation], [100, 500], 'r', 'LineWidth',2);
    hold on
    plot([Tid_real(boundery), Tid_real(boundery)], [100, 500], 'r', 'LineWidth',2);
    hold on
    plot([Tid_real(boundery), Tid_real(end)], [mean_CI, mean_CI], 'r', 'LineWidth',2);
    %tids_array(end)
    hold off
    xlabel('Tid (s)');
    ylabel('Koncentration (nM)');
    legend('CI_m','CI_d','MOR_m','CI:MOR', 'CI total', 'MOR total')
    % set (gca, 'XTick', [1:delta_t:Tid*60*delta_t]);
    %axis([1 Tid 0 10^3])
    title('Protein koncentrationer')
end

end %End Function