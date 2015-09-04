function TP901_running %[win_ratio_pl]=
% Denne funktion udf�rer beregninger af promotoraktiviteten for PR og PL i
% TP901-1 bakteriofagen. 
% Der arbejdes med valgte differentialligninger, via iterativ metode. 
% Konstanter vurderet ud fra A. K. Alsings arbejde. se side 43 og 48.
% NEW: Denne function er kopieret fra "TP901_work d. 3. april - CLEAN god version.m" d. 4. september 2015. 

T = 298.15; % [K]
R = 8.3144621/(4.184*10^3); % [kcal/(K*mol)]
kbT = R*T; % [J/mol] DETTE ER FAKTISK R*T (s� d�rligt variabelnavn)

% BINDINGSKONSTANTER (bestemt eksperimentelt) MOL�R
K_or = 2000*10^-9; %[M] NORMALT: 2000
K_ol = 28*10^-9; %[M]
K_od = 28*10^-9; %[M]
K_om = 10*10^-9; %[M] DENNE ER IKKE VIST EKSP.

% BINDINGSENERGIER (i JOULE)
G_or = -kbT*log(1/K_or); % [kcal/mol]
G_ol = -kbT*log(1/K_ol);
G_od = -kbT*log(1/K_od);
G_om = -kbT*log(1/K_om);

% En faktor 10 �ndring i G_kooperativtet giver en faktor ca. 1200 i v�gten
G_pair = -10000/(4.184*10^3); % [kcal/mol] Udregnet til: 20000
G_trip = -15000/(4.184*10^3); % [kcal/mol] Udregnet til: 30000

K_d = (20*10^-9); %[M]
K_CIMOR=(10^-1*10^-9);  % Bindingskonstant mellem CI og MOR

% Valg af tidsintervallet for k�rsel (default = 20 min)
% HER HVILER delta_t
Tid=1000; % tidstrin (3600 tidstrin ~ 3 timer
%tids_array=[1:1:Tid*60/delta_t];%/delta_t]; %lidt m�rkeligt system, sorry...
% Omregning af "tid" til sekunder. Derefter divison med vores tidsskridt
% alts� delta_t. Dette giver antal opdateringer indtil �nskede "tid"

%Estimeret antal tidstrin
tids_array=[1:1:Tid];% lidt m�rkeligt system, sorry...
Tid_real=zeros(1,length(tids_array));
% Gr�nse v�rdi for CI_tot(1) og PL dominans. 19.36-19.37
% Begyndelsesbetingelser
CI_total(1)=0;
MOR_total(1)=0;

% VALG AF start CI-koncentrationer (overskriver vist nok
% "begyndelsesbetingelser" for CI_total)
CI_koncentrationer=0; %[0:1:10];

% Antal trails: antal k�rsler for hver CI_startkoncentration
runs = 10;

%Initiering af arrays
CI_total=zeros(1,length(tids_array));
MOR_total=zeros(1,length(tids_array));
CI_m=zeros(1,length(tids_array));
CI_d=zeros(1,length(tids_array));
MOR_m=zeros(1,length(tids_array));
CIMOR=zeros(1,length(tids_array));
PR=zeros(1,length(tids_array));
PL=zeros(1,length(tids_array));
FORHOLD_MOR_CI=zeros(1,length(tids_array)); % Til plotting af forhold
PL_win=zeros(1,length(CI_koncentrationer));
PR_win=zeros(1,length(CI_koncentrationer));

% Produktionskonstanter - Gillespie fra A.A
T_generation=3600/2; % 3600= en halv time i sekunder
%PARAMETRE FOR CI-PRODUKTION
r_CItrsc=0.22;% Transcriptionsrate, Egan(2004). Bruges til at udregne raten, a_1=a_CIprod
n_1=30; %n_CIprod{30;300}
N_CI_total=200; % Steady state concentrations
r_CI=N_CI_total/(T_generation*n_1); % Skalering r_CItrsc/F for at opfylde lignignen N_total/T_generation ca lig med n_CIprod*r_CI
tau_CI=1/3000; % [s^-1] rate for nedbrydning DENNE V�RDI ANGIVER 
n_2=1; % �ndring i CI_total pr. nedbrydningsevent

%PARAMETRE FOR MOR-PRODUKTION
r_MORtrsc=0.22;% Transcriptionsrate, {0.1/s;1/s} Bruges til at udregne raten, a_3=a_MORprod
n_3=30; %n_CIprod{30;300}
N_MOR_total=800; % steady state
r_MOR=N_MOR_total/(T_generation*n_3); % Skalering r_CItrsc/F for at opfylde lignignen N_total/T_generation ca lig med n_CIprod*r_CI
tau_MOR=1/3000;
n_4=1 ;% �ndring i MOR_total pr. nedbrydningsevent

%VIL DU HAVE VARIANSBESTEMMELSE?
varians_bestemmelse = 1;

% M�ling af k�rselstid
tic_value1=tic;

tic
for j=[1:1:length(CI_koncentrationer)]
    CI_total(1)=CI_koncentrationer(j);
    tic
    % Antal trials (antal k�rsler for hver CI_startkoncentration)
    for trial=1:runs % For hver 10. trial tager det 3-4 sekunder for en CI-koncentration.
        % INITIALYZING FLAGS:
        % Idea is that initial value is 1 and when found then setting value to 0
        been_here_before = 1; % this is for udregning af "VINDER"
        flag_kritisk_tid(trial) = 1; % this if for saving critic times
        
        for t=tids_array % Hver trial tager ca. 0.26sek
            
            % UDREGNING AF CI, MOR OG CIMOR KONCENTRATIONER
            %-------------------------------------------------------%
            %MOR isoleret af lign2 og indsat i lign1. Derefter solving for CI_m.
            % OK L�SNING, men giver komplekse tal (ikke noget problem)
            CI_m(t)=(1/6)*(-9*CI_total(t)*K_d^2+36*CI_total(t)*K_d*K_CIMOR+3*K_d^2*K_CIMOR+6*K_d*K_CIMOR^2+9*MOR_total(t)*K_d^2+18*MOR_total(t)*K_d*K_CIMOR-K_d^3-8*K_CIMOR^3+3*sqrt(-12*CI_total(t)*K_d^3*K_CIMOR*MOR_total(t)+240*CI_total(t)*K_d^2*MOR_total(t)*K_CIMOR^2-24*CI_total(t)^3*K_d^3-3*CI_total(t)^2*K_d^4+12*K_d^3*K_CIMOR^3-3*K_d^4*K_CIMOR^2-12*K_d^2*K_CIMOR^4+24*MOR_total(t)^3*K_d^3-3*MOR_total(t)^2*K_d^4+96*CI_total(t)^2*K_d^2*K_CIMOR^2-12*CI_total(t)*K_d^3*K_CIMOR^2-96*CI_total(t)*K_d*K_CIMOR^4-6*CI_total(t)*K_d^4*K_CIMOR+96*CI_total(t)*K_d^2*K_CIMOR^3-72*CI_total(t)*K_d^3*MOR_total(t)^2+6*CI_total(t)*K_d^4*MOR_total(t)-48*CI_total(t)^2*K_d^3*K_CIMOR+72*CI_total(t)^2*K_d^3*MOR_total(t)+60*K_d^3*K_CIMOR*MOR_total(t)^2-6*K_d^4*K_CIMOR*MOR_total(t)+48*K_d^3*K_CIMOR^2*MOR_total(t)-24*K_d^2*K_CIMOR^3*MOR_total(t)-12*MOR_total(t)^2*K_d^2*K_CIMOR^2))^(1/3)-(6*(-(1/6)*CI_total(t)*K_d+(1/18)*K_d*K_CIMOR+(1/6)*MOR_total(t)*K_d-(1/36)*K_d^2-(1/9)*K_CIMOR^2))/(-9*CI_total(t)*K_d^2+36*CI_total(t)*K_d*K_CIMOR+3*K_d^2*K_CIMOR+6*K_d*K_CIMOR^2+9*MOR_total(t)*K_d^2+18*MOR_total(t)*K_d*K_CIMOR-K_d^3-8*K_CIMOR^3+3*sqrt(-12*CI_total(t)*K_d^3*K_CIMOR*MOR_total(t)+240*CI_total(t)*K_d^2*MOR_total(t)*K_CIMOR^2-24*CI_total(t)^3*K_d^3-3*CI_total(t)^2*K_d^4+12*K_d^3*K_CIMOR^3-3*K_d^4*K_CIMOR^2-12*K_d^2*K_CIMOR^4+24*MOR_total(t)^3*K_d^3-3*MOR_total(t)^2*K_d^4+96*CI_total(t)^2*K_d^2*K_CIMOR^2-12*CI_total(t)*K_d^3*K_CIMOR^2-96*CI_total(t)*K_d*K_CIMOR^4-6*CI_total(t)*K_d^4*K_CIMOR+96*CI_total(t)*K_d^2*K_CIMOR^3-72*CI_total(t)*K_d^3*MOR_total(t)^2+6*CI_total(t)*K_d^4*MOR_total(t)-48*CI_total(t)^2*K_d^3*K_CIMOR+72*CI_total(t)^2*K_d^3*MOR_total(t)+60*K_d^3*K_CIMOR*MOR_total(t)^2-6*K_d^4*K_CIMOR*MOR_total(t)+48*K_d^3*K_CIMOR^2*MOR_total(t)-24*K_d^2*K_CIMOR^3*MOR_total(t)-12*MOR_total(t)^2*K_d^2*K_CIMOR^2))^(1/3)-(1/6)*K_d-(1/3)*K_CIMOR;
            CI_m(t)=real(CI_m(t)); % Dont worry about this :)
            
            % De f�lgende to m�de at regne MOR_m ud p�, giver tilsyneladende sammme plot
            %MOR_m=-(CI_m(1,j)*K_d+2*CI_m(1,j)^2-CI_total*K_d)*K_CIMOR/(K_d*CI_m(1,j)); % isolering af MOR_m fra lign1
            CI_d(t)=(CI_m(t)^2)/K_d;
            MOR_m(t)=MOR_total(t)*K_CIMOR/(K_CIMOR+CI_m(t)); % BRUG DENNE, da dette udtryk er indsat i monosolveren. Isolering af MOR_m fra lign2
            CIMOR(t)=CI_m(t)*MOR_m(t)/K_CIMOR;
            %-------------------------------------------------------%
            
            %V�GTE
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
            %     w(16) = CIMOR(t)*exp(-G_om/kbT);
            %     w(17) = w(2)*w(16);
            %     w(18) = w(3)*w(16);
            %     w(19) = w(4)*w(16);
            %     w(20) = w(5)*w(16);
            %     w(21) = w(6)*w(16);
            %     w(22) = w(7)*w(16);
            %     w(23) = w(8)*w(16);
            %     w(24) = w(9)*w(16);
            %     w(25) = w(10)*w(16);
            %     w(26) = w(11)*w(16);
            %     w(27) = w(12)*w(16);
            %     w(28) = w(13)*w(16);
            %     w(29) = w(14)*w(16);
            %     w(30) = w(15)*w(16);
            w_length=length(w);
            
            % Making a weight table?
            pr=zeros(1,15);
            pr_on=[1 3 4 7 11];
            pr(1,pr_on)=1;
            pl=zeros(1,15);
            pl_on=[1 4];
            pl(1,pl_on)=1;
            
            % Finding the weights
            w_total=[0;0]; % First row is for PL and second is for PR
            for k1=[1:w_length];
                w_total(1)=pl(k1)*w(k1)+w_total(1);
                w_total(2)=pr(k1)*w(k1)+w_total(2);
            end
            % Udregning af promoteraktiviteten
            Z=sum(w(:));
            PL(t)=(w_total(1)/Z); %*beta_PL
            PR(t)=(w_total(2)/Z); %*beta_PR
            
            if t+1<=length(tids_array)
                % Gillespie algorimte
                a_1=PR(t)*r_CI;
                a_2=CI_total(t)/(0.831*10^-9)*tau_CI;
                a_3=PL(t)*r_MOR;
                a_4=MOR_total(t)/(0.831*10^-9)*tau_MOR;
                
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
                % husk at antal proteiner produceret konverteres til koncentration!
                % For E. coli: molekyler * 0.831 = koncentration (nM)
                CI_total(t+1)=CI_total(t)+CI_update*0.831*10^-9;
                MOR_total(t+1)=MOR_total(t)+MOR_update*0.831*10^-9;
                
            end % end if tjek af tidsarray
            
            
            % GEMMER CI_total og MOR_total til plotting - NYT d. 15/03/2012 - Pascal
            % Tiderne er r�kker, trials er columns
            CI_total_array(t,trial)=CI_total(t);
            MOR_total_array(t,trial)=MOR_total(t);
            
            % Gemmer kritiske v�rdier
            % if CI_total_var(t,trial) >= 295 && CI_total_var(t,trial) <= 305
            %     kritisk_tid_trin(trial) = t;
            %     kritisk_tid_real(trial) = Tid_real(t);
            % end
            
            % GEMMER tiden for tid_array for real tider omkring 15 min
            if Tid_real(t) >= 60*15-30 && flag_kritisk_tid(trial) == 1
                kritisk_tid_trin(trial) = t;
                flag_kritisk_tid(trial) = 0;
            end
            
            % Udregning af FORHOLD_MOR_CI forhold (TIL PLOTTING)
            FORHOLD_MOR_CI(t)=MOR_total(t)/CI_total(t);
            
            %UDREGNING AF "VINDER" Lytisk eller lysogen
            % PL_win(j): hvor j er for en given start CI koncentration.
            % Der overskrives IKKE selvom der er flere hits i intervallet for Tid_real
            if Tid_real(t) >= 60*15-15 && been_here_before == 1 %&& Tid_real(t) <= 60*15+15
                been_here_before = 0;
                FORHOLD_MOR_CI(t)=MOR_total(t)/CI_total(t);
                if FORHOLD_MOR_CI(t) >= 0.8
                    PL_win(j)=PL_win(j)+1;
                    
                else %FORHOLD_MOR_CI(t) < 0.8
                    PR_win(j)=PR_win(j)+1;
                end
                
            end % end for bestemmelse af VINDER
            
        end % End for-loop t=tidsarray.
    end% End for-loop (trials).
    toc
end % For CI koncentrationer (j)
toc

%-----------------------------END FOR SIMULERING----------------------------
%Her f�lger plots og analyse af simulering
% - Varians bestemmelse
% - Plots
% - Fordeling

CI_total_array(:,:)=CI_total_array(:,:)*10^9;
MOR_total_array(:,:)=MOR_total_array(:,:)*10^9;
CI_m(:)=CI_m(:)*10^9;
CI_d(:)=CI_d(:)*10^9;
MOR_m(:)=MOR_m(:)*10^9;
CIMOR(:)=CIMOR(:)*10^9;
CI_total(:)=CI_total(:)*10^9;
MOR_total(:)=MOR_total(:)*10^9;


% Bestemmelse af standard afvigelse d. 2. april - Pascal:
if varians_bestemmelse == 1
    boundery = 0.5*Tid; %2*10^4;
    std_CI = std(CI_total_array([boundery:end],1),0,1);
    std_MOR = std(MOR_total_array([boundery:end],1),0,1);
    mean_CI = mean(CI_total_array([boundery:end],1));
    mean_MOR = mean(MOR_total_array([boundery:end],1));

    % Pascal's implementation of calculating the standard deviation
    % NB: the results are identical.
%     SS = 0;
%     for n=[boundery:tids_array(end)]
%         SS = SS+sum((CI_total_array(n,1)-mean_CI)^2);
%     end
%     std_own = sqrt(SS/(length(CI_total_array([boundery:end],1))-1));
%     display(['Standard deviation for CI_total_OWN was: ', num2str(std_own), ...
%             ' (', num2str(std_own/mean_CI*100), '% of the mean_CI)'])

    display(['Standard deviation for CI_total was: ', num2str(std_CI), ...
            ' (', num2str(std_CI/mean_CI*100), '% of the mean_CI)'])
    display(['Standard deviation for MOR_total was: ', num2str(std_MOR), ...
            ' (', num2str(std_MOR/mean_MOR*100), '% of the mean_MOR)'])
end % end if for bestemmelse af varians

% DISPLAYER KRITISKE TIDER:
% display (['Kritisk_tid_trin:', num2str(kritisk_tid_trin(:))])
% display (['Kritisk_tid_real:', num2str(kritisk_tid_real(:))])

%% PLOT
% PLOT - NYT d. 15/03/2012 - Pascal
figure('Name','CI[total] vs MOR[total] plot','NumberTitle','off')
cmap = hsv(runs); % loading colormap
for l=1:runs
    plot(CI_total_array(:,l),MOR_total_array(:,l), 'Color', cmap(l,:))
    hold on
    %NEDENST�ENDE KAN SLETTES. PRODUCERET d. 2. april - Pacal
    %plot(CI_total_var(kritisk_tid_trin(l),l),MOR_total_var(kritisk_tid_trin(l),l), '*')
    %hold on
    s{l}=sprintf('%s%g','Run=',l); % defining title for legend
    legend(s) % setting legend
    
    %Plotter punkt for kritisk tid (15 min Tid_real)
    %if point_maker_exist(l) = 1
    %end
end
for l=1:runs
    plot(CI_total_array(kritisk_tid_trin(l),l),MOR_total_array(kritisk_tid_trin(l),l), '*', 'MarkerSize', 10)
    hold on
end
% Plot af ret linje til indikation af forhold mellem CI og MOR
plot ([1:1:1.5*N_MOR_total],[1:1:1.5*N_MOR_total], 'Linewidth', 2);
hold off;
% Pr�v med scaling constrained, og plot en ret linje
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
%plot([T_generation, T_generation], [100, 500], 'r', 'LineWidth',2);
%hold on
%plot([Tid_real(boundery), Tid_real(boundery)], [100, 500], 'r', 'LineWidth',2);
%hold on
%plot([Tid_real(boundery), Tid_real(end)], [mean_CI, mean_CI], 'r', 'LineWidth',2);
%tids_array(end)
hold off
xlabel('Tid (s)');
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

%% PLOT
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
{
figure('Name','Forhold mellem MOR og CI','NumberTitle','off')
plot(Tid_real(:),FORHOLD_MOR_CI(1,:),'m','LineWidth',2)
xlabel('Tid (s)');
ylabel('Forhold mellem MOR_t og CI_t');
}

% Lytisk/lysogen main plot
win_ratio_pl=PL_win(:)./(PL_win(:)+PR_win(:));
%{
figure('Name','MAINPLOT: %lytiske cykler','NumberTitle','off')
plot(CI_koncentrationer(:),win_ratio_pl(:),'--bs','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','b',...
                'MarkerSize',5);
xlabel('CI startkoncentration (nM)');
ylabel('% lytiske cykler');
axis([0 length(CI_koncentrationer) 0 1]);
%}

disp(['The total time was approximately:' num2str(toc(tic_value1))]);
disp(['PR er lig: ' num2str(PR_win) '. PL er lig: ' num2str(PL_win) '']);
disp(['Fordelling er lig: ' num2str(win_ratio_pl) '']);
end %End Function