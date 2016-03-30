%% SYNOPSIS
% This function calculates the promoter activities for PR and PL in TP901-1 bacteriophage. 
% The models differential equations was established under the supervision of Kim Sneppen and Anne Alsing (NBI).
% The models differential equations was solved in Mapple. 
% The solutions to the equations were copied into MATALB and the expressions are updated iteratively during the simulation time.
% The biochemical constants are taken from the work of A. K. Alsings (Master thesis page 43 and 48). See also her papers.

%% Matlab Function in Script
% SEE: http://stackoverflow.com/questions/5363397/in-matlab-can-i-have-a-script-and-a-function-definition-in-the-same-file


function [time_steps, time_real, ...
    CI_total, MOR_total, ...
    CI_m, CI_d, MOR_m, CIMOR, ...
    RATIO_MORtoCI, ...
    PR_activity, PL_activity, ...
    decision_time_step, decision_time_real, ...
    PL_win, PR_win] = simulatePhageLifeCycle_GillespieAlgo(n_time_steps,...
                                      threshold_time_switch_decision, ...
                                      threshold_winning_ratio_MORtoCI, ...
                                      CI_initial_concentration, ...
                                      MOR_initial_concentration, ...
                                      PSR_PL2PR, ...
                                      p_include_O_M_site)


%% Explanation of arguments
% threshold_time_switch_decision: [in seconds]. threshold for when the decision is made and the lytic/lysogen cyclus is determined
% ... the other arguments are self explanatory.
          
%% GLOBAL VARIABLES

global kbT;
global G_or;
global G_ol;
global G_od;
global G_om;
global G_pair;
global G_trip;
global K_d;
global K_CIMOR;
global n_1;
global r_CI;
global tau_CI;
global n_2;
global n_3;
global r_MOR;
global tau_MOR;
global n_4;

%% Initiering af arrays
% Antal tidstrin
time_steps=[1:1:n_time_steps]; % this array contains the time steps (sequential, increments by 1)
time_real=zeros(1,length(time_steps));

% Protein concentrations
CI_total=zeros(1,length(time_steps)); 
MOR_total=zeros(1,length(time_steps));
CI_m=zeros(1,length(time_steps));
CI_d=zeros(1,length(time_steps));
MOR_m=zeros(1,length(time_steps));
CIMOR=zeros(1,length(time_steps));
RATIO_MORtoCI=zeros(1,length(time_steps)); % plotting of ratio

% Promoter activity levels
PR_activity=zeros(1,length(time_steps));
PL_activity=zeros(1,length(time_steps));

% *Initial conditions*
CI_total(1) = CI_initial_concentration;
MOR_total(1) = MOR_initial_concentration;


PL_win = 0;
PR_win = 0;

% INITIALYZING FLAGS:
% Idea is that initial value is 0 and when found then setting value to 1
flag_decision_made = 0; % this flag is "raised" (i.e. a value of 1) once the decision is made (e.g. time exceeding the threshold in threshold_time_switch_decision)


for t=time_steps % looping over time steps
    % UDREGNING AF CI, MOR OG CIMOR KONCENTRATIONER
    %-------------------------------------------------------%
    % MOR isoleret af lign2 og indsat i lign1. Derefter solving for CI_m.
    % OK LØSNING, men giver komplekse tal (ikke noget problem)
    CI_m(t)=(1/6)*(-9*CI_total(t)*K_d^2+36*CI_total(t)*K_d*K_CIMOR+3*K_d^2*K_CIMOR+6*K_d*K_CIMOR^2+9*MOR_total(t)*K_d^2+18*MOR_total(t)*K_d*K_CIMOR-K_d^3-8*K_CIMOR^3+3*sqrt(-12*CI_total(t)*K_d^3*K_CIMOR*MOR_total(t)+240*CI_total(t)*K_d^2*MOR_total(t)*K_CIMOR^2-24*CI_total(t)^3*K_d^3-3*CI_total(t)^2*K_d^4+12*K_d^3*K_CIMOR^3-3*K_d^4*K_CIMOR^2-12*K_d^2*K_CIMOR^4+24*MOR_total(t)^3*K_d^3-3*MOR_total(t)^2*K_d^4+96*CI_total(t)^2*K_d^2*K_CIMOR^2-12*CI_total(t)*K_d^3*K_CIMOR^2-96*CI_total(t)*K_d*K_CIMOR^4-6*CI_total(t)*K_d^4*K_CIMOR+96*CI_total(t)*K_d^2*K_CIMOR^3-72*CI_total(t)*K_d^3*MOR_total(t)^2+6*CI_total(t)*K_d^4*MOR_total(t)-48*CI_total(t)^2*K_d^3*K_CIMOR+72*CI_total(t)^2*K_d^3*MOR_total(t)+60*K_d^3*K_CIMOR*MOR_total(t)^2-6*K_d^4*K_CIMOR*MOR_total(t)+48*K_d^3*K_CIMOR^2*MOR_total(t)-24*K_d^2*K_CIMOR^3*MOR_total(t)-12*MOR_total(t)^2*K_d^2*K_CIMOR^2))^(1/3)-(6*(-(1/6)*CI_total(t)*K_d+(1/18)*K_d*K_CIMOR+(1/6)*MOR_total(t)*K_d-(1/36)*K_d^2-(1/9)*K_CIMOR^2))/(-9*CI_total(t)*K_d^2+36*CI_total(t)*K_d*K_CIMOR+3*K_d^2*K_CIMOR+6*K_d*K_CIMOR^2+9*MOR_total(t)*K_d^2+18*MOR_total(t)*K_d*K_CIMOR-K_d^3-8*K_CIMOR^3+3*sqrt(-12*CI_total(t)*K_d^3*K_CIMOR*MOR_total(t)+240*CI_total(t)*K_d^2*MOR_total(t)*K_CIMOR^2-24*CI_total(t)^3*K_d^3-3*CI_total(t)^2*K_d^4+12*K_d^3*K_CIMOR^3-3*K_d^4*K_CIMOR^2-12*K_d^2*K_CIMOR^4+24*MOR_total(t)^3*K_d^3-3*MOR_total(t)^2*K_d^4+96*CI_total(t)^2*K_d^2*K_CIMOR^2-12*CI_total(t)*K_d^3*K_CIMOR^2-96*CI_total(t)*K_d*K_CIMOR^4-6*CI_total(t)*K_d^4*K_CIMOR+96*CI_total(t)*K_d^2*K_CIMOR^3-72*CI_total(t)*K_d^3*MOR_total(t)^2+6*CI_total(t)*K_d^4*MOR_total(t)-48*CI_total(t)^2*K_d^3*K_CIMOR+72*CI_total(t)^2*K_d^3*MOR_total(t)+60*K_d^3*K_CIMOR*MOR_total(t)^2-6*K_d^4*K_CIMOR*MOR_total(t)+48*K_d^3*K_CIMOR^2*MOR_total(t)-24*K_d^2*K_CIMOR^3*MOR_total(t)-12*MOR_total(t)^2*K_d^2*K_CIMOR^2))^(1/3)-(1/6)*K_d-(1/3)*K_CIMOR;
    CI_m(t)=real(CI_m(t)); % Dont worry about this... K. Sneppen says it is okay.
    
    % De følgende to måde at regne MOR_m ud på, giver tilsyneladende sammme plot
    % MOR_m=-(CI_m(1,j)*K_d+2*CI_m(1,j)^2-CI_total*K_d)*K_CIMOR/(K_d*CI_m(1,j)); % isolering af MOR_m fra lign1
    CI_d(t)=(CI_m(t)^2)/K_d;
    MOR_m(t)=MOR_total(t)*K_CIMOR/(K_CIMOR+CI_m(t)); % BRUG DENNE, da dette udtryk er indsat i monosolveren. Isolering af MOR_m fra lign2
    CIMOR(t)=CI_m(t)*MOR_m(t)/K_CIMOR;
    %-------------------------------------------------------%
    
    % WEIGHTS
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
        w_length=30; % set to 30 ==> WITH O_M
    else
        w_length=15; % set to 15 ==> WITHOUT O_M
    end
    
    % Making a weight table
    pr=zeros(1,30);
    pr_on=[1 3 4 7 11];
    pr(1,pr_on)=1;
    
    pl=zeros(1,30);
    pl_on=[1 4 16 19];
    pl(1,pl_on)=1;
    
    % Finding the weights
    % TODO: use vectorized operation instead of for-loop.
    w_total=[0;0]; % First row is for PL and second is for PR
    for k1=[1:w_length];
        w_total(1)=pl(k1)*w(k1)+w_total(1);
        w_total(2)=pr(k1)*w(k1)+w_total(2);
    end
    % Calculating the promoter activity.
    Z=sum(w(1:w_length));
    % PL_activity(t)=(w_total(1)/Z)*beta_PL; % old way of controling the promoter strengths
    % PR_activity(t)=(w_total(2)/Z)*beta_PR; % old way of controling the promoter strengths
    PL_activity(t)=(w_total(1)/Z)*PSR_PL2PR; % *IMPORTANT* scaling PL_activity according to PSR_PL2PR
    PR_activity(t)=(w_total(2)/Z);
    
    if t+1<=length(time_steps) % avoiding out of bounds errors? check if we can update once more
        % Gillespie algorithm
        a_1=PR_activity(t)*r_CI;
        a_2=CI_total(t)/(0.831*10^-9)*tau_CI;  % For E. coli: N_molecules = concentration(nM)/0.831
        a_3=PL_activity(t)*r_MOR;
        a_4=MOR_total(t)/(0.831*10^-9)*tau_MOR; % For E. coli: N_molecules = concentration(nM)/0.831
        
        a_0=a_1+a_2+a_3+a_4;
        r1=random('unif', 0, 1); % random number from uniform distrib | Determines the time to next event.
        r2=random('unif', 0, 1); % random number from uniform distrib. | Determines the next event.
        time_step_gillespie=(1/a_0)*log(1/r1); % determining time to next event
        time_real(t+1)=time_real(t)+time_step_gillespie; % setting time for next event
        
        % Hilberts way of determining which j-value to use
        if r2*a_0<=a_1
            event=1; % CI production
        elseif r2*a_0<=a_1+a_2
            event=2; % CI degradation
        elseif r2*a_0<=a_1+a_2+a_3
            event=3; % MOR production
        elseif r2*a_0<=a_1+a_2+a_3+a_4
            event=4; % MOR degradation
        else
            error 'ERROR: unexpected conditions in if statements'
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
        
        % Opdating concentrations:
        % *REMEMBER*: the number of proteins produced is converted to concentrations via "N_molecules = concentration(nM)/0.831" (E. coli)
        CI_total(t+1)=CI_total(t)+CI_update*0.831*10^-9;
        MOR_total(t+1)=MOR_total(t)+MOR_update*0.831*10^-9;
        
    end % end if tjek af tidsarray
    
    
    % Saving protein levels
    CI_total_array(t)=CI_total(t);
    MOR_total_array(t)=MOR_total(t);
    RATIO_MORtoCI(t)=MOR_total(t)/CI_total(t);
    

    %% DECISION BETWEEN LYTIC AND LYSOGENIC life cycle
    % Another way of "making the decision" would be to use a threshold for the protein levels.
    % E.g. 295 <= CI_total <= 305
    % flag_decision_made ensures that we only enter this statement ONCE! 
    %  --> i.e. the decision between lytic and lysogen cyclus is made only ONCE
    if time_real(t) >= threshold_time_switch_decision && flag_decision_made == 0
        flag_decision_made = 1;
        decision_time_step = t;
        decision_time_real = time_real(t); % Not really needed, but fun to save.
        if RATIO_MORtoCI(t) >= threshold_winning_ratio_MORtoCI
            PL_win = 1;
        else %RATIO_MORtoCI(t) < threshold_winning_ratio_MORtoCI
            PR_win = 1;
        end
    end

    
end % End for-loop

end % end function