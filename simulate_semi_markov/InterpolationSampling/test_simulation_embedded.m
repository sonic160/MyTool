% This is a test of the simulation algorithm for inhomogeneous semi-Markov process.
% Invere cdf based on interpolation is used to sample the holding time distribution.
% Benchmark: A four state system, where state transitions are caused by
% arrival of disruptive events, performance of safety barriers, and
% recovery time distributions.
%            State 0: Perfect
%            State 1: Dummy perfect state.
%            State 2: Performance degradation.
%            State 3: Total failure (absorption state).
% Last edit: 20191223 - Created by ZZ
clear; clc;
%% Parameter definition
% First arrival time of the disruptive event
lambda_D = 1e-3;
pd_disp = makedist('exponential','mu',1/lambda_D);
% Safety barriers
p_large = 1.6e-3; % Probability of large damage, leading to state 3
eta_s = 100; beta_s = 2;
pd_t_f_s = makedist('weibull','a',eta_s,'b',beta_s); % Time to failure distribution of the safety barrier
% Time to repair distribution
eta_r = 1; beta_r = 2;
pd_t_r = makedist('weibull','a',eta_r,'b',beta_r); 

% Simulation to calculate p(t) for each state, for 0<t<T_max
evaluation_horizon = 1e3;
t = linspace(.1,evaluation_horizon,10); % Evaluate time

ns = 1e5; % Sample size

%% Simulate using holding time distribution and embedded chain. Sampling with rejection method.
load('inv_cdf.mat','inv_cdf_f_02', 'inv_cdf_f_20');

% Define the embedded chain and the holding time dist.
pi_0 = [1,0,0,0]; % Initial distribution

p_01 = @(tau) (-0.887024E-1).*exp(1).^((1/1000).*tau).*...
    ((-0.1E1)+erf((1/100).*(5+tau)));
p_02 = @(tau) 0.9984E0+0.887024E-1.*exp(1).^((1/1000).*tau).*((-0.1E1)+...
    erf(0.5E-1+ 0.1E-1.*tau));
% p_03 = @(tau) .16e-2*ones(size(tau));
p_03 = @(tau) 1 - p_01(tau) - p_02(tau);
p_10 = p_01;
p_12 = p_02;
p_13 = p_03;
p_20 = @(tau) .999114*ones(size(tau));
p_23 = @(tau) 1 - p_20(tau);

% Transition probability for state 1
cal_p_tr_0 = @(tau) [zeros(size(tau)), p_01(tau), p_02(tau), p_03(tau)];
cal_p_tr_1 = @(tau) [p_01(tau), zeros(size(tau)), p_02(tau), p_03(tau)];
cal_p_tr_2 = @(tau) [p_20(tau), zeros(size(tau)), zeros(size(tau)), p_23(tau)];
cal_p_tr_3 = @(tau) [zeros(size(tau)), zeros(size(tau)),...
    zeros(size(tau)), ones(size(tau))];
cal_p = {cal_p_tr_0, cal_p_tr_1, cal_p_tr_2, cal_p_tr_3};

% Inverse functions of the CDF.
inv_01 = @(p,tau) -5 - tau + 100*erfinv(p + (1-p).*erf(.05+.01*tau));
inv_23 = @(p,tau) (-1 + 2000*erfinv(p + erf(1/2000) - p.*erf(1/2000)))/2000;
inv_03 = @(p,tau) -1e3*log(1-p);

inv_02 = @(p,tau) inv_cdf_f_02([p,tau]);
inv_20 = @(p,tau) inv_cdf_f_20(p);

handle_zero = @(p,tau) 0;
inv_F_t = {handle_zero, inv_01, inv_02, inv_03;...
    inv_01, handle_zero, inv_02, inv_03;...
    inv_20, handle_zero, handle_zero, inv_23;...
    handle_zero, handle_zero, handle_zero,handle_zero};

% Calculate the same probability distribution
fprintf('\n Embedded chain approach\n')
n_state = 4;
p_t_sim = zeros(n_state,length(t)); % p(t): 4*1, each column corresponds to one state

n_jump = 15;

tic;

[y, tau] = simulate_semi_nhsmp_embedded(ns, n_jump, evaluation_horizon, ...
    pi_0, cal_p, inv_F_t);

for i = 1:length(t)
    tt = t(i);
    p_tt_j = zeros(n_state,1);
    
    % n(t) = 0.
    j = 1;
    index_tt = find(tau(:,j+1)>=tt);
    tau_tt = tau(index_tt,j);
    y_tt = y(index_tt,j);
    bin_edge = 1:n_state+1;
    p_tt_j = p_tt_j + transpose(histcounts(y_tt, bin_edge))/n_jump/ns;
    
    for j = 2:n_jump+1
        index_tt = find(tau(:,j)<tt);
        if isempty(index_tt)
            break;
        else
            tau_tt = tau(index_tt,j);
            y_tt = y(index_tt,j);
            bin_edge = 1:n_state+1;
            p_tt_j = p_tt_j + transpose(histcounts(y_tt, bin_edge))/n_jump/ns;        
        end
    end
    p_t_sim(:,i) = p_tt_j;
end

toc


 
    


