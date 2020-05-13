clear; clc;
addpath('../simulate_semi_markov/InterpolationSampling')
addpath('../simulate_semi_markov/')
load('result_ref.mat'); % Load the benchmark.

%% Defining the renewal kernals.
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
t = linspace(.1,evaluation_horizon,20); % Evaluate time

r = [4,3,2,1];
y_th = 3000;
F_nuit = zeros(1,length(t));

ns = 1e6; % Sample size

% Define the embedded chain and the holding time dist.
% Elements in the transition probability matrix
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

% % Handles to the transition probability matrix for each state: Returns a
% % n_t*4 matrix, where each row is the transition probability vector
% % p_{i,j}, j = 0,1,2,3, and the column represents the current time tau.
% cal_p_tr_0 = @(tau) [zeros(size(tau)), p_01(tau), p_02(tau), p_03(tau)];
% cal_p_tr_1 = @(tau) [p_01(tau), zeros(size(tau)), p_02(tau), p_03(tau)];
% cal_p_tr_2 = @(tau) [p_20(tau), zeros(size(tau)), zeros(size(tau)), p_23(tau)];
% cal_p_tr_3 = @(tau) [zeros(size(tau)), zeros(size(tau)),...
%     zeros(size(tau)), ones(size(tau))];
% cal_p = {cal_p_tr_0, cal_p_tr_1, cal_p_tr_2, cal_p_tr_3};

% Elements of f_ij
f_01 = @(x,tau) -lambda_D.*exp(-lambda_D.*x).*exp(-((tau + x)/eta_s).^beta_s).*(p_large - 1);
f_02 = @(x,tau) lambda_D.*exp(-lambda_D.*x).*(exp(-((tau + x)/eta_s).^beta_s) - 1).*(p_large - 1);
f_03 = @(x,tau) lambda_D.*p_large.*exp(-lambda_D.*x);
f_20 = @(x,tau) (beta_r.*exp(-lambda_D.*x).*exp(-((tau + x)/eta_r).^beta_r).*((tau + x)/eta_r).^(beta_r - 1))/eta_r;
f_23 = @(x,tau) lambda_D.*exp(-lambda_D.*x).*exp(-((tau + x)/eta_r).^beta_r);

cal_f_0 = @(x,tau) [zeros(size(tau)), f_01(x,tau), f_02(x,tau), f_03(x,tau)];
cal_f_1 = @(x,tau) [f_01(x,tau), zeros(size(tau)), f_02(x,tau), f_03(x,tau)];
cal_f_2 = @(x,tau) [f_20(x,tau), zeros(size(tau)), zeros(size(tau)), f_23(x,tau)];
cal_f_3 = @(x,tau) [zeros(size(tau)), zeros(size(tau)),...
    zeros(size(tau)), ones(size(tau))];
cal_f = {cal_f_0, cal_f_1, cal_f_2, cal_f_3};

%% Do the numerical integration.
% Write in para.
para.pi_0 = [1,0,0,0]; % Initial distribution
para.S = [1,2,3,4]; % State space
para.r = r; % Rewards per unit time of sojourn in each state
para.cal_p = cal_p;
para.cal_f = cal_f;

Delta = 1;

% Run the simulation.
tic;
f_appr = cal_reward_numerical_int(evaluation_horizon,y_th,Delta,para);
%% Post-processing.
n_t = 20;
t_int = Delta:Delta:evaluation_horizon;
for i = 2:n_t
    temp_t = t(i);
    index = find(t_int<temp_t,1,'last');
    F_nuit(i) = sum(sum(f_appr(:,:,index)*Delta));
end
elasped_time_nuit = toc;