clear; clc;
addpath('../simulate_semi_markov/InterpolationSampling')
addpath('../simulate_semi_markov/')

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
t = linspace(.1,evaluation_horizon,10); % Evaluate time

r = [4,3,2,1];
y_th = 3000;
F_nuit = zeros(1,length(t));

ns = 1e6; % Sample size

% Elements of f_ij
f_01 = @(x,tau) -lambda_D.*exp(-lambda_D.*x).*exp(-((tau + x)/eta_s).^beta_s).*(p_large - 1);
f_02 = @(x,tau) lambda_D.*exp(-lambda_D.*x).*(exp(-((tau + x)/eta_s).^beta_s) - 1).*(p_large - 1);
f_03 = @(x,tau) lambda_D.*p_large.*exp(-lambda_D.*x);
f_20 = @(x,tau) (beta_r.*exp(-(x./eta_r).^beta_r).*exp(-lambda_D.*x).*(x./eta_r).^(beta_r - 1))./eta_r;
f_23 = @(x,tau) lambda_D.*exp(-(x./eta_r).^beta_r).*exp(-lambda_D.*x);

cal_f_0 = @(x,tau) [zeros(size(tau)), f_01(x,tau), f_02(x,tau), f_03(x,tau)];
cal_f_1 = @(x,tau) [f_01(x,tau), zeros(size(tau)), f_02(x,tau), f_03(x,tau)];
cal_f_2 = @(x,tau) [f_20(x,tau), zeros(size(tau)), zeros(size(tau)), f_23(x,tau)];
cal_f_3 = @(x,tau) [zeros(size(tau)), zeros(size(tau)),...
    zeros(size(tau)), zeros(size(tau))];
cal_d_Q = @(x,tau) [cal_f_0(x,tau); cal_f_1(x,tau); cal_f_2(x,tau); cal_f_3(x,tau)];

% D
cal_D = @(x,tau) diag([exp(-lambda_D.*x),...
    exp(-lambda_D.*x),...
    exp(-(x./eta_r).^beta_r).*exp(-lambda_D.*x),...
    1]);

%% Do the numerical integration.
% Write in para.
para.pi_0 = [1,0,0,0]; % Initial distribution
para.cal_d_Q = cal_d_Q;
para.cal_D = cal_D;
Delta = .1;

% Run the simulation.
state_probability = cal_state_p_numerical_int(evaluation_horizon,Delta,para)

load('result_ref.mat');

% Display results
p_t_ref(:,end)