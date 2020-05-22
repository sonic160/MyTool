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

% Elements of f_ij
f_01 = @(x,tau) -lambda_D.*exp(-lambda_D.*x).*exp(-((tau + x)/eta_s).^beta_s).*(p_large - 1);
f_02 = @(x,tau) lambda_D.*exp(-lambda_D.*x).*(exp(-((tau + x)/eta_s).^beta_s) - 1).*(p_large - 1);
f_03 = @(x,tau) lambda_D.*p_large.*exp(-lambda_D.*x);
f_20 = @(x,tau) (beta_r.*exp(-(x./eta_r).^beta_r).*exp(-lambda_D.*x).*(x./eta_r).^(beta_r - 1))./eta_r;
f_23 = @(x,tau) lambda_D.*exp(-(x./eta_r).^beta_r).*exp(-lambda_D.*x);

% cal_Q_01 = @(x,tau) integral(@(t)f_01(t,tau),0,x);
% cal_Q_02 = @(x,tau) integral(@(t)f_02(t,tau),0,x);
% cal_Q_03 = @(x,tau) integral(@(t)f_03(t,tau),0,x);
% cal_Q_20 = @(x,tau) integral(@(t)f_20(t,tau),0,x);
% cal_Q_23 = @(x,tau) integral(@(t)f_23(t,tau),0,x);

% cal_Q = @(x,tau) [0,cal_Q_01(x,tau),cal_Q_02(x,tau),cal_Q_03(x,tau);...
%     cal_Q_01(x,tau),0,cal_Q_02(x,tau),cal_Q_03(x,tau);...
%     cal_Q_20(x,tau),0,0,cal_Q_23(x,tau);...
%     0,0,0,1];

% D
cal_D = @(x,tau) diag([exp(-lambda_D.*x),...
    exp(-lambda_D.*x),...
    exp(-(x./eta_r).^beta_r).*exp(-lambda_D.*x),...
    1]);

%% Do the numerical integration.
% Write in para.
para.pi_0 = [1,0,0,0]; % Initial distribution
para.cal_D = cal_D;
Delta = .5;

cal_Q_01 = @(x_0,x_1,tau) integral(@(t)f_01(t,tau),x_0,x_1);
cal_Q_02 = @(x_0,x_1,tau) integral(@(t)f_02(t,tau),x_0,x_1);
cal_Q_03 = @(x_0,x_1,tau) integral(@(t)f_03(t,tau),x_0,x_1);
cal_Q_20 = @(x_0,x_1,tau) integral(@(t)f_20(t,tau),x_0,x_1);
cal_Q_23 = @(x_0,x_1,tau) integral(@(t)f_23(t,tau),x_0,x_1);

cal_Q = @(x_0,x_1,tau) [0,f_01(x_1,tau)*Delta,f_02(x_1,tau)*Delta,f_03(x_1,tau)*Delta;...
    f_01(x_1,tau)*Delta,0,f_02(x_1,tau)*Delta,f_03(x_1,tau)*Delta;...
    cal_Q_20(x_0,x_1,tau),0,0,cal_Q_23(x_0,x_1,tau);...
    0,0,0,0];
para.cal_Q = cal_Q;

% Run the simulation.
tic;
state_probability = cal_state_p_numerical_int_use_Q(evaluation_horizon,Delta,para)
toc

load('result_ref.mat');

% Display results
p_t_ref(:,end)