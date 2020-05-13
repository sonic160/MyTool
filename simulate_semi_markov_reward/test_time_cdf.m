% This is to compare use default integral funcion in a loop, and write a
% vectorize function, to evaluation the CDF.

clear; clc;
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

%% Test set-up
n_test = 1e4;
evaluation_horizon = 1e3;
test_x = linspace(0,evaluation_horizon,n_test); % Test points.
cdf_default = zeros(1,n_test); % Initial values for the result: Default method.
% cdf_test = zeros(1,n_test); % Initial values for the result: Tested method.
tau = 0;

cal_Q_01 = @(x,tau) integral(@(t)f_01(t,tau),0,x);
cal_Q_02 = @(x,tau) integral(@(t)f_02(t,tau),0,x);
cal_Q_03 = @(x,tau) integral(@(t)f_03(t,tau),0,x);
cal_Q_20 = @(x,tau) integral(@(t)f_20(t,tau),0,x);
cal_Q_23 = @(x,tau) integral(@(t)f_23(t,tau),0,x);

test_Q = cal_Q_20;
test_f = f_20;

%% Default method
tic;
for i = 1:n_test
    cdf_default(i) = test_Q(test_x(i),tau);
end
toc
figure
plot(test_x,cdf_default,'-k');
hold on;

%% Cumulative integration
tic;
cdf_test = int_cum(test_f,test_x',tau*ones(size(test_x')));
toc
plot(test_x,cdf_test,'--r')

[max_relative_error, position] = max((cdf_test'-cdf_default)./cdf_default)

%% Cumulative integration
function cdf = int_cum(cal_pdf,x,tau)
    pdf = cal_pdf(x,tau);
    cdf = cumtrapz(x,pdf);
end