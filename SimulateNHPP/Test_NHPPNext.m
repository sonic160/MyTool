% This is to test NHPPNext at different values of T_max
clear; clc;
%% Parameter settings
lambda_0 = 2.5e-5;
gamma_a = 1e-5;
beta_s = 1e-4;
phi = 0;
handle_lambda = @(t) (lambda_0 + gamma_a*(phi + beta_s*t)); % intensity
handle_lambda_int = @(t) (lambda_0*t + gamma_a*(phi*t + .5*beta_s*t.^2)); % intensity
%% Test using the distribution of first arrival time
% Using the fact that CDF of waiting time = t is $F_s(t) =
% 1-exp(-\int_s^{s+t}\lambda (u) du)$
NS = 1e5;
t_next_sim = zeros(NS,1);
for i = 1:NS
    disp([num2str(i) '/' num2str(NS)])
    t_next_sim(i) = NHPPNext_thinning(handle_lambda,t_0,T_max);
end
P = (1:NS)/NS;
t_next_sim = sort(t_next_sim);
CDF_true = 1-exp(-1*handle_lambda_int(t)+handle_lambda_int(t_0));
figure
plot(t,CDF_true,'-k')
hold on
plot(t_next_sim,P,'r--')
