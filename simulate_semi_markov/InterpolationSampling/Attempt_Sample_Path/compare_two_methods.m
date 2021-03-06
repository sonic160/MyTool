% This is to test the performance of interpolation-based sampling by comparing it to
% the rejection method.

clear; clc;
load('inv_cdf.mat')

%% Test rejection method
ns = 1e6;
sample_f_20 = zeros(1,ns); % Preallocate memory
tau = 0;

% Define parameters
lambda_D = 1e-3;
pd_disp = makedist('exponential','mu',1/lambda_D);

eta_r = 1; beta_r = 2;
pd_t_r = makedist('weibull','a',eta_r,'b',beta_r); 

F_02 = @(x,tau) (.9984 - .9984*exp(-1*x/1000) + .0887024*exp(tau/1000).*(...
    erf(.05+.01*tau) - erf(.01*(5+x+tau))))./...
    (.9984 + .0887024*exp(tau/1000).*(-1 + erf(.05 + .01*tau)));
F_20 = @(t) (1 - exp(-t/1000 - t.^2) + exp(1/4000000)*sqrt(pi)*erf(1/2000)/2000 -...
    exp(1/4000000)*sqrt(pi)*erf(1/2000+t)/2000)./...
    (1 + (exp(1/4000000)*sqrt(pi)*(-1+erf(1/2000)))/2000);

f_02 = @(x,tau) 0.9984E-3.*exp(1).^((-1/1000).*x).*(1+(-1).*exp(1).^((-1/10000).*...
    (x+tau).^2)).*(0.9984E0+0.887024E-1.*exp(1).^((1/1000).*tau).*...
    ((-0.1E1)+erf(0.5E-1+0.1E-1.*tau))).^(-1);
f_20 = @(x,tau) 0.200177E1.*exp(1).^((-1/1000).*x+(-1).*x.^2).*x;

g_02 = @(t,tau) pd_disp.pdf(t);
g_20 = @(t,tau) pd_t_r.pdf(t);

g_rng_02 = @(tau) pd_disp.random();
g_rng_20 = @(tau) pd_t_r.random();

c_02 = @(tau) 9.9984e-4./(.9984 + .0887024*exp(tau/1000).*(-1 + erf(.05 + .01*tau)));
c_20 = @(tau) 2.00177/2;

f = @(x) f_20(x,tau);
g = @(x) g_20(x,tau);
g_rng = @() g_rng_20(tau);
c = c_20(tau);

rejected_f_20 = 0;
tic;

for i = 1:ns
    accept = false;
    while accept == false
        u = rand();
        v = g_rng();
        if c*u <= f(v)/g(v)
           sample_f_20(i) = v;
           accept = true;
        else
           rejected_f_20 = rejected_f_20 + 1;
        end
    end
end
time_rejection = toc

%% Test interpolation method
cal_inv_cdf_f_20 = @(p) inv_cdf_f_20(p);

tic;
u = rand(ns,1);
sample_f_20_int = cal_inv_cdf_f_20(u);
time_int = toc

% % Test if run at one-at-a-time, what's the performance.
% sample_f_20_int = zeros(1,ns);
% tic
% for i = 1:ns
%     u = rand;
%     sample_f_20_int(i) = cal_inv_cdf_f_20(u);
% end
% time_int = toc

% Compare results
histogram(sample_f_20)
hold on;
histogram(sample_f_20_int)
legend('Rejection method','Interpolation-based sampling')

%% Sampling f_02: rejection method: tau = 0
sample_f_02_0 = zeros(1,ns); % Preallocate memory
tau = 0;
f = @(x) f_02(x,tau);
g = @(x) g_02(x,tau);
g_rng = @() g_rng_02(tau);
c = c_02(tau);

rejected_f_02 = 0;
tic
for i = 1:ns
    accept = false;
    while accept == false
        u = rand();
        v = g_rng();
        if c*u <= f(v)/g(v)
           sample_f_02_0(i) = v;
           accept = true;
        else
           rejected_f_02 = rejected_f_02 + 1;
        end
    end
end
time_rejection_f_02 = toc

%% Test NN
cal_inv_cdf_f_02 = @(p,t_tau) inv_cdf_f_02([p,t_tau]);

tic;
u = rand(ns,1);
t_tau = tau*ones(ns,1);
sample_f_20_0_int = cal_inv_cdf_f_02(u,t_tau);
time_int_f_02 = toc

% Compare results
figure
histogram(sample_f_02_0)
hold on;
histogram(sample_f_20_0_int)
legend('Rejection method','Interpolation-based sampling')

%% Sampling f_02: rejection method: tau = 0
sample_f_02_1000 = zeros(1,ns); % Preallocate memory
tau = 500;
f = @(x) f_02(x,tau);
g = @(x) g_02(x,tau);
g_rng = @() g_rng_02(tau);
c = c_02(tau);

rejected_f_02 = 0;
for i = 1:ns
    accept = false;
    while accept == false
        u = rand();
        v = g_rng();
        if c*u <= f(v)/g(v)
           sample_f_02_1000(i) = v;
           accept = true;
        else
           rejected_f_02 = rejected_f_02 + 1;
        end
    end
end

%% Test NN
u = rand(ns,1);
t_tau = tau*ones(ns,1);
sample_f_02_1000_int = cal_inv_cdf_f_02(u,t_tau);

% Compare results
figure
histogram(sample_f_02_1000)
hold on;
histogram(sample_f_02_1000_int)
legend('Rejection method','NN-based sampling')