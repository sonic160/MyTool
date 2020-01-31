% This is to compare three ways of generating samples: rejection method,
% NN-based, and interpolation-based method.
% Benchmark: A four state system, where state transitions are caused by
% arrival of disruptive events, performance of safety barriers, and
% recovery time distributions.
%            State 0: Perfect
%            State 1: Dummy perfect state.
%            State 2: Performance degradation.
%            State 3: Total failure (absorption state).
% Last edit: 20200129 - Created by ZZ

clear; clc;

%% Parameter definitions
ns = 1e7; % Sample size.
tau_lim = 1e3;
tau = rand(1,ns)*tau_lim; % Current time.

% Parameters for the benchmark case study.
lambda_D = 1e-3;
pd_disp = makedist('exponential','mu',1/lambda_D);

eta_r = 1; beta_r = 2;
pd_t_r = makedist('weibull','a',eta_r,'b',beta_r); 

F_02 = @(x,tau) (.9984 - .9984*exp(-1*x/1000) + .0887024*exp(tau/1000).*(...
    erf(.05+.01*tau) - erf(.01*(5+x+tau))))./...
    (.9984 + .0887024*exp(tau/1000).*(-1 + erf(.05 + .01*tau)));
f_02 = @(x,tau) 0.9984E-3.*exp(1).^((-1/1000).*x).*(1+(-1).*exp(1).^((-1/10000).*...
    (x+tau).^2)).*(0.9984E0+0.887024E-1.*exp(1).^((1/1000).*tau).*...
    ((-0.1E1)+erf(0.5E-1+0.1E-1.*tau))).^(-1);

%% Rejection method
% Parameters
g_02 = @(t) pd_disp.pdf(t);
g_rng_02 = @(m,n) pd_disp.random(m,n);
c_02 = @(tau_tau) 9.9984e-4./(.9984 + .0887024*exp(tau_tau/1000).*(-1 + erf(.05 + .01*tau_tau)));

% This is a function to do sampling using rejection method.
% Inputs: f - pdf.
%         g - The function that bounds f.
%         g_rng - Random number generator for g.
%         c - Constant that satisfies f <= c*g.
%         ns - Number of samples to be generated.
% Outputs: X - The random numbers.

f = @(x) f_02(x,tau);
g = g_02;
g_rng = g_rng_02;
c = c_02(tau);

tic;
sample_rej = zeros(1,ns); % Preallocate memory

index = 0;
counter = 0;
while index < ns && counter < 1000
    u = rand(1,ns); % Random number from U[0,1].
    x = g_rng(1,ns); % Candidates.
    x = x(c.*u <= f(x)./g(x)); % Only keep the accepted ones.
    sample_rej(index+1:min([index+length(x),ns])) =...
        x(1:min([length(x), ns-index]));
    index = index + length(x);
    counter = counter+1;
end
if counter == 1000
    error('No sample generated.')
end
time_rej = toc

%% NN training
%% Training data generator
t_lim = 2e4;
% t_f_02 = linspace(0,t_lim,400);
t_f_02 = [linspace(0,200,400),linspace(200,t_lim,200)];
tau_f_02 = linspace(0,tau_lim,100);
% Generate training data
[mesh_t_f_02, mesh_tau_f_02] = meshgrid(t_f_02,tau_f_02);
mesh_p_f_02 = F_02(mesh_t_f_02, mesh_tau_f_02);

% Divide the training data.
p_th_f_02_0 = 5e-3;
p_th_f_02_1 = 1e-1;
p_th_f_02_2 = .999;
% p_th_f_02_0 = .999;
% p_th_f_02_1 = 1;
% % Sort the p of the training samples.
training_t_f_02 = transpose(mesh_t_f_02(:));
training_p_f_02 = transpose(mesh_p_f_02(:));
training_tau_f_02 = transpose(mesh_tau_f_02(:));

%% Training data for the first range
index_0 = training_p_f_02<p_th_f_02_0;
training_data_0 = [training_p_f_02(index_0);...
    training_tau_f_02(index_0)];
training_t_0 = training_t_f_02(index_0);
% Training data for the lower range
index_1 = ~index_0 & training_p_f_02<p_th_f_02_1;
training_data_1 = [training_p_f_02(index_1);...
    training_tau_f_02(index_1)];
training_t_1 = training_t_f_02(index_1);
% Training data for the middle range
index_2 = ~(index_0 | index_1) & training_p_f_02<p_th_f_02_2;
training_data_2 = [training_p_f_02(index_2);...
    training_tau_f_02(index_2)];
training_t_2 = training_t_f_02(index_2);
% Training data for the upper range
index_3 = ~(index_0 | index_1 | index_2);
training_data_3 = [training_p_f_02(index_3);...
    training_tau_f_02(index_3)];
training_t_3 = training_t_f_02(index_3);

% Training for the first region
net_f_02_0 = fitnet(10);
net_f_02_0 = train(net_f_02_0, training_data_0, training_t_0);
% Training lower reange
net_f_02_1 = fitnet(10);
net_f_02_1 = train(net_f_02_1, training_data_1, training_t_1);
% Training lower reange
net_f_02_2 = fitnet(10);
net_f_02_2 = train(net_f_02_2, training_data_2, training_t_2);
% Training lower reange
net_f_02_3 = fitnet(10);
net_f_02_3 = train(net_f_02_3, training_data_3, training_t_3);

%% NN testing
tic;
u = rand(1,ns);
u = sort(u);
index_0 = u<p_th_f_02_0;
index_1 = ~index_0 & u<p_th_f_02_1;
index_2 = ~(index_0 | index_1) & u<p_th_f_02_2;
index_3 = ~(index_0 | index_1 | index_2);
sample_NN = zeros(1,ns);
sample_NN(index_0) = net_f_02_0([u(index_0); tau(index_0)]);
sample_NN(index_1) = net_f_02_1([u(index_1); tau(index_1)]);
sample_NN(index_2) = net_f_02_2([u(index_2); tau(index_2)]);
sample_NN(index_3) = net_f_02_3([u(index_3); tau(index_3)]);
time_NN = toc

%% Training for interpolation-based methods.
% Train a network using all the data.
training_data_f_02 = [mesh_p_f_02(:), mesh_tau_f_02(:)];
training_t_f_02 = mesh_t_f_02(:);

% Delete repeted samples
[training_data_f_02, i_unique_p] = unique(training_data_f_02,'rows');
training_t_f_02 = training_t_f_02(i_unique_p);

% Training
[inv_cdf_f_02,gof_f_02] = fit(training_data_f_02,training_t_f_02,'linearinterp');

%% Interpolation testing
tic;
u = rand(1,ns);
u = sort(u);
sample_Int = feval(inv_cdf_f_02,[u', tau']);
time_Int = toc

%% Compare results
figure
histogram(sample_rej)
hold on;
histogram(sample_NN)
hold on;
histogram(sample_Int)
legend('Rejection method','NN-based sampling','Interpolation-based sampling')

mean(sample_rej)
mean(sample_NN)
mean(sample_Int(~isnan(sample_Int)))

std(sample_rej)
std(sample_NN)
std(sample_Int(~isnan(sample_Int)))