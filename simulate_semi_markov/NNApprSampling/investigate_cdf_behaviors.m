clear; clc;

%% Generate samples using rejection method.
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

% Rejection sampling for g_02
ns = 1e4;
sample_f_02 = zeros(1,ns); % Preallocate memory

tau = 0;
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
           sample_f_02(i) = v;
           accept = true;
        else
           rejected_f_02 = rejected_f_02 + 1;
        end
    end
end

% Rejection sampling for g_20
ns = 1e4;
sample_f_20 = zeros(1,ns); % Preallocate memory

tau = 1e3;
f = @(x) f_20(x,tau);
g = @(x) g_20(x,tau);
g_rng = @() g_rng_20(tau);
c = c_20(tau);

rejected_f_20 = 0;
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

%% Training inverse CDF for f_20
% Generating training samples in the tail region.
t_th = 3;
training_t_tail = linspace(t_th,4,5e3);
training_p_tail = F_20(training_t_tail);
p_th = training_p_tail(1);

% Training the tail network
net_f_20_tail = fitnet(10);
net_f_20_tail = train(net_f_20_tail,training_p_tail,training_t_tail);

% Generating samples in the upper region
training_t_upper = linspace(0,t_th,5e3);
training_p_upper = F_20(training_t_upper);

% Train the upper network
net_f_20_upper = fitnet(10);
net_f_20_upper = train(net_f_20_upper,training_p_upper,training_t_upper);

%% Testing
% Testing samples
testing_t = sort(sample_f_20);
testing_p = F_20(testing_t);

% Calculate testing values
testing_t_pred = zeros(size(testing_t));
for i = 1:length(testing_t_pred)
    if testing_p(i) < p_th
        testing_t_pred(i) = net_f_20_upper(testing_p(i));
    else
        testing_t_pred(i) = net_f_20_tail(testing_p(i));
    end
end

plot(testing_p,testing_t,'-k',testing_p,testing_t_pred,'or')
perform(net_f_20_upper,testing_t,testing_t_pred)

%% Training inverst CDF for f_02
% Visualize the CDF
t_f_02 = linspace(0,1e4,1e3);
tau_f_02 = linspace(0,1e3,1e3);
[mesh_t_f_02, mesh_tau_f_02] = meshgrid(t_f_02,tau_f_02);
mesh_p_f_02 = F_02(mesh_t_f_02, mesh_tau_f_02);
surf(mesh_p_f_02,mesh_tau_f_02,mesh_t_f_02);
xlabel('p')
ylabel('\tau')
zlabel('t')

% % Train a network using all the data.
% training_t_f_02 = transpose(mesh_t_f_02(:));
% training_p_f_02 = transpose(mesh_p_f_02(:));
% training_tau_f_02 = transpose(mesh_tau_f_02(:));
% 
% net_f_02 = fitnet(10);
% net_f_02 = train(net_f_02,[training_p_f_02; training_tau_f_02], training_t_f_02);

% Train the network in region upper: p<p_th_f_02
% Sort the p of the training samples.
training_t_f_02 = transpose(mesh_t_f_02(:));
training_p_f_02 = transpose(mesh_p_f_02(:));
training_tau_f_02 = transpose(mesh_tau_f_02(:));
[training_p_f_02, sort_p_index] = sort(training_p_f_02);
training_t_f_02 = training_t_f_02(sort_p_index);
training_tau_f_02 = training_tau_f_02(sort_p_index);

% Get the training sample for the upper region.
p_th_f_02 = .999;
index_p_l_p_th = find(training_p_f_02<p_th_f_02, 1, 'last');
training_data_f_02_uppper = [training_p_f_02(1:index_p_l_p_th); ...
    training_tau_f_02(1:index_p_l_p_th)];
training_t_f_02_uppper = training_t_f_02(1:index_p_l_p_th);

% Training the first network
net_f_02_upper = fitnet(10);
net_f_02_upper = train(net_f_02_upper, training_data_f_02_uppper, training_t_f_02_uppper);

% Train the network in region tail: p>p_th_02
training_data_f_02_tail = [training_p_f_02(index_p_l_p_th+1:end); ...
    training_tau_f_02(index_p_l_p_th+1:end)];
training_t_f_02_tail = training_t_f_02(index_p_l_p_th+1:end);

% Training the first network
net_f_02_tail = fitnet(10);
net_f_02_tail = train(net_f_02_tail, training_data_f_02_tail, training_t_f_02_tail);
