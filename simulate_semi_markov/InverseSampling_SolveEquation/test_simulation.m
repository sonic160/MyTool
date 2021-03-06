% This is a test of the simulation algorithm.
% Benchmark: A two state system.
% Time to failure: Weibull, eta = 10, beta = 1.1
% Time to repair: Lognormal, mu = 1, sigma = 0.1
% Last edit: 20191223 - Update the rejection sampling method.
%                     - Now use proposal density, rather than a rectangle
%                     area.
%            20191221 - Created by ZZ

clear; clc;
%% Create benchmark
a = 10; b = 1.1;
pd_t_f = makedist('weibull','a',a,'b',b); % Time to failure distribution
mu = 1; sigma = .1;
pd_t_r = makedist('lognormal','mu',mu,'sigma',sigma); % Time to repair distribution
T_max = 20;
t = linspace(.1,T_max,10); % Evaluate time
ns = 1e3; % Sample size
A_t_ref = zeros(1,length(t)); % Availability
for i = 1:length(t)
    fprintf('%d/%d\n',i,length(t));
    count = 0; % Counter
    for j = 1:ns
        t_cur = 0; % Initial value of time
        state = 0; % Start from working
        while 1
            if state == 0 % if current state is working
                t_next = pd_t_f.random(); % Generate time to failure
                t_cur = t_cur + t_next; % Update current time
                if t_cur > t(i) % If beyond simulation horizon
                    count = count + 1;
                    break;
                else
                    state = 1; % System fails
                end
            else
                t_next = pd_t_r.random(); % Generate time to failure
                t_cur = t_cur + t_next; % Update current time
                if t_cur > t(i) % If beyond simulation horizon
                    break;
                else
                    state = 0; % System fails
                end
            end
        end       
    end
    A_t_ref(i) = count/ns;
end
figure
plot(t,A_t_ref,'-o')
hold on
%% Run the simulation using the algorithm being tested - CDF approach
pi = [1,0];
Q_01 = @(t) wblcdf(t,a,b);
Q_10 = @(t) logncdf(t,mu,sigma);
Q_t = @(t) [0, Q_01(t);...
    Q_10(t),0];
P = [0, 1;...
    1, 0];
% Calculate the same availability
fprintf('\n CDF approach\n')
A_t_sim = zeros(1,length(t)); % Availability
for i = 1:length(t)
    fprintf('%d/%d\n',i,length(t));
    count = 0; % Counter
    for j = 1:ns
        [temp_t,y] = simulate_semi_markov(Q_t,P,pi,t(i));
        if y(end) == 1
            count = count + 1;
        end
    end
    A_t_sim(i) = count/ns;
end
plot(t,A_t_sim,'-dr')
%% PDF approach
% Defining the pdf of the holding time.
pi = [1,0];
P = [0, 1;...
     1, 0];
% Pdf of the holding time dist.
f_01 = @(t) wblpdf(t,a,b);
f_10 = @(t) lognpdf(t,mu,sigma);
f_t = {0, f_01;...
    f_10,0};
% Proposal density
g_01 = @(t) exppdf(t,10);
g_10 = @(t) normpdf(t,2.7,.35);
g = {0,g_01;...
    g_10,0};
% Random number generator for the proposal density
g_rng_01 = @(t) exprnd(10);
g_rng_10 = @(t) normrnd(2.7,.35);
g_rng = {0,g_rng_01;...
    g_rng_10,0};
% Adjustment factor
c = [0,1.15; 1.35,0];
% Calculate the same availability
fprintf('\n PDF approach\n')
A_t_sim_pdf = zeros(1,length(t)); % Availability
for i = 1:length(t)
    fprintf('%d/%d\n',i,length(t));
    count = 0; % Counter
    for j = 1:ns
        [temp_t,y] = simulate_semi_markov_rejection(f_t,g,g_rng,c,...
            P,pi,t(i));
        if y(end) == 1
            count = count + 1;
        end
    end
    A_t_sim_pdf(i) = count/ns;
end
plot(t,A_t_sim_pdf,'-hg')
