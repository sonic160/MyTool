% This is to calculate the cdf at a vector of t.

function [cdf, ci, elasped_time_MC] = create_benchmark_t(t,n_s)
T_max = t(end); % maximum time
x_th = 110;

pi = [1,0,0,0,0]; % initial distribution
% Define the Q-matrix
Q = [0,.1,.5,.25,0;...
    1,0,.5,0,.5;...
    0,.25,0,1,.5;...
    1,0,0,0,.5;...
    1,.25,0,1,0];
for i = 1:5
    Q(i,i) = -1*sum(Q(i,:));
end
r = [5,2,4,5,4]; % Reward matrix

tic;
% Simulation to generate the sample paths.
[sample_jump_time,sample_path] = simulateCTMCQMatrix(Q,pi,T_max,n_s);
% Post-processing.
% Confidence interval
alpha = .05;
cal_sigma_hat = @(N,p) sqrt(1/(N-1).*(p.*(1-p)));
cal_ci = @(N,p,alpha) [p-norminv(1-alpha/2)*cal_sigma_hat(N,p),...
    p+norminv(1-alpha/2)*cal_sigma_hat(N,p)];

cdf = zeros(1,length(t));
ci = zeros(2,length(t));
for i = 1:length(t)
    tt = t(i);
    %     fprintf('%d/%d\n\n',i,length(t));
    count = 0;
    for j = 1:n_s
        % Get the current sample path and jump time.
        temp_temp_sample_path = sample_path{j};
        temp_temp_sample_jump_time = sample_jump_time{j};
        % Get the part ending before tt.
        index = temp_temp_sample_jump_time<tt;
        temp_temp_sample_path = temp_temp_sample_path(index);
        temp_temp_sample_jump_time = temp_temp_sample_jump_time(index);
        temp_sample_path = [temp_temp_sample_path,temp_temp_sample_path(end)];
        temp_sample_jump_time = [temp_temp_sample_jump_time,tt];
        % Get time interval.
        time_interval = temp_sample_jump_time(2:end)-temp_sample_jump_time(1:end-1); 
        % Get the current state.
        state_cur = temp_sample_path(1:end-1);
        % Calculate the accumulated reward.
        count = count + (sum(r(state_cur).*time_interval) <= x_th);
    end
    cdf(i) = count/n_s;
    ci(:,i) = transpose(cal_ci(n_s,cdf(i),alpha));        
end
elasped_time_MC = toc;
end

% This is to simulate a trajectory of a CTMC using Q matrix.
% Input: Q - Transition rate matrix of the embedded chain: n*n
%        pi - Initial distribution, n*1
%        T_max - End time, simulation is performed in (0,T_max)
%        NS - Number of sample paths: default value: NS = 1
% Output: t - times of the jumps, NS*1 cell array, each element is a
%         vector, containing jump times of a particular sample path 
%         y - the state after the jumps NS*1 cell array, each element is a
%         vector, containing jump times of a particular sample path
% History:  20171208 - Created by ZZ
function [t,y] = simulateCTMCQMatrix(Q,pi,T_max,NS)
if nargin == 3 % If NS is not input
    NS = 1; % Default value is NS = 1
end
% Initial values for t and y
t = cell(NS,1);
y = cell(NS,1);
[a,P] = CalculatePEmbedded(Q); % Contruct the embedded chain
% Do the sampling
for i = 1:NS 
    [temp_t,temp_y] = SimulateCTMCEmbedded(a,P,pi,T_max); % Simulate a sample path using the embedded chain
    % Format transform: from numerical to cell
    t(i) = {temp_t};
    y(i) = {temp_y};
end
end

% This is to simulate a trajectory of a CTMC using embedded chain.
% Input: a - Rates of the holding time distributions: n*1
%        P - Transition matrix of the embedded chain: n*n
%        pi - Initial distribution, n*1
%        T_max - End time, simulation is performed in (0,T_max)
% Output: t - times of the jumps
%         y - the state after the jumps
% History:  20171208 - Modified by ZZ
%                      Fix errors in the annotations in the beginning of
%                      the code
%           20170205 - Modified by ZZ
%                      Fix if t_arrive > T_max
%           20170203 - Created by ZZ
function [t,y] = SimulateCTMCEmbedded(a,P,pi,T_max)
% Simulate the embedded chain and holding time
t(1) = 0; % start times at 0
y(1) = rando(pi); % generate the initial y
i = 1;
while 1
    t(i+1) = t(i) + exprnd(1/a(y(i))); % Next transitions
    if t(i+1) < T_max
        y(i+1) = rando(P(y(i),:)); % Simulate next state using the transitin matrix of the embedded chain
        i=i+1;
    else
        t(i+1) = T_max;
        y(i+1) = y(i);
        break;
    end
end
end

%  rando.m generates a random variable in 1, 2, ..., n given a distribution 
%  vector. 
function index = rando(p)
u = rand;
i = 1;
s = p(1);
while ((u > s) && (i < length(p)))
    i=i+1;
    s=s+p(i);
end
index=i;
end

% This is to generate the embedded chain from Q matrix
% Output: a - Rates of the holding time distributions
%         P - Transition matrix of the embedded chain
% History: 20170203, Created by ZZ
function [a,P] = CalculatePEmbedded(Q)
a = -1*diag(Q); % Calculate the rates for the holding time distributions
n_state = length(a); % Number of states
% Construct the transition matrix for the embedded chain
P = Q./repmat(a,1,n_state);
for i = 1:n_state
    P(i,i) = 0;
end
end