% This is to 
clear; clc;
S = 1:5; % state space
T_max = 25; % maximum time
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
NS = 1e7;
time = T_max;
P_O_time = zeros(length(time),NS);
tic;
for i = 1:length(time)    
    count_success = 0;
    [t,y] = simulateCTMCQMatrix(Q,pi,time(i),NS);
    temp_A_R = 0;
    for j = 1:NS
        temp_t = t{j};
        temp_y = y{j};
        Delta_t = temp_t(2:end) - temp_t(1:end-1);
        temp_A_R = sum(Delta_t.*r(temp_y(1:end-1)));
        P_O_time(i,j) = temp_A_R;
    end
end
elasped_time_MC = toc;
save('result_benchmark.mat','P_O_time','elasped_time_MC');

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
while ((u > s) && (i < length(p))),
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