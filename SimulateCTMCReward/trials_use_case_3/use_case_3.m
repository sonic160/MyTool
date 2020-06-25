% Use case from Cloth (2006).
clear; clc;
addpath('../compare_different_codings')
%% Parameter definition
S = 1:4; % state space
pi = [1,0,0,0]; % initial distribution
% Define the Q-matrix
% NOTE: the paper has an error in Q(2,1)!!!
Q = [0, 3, 6, 1;...
    1, 0, 0, 0;...
    8, 0, 0, 0;...
    1, 0, 0, 0];
for i = 1:4
    Q(i,i) = -1*sum(Q(i,:));
end
r = [50,20,100,0]; % Reward matrix
% Defining the structure para
para.S = S;
para.pi = pi;
para.Q = Q;
para.r = r;
x_max = 5; % Limit x values
T_max = .2; % Time limit

% r = r/2;
% x_max = x_max/2;

Delta = 1e-4;
tic;
result = accumulated_reward_numit_ext_markov(T_max,x_max,Delta,para)
toc

% %% Consider tolerance
% Delta = 1e-3.*(1/2).^(0:8); % Search space.
% tol = 1e-5;
% 
% tic;
% [cdf, epsilon] = accumulated_reward_numit_ext_markov(T_max,x_max,Delta,para,tol)
% toc
% save('result');