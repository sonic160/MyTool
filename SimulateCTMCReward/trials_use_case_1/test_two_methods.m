% This is to implement the fast algorithm in Tijms and Veldman (2000)
clear; clc;
%% Parameter definition
S = 1:5; % state space
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
x_max = 110; % Limit x values
T_max = 25; % Time limit
offset = -2;
r = r+offset;
x_max = x_max + offset*T_max;
% Defining the structure para
para.S = S;
para.pi = pi;
para.Q = Q;
para.r = r;
addpath('..\');
Delta = 1/4*(1/2).^(0:10);
tic;
result = 1 - accumulated_reward_numit(T_max,x_max,Delta(9),para)
toc
% tic;
% result_direct = 1 - accumulated_reward_numit_direct(T_max,x_max,Delta(9),para)
% toc
tic;
result_markov = 1 - accumulated_reward_numit_ext_markov(T_max,x_max,Delta(9),para)
toc