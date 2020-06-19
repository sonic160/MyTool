clear; clc;
%% Parameter definition
S = 1:5; % state space
pi = [1,0,0,0,0]; % initial distribution
% Define the Q-matrix
Q = [0, 1.2298e-4, 1.2298e-4, 2.5939e-4, 3.0095e-4;...
    .76, 0, 0, 0, 0;...
    .76, 0, 0, 0, 0;...
    0, .76, 0, 0, 0;...
    0, 0, 0, 0, 0];
r = [0, 7.24, 18.2, 25.4, 25.4]*100/4; % Reward matrix
offset = 0;
r = r+offset;
% pi = [0,0,0,0,1]; % initial distribution
% % Define the Q-matrix
% Q = [0, 0, 0, 0, 0;...
%     0, 0, 0, .76, 0;...
%     0, 0, 0, 0, .76;...
%     0, 0, 0, 0, .76;...
%     3.0095e-4, 2.5939e-4, 1.2298e-4, 1.2298e-4, 0];
% r = [25.4, 25.4, 18.2, 7.24, 0]*100/4; % Reward matrix
for i = 1:5
    Q(i,i) = -1*sum(Q(i,:));
end
% Defining the structure para
para.S = S;
para.pi = pi;
para.Q = Q;
para.r = r;
addpath('..\');
x_max = 25.4*100/4; % Limit x values
T_max = 40; % Time limit
x_max = x_max+offset*T_max;
Delta = 1/4*(1/2).^(0:7);

tic;
result_markov = accumulated_reward_numit_ext_markov(T_max,x_max,Delta(6),para)
toc

tic;
result = accumulated_reward_numit(T_max,x_max,Delta(6),para)
toc

% tic;
% result_direct_sp_sum = accumulated_reward_numit_dir_sp_sum(T_max,x_max,Delta(6),para)
% toc
% 
tic;
result_direct = accumulated_reward_numit_direct(T_max,x_max,Delta(6),para)
toc

% tic;
% [result_trun, cdf_trun_error, epsilon] = accumulated_reward_numit_trun(T_max,x_max,Delta(8),para)
% toc;
