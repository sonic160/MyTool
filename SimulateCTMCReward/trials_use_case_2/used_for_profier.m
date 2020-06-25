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
for i = 1:5
    Q(i,i) = -1*sum(Q(i,:));
end
r = [0, 7.24, 18.2, 25.4, 25.4]*100/4; % Reward matrix
% Defining the structure para
para.S = S;
para.pi = pi;
para.Q = Q;
para.r = r;
addpath('..\');
x_max = 25.4*100/4; % Limit x values
T_max = 40; % Time limit
% Algorithm parameters
Delta = .05*1/2;
result = zeros(1,length(Delta));
elasped_time = zeros(1,length(Delta));
for i = 1:length(Delta)
    tic;
    result(i) = accumulated_reward_numit_ext_markov(T_max,x_max,Delta(i),para);
    elasped_time(i) = toc;
    % Output the results.
    fprintf('%d/%d\n',i,length(Delta));
    fprintf('Delta = %f\n',Delta(i));
    disp(result(i));
    disp(elasped_time(i));
end