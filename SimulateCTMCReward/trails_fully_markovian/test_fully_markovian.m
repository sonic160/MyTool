% This is to test the fully Markovian approximation using the use case 1 as
% the other files.

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

Delta = 1/4*(1/2).^(0:9);
result = zeros(1,length(Delta));
elasped_time = zeros(1,length(Delta));
% Run the test.
epsilon = 1e-3;
format longE;
for i = 1:length(Delta)
    % First method.
    tic;
    result(1,i) = 1 - accumulated_reward_fully_markovian(T_max,x_max,Delta(i),para,epsilon);
    elasped_time(1,i) = toc;

    % Output the results.
    fprintf('%d/%d\n',i,length(Delta));
    fprintf('Delta = %f\n',Delta(i));
    disp(result(:,i));
    disp(elasped_time(:,i));
end