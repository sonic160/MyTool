% This is to implement the fast algorithm in Tijms and Veldman (2000)
clear; clc;
%% Parameter definition
S = [1:5]; % state space
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
% Algorithm parameters
Delta = 1/64; % Step size
x_max = 110; % Limit x values
T_max = 25; % Time limit
% Defining the structure para
para.S = S;
para.pi = pi;
para.Q = Q;
para.r = r;
%% Calculate P(O(t)>x)
f_appr = CalcARFastPDF(T_max,x_max,Delta,para);
ICDF_O_t = 1-sum(sum(f_appr(:,:,end)*Delta))