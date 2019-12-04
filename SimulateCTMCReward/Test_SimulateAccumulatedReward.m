% This is to 
clear; clc;
S = [1:5]; % state space
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
time = 25;
x = 110;
P_O_time = zeros(length(time),NS);
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
x = 115;
P_O_time_x = length(find(P_O_time>x))/NS
