% This is to test SimulateCTMCEmbedded.m
clear; clc;
S      = [1 2 3];                               % state space
T_max   = 80;                                    % maximum time
pi     = [0 0 1];                               % initial distribution
lambda = [0.1 0.08 0.5];                        % sojourn parameters
P      = [[0 .8 .2]; [.7 0 .3]; [.9 .1 0]];     % jump matrix of the embedded chain
[t,y] = SimulateCTMCEmbedded(lambda,P,pi,T_max);
% Plot the trajectory
stairs(t,y);                 % sometime this is more appropriate
axis([0 T_max 0 (length(pi)+1)]);
xlabel('Time');