% This is to simulate a trajectory of a CTMC using Q matrix.
% Input: lambda - Rates of the holding time distributions: n*1
%        P - Transition matrix of the embedded chain: n*n
%        pi - Initial distribution, n*1
%        T_max - End time, simulation is performed in (0,T_max)
% Output: t - times of the jumps
%         y - the state after the jumps
% History:  20170205 - Modified by ZZ
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