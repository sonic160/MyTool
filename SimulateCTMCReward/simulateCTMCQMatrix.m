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