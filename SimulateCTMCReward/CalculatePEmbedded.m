% This is to generate the embedded chain from Q matrix
% Output: a - Rates of the holding time distributions
%         P - Transition matrix of the embedded chain
% History: 20170203, Created by ZZ
function [a,P] = CalculatePEmbedded(Q)
a = -1*diag(Q); % Calculate the rates for the holding time distributions
n_state = length(a); % Number of states
% Construct the transition matrix for the embedded chain
P = Q./repmat(a,1,n_state);
for i = 1:n_state
    P(i,i) = 0;
end