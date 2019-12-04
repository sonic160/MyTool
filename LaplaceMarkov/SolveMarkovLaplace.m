% This is to find the state proability vector p_state, given a Q-matrix, using
% C-K equation and Laplace transform.
% Input: Q - Q-matrix n*n
%        pi - Initial distribution, n*1
% Output: p_state - Vector of state probabilities, n*1, symbolic
% History:  20170203 - Created by ZZ
function p_state = SolveMarkovLaplace(Q,pi)
n_state = size(Q,1); % Number of states in the Markov chain.
% Make sure pi is a row vector
if iscolumn(pi) == 1 % If it is a column vector, transpose it
    pi = pi';
end
syms s t; % s-Lapace domain, t-time domain
% Solve p_state in frequency domain, using CK equation
p_state_s = pi/(s*eye(n_state)-Q); % This is the solutiono in frequency domain
p_state = ilaplace(p_state_s,s,t); % Inverse lapace transform
p_state = simplify(p_state); % Simplify the result