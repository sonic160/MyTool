% This is to simulate a trajectory of a homogeneous semi-Markov process.
% Input: Q - A handle to a function that generate values of Q_{i,j}(t).
%        P - Transition probability matrix of the embedded Markov chain.
%        pi - Initial distribution.
%        T_max - Evaluation horizon.
% Output: t - times of the jumps.
%         y - the state after the jumps.
% History:  20191205 - Created by ZZ

function [t,y] = simulate_semi_markov(Q_t,P,pi,T_max)
t(1) = 0; % start times at 0
y(1) = rando(pi); % generate the initial y

i = 1;
while 1
    % Generate next state
    y_next = rando(P(y(i),:)); % Simulate next state using the transitin matrix of the embedded chain
    % Generate holding time
    t_next = generate_t_next(t,Q_t,P,y(i),y_next);
    t_cur = t(i) + t_next; % Update time
    if t_cur < T_max % If this transition happens in the evaluation horizon
        % Update the parameters
        t(i+1) = t_cur;
        y(i+1) = y_next;
        i = i+1;
    else
        break;
    end
end

%  rando.m generates a random variable in 1, 2, ..., n given a distribution 
%  vector. 
function index = rando(p)
u = rand;
i = 1;
s = p(1);
while ((u > s) && (i < length(p)))
    i=i+1;
    s=s+p(i);
end
index=i;

% This is to generate the next arrival time
function t_next = generate_t_next(t,Q_t,P,y_prev,y_next)
r = rand; % A random number from uniform dist
eq = @(t) cal_holding_time_cdf(t,Q_t,P,y_prev,y_next) - r; % Equation for inverse function of holding time cdf
t_next = fzero(eq,0); % Inverse function method to generate random number


% This is to calculate the holding time distribution
function cdf_ht_t = cal_holding_time_cdf(t,Q_t,P,y_prev,y_next)
F = Q_t(t)./P;
cdf_ht_t = F(y_prev,y_next);