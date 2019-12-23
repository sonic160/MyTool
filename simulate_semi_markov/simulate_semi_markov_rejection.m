% This is to simulate a trajectory of a homogeneous semi-Markov process.
% Rejection method is used for generating random samples
% Input: f_t - A handle to a function that generate values of f_{i,j}(t).
%        f_{i,j}(t) is the pdf of the corresponding holding time
%        distributions
%        max_value - Maximum of the pdf in its support. A n*n matrix.
%        support - A vector [lower,upper] that defines the support of the
%        pdf. A n*n cell matrix.
%        P - Transition probability matrix of the embedded Markov chain.
%        pi - Initial distribution.
%        T_max - Evaluation horizon.
% Output: t - times of the jumps.
%         y - the state after the jumps.
% History: 20191223 - Update the rejection sampling subfunction.
%                     Extend to be applicable using proposal density.
%          20191221 - Created by ZZ

function [t,y] = simulate_semi_markov_rejection(f_t,g_t,g_rng_t,c_t,P,pi,T_max)
    t(1) = 0; % start times at 0
    y(1) = rando(pi); % generate the initial y

    i = 1; % Index of the current state
    while 1
        % Generate next state
        y_next = rando(P(y(i),:)); % Simulate next state using the transitin matrix of the embedded chain
        % Generate holding time
        pdf_holding_time = f_t{y(i),y_next}; % handle to the holding time distribution pdf
        prop_den = g_t{y(i),y_next}; % Handle to the proposal density.
        prop_rng = g_rng_t{y(i),y_next}; % Handle for generating samples from the proposal density.
        c = c_t(y(i),y_next);
        t_next = accrejrnd(pdf_holding_time,prop_den,prop_rng,c,1,1); % Generate random number from pdf.
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
end
    
% This is to generate a random number given a pdf using rejection method.
% f - target pdf.
% g - proposal density.
% grnd - random number generator for g
% c - coefficient that make f(x) <= c g(x), for all x
% m, n, dimention of the samples
% X - matrix containing all the samples
function X = accrejrnd(f,g,grnd,c,m,n)
    X = zeros(m,n); % Preallocate memory
    for i = 1:m*n
        accept = false;
        while accept == false
            u = rand();
            v = grnd();
            if c*u <= f(v)/g(v)
               X(i) = v;
               accept = true;
            end
        end
    end
end