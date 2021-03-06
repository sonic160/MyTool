% This is to simulate a trajectory of a nonhomogeneous semi-Markov process.
% Input: inv_F_t - A matrix that contains handles to a function that 
%        generate values of the inverse function of F_{i,j}(t,tau).
%        F_{i,j}(t,tau) is the PDF of the corresponding holding time
%        distributions, where t is the holding time, tau is the age.
%        P - Function handle matrix to transition probability matrix of the
%        embedded Markov chain. Also a function of t and tau.
%        pi_0 - Initial distribution.
%        T_max - Evaluation horizon.
% Output: t - times of the jumps.
%         y - the state after the jumps.
% History:  20191231: Created by ZZ

function [t,y] = simulate_semi_markov_NNinv(inv_F_t,P,pi_0,T_max)
    t(1) = 0; % start times at 0
    y(1) = rando(pi_0); % generate the initial y

    i = 1;
    while 1
        tau = t(i);

        % Generate next state
        n_state = length(P);
        p_next_state = zeros(1,n_state);
        for j = 1:n_state
            handle_p = P{y(i),j};
            p_next_state(j) = handle_p(tau);
        end        
        y_next = rando(p_next_state);  % Simulate next state using the transitin matrix of the embedded chain
        
        % Generate holding time
        handle_inv_F_t = inv_F_t{y(i),y_next};
        r = rand; % A random number from uniform dist
        theta = handle_inv_F_t(r,tau);
        
        t_cur = t(i) + theta; % Update time
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


