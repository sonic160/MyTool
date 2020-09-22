% This is to implement the fully Markovian approximation method in Cloth 
% (2006). Cloth, L. and B. R. Haverkort (2006). Five performability algorithms:
% A comparison. Markov Anniversary Meeting, Charleston, SC, USA.

% Input parameters: x_max - x limit for the evaluation
%                   Delta - step size for the approximation
%                   para - parameter structure: 
%                          para.S - state space, a row vector that starts from 1
%                          para.pi - Initial probability distribution
%                          para.Q - Q-matrix for Markov chain
%                          para.r - Reward vector, a row vector that takes
%                                   non-negtive integer values (IMPORTANT,
%                                   must be integer).
%                   t - A vector of different time points to evaluate cdf.
%                   epsilon - error limit for the uniformization.
% Output parameter: cdf - The cdf F(t,y<x), at t.
% Version history: 18/09/2020: Created.

function cdf = accumulated_reward_fully_markovian(t,x_max,Delta,para,epsilon)
    %% Initialize parameters
    S = para.S;
    pi_0 = para.pi;
    Q = para.Q;
    r = para.r;
    % Initial values for cdf.
    cdf = zeros(1,length(t));
    
    % Validity check: Row sum of Q should be 1.
    if sum(~(abs(sum(Q,2)) <= 100*eps*ones(length(S),1))) ~= 0
        error('Row sum of Q should be 0!');
    end

    %% Algorithm begin: Initialization
    n_x = floor(x_max/Delta); % Number of steps in x.    
    n_state = length(S); % Number of states.
    
    % Derive the transition rate matrix of the CTMC.
    D = diag(r); % Create a diagnose matrix where the main diagnose elements are the reward rates.
    % Initialize the transition rate matrix.
    % It is a matrix of the size n_state*n_x*n_state*n_x.
    A = zeros(n_state*n_x,n_state*n_x); 
    % Fill the blocks of A.
    for i = 0:n_x-2
        A(i*n_state+(1:n_state), i*n_state+(1:n_state)) = Q;
        A(i*n_state+(1:n_state), (i+1)*n_state+(1:n_state)) = D/Delta;
        % Recalculate the main diagnose elements.
        for j = 1:n_state
            A(i*n_state+j, i*n_state+j) = -sum(A(i*n_state+j,:)) + A(i*n_state+j, i*n_state+j);
        end                
    end
    % The last block.
    i = n_x-1;
    A(i*n_state+(1:n_state), i*n_state+(1:n_state)) = Q;
    % Recalculate the main diagnose elements.
    for j = 1:n_state
        A(i*n_state+j, i*n_state+j) = Q(j,j) - sum(D(j,:))/Delta;
    end
    
    %% Uniformization.
    % Get $\lambda$: the maximized transition rate.
    lambda = max(-1*diag(A));
    % Calculate the matrix U.
    U = eye(size(A)) + 1/lambda*A;
    % Create a sparse version of U.
    U = sparse(U);
    % Get initial state matrix of the uniform chain.
    pi_0_U = zeros(1,n_state*n_x);
    pi_0_U(pi_0>0) = pi_0(pi_0>0);
    
    % Consider all the values in t.
    for i = 1:length(t)
        % Calculate the state probability vector.
        tt = t(i);        
        n_t = 0; % Needed iterations.
        e_c = 1 - cal_PP(lambda,tt,n_t); % Current error.
        
        pi_U_t = 0; % State probablity at the current t step.
        pi_temp = pi_0_U; % Second part in the summation equation.
                
        while e_c >= epsilon
            n_t = n_t + 1;
            PP_n_t = cal_PP(lambda,tt,n_t);
            e_c = e_c - PP_n_t; % Current error.
            pi_temp = pi_temp*U;
            pi_U_t = pi_U_t + pi_temp*PP_n_t;
        end
        % Output the desired cdf.
        cdf(i) = sum(pi_U_t);
    end    
end

% This is the subfunction to calculate PP(t): the probability of having $n$
% arrivals according to the Poisson process.
function PP_n = cal_PP(lambda,t,n)
    PP_n = poisspdf(n,lambda*t);
end