% This is to extend the fast algorithm in Tijms and Veldman (2000) by
% changing to trapezoid integration method.
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
% Output parameter: cdf - The cdf F(t,y<x), at t.
% Version history: 03/07/2020: Change name from .._numit to .._numit_tijm.
%                               Implement strictly Tijms' method: without
%                               considering elenimating non-zero elements.
%                  19/06/2020: Fix bug on the truncation by aggregation.
%                  10/06/2020: Improve the efficiency of resetting all_reward_cur.
%                  08/06/2020: - Avoid create zeros(n_x+1,n_state) in each cycle.
%                              - Change index_set from cell array to logical array.
%                              - Delete a redundant calculation of all_reward_prev.
%                  07/06/2020: Branch based on t_cur: Below it, all the states
%                              are possible. If t_cur_th < 2, then, automatically, 
%                              the first branch will not be executed.
%                  04/06/2020: Replace the inner-most loop by vector.
%                              Replace row indexing by column indexing.
%                  03/06/2020: Remove n_eval. 
%                              Only consider the next states with non-zero 
%                              transition probabilities.
%                              Optimize performance by switching order of
%                              iteration.
%                  01/06/2020: Add truncatiopn by state aggregation.
%                  21/05/2020: Created.

function cdf = accumulated_reward_numit_tijm(t,x_max,Delta,para)
%% Initialize parameters
S = para.S;
pi_0 = para.pi;
Q = para.Q;
r = para.r;

%% Algorithm begin: Initialization
n_x = floor(x_max/Delta); % Number of steps in x
T_max = t(end); % Get the evaluation horizon.
n_t = floor(T_max/Delta); % Number of steps in t
n_state = length(S); % Number of states
% A matrix that stores the previous density f_i((k-1)\Delta,y).
% The rows of this matrix represents each discretized point on y.
density_matrix_prev = zeros(n_x+1,n_state); 
density_matrix_cur_0 = zeros(n_x+1,n_state); 
density_matrix_cur = density_matrix_cur_0;
% Calculate Q*Delta: This is a matrix whose element is what to be
% multiplied in the numerical integration.
Q_times_Delta = Q*Delta;
% The main diagnose: The probablity of remaining at state i.
for i = 1:n_state
    Q_times_Delta(i,i) = 1+Q(i,i)*Delta;
end

n_time_point = length(t); % Get the number of time points that need to be evaluated.
cdf = zeros(1,n_time_point); % Initial values for cdf.
index_t = floor(t/Delta); % Transform t into the index considering Delta.
current_position_t = 1; % A flag variable indicating the current t to be evaluated.
% If the first point is zero, set cdf(1) to zero directly.
if index_t(current_position_t) == 0
    cdf(current_position_t) = (0<x_max);
    current_position_t = current_position_t + 1; % Go to next time point.
end

%% Get the index of non-zero transition rates.
% index_set = Q_times_Delta>0;

%% Step zero: The intial states and initial distributions.
for state_prev = 1:n_state
    reward_cur = r(state_prev); % Update current reward value.
    % Judge if current reward value already beyond evaluation
    % threshold.
    if reward_cur >= n_x
        continue;
    else
        % Consider all the possible outbound states.
        for state_cur = 1:n_state
            density_matrix_prev(reward_cur+1,state_cur) = density_matrix_prev(reward_cur+1,state_cur) +...
                pi_0(state_prev)/Delta*Q_times_Delta(state_prev,state_cur);
        end
    end
end    
% If the cdf at this time step needs to be saved.
if index_t(current_position_t) == 1
    cdf(current_position_t) = sum(sum(density_matrix_prev(1:(n_x),:)*Delta) + ...
        density_matrix_prev(n_x+1,:)*Delta/2); % Calculate the cdf at this time point using trapedoid method.
    current_position_t = current_position_t + 1; % Go to next time point.
end

%% Begin iterative processes.
for t_cur = 2:n_t
    % For each row in density_matrix_prev    
    for state_prev = 1:n_state                              
        % Consider only the outbounding states with non-zero density
        % and update the current density.        
%         density_matrix_cur(r(state_prev)+1:end,index_set(state_prev,:)) = ...
%             density_matrix_cur(r(state_prev)+1:end,index_set(state_prev,:)) + ...
%             density_matrix_prev(1:(n_x-r(state_prev)+1),state_prev).*...
%             Q_times_Delta(state_prev,index_set(state_prev,:));
        density_matrix_cur(r(state_prev)+1:end,:) = ...
            density_matrix_cur(r(state_prev)+1:end,:) + ...
            density_matrix_prev(1:(n_x-r(state_prev)+1),state_prev).*...
            Q_times_Delta(state_prev,:);
    end
    density_matrix_prev = density_matrix_cur;
    density_matrix_cur = density_matrix_cur_0;
    
    % If the cdf at this time step needs to be saved.
    if index_t(current_position_t) == t_cur
        cdf(current_position_t) = sum(sum(density_matrix_prev(1:(n_x),:)*Delta) + ...
            density_matrix_prev(n_x+1,:)*Delta/2); % Calculate the cdf at this time point using trapedoid method.
        current_position_t = current_position_t + 1; % Go to next time point.
    end
end

end
