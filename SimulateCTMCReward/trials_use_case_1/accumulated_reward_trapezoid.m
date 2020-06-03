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
% Version history: 01/06/2020: Add truncatiopn by state aggregation.
%                  21/05/2020: Created.

function [cdf,n_eval] = accumulated_reward_trapezoid(t,x_max,Delta,para)
%% Initialize parameters
S = para.S;
pi_0 = para.pi;
Q = para.Q;
r = para.r;
n_eval = 0; % Counter for evaluation numbers.

%% Algorithm begin: Initialization
n_x = floor(x_max/Delta); % Number of steps in x
T_max = t(end); % Get the evaluation horizon.
n_t = floor(T_max/Delta); % Number of steps in t
n_state = length(S); % Number of states
% A matrix that stores the previous density f_i((k-1)\Delta,y).
% The columns of this matrix represents each discretized point on y.
density_matrix_prev = zeros(n_state,n_x+1); 
density_matrix_cur = zeros(n_state,n_x+1); 
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
            density_matrix_prev(state_cur,reward_cur+1) = density_matrix_prev(state_cur,reward_cur+1) +...
                pi_0(state_prev)/Delta*Q_times_Delta(state_prev,state_cur);
            n_eval = n_eval + 1; % Update the counter.
        end
    end
end    
% If the cdf at this time step needs to be saved.
if index_t(current_position_t) == 1
    cdf(current_position_t) = sum(sum(density_matrix_prev(:,1:(n_x))*Delta,2) + ...
        density_matrix_prev(:,n_x+1)*Delta/2); % Calculate the cdf at this time point using trapedoid method.
    current_position_t = current_position_t + 1; % Go to next time point.
end

%% Begin iterative processes.
additional_term = 0; % This term collects the density that is not possible to exceed n_x.
additional_error = 0; % Error term in Kahan summation algorithm.
for t_cur = 2:n_t
    % For each row in density_matrix_prev
    for state_prev = 1:n_state
        % Consider only the states with non-zeros density at the previous
        % states.
        all_reward_prev = find(density_matrix_prev(state_prev,:)) - 1; % Take the corresponding row.
        % Judge if it is possible to exceed n_x.
        t_remain = n_t - t_cur + 1; % Remaining time.
        % Get the index that cannot exceed.
        exceed_index = (all_reward_prev + max(r)*t_remain) <= n_x;
        % Store the density of these points. 
        % Sum using Kahan summation algorithm to avoid round-off errors.
        [additional_term, additional_error] = ...
            kahan_summation(additional_term,...
            sum(density_matrix_prev(state_prev, all_reward_prev(exceed_index)+1)),...
            additional_error);
%         additional_term = additional_term + sum(density_matrix_prev(state_prev, all_reward_prev(exceed_index)+1));
        % Remove the impossible density and continue.
        all_reward_prev = all_reward_prev(~exceed_index);
        for reward_prev = all_reward_prev
            reward_cur = reward_prev + r(state_prev); % Update current reward value. Delta is scaled out. No need to consider.
            % Judge if current reward value already beyond evaluation
            % threshold.
            if reward_cur > n_x
                continue;
            else
                % Consider all the outbounding states and update the current
                % density.
                for state_cur = 1:n_state
                    p = density_matrix_prev(state_prev,reward_prev+1)*Q_times_Delta(state_prev,state_cur);
                    density_matrix_cur(state_cur,reward_cur+1) = ...
                        density_matrix_cur(state_cur,reward_cur+1) + p;
                    n_eval = n_eval + 1; % Update the counter.
                end                
            end
        end
    end
    density_matrix_prev = density_matrix_cur;
    density_matrix_cur = zeros(n_state,n_x+1);
    
    % If the cdf at this time step needs to be saved.
    if index_t(current_position_t) == t_cur
        cdf(current_position_t) = sum(sum(density_matrix_prev(:,1:(n_x))*Delta,2) + ...
            density_matrix_prev(:,n_x+1)*Delta/2) + ...
            additional_term*Delta; % Calculate the cdf at this time point using trapedoid method.
        current_position_t = current_position_t + 1; % Go to next time point.
    end
end

end

% Kahan summation algorithm: Return sum_result = sum_result + input.
% Avoid round-off errors.
% Input: c - Accumulated errors.
function [sum_result, c] = kahan_summation(sum_result,input,c)
    % Correct the accumulated error.
    y = input - c;
    % sum_result is big, y small, so low-order digits of y are lost.
    t = sum_result + y;
    % (t - sum) cancels the high-order part of y; subtracting y recovers negative (low part of y)
    c = (t - sum_result) - y;
    sum_result = t;       
end