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
%                   trun_prob - When provided, use this value to trucation
%                               the iteration.
%                               If not provided, do not use truncation.
% Output parameter: cdf - The cdf F(t,y<x), at t.
% Version history: 01/06/2020: Add truncatiopn by probability.
%                  21/05/2020: Created.

function [cdf,n_eval] = accumulated_reward_trapezoid_trial(t,x_max,Delta,para,trun_prob)
    % Judge if truncation is activated.
    if nargin <= 4 % No truncation
        [cdf,n_eval] = accumulated_reward_trapezoid_no_trun(t,x_max,Delta,para);
    else % Truncation
        [cdf,n_eval] = accumulated_reward_trapezoid_trun(t,x_max,Delta,para,trun_prob);
    end       
end

% This is the subfunction for not using truncation.
function [cdf,n_eval] = accumulated_reward_trapezoid_no_trun(t,x_max,Delta,para)
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
for t_cur = 2:n_t
    % For each row in density_matrix_prev
    for state_prev = 1:n_state
        % Consider only the states with non-zeros density at the previous
        % states.
        all_reward_prev = find(density_matrix_prev(state_prev,:)) - 1; % Take the corresponding row.
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
            density_matrix_prev(:,n_x+1)*Delta/2); % Calculate the cdf at this time point using trapedoid method.
        current_position_t = current_position_t + 1; % Go to next time point.
    end
end

end

% This is the subfunction for using truncation.
function [cdf,n_eval] = accumulated_reward_trapezoid_trun(t,x_max,Delta,para,trun_prob)
%% Initialize parameters
S = para.S;
pi_0 = para.pi;
Q = para.Q;
r = para.r;
n_eval = 0; % Counter for evaluation numbers.
% trun_den = trun_prob/Delta; % Change from truncation probabiltiy to truncation density.

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
for t_cur = 2:n_t
    % For each row in density_matrix_prev
    for state_prev = 1:n_state
        % Consider only the states with non-zeros density at the previous
        % states.
        all_reward_prev = find(density_matrix_prev(state_prev,:)) - 1; % Take the corresponding row.       
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
                    % Judge if the current transition can be ignored.
                    if p <= trun_prob
                        continue;
                    else
                        density_matrix_cur(state_cur,reward_cur+1) = ...
                            density_matrix_cur(state_cur,reward_cur+1) + p;
                    end
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
            density_matrix_prev(:,n_x+1)*Delta/2); % Calculate the cdf at this time point using trapedoid method.
        current_position_t = current_position_t + 1; % Go to next time point.
    end
end

end