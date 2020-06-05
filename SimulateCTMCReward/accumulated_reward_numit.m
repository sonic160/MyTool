% This function calculate the cdf of the accumulated reward of a Markov
% reward process at t: P(AR(T)max) <= x_max). The algorithm used is based
% on the one in Tijms and Veldman (2000), but an improved version:
% - Trapedoid, rather than rectangle integration is used.
% - Aitken's delta-2 process is used to accelerate the convergence.
% The integration is continued until the difference between two adjacent
% simulations is less than tol.

% Input parameters: T_max - time limit for the evalutation
%                   x_max - x limit for the evaluation
%                   Delta - step size for the approximation
%                         - A row vector that defines the search space. The
%                         algorithm ends when tol is reached, or this
%                         search space is empty.
%                   para - parameter structure: 
%                          para.S - state space, a row vector that starts from 1
%                          para.pi - Initial probability distribution
%                          para.Q - Q-matrix for Markov chain
%                          para.r - Reward vector, a row vector that takes
%                                   non-negtive integer values (IMPORTANT,
%                                   must be integer).
%                   tol - Error tolerance. If not provided, do not
%                         accelerate, but calculate all the Delta.
%                   t   - A vector of time points for cdf evaluation.
% Output parameter: cdf - The cdf F(t,y<x), at T_max.
%                   epsilon - The maximal error at the last time point.
% Version history:  03/06/2020: Add input switch.
%                   01/06/2020: Add truncatiopn by state aggregation.
%                   25/05/2020: Add cdf at different time points.
%                   23/05/2020: Created.

function [cdf, epsilon] = accumulated_reward_numit(t,x_max,Delta,para,tol)
    if nargin == 5 % If tol is provided.        
        % Initial values.       
        result = zeros(length(Delta),length(t)); % A vector of cdfs for each point in Delta.
        AX = zeros(length(Delta)-2,length(t)); % Result of the Aitken's series.
        diff = 1; % Difference between two adjacent simulations.

        % The first point.
        i = 1;
        % Calculate the first three points in Delta. This is because to
        % calculate AX, it requires three consecutive points.
        result(i,:) = accumulated_reward_trapezoid(t,x_max,Delta(i),para);
        result(i+1,:) = accumulated_reward_trapezoid(t,x_max,Delta(i+1),para);
        result(i+2,:) = accumulated_reward_trapezoid(t,x_max,Delta(i+2),para);
        % Calculate the Aitken's series.
        AX(i,:) = aitken_extr(result(i,:),result(i+1,:),result(i+2,:));
        % Go to next point.
        i = i+1;

        % Continue search until Delta is empty or the desired accuracy is
        % reached.
        while (i <= length(Delta)-2) && (abs(diff(end)) >= tol)
            % Calculate the original series.
            result(i+2,:) = accumulated_reward_trapezoid(t,x_max,Delta(i+2),para);
            % Calculate the Aitken's series.
            AX(i,:) = aitken_extr(result(i,:),result(i+1,:),result(i+2,:));      
            % Update the difference between two adjacent points.
            diff = AX(i,:) - AX(i-1,:);
            % Save the result.
            cdf = AX(i,:);
            epsilon = diff;
            % Go to next point.
            i = i + 1;
        end

        % Consider different exit conditions.  
        if i > length(Delta)-2
            fprintf('Warning: Itegration terminated because all the points in Delta has been tested but the required accuracy not reached!\n');
        end
    else
        if nargin == 4 % If tol not provided.
            cdf = zeros(length(Delta),length(t)); % A vector of cdfs for each point in Delta.
            % Do a loop to calculate each cdf.
            for i = 1:length(Delta)
                cdf(i,:) = accumulated_reward_trapezoid(t,x_max,Delta(i),para);
            end
            epsilon = 'Not applicable';            
        else
            error('Wrong inputs!');
        end
    end
end

% Calculate the accelerated series AX.
function AX = aitken_extr(x_i,x_ip1,x_ip2)   
    Delta_x_n = x_ip1 - x_i;
    Delta_x_np1 = x_ip2 - x_ip1;
    Delta_square_x_n = Delta_x_np1 - Delta_x_n;
    AX = x_i - Delta_x_n.^2./Delta_square_x_n;
end

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
% Version history: 04/06/2020: Replace the inner-most loop by vector.
%                              Replace row indexing by column indexing.
%                  03/06/2020: Remove n_eval. 
%                              Only consider the next states with non-zero 
%                              transition probabilities.
%                              Optimize performance by switching order of
%                              iteration.
%                  01/06/2020: Add truncatiopn by state aggregation.
%                  21/05/2020: Created.

function cdf = accumulated_reward_trapezoid(t,x_max,Delta,para)
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
density_matrix_cur = zeros(n_x+1,n_state); 
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
index_set = cell(n_state,1);
for i = 1:n_state
    index_set{i} = find(Q_times_Delta(i,:));
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
additional_term = 0; % This term collects the density that is not possible to exceed n_x.
additional_error = 0; % Error term in Kahan summation algorithm.
for t_cur = 2:n_t
    % For each row in density_matrix_prev
    for state_prev = 1:n_state
        % Consider only the states with non-zeros density at the previous
        % states.
        all_reward_prev = find(density_matrix_prev(:,state_prev)) - 1; 
        % Judge if it is possible to exceed n_x.
        t_remain = n_t - t_cur + 1; % Remaining time.
        % Get the index that cannot exceed.
        exceed_index = (all_reward_prev + max(r)*t_remain) <= n_x;
        % Store the density of these points. 
        % Sum using Kahan summation algorithm to avoid round-off errors.
        [additional_term, additional_error] = ...
            kahan_summation(additional_term,...
            sum(density_matrix_prev(all_reward_prev(exceed_index)+1,state_prev)),...
            additional_error);
%         additional_term = additional_term + sum(density_matrix_prev(state_prev, all_reward_prev(exceed_index)+1));
        % Remove the impossible density and continue.
        all_reward_prev = all_reward_prev(~exceed_index);        
        % Update current reward value. Delta is scaled out. No need to consider.
        all_reward_cur = all_reward_prev + r(state_prev); 
        % Judge if current reward value already beyond evaluation
        % threshold.
        flag = (all_reward_cur <= n_x);
        all_reward_cur = all_reward_cur(flag);
        all_reward_prev = all_reward_prev(flag);
        % Consider only the outbounding states with non-zero density
        % and update the current density.
        index = index_set{state_prev};
        density_matrix_cur(all_reward_cur+1,index) = ...
            density_matrix_cur(all_reward_cur+1,index) + ...
            density_matrix_prev(all_reward_prev+1,state_prev).*...
            Q_times_Delta(state_prev,index);
    end
    density_matrix_prev = density_matrix_cur;
    density_matrix_cur = zeros(n_x+1,n_state);
    
    % If the cdf at this time step needs to be saved.
    if index_t(current_position_t) == t_cur
        cdf(current_position_t) = sum(sum(density_matrix_prev(1:(n_x),:)*Delta) + ...
            density_matrix_prev(n_x+1,:)*Delta/2) + ...
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
