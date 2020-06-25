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

function [cdf, epsilon] = accumulated_reward_numit_ext_markov(t,x_max,Delta,para,tol)
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
% The algorithm is implemented by constructing an extended discrete time
% discrete state Markov chain {r,X}, whose transition probability matrix is
% P_ext.
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
% Version history: 20/06/2020: - Implement state aggregation.
%                              - Change the definition of states in P_ext.
%                  19/06/2020: Created.

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
% Calculate Q*Delta: This is a matrix whose element is what to be
% multiplied in the numerical integration.
Q_times_Delta = Q*Delta;
% The main diagnose: The probablity of remaining at state i.
for i = 1:n_state
    Q_times_Delta(i,i) = 1+Q(i,i)*Delta;
end
% Get the index of non-zero transition rates.
index_set = cell(n_state,1);
for i = 1:n_state
    index_set{i} = find(Q_times_Delta(i,:));
end
% Calculate total number of non-zero elements in Q_times_Delta.
n_nz_Q = sum(sum(Q_times_Delta>0));

%% Construct the extended Markov chain.
% The transition matrix, P_ext is represended as a sparse matrix.
% Initialization.
idx_i = zeros((n_x+1)*n_nz_Q,1); % Row index of the non-zeros in P_ext.
idx_j = zeros((n_x+1)*n_nz_Q,1); % Column index of the non-zeros in P_ext.
p_i_j = zeros((n_x+1)*n_nz_Q,1); % Values of the transition probabilities.
% Fill the three vectors.
start_pos = 0; % Starting positions for the index, updated after each cycle.
% Consider all the possible previous states, but only the current states
% which are possible to move to (with non-zero probabilities).
for state_prev = 1:n_state
    for state_cur = index_set{state_prev}
        end_pos = start_pos+n_x+1-r(state_prev);
        % Only consider the range where transitions will not exceed n_x+1.
        idx_i(start_pos+1:end_pos) =...
            state_prev + n_state*(0:n_x-r(state_prev));
        % The current reward values are shifted by r(state_prev).
        idx_j(start_pos+1:end_pos) =...
            state_cur + n_state*(r(state_prev):n_x);
        p_i_j(start_pos+1:end_pos) =...
            Q_times_Delta(state_prev,state_cur);
        % Update the current starting position.
        start_pos = end_pos;
    end
end
% Create a sparse matrix Q_ext.
% Only keep the range with non-zero probabilities.
eff_range = idx_i>0;
idx_i = idx_i(eff_range);
idx_j = idx_j(eff_range);
p_i_j = p_i_j(eff_range);
P_ext = sparse(idx_i,idx_j,p_i_j,n_state*(n_x+1),n_state*(n_x+1));
% Initial probability.
pi_0_ext = zeros(1,n_state*(n_x+1));
pi_0_ext(pi_0>0) = pi_0(pi_0>0)/Delta;

%% Calcualte the density iteratively.
% Initialize variables.
n_time_point = length(t); % Get the number of time points that need to be evaluated.
cdf = zeros(1,n_time_point); % Initial values for cdf.
index_t = floor(t/Delta); % Transform t into the index considering Delta.
current_position_t = 1; % A flag variable indicating the current t to be evaluated.
additional_term = 0; % This term collects the density that is not possible to exceed n_x-1.
% additional_error = 0; % Error term in Kahan summation algorithm.
% Iterations begins.
% First branch: t_cur = 2:t_cur_th.
% Do not need to consider truncation.
% Threshold t_cur: Below it, all the states are possible.
t_cur_th = min(n_t,floor(n_t+1-(n_x-1)/max(r)));
for t_cur = 1:t_cur_th
    % Kolmogorov equation.
    pi_0_ext = pi_0_ext*P_ext;
    % If the cdf at this time step needs to be saved.
    if index_t(current_position_t) == t_cur
        % Calculate the cdf at this time point using trapedoid method.
        cdf(current_position_t) = ...
            (sum(pi_0_ext) - 1/2*sum(pi_0_ext(end-n_state+1:end)))*Delta;
        current_position_t = current_position_t + 1; % Go to next time point.
    end    
end

% Second branch: t_cur = t_cur_th+1:n_t.
% When next states are surely below n_x, do not further continue on this 
% branch. Accumulate the density of such branches in additional_term.
for t_cur = t_cur_th+1:n_t
    % Judge if it is possible to exceed n_x-1.
    % Get the index that cannot exceed.
    exceed_index = ...
        (((n_x-max(r)*(n_t-t_cur+2))+1:(n_x-max(r)*(n_t-t_cur+1)))-1)*n_state + ...
        (1:n_state)';
%     exceed_index = ((1:n_state)-1)*(1+n_x) + ...
%         ((n_x-max(r)*(n_t-t_cur+2))+1:(n_x-max(r)*(n_t-t_cur+1)))';
    % Store the density of these points. 
    % Sum using Kahan summation algorithm to avoid round-off errors.
%     [additional_term, additional_error] = ...
%         kahan_summation(additional_term,...
%         sum(pi_0_ext(exceed_index)), additional_error);
    additional_term = additional_term + sum(pi_0_ext(exceed_index(:)));
    % Remove the impossible density and continue.
    pi_0_ext(exceed_index(:)) = 0;    
    P_ext(exceed_index(:),exceed_index(:)) = 0;
       
    % Kolmogorov equation.
    pi_0_ext = pi_0_ext*P_ext;
%     pi_0_ext(exceed_index(end)+1:end) =...
%         pi_0_ext(exceed_index(end)+1:end)*...
%         P_ext((exceed_index(end)+1:end),(exceed_index(end)+1:end));

    % If the cdf at this time step needs to be saved.
    if index_t(current_position_t) == t_cur
        % Calculate the cdf at this time point using trapedoid method.
        cdf(current_position_t) = ...
            (sum(pi_0_ext) - 1/2*sum(pi_0_ext(end-n_state+1:end)) ...
            + additional_term)*Delta;
        current_position_t = current_position_t + 1; % Go to next time point.
    end
end
end

% % Kahan summation algorithm: Return sum_result = sum_result + input.
% % Avoid round-off errors.
% % Input: c - Accumulated errors.
% function [sum_result, c] = kahan_summation(sum_result,input,c)
%     % Correct the accumulated error.
%     y = input - c;
%     % sum_result is big, y small, so low-order digits of y are lost.
%     t = sum_result + y;
%     % (t - sum) cancels the high-order part of y; subtracting y recovers negative (low part of y)
%     c = (t - sum_result) - y;
%     sum_result = t;       
% end