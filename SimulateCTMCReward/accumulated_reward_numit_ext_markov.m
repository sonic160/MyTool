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
% Version history:  26/06/2020: Change the acceleration series from
%                               Aitken's Delta 2 to adjacent error
%                               estimations.
%                   03/06/2020: Add input switch.
%                   01/06/2020: Add truncatiopn by state aggregation.
%                   25/05/2020: Add cdf at different time points.
%                   23/05/2020: Created.

function [cdf, epsilon] = accumulated_reward_numit_ext_markov(t,x_max,Delta,para,tol)
    if nargin == 5 % If tol is provided.        
        % Initial values.       
        result = zeros(length(Delta),length(t)); % A vector of cdfs for each point in Delta.
        result_cor = zeros(length(Delta)-1,length(t)); % Result of the Aitken's series.
        diff = 1; % Difference between two adjacent simulations.

        % The first point.
        i = 1;
        % Calculate the first three points in Delta. This is because to
        % calculate AX, it requires three consecutive points.
        result(i,:) = accumulated_reward_trapezoid(t,x_max,Delta(i),para);
        result(i+1,:) = accumulated_reward_trapezoid(t,x_max,Delta(i+1),para);
        % Calculate the Aitken's series.
        result_cor(i,:) = 2*result(i+1,:)-result(i,:);
        % Go to next point.
        i = i+1;

        % Continue search until Delta is empty or the desired accuracy is
        % reached.
        while (i <= length(Delta)-1) && (abs(diff(end)) >= tol)
            % Calculate the original series.
            result(i+1,:) = accumulated_reward_trapezoid(t,x_max,Delta(i+1),para);
            % Calculate the Aitken's series.
            result_cor(i,:) = 2*result(i+1,:)-result(i,:);
            % Update the difference between two adjacent points.
            diff = result_cor(i,:) - result_cor(i-1,:);
            % Save the result.
            cdf = result_cor(i,:);
            epsilon = diff;
            % Go to next point.
            i = i + 1;
        end

        % Consider different exit conditions.  
        if i > length(Delta)-1
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
% Version history: 01/07/2020: Add a validity check of Q: row sum should be
%                              0.
%                  19/06/2020: Created.

function cdf = accumulated_reward_trapezoid(t,x_max,Delta,para)
%% Initialize parameters
S = para.S;
pi_0 = para.pi;
Q = para.Q;
r = para.r;

% Validity check: Row sum of Q should be 1.
if sum(~(abs(sum(Q,2)) <= 100*eps*ones(length(S),1))) ~= 0
    error('Row sum of Q should be 0!');
end

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
        % Only consider the range where transitions will not exceed n_x+1.
        updated_range = start_pos + (1:n_x+1-r(state_prev));
        idx_i(updated_range) = (state_prev-1)*(1+n_x) + (1:n_x+1-r(state_prev));
        % The current reward values are shifted by r(state_prev).
        idx_j(updated_range) = (state_cur-1)*(n_x+1) +...
            ((r(state_prev)+1):(n_x+1));
        p_i_j(updated_range) = Q_times_Delta(state_prev,state_cur);
        % Update the current starting position.
        start_pos = start_pos + n_x + 1;
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
pi_0_ext((find(pi_0)-1)*(n_x+1)+1) = pi_0(pi_0>0)/Delta;

%% Calcualte the density iteratively.
% Initialize variables.
n_time_point = length(t); % Get the number of time points that need to be evaluated.
cdf = zeros(1,n_time_point); % Initial values for cdf.
index_t = floor(t/Delta); % Transform t into the index considering Delta.
current_position_t = 1; % A flag variable indicating the current t to be evaluated.
if index_t(current_position_t) == 0
    cdf(current_position_t) = 1;
    current_position_t = current_position_t + 1;
end
% Iterations.
for t_cur = 1:n_t
    % Kolmogorov equation.
    pi_0_ext = pi_0_ext*P_ext;
    % If the cdf at this time step needs to be saved.
    if index_t(current_position_t) == t_cur
        % Calculate the cdf at this time point using trapedoid method.
        cdf(current_position_t) = ...
            (sum(pi_0_ext) - 1/2*sum(pi_0_ext((1:n_state)*(n_x+1))))*Delta;
        current_position_t = current_position_t + 1; % Go to next time point.
    end    
end

end