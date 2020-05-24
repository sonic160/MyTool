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
%                   tol - Error tolerance
% Output parameter: cdf - The cdf F(t,y<x), at T_max.
% Version history: 23/05/2020: Created.

function cdf = accumulated_reward_numit(T_max,x_max,Delta,para,tol)
    % Initial values.       
    result = zeros(1,length(Delta)); % A vector of cdfs for each point in Delta.
    AX = zeros(1,length(Delta)-2); % Result of the Aitken's series.
    diff = 1; % Difference between two adjacent simulations.
    
    % The first point.
    i = 1;
    % Calculate the first three points in Delta. This is because to
    % calculate AX, it requires three consecutive points.
    result(i) = accumulated_reward_trapezoid(T_max,x_max,Delta(i),para);
    result(i+1) = accumulated_reward_trapezoid(T_max,x_max,Delta(i+1),para);
    result(i+2) = accumulated_reward_trapezoid(T_max,x_max,Delta(i+2),para);
    % Calculate the Aitken's series.
    AX(i) = aitken_extr(result(i),result(i+1),result(i+2));
    % Go to next point.
    i = i+1;
    
    % Continue search until Delta is empty or the desired accuracy is
    % reached.
    while (i <= length(Delta)-2) && (diff >= tol)
        % Calculate the original series.
        result(i+2) = accumulated_reward_trapezoid(T_max,x_max,Delta(i+2),para);
        % Calculate the Aitken's series.
        AX(i) = aitken_extr(result(i),result(i+1),result(i+2));      
        % Update the difference between two adjacent points.
        diff = abs(AX(i) - AX(i-1));
        % Save the result.
        cdf = AX(i);
        % Go to next point.
        i = i + 1;
    end
    
    % Consider different exit conditions.  
    if i > length(Delta)-2
        fprintf('Warning: Itegration terminated because all the points in Delta has been tested but the required accuracy not reached!\n');
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
% Input parameters: T_max - time limit for the evalutation
%                   x_max - x limit for the evaluation
%                   Delta - step size for the approximation
%                   para - parameter structure: 
%                          para.S - state space, a row vector that starts from 1
%                          para.pi - Initial probability distribution
%                          para.Q - Q-matrix for Markov chain
%                          para.r - Reward vector, a row vector that takes
%                                   non-negtive integer values (IMPORTANT,
%                                   must be integer).
% Output parameter: cdf - The cdf F(t,y<x), at T_max.
% Version history: 21/05/2020: Created.

function cdf = accumulated_reward_trapezoid(T_max,x_max,Delta,para)
%% Initialize parameters
S = para.S;
pi_0 = para.pi;
Q = para.Q;
r = para.r;

%% Algorithm begin: Initialization
n_x = floor(x_max/Delta); % Number of steps in x
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
        end
    end
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
                end                
            end
        end
    end
    density_matrix_prev = density_matrix_cur;
    density_matrix_cur = zeros(n_state,n_x+1);
end
cdf = sum(sum(density_matrix_prev(:,1:(n_x))*Delta,2) + ...
    density_matrix_prev(:,n_x+1)*Delta/2);

end
