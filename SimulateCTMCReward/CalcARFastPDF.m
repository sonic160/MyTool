% This is to implement the fast algorithm in Tijms and Veldman (2000)
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
%                   direction - Determine whether the cdf or icdf is
%                          calculated. 'forward': calculate cdf by
%                          integrate from 0 to x_max. 'inverse': calculate
%                          icdf by integrate from x_max to max(r)*T_max.
%                          Default value: 'forward'.
% Output parameter: cdf - The cdf F(t,y<x), at T_max.
% Version history: 20/05/2020: Optimize performance.
%                  15/01/2019: Created by ZZ.

function cdf = CalcARFastPDF(T_max,x_max,Delta,para,direction)
if nargin == 4 % Default value for direction is 'forward'
    direction = 'forward';
else
    if nargin ~= 5
        error('Wrong input parameters!')
    end
end

switch direction
    case 'forward' % Forward integration: from 0 to x
        cdf = CalcARFastPDF_forward(T_max,x_max,Delta,para);
    case 'inverse' % Integrate from x to y_max, to calculate icdf.
        icdf = CalcARFastPDF_inverse(T_max,x_max,Delta,para);
        cdf = 1 - icdf;
    otherwise
        error('Unrecgonized direction input!')
end

end

% Forward integration: from 0 to x to calculate cdf.
function cdf = CalcARFastPDF_forward(T_max,x_max,Delta,para)
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
density_matrix_prev = zeros(n_state,n_x); 
density_matrix_cur = zeros(n_state,n_x); 
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
            density_matrix_prev(state_cur,reward_cur) = density_matrix_prev(state_cur,reward_cur) +...
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
        all_reward_prev = find(density_matrix_prev(state_prev,:)); % Take the corresponding row.       
        for reward_prev = all_reward_prev
            reward_cur = reward_prev + r(state_prev); % Update current reward value. Delta is scaled out. No need to consider.
            % Judge if current reward value already beyond evaluation
            % threshold.
            if reward_cur >= n_x
                continue;
            else
                % Consider all the outbounding states and update the current
                % density.
                for state_cur = 1:n_state
                    p = density_matrix_prev(state_prev,reward_prev)*Q_times_Delta(state_prev,state_cur);
                    density_matrix_cur(state_cur,reward_cur) = ...
                        density_matrix_cur(state_cur,reward_cur) + p;
                end                
            end
        end
    end
    density_matrix_prev = density_matrix_cur;
    density_matrix_cur = zeros(n_state,n_x);
end
cdf = sum(sum(density_matrix_prev*Delta,2));
end

% Inverse integration: from x to y_max to calculate icdf.
% Difference compared to the previous function is that, here we calculate
% all the values for y, all the way to y_max. But in the previous function,
% once y is larger than x_max, no further calculation is needed (density
% treated as zeros).
function icdf = CalcARFastPDF_inverse(T_max,x_max,Delta,para)
%% Initialize parameters
S = para.S;
pi_0 = para.pi;
Q = para.Q;
r = para.r;
% Calculate y_max
y_max = max(r)*T_max;

%% Algorithm begin: Initialization
n_y = floor(y_max/Delta); % Number of steps in x
n_th = floor(x_max/Delta); % Lower limit of integration
n_t = floor(T_max/Delta); % Number of steps in t
n_state = length(S); % Number of states
% A matrix that stores the previous density f_i((k-1)\Delta,y).
% The columns of this matrix represents each discretized point on y.
density_matrix_prev = zeros(n_state,n_y); 
density_matrix_cur = zeros(n_state,n_y); 
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
    % Consider all the possible outbound states.
    for state_cur = 1:n_state
        density_matrix_prev(state_cur,reward_cur) = density_matrix_prev(state_cur,reward_cur) +...
            pi_0(state_prev)/Delta*Q_times_Delta(state_prev,state_cur);
    end
end    

%% Begin iterative processes.
for t_cur = 2:n_t
    % For each row in density_matrix_prev
    for state_prev = 1:n_state
        % Consider only the states with non-zeros density at the previous
        % states.
        all_reward_prev = find(density_matrix_prev(state_prev,:)); % Take the corresponding row.       
        for reward_prev = all_reward_prev
            reward_cur = reward_prev + r(state_prev); % Update current reward value. Delta is scaled out. No need to consider.
            % Consider all the outbounding states and update the current
            % density.
            for state_cur = 1:n_state
                p = density_matrix_prev(state_prev,reward_prev)*Q_times_Delta(state_prev,state_cur);
                density_matrix_cur(state_cur,reward_cur) = ...
                    density_matrix_cur(state_cur,reward_cur) + p;
            end                
        end
    end
    density_matrix_prev = density_matrix_cur;
    density_matrix_cur = zeros(n_state,n_y);
end
icdf = sum(sum(density_matrix_prev(:,n_th:end)*Delta,2));
end