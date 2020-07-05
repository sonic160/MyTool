%% This is to implement the vectorized simulation.
% Simulate the semi-markov resilience case study using interpolation
% sampling.
function [result_vector, et_vector] = ...
    method_vec_mc(n_s,Q,pi_0,r,evaluation_horizon,x_max)
    addpath(genpath('C:/Users/Zhiguo/OneDrive - CentraleSupelec/Code/MATLAB/MyTool'));

    %% Parameters related to simulation.
    % Transition probability matrix for the embedded chain.
    [a,P] = CalculatePEmbedded(Q);

    %% Simulation
    tic;
    % y and tau are matrix with dimension ns*(n_jump+1).
    % Row of y represents the state after each jump, from 0 (initial values) to n_jump.
    % Row of tau represents the occurrence time after each jump.
    [y_int, tau_int] = simulate_hmp_embedded(n_s, evaluation_horizon, ...
        pi_0, P, a);

    %% Post-processing for interpolation method
    % Duplicate the last state.
    y_int = [y_int,y_int(:,end)];
    tau_int = [tau_int,evaluation_horizon*ones(size(tau_int,1),1)];
    Delta_t = tau_int(:,2:end) - tau_int(:,1:end-1);
    r(end+1) = 0;
    P_O_time = sum(Delta_t.*r(y_int(:,1:end-1)),2);
    result_vector = sum(P_O_time<=x_max)/n_s;
    et_vector = toc;

end

% This is to generate the embedded chain from Q matrix
% Output: a - Rates of the holding time distributions
%         P - Transition matrix of the embedded chain
% History: 20170203, Created by ZZ
function [a,P] = CalculatePEmbedded(Q)
    a = -1*diag(Q); % Calculate the rates for the holding time distributions
    n_state = length(a); % Number of states
    % Construct the transition matrix for the embedded chain
    P = Q./repmat(a,1,n_state);
    for i = 1:n_state
        P(i,i) = 0;
    end
end

% This is to simulate a trajectory of a homogeneous Markov process
% using embedded chain, conditioned on the number of occurrence and
% vectorized programming.
% Input: ns - sample size.
%        evaluation_horizon - Time limit of the analysis.
%        pi_0_old - Initial distribution.
%        P_old - Transition probability matrix of the embedded chain.
%        a_old - Exit rates of each state. n_state*1.
%        jump_lim: Optional. jump_lim defines the
%          maximal number of jumps allowed. The simulation is continued until
%          all the sample paths reached evaluation horizon, or jump_lim is
%          reached. Default value: 1e3.
% Outputs: y - A ns*(n_jump+1) matrix. Each row is one sample path. Each column
%              is a sample.
%          tau - The same dimension as y, containing the time for the jump.
%          Note: y has one extra state, representing that the current jump
%                is beyond the evalution horizon. The corresponding tau is
%                evaluation horizon.

function [y, tau] = simulate_hmp_embedded(ns, evaluation_horizon,...
    pi_0_old, P_old, a_old, jump_lim)
        non_exceed_rate = 0; % This parameter controls how many sample paths that do not exceed evaluation_horizon are allowed.
        % Reset maximal number of jumps: to avoid dead lock.
        if nargin < 6 % If jump_lim is not given.
            jump_lim = 1e3; % Use its default value.
        end
        % Run simulation: All paths need to reach the end
        [y, tau] = all_reach_end(ns, jump_lim, evaluation_horizon,...
                    pi_0_old, P_old, a_old, non_exceed_rate);    
end

% This is the case where max_n_jump<0: Run the simulation, until all the
% sample paths reach evaluation_horizon,  or the number of jumps reaches
% jump_lim.
% Inputs and outputs have the same meaning as the main function.
function [y, tau] = all_reach_end(ns, jump_lim, evaluation_horizon,...
    pi_0_old, P_old, a_old, non_exceed_rate)
    max_n_jump = jump_lim; % Rest the max_n_jump. Now it is a upper limit.
    index_non_exceed = transpose(1:ns); % Set of indexes of the sample paths that do not exceed evaluation_horizon.
    
    %% Add a virtual state (the last one) to represent right censored data when t>evaluation horizon.
    n_state = length(pi_0_old);
    pi_0 = [pi_0_old, 0];

    % Update the transition probability matrix: Make sure that once it
    % enters the extra right censored state, it remains there.
    P = zeros(1+n_state, 1+n_state);
    P(1:n_state, 1:n_state) = P_old;
    P(1+n_state, 1+n_state) = 1;
    % When it enters the censored state, the holding time will be zero.
    a = zeros(1+n_state,1);
    a(1:n_state) = a_old;
    
    %% Inital states of the embedded chain.
    tau = zeros(ns, max_n_jump+1); % Cumulative times at transitions.
    y = zeros(ns, max_n_jump+1);

    % Generate initial states.
    pd_pi_0 = makedist('Multinomial','probabilities',pi_0);
    y_initial = random(pd_pi_0,ns,1);
    tau_initial = zeros(ns,1);

    i_jump = 1;
    y(:,i_jump) = y_initial;

    %% Simulate the upcoming jumps.
    y_cur = y_initial;
    tau_cur = tau_initial;

    for i_jump = 2:max_n_jump+1
        % Get unique current states.  
        unique_state = unique(y_cur);
        y_next = y_cur;
        
        % Do a loop for each unique state:
        for i = 1:length(unique_state)
            % Find the tau associated with each unique state.
            index = (y_cur==unique_state(i));
                        
            % Calcualte the transition probabilty matrix.
            p_tran = P(unique_state(i),:);
            
            % Generate next state.
            temp_y_next = generate_y_next(p_tran,sum(index));

            % Handle unexceptions: temp_y_next might be 0, which is caused
            % by the fact that the tau is beyond the range of the trained
            % surrogate model for the inverse CDF. Solution: replace these
            % points with the censoring state. The difference between these
            % states and the true cencoring states is that in this case,
            % the tau is not evaluation horizon.
            if min(temp_y_next) == 0                
%                 temp_y_next(temp_y_next==0) = n_state+1;
                error('Invalid value of y_next generated!');
            end
            
            % Simulate holding time.
            temp_theta = exprnd(1/a(unique_state(i)), sum(index), 1);            

            % Update the current state and time
            y_next(index) = temp_y_next;
            tau_cur(index) = ...
                tau_cur(index) + temp_theta;
        end
        % Update the current state and time
        y_cur = y_next;
        
        % Check the termination time.
        ind_out = (tau_cur>evaluation_horizon);
        tau_cur(ind_out) = evaluation_horizon;
        y_cur(ind_out) = n_state+1;
            
        % Save the current state and current time
        y(:,i_jump) = y_cur;
        tau(:,i_jump) = tau_cur;
        
        % Judge if all the paths exceed evaluation_horizon
        index_non_exceed = index_non_exceed((tau(index_non_exceed,i_jump)-tau(index_non_exceed,i_jump-1))~=0);        
        if length(index_non_exceed) <= ns*non_exceed_rate
            y = y(:,1:i_jump);
            tau = tau(:,1:i_jump);
            return;
        end            
    end
    % If jump_lim is reached: Ask to increase jump_lim
    error('Not all sample paths reached the evaluation horizon. Jump_lim too small!')
end

function temp_y_next = generate_y_next(p_tran,dim)
    temp_y_next = transpose(datasample(1:length(p_tran), dim, ...
        'Weights', p_tran));
end

