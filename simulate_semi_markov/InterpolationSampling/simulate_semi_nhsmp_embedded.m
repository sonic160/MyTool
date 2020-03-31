% This is to simulate a trajectory of a non-homogeneous semi-Markov process
% using embedded chain, conditioned on the number of occurrence.
% Input: ns - sample size.
%        max_n_jump - Maximum number of jumps simulated.
%        evaluation_horizon - Time limit of the analysis.
%        pi_0_old - Initial distribution.
%        cal_p_old - A cell matrix that contains handles for caluculating
%                    P_ij(tau). The handle takes only one input of tau.
%        rnd_matrix_theta_old - A cell matrix that contains handles for
%        generating random numbers from the holding time distribution f_ij.
%        The random number generator takes only one input, which is tau.
% Outputs: y - A ns*(n_jump+1) matrix. Each row is one sample path. Each column
%              is a sample.
%          tau - The same dimension as y, containing the time for the jump.
%          Note: y has one extra state, representing that the current jump
%                is beyond the evalution horizon. The corresponding tau is
%                evaluation horizon.

function [y, tau] = simulate_semi_nhsmp_embedded(ns, max_n_jump, evaluation_horizon,...
    pi_0_old, cal_p_old, rnd_matrix_theta_old)
    %% Add a virtual state (the last one) to represent right censored data when t>evaluation horizon.
    n_state = length(pi_0_old);
    pi_0 = [pi_0_old, 0];

    % Update the transition probability matrix: Make sure that once it
    % enters the extra right censored state, it remains there.
    cal_p = cell(1, n_state+1);
    for i = 1:n_state
        cal_p_state_i_old = cal_p_old{i};
        cal_p_state_i = @(tau) [cal_p_state_i_old(tau), zeros(size(tau))];
        cal_p{i} = cal_p_state_i;
    end
    cal_p{n_state+1} = @(tau) [zeros(length(tau),n_state), ones(size(tau))];

    % When it enters the censored state, the holding time will be zero.
    rnd_matrix_theta = cell(n_state+1, n_state+1);
    for i = 1:n_state
        for j = 1:n_state
            rnd_matrix_theta{i,j} = rnd_matrix_theta_old{i,j};
        end
    end
    handle_zero = @(tau) 0;
    for j = 1:n_state
        rnd_matrix_theta{j, n_state+1} = handle_zero;
        rnd_matrix_theta{n_state+1, j} = handle_zero;
    end
    rnd_matrix_theta{n_state+1, n_state+1} = handle_zero;

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
            temp_tau_cur = tau_cur(index);
            
            % Calcualte the transition probabilty matrix.
            cal_p_tr = cal_p{unique_state(i)};    
            p_tran = cal_p_tr(temp_tau_cur);
            
            % Generate next state.
            temp_y_next = generate_y_next(p_tran);

            % Handle unexceptions: temp_y_next might be 0, which is caused
            % by the fact that the tau is beyond the range of the trained
            % surrogate model for the inverse CDF. Solution: replace these
            % points with the censoring state. The difference between these
            % states and the true cencoring states is that in this case,
            % the tau is not evaluation horizon.
            if min(temp_y_next) == 0                
                temp_y_next(temp_y_next==0) = n_state+1;
            end
            
            % Simulate holding time.
            temp_theta = ...
                generate_theta(temp_y_next, temp_tau_cur, rnd_matrix_theta, unique_state(i)); 

            % Update the current state and time
            y_next(index) = temp_y_next;
            tau_cur(index) = ...
                temp_tau_cur + temp_theta;
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
    end
end

% This subfunction generates the holding time.
% Input: rnd_matrix_theta is a handle matrix, whose rnd_theta(i,j) gives a hanlde
%        to a random number generate of a holiding time from i to j.
%        Calling format: rnd_theta = rnd_matrix_theta{i,j}; 
%                        holing_time = rnd_theta(tau), returns a vector the same size as
%                        tau, which contains the holding time from i to j.
function temp_theta = generate_theta(temp_y_next, temp_tau_cur, rnd_matrix_theta, current_state)
    n_range = length(temp_y_next); % Get the number of samples in this block.    
    % Get unique states in the next states.
    us_temp_y_next = unique(temp_y_next);
    
    % Generate the random number.
    temp_theta = zeros(n_range,1);
    for j = 1:length(us_temp_y_next)             
        % Get the inverse cdf handle.
        handle_rnd_theta = rnd_matrix_theta{current_state, us_temp_y_next(j)};           
        % Generate samples.
        index = (temp_y_next==us_temp_y_next(j));
        temp_theta(index) = handle_rnd_theta(temp_tau_cur(index));
    end    
end


function temp_y_next = generate_y_next(p_tran)
    n_range = size(p_tran,1); % Get the number of samples in this block.
    
    % Generate next state
    cum_p = cumsum(p_tran,2);
    u = rand(n_range,1);
    diff = cum_p - u;
    temp_y_next = findfirst(diff>0,2,1);    
end

