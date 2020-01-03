function [y, tau] = simulate_semi_nhsmp_embedded(ns, n_jump, evaluation_horizon, ...
    pi_0_old, cal_p_old, inv_F_t_old)
    %% Add a virtual state (the last one) to represent right censored data when t>evaluation horizon.
    n_state = length(pi_0_old);
    pi_0 = [pi_0_old, 0];

    cal_p = cell(1, n_state+1);
    for i = 1:n_state
        cal_p_state_i_old = cal_p_old{i};
        cal_p_state_i = @(tau) [cal_p_state_i_old(tau), zeros(size(tau))];
        cal_p{i} = cal_p_state_i;
    end
    cal_p{n_state+1} = @(tau) [zeros(length(tau),n_state), ones(size(tau))];

    inv_F_t = cell(n_state+1, n_state+1);
    for i = 1:n_state
        for j = 1:n_state
            inv_F_t{i,j} = inv_F_t_old{i,j};
        end
    end
    handle_zero = @(p,tau) 0;
    for j = 1:n_state
        inv_F_t{j, n_state+1} = handle_zero;
        inv_F_t{n_state+1, j} = handle_zero;
    end
    inv_F_t{n_state+1, n_state+1} = handle_zero;

    %% Inital states of the embedded chain.
    tau = zeros(ns, n_jump+1); % Cumulative times at transitions.
    y = zeros(ns, n_jump+1);

    % Generate initial states.
    pd_pi_0 = makedist('Multinomial','probabilities',pi_0);
    y_initial = random(pd_pi_0,ns,1);
    tau_initial = zeros(ns,1);

    i_jump = 1;
    y(:,i_jump) = y_initial;

    %% Simulate the upcoming jumps.
    y_cur = y_initial;
    tau_cur = tau_initial;

    for i_jump = 2:n_jump+1
        % Rearrange y_cur and tau_cur. Gather the same states together.
        [y_cur, index_y] = sort(y_cur);
        tau_cur = tau_cur(index_y);    
        [unique_state, ia, ~] = unique(y_cur);
        loc_labels = [ia; ns+1];

        % Do a loop for each state:
        for i = 1:length(unique_state)
            % Get the parts associated with each state.
            range = loc_labels(i):(loc_labels(i+1)-1);
            n_range = length(range);
            temp_tau_cur = tau_cur(range);

            % Calcualte the transition probabilty matrix.
            cal_p_tr = cal_p{unique_state(i)};    
            p_tran = cal_p_tr(temp_tau_cur);

            % Generate next state.
            [temp_y_next,temp_tau_cur] = ...
                generate_y_next(n_range, p_tran, temp_tau_cur);

            % Simulate holding time.
            temp_theta = ...
                generate_theta(temp_y_next, temp_tau_cur, n_range, inv_F_t, unique_state(i)); 

            % Update the current state and time
            y_cur(loc_labels(i):loc_labels(i)+length(temp_y_next)-1) = temp_y_next;
            tau_cur(loc_labels(i):loc_labels(i)+length(temp_y_next)-1) = ...
                temp_tau_cur + temp_theta;

            % Check the termination time.
            ind_out = find(tau_cur>evaluation_horizon);
            tau_cur(ind_out) = evaluation_horizon;
            y_cur(ind_out) = n_state+1;
        end

        % Save the current state and current time
        y(:,i_jump) = y_cur;
        tau(:,i_jump) = tau_cur;    
    end
end


function temp_theta = generate_theta(temp_y_next, temp_tau_cur, n_range, inv_F_t, current_state)
    [us_temp_y_next, ia_temp_temp_y_next, ~] = unique(temp_y_next);
    ia_temp_y_next = [ia_temp_temp_y_next; n_range+1];
    
    temp_theta = zeros(n_range,1);
    for j = 1:length(us_temp_y_next)
        range_local = ia_temp_y_next(j):(ia_temp_y_next(j+1)-1);
        n_range_local = length(range_local);
        
        % Get the inverse cdf handle.
        handle_inv_F_t = inv_F_t{current_state, us_temp_y_next(j)};

        % Generate samples.
        u = rand(n_range_local,1);
        temp_theta(range_local) = handle_inv_F_t(u,temp_tau_cur(range_local));
    end
end


function [temp_y_next,temp_tau_cur] = generate_y_next(n_range, p_tran, temp_tau_cur)
%     % Generate next state
%     temp_y_next = zeros(n_range,1);
%     for j = 1:n_range
%         pd_p_tr = makedist('Multinomial','probabilities',p_tran(j,:));
%         temp_y_next(j) = random(pd_p_tr, 1, 1);
%     end
    
    % Generate next state
    cum_p = cumsum(p_tran,2);
    u = rand(n_range,1);
    diff = cum_p - u;
    temp_y_next = findfirst(diff>0,2,1);
    
    [temp_y_next, ind_sort_temp_y_next] = sort(temp_y_next);
    temp_tau_cur = temp_tau_cur(ind_sort_temp_y_next);
end

