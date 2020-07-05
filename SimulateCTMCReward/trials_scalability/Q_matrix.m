% This subfunction derive the Q_matrix of the case study, given input m,
% which is the dimension of the problem.
% Output: Q_sp - Q-matrix in sparse matrix form;
%         Q_full - Q-matrix in full matrix form.

function [Q_sp, Q_full, r, pi_0, nnz] = Q_matrix(m, para, x_0)
    lambda = para.lambda; % Failure rate for the processor.
    gamma = para.gamma; % Failure rate for the memory.
    delta = para.delta; % Failure rate for the switch.
    nu = para.nu; % Local repair rate for the processor.
    eta = para.eta; % Local repair rate for the memory.
    epsilon = para.epsilon; % Local repair rate for the switch.
    mu = para.mu; % Global repair rate.
    c = para.c; % Coverage probability.

    % Initialize the parameters.
    end_state = (m-2)^2+(m-3)^2;
    nnz = 4*(m-3)^2 + 2*(m-4)*(m-3) + ...
        2*(1 + 2*(m-4)) + ...
        1 + 2*(m-4) + (m-4)*(2+3*(m-4)) + ...
        1 + ...
        end_state; % Last row: Diagnose elements.
    state_cur = zeros(nnz,1);
    state_next = zeros(nnz,1);
    transition_rate = zeros(nnz,1);
    r = zeros(1,end_state);

    % Do a loop to fill the elements in the Q-matrix.
    current = 1;
    previous = 1;
    % First consider k == 1.
    k = 1;
    % From any working state, 4 possible failure transitions + 2 possible repairs.
    for i = m:-1:4
        for j = m:-1:4
            % Failures.
            state_cur(current:current+3) = to_state_idx(i,j,k,m);
            state_next(current:current+3) = [to_state_idx(i-1,j,k,m);...
                to_state_idx(i,j-1,k,m); to_state_idx(i,j,k-1,m); end_state];
            transition_rate(current:current+3) = [i*c*lambda; j*c*gamma; delta; ...
                (1-c)*(i*lambda+j*gamma)];
            current = current + 4;
            % If recoverable in i.
            if i < m
                state_cur(current) = to_state_idx(i,j,k,m);
                state_next(current) = to_state_idx(i+1,j,k,m);
                transition_rate(current) = nu;
                current = current + 1;
            end
            % If recoverable in j.
            if j < m
                state_cur(current) = to_state_idx(i,j,k,m);
                state_next(current) = to_state_idx(i,j+1,k,m);
                transition_rate(current) = eta;
                current = current + 1;
            end

            % Diagnose elements.
            state_cur(current) = to_state_idx(i,j,k,m);
            state_next(current) = to_state_idx(i,j,k,m);
            transition_rate(current) = -1*sum(transition_rate(previous:current-1));
            current = current + 1;
            previous = current;

            % Fill the reward rates.
            r(to_state_idx(i,j,k,m)) = max([i,j])*(1-(1-1/max([i,j]))^min([i,j]));        
        end
    end
    % Failure states. Only recovery is possible.
    i = 3;
    for j = m:-1:4
        if j == m % j is not repairable
            state_cur(current) = to_state_idx(i,j,k,m);
            state_next(current) = to_state_idx(i+1,j,k,m);
            transition_rate(current) = nu;
            current = current + 1;
        else
            state_cur(current:current+1) = to_state_idx(i,j,k,m);
            state_next(current:current+1) = [to_state_idx(i+1,j,k,m);...
                to_state_idx(i,j+1,k,m)];
            transition_rate(current:current+1) = [nu; eta];
            current = current + 2;
        end

        % Diagnose elements.
        state_cur(current) = to_state_idx(i,j,k,m);
        state_next(current) = to_state_idx(i,j,k,m);
        transition_rate(current) = -1*sum(transition_rate(previous:current-1));
        current = current + 1;
        previous = current;
    end
    j = 3;
    for i = m:-1:4
        if i == m % i is not repairable
            state_cur(current) = to_state_idx(i,j,k,m);
            state_next(current) = to_state_idx(i,j+1,k,m);
            transition_rate(current) = eta;
            current = current + 1;
        else
            state_cur(current:current+1) = to_state_idx(i,j,k,m);
            state_next(current:current+1) = [to_state_idx(i+1,j,k,m);...
                to_state_idx(i,j+1,k,m)];
            transition_rate(current:current+1) = [nu; eta];
            current = current + 2;
        end

        % Diagnose elements.
        state_cur(current) = to_state_idx(i,j,k,m);
        state_next(current) = to_state_idx(i,j,k,m);
        transition_rate(current) = -1*sum(transition_rate(previous:current-1));
        current = current + 1;
        previous = current;
    end
    % Then, k == 0.
    k = 0;
    for i = m:-1:4
        if i == m % i is not repairable
            for j = m:-1:4
                if j == m % only k is repairable.
                    state_cur(current) = to_state_idx(i,j,k,m);
                    state_next(current) = to_state_idx(i,j,k+1,m);
                    transition_rate(current) = epsilon;
                    current = current + 1;
                else % j and k are repairable.
                    state_cur(current:current+1) = to_state_idx(i,j,k,m);
                    state_next(current:current+1) = [to_state_idx(i,j+1,k,m);...
                        to_state_idx(i,j,k+1,m)];
                    transition_rate(current:current+1) = [eta; epsilon];
                    current = current + 2;
                end

                % Diagnose elements.
                state_cur(current) = to_state_idx(i,j,k,m);
                state_next(current) = to_state_idx(i,j,k,m);
                transition_rate(current) = -1*sum(transition_rate(previous:current-1));
                current = current + 1;
                previous = current;
            end
        else % i is repairable.
            for j = m:-1:4
                if j == m % only i and k is repairable.
                    state_cur(current:current+1) = to_state_idx(i,j,k,m);
                    state_next(current:current+1) = [to_state_idx(i+1,j,k,m);...
                        to_state_idx(i,j,k+1,m)];
                    transition_rate(current:current+1) = [nu; epsilon];
                    current = current + 2;
                else % i, j and k are repairable.
                    state_cur(current:current+2) = to_state_idx(i,j,k,m);
                    state_next(current:current+2) = [to_state_idx(i+1,j,k,m);...
                        to_state_idx(i,j+1,k,m); to_state_idx(i,j,k+1,m)];
                    transition_rate(current:current+2) = [nu; eta; epsilon];
                    current = current + 3;
                end

                % Diagnose elements.
                state_cur(current) = to_state_idx(i,j,k,m);
                state_next(current) = to_state_idx(i,j,k,m);
                transition_rate(current) = -1*sum(transition_rate(previous:current-1));
                current = current + 1;
                previous = current;
            end
        end    
    end
    % Transition from the global failure state.
    state_cur(current) = end_state;
    state_next(current) = to_state_idx(m,m,1,m);
    transition_rate(current) = mu;
    % Diagnose elements.
    state_cur(current+1) = end_state;
    state_next(current+1) = end_state;
    transition_rate(current+1) = -mu;

    % Construct the Q-matrix
    n_state = (m-2)^2 + (m-3)^2;
    Q_sp = sparse(state_cur, state_next, transition_rate, n_state, n_state);
    Q_full = full(Q_sp);

    % Calculate the initial distribution.
    pi_0 = zeros(1,end_state);
    pi_0(to_state_idx(x_0.i,x_0.j,x_0.k,m)) = 1;
end

% This subfunction translate a state combination of (i,j,k) to the
% corresponding state index. It applies to the case study in Smith, R., et 
% al. (1988). "Performability analysis: measures, an algorithm, and a case study." IEEE Transactions on Computers 37(4): 406-417.

% i, j, k: numbers of working processors, memories and switches.
% m: total number of processors and memories (both assumed to be m).

function state_idx = to_state_idx(i,j,k,m)
    switch k % If the switch is up or down.
        % Switch is down: all the states are failed. So we don't need to 
        % consider i or j equals to 3, as no failures in down states.
        case 0
            % Verify if the inputs are valid or not.
            if i<4 || j<4
                error('i or j must be larger than 3, when k == 0!');
            else
                % Index: from 1 to (m-3)^2.
                state_idx = (m-i)*(m-3) + (m-j+1); 
            end
        % If the switch is on: then i or j equals to 3 is possible, but not
        % both.
        case 1
            % Verify if the inputs are valid or not.
            if i<3 || j<3
                error('i or j must be larger than 2, when k == 1!');
            end
            if i==3 && j==3
                error('i and j cannot be all 3 at the same time.');
            end
            % Index: from (m-3)^2 to (m-3)^2 + (m-2)^2.
            state_idx = (m-3)^2 + (m-i)*(m-2) + (m-j+1);
        otherwise
            error('k is not a valid state for the switch!');
    end
end