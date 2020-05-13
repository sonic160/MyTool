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
% Output parameter: f_appr - A n_x*n_S*n_t vector, which contains the
%                            approximated PDF at time prior to T_max (the $f_i^\Delta(u,y)$ in Tijms and Veldman (2000).
% Use of output: sum(sum(f_appr(:,:,end)*Delta)) gives P(O(t)<=x) (See Eq.
% (1) of Tijms and Veldman (2000)
% Version history: 15/01/2019: Created by ZZ
function f_appr = CalcARFastPDF(T_max,x_max,Delta,para)
%% Initialize parameters
S = para.S;
pi = para.pi;
Q = para.Q;
r = para.r;
%% Algorithm begin: Initialization
n_x = floor(x_max/Delta); % Number of steps in x
n_t = floor(T_max/Delta); % Number of steps in t
n_state = length(S); % Number of states
f_appr = zeros(n_x,n_state,n_t); % Initial value for the approximated pdf
Is_eff_y = zeros(1,n_x); % This is a flag vector: 1 means the corresponding y has non zero pdf.
%% The first time step: t = Delta
t_cur = 1; % Current time step
for temp_i = 1:n_x
    for S_prev = 1:n_state
        if temp_i == r(S_prev) % If the current y is jumped from S_prev
            f_appr(temp_i,S_prev,t_cur) = pi(S_prev)/Delta; % Set the pdf
            if pi(S_prev) > 0
                Is_eff_y(temp_i) = 1; % Label the non-zero pdf
            end
        else
            continue;
        end
    end
end
%% Begin iterative processes
for t_cur = 2:n_t
    t_prev = t_cur - 1; % Previous time
    % Find the possible value of y, in the previous time step
    ind_eff_y = find(Is_eff_y); % Index for the possible y
%     n_eff_y = length(ind_eff_y); % Number of possible y
    Is_eff_y = zeros(1,n_x); % Set the flag vector to zero
    for ind_y_prev = ind_eff_y % For each value of y_prev, do the one-step transition
        ind_S_prev = find(f_appr(ind_y_prev,:,t_prev)); % Locate the state with non zero pdf
        for S_prev = ind_S_prev % Consider the states one by one
            for S_cur = 1:n_state % It can jump to n_state possible states
                ind_y_cur = ind_y_prev + r(S_prev); % Update the reward
                if S_cur == S_prev % If the transition remains in the previous state
                    if ind_y_cur > n_x
                        continue;
                    else
                        Is_eff_y(ind_y_cur) = 1;
                        p = f_appr(ind_y_prev,S_prev,t_prev)*(1+Q(S_prev,S_prev)*Delta);
                        f_appr(ind_y_cur,S_cur,t_cur) = f_appr(ind_y_cur,S_cur,t_cur) + p;
                    end
                else % If it jumps to a new state
                    if ind_y_cur > n_x
                        continue;
                    else
                        Is_eff_y(ind_y_cur) = 1;
                        p = f_appr(ind_y_prev,S_prev,t_prev)*Q(S_prev,S_cur)*Delta;
                        f_appr(ind_y_cur,S_cur,t_cur) = f_appr(ind_y_cur,S_cur,t_cur) + p;
                    end                    
                end
            end
        end
    end    
end