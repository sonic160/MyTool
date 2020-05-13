% This is a test of the simulation algorithm for inhomogeneous semi-Markov process.
% Invere cdf based on interpolation is used to sample the holding time distribution.
% Benchmark: A four state system, where state transitions are caused by
% arrival of disruptive events, performance of safety barriers, and
% recovery time distributions.
%            State 0: Perfect
%            State 1: Dummy perfect state.
%            State 2: Performance degradation.
%            State 3: Total failure (absorption state).
% Last edit: 20200128 - Add comments.
%            20191223 - Created by ZZ
clear; clc;
%% Parameter definition
% First arrival time of the disruptive event
lambda_D = 1e-3;
pd_disp = makedist('exponential','mu',1/lambda_D);
% Safety barriers
p_large = 1.6e-3; % Probability of large damage, leading to state 3
eta_s = 100; beta_s = 2;
pd_t_f_s = makedist('weibull','a',eta_s,'b',beta_s); % Time to failure distribution of the safety barrier
% Time to repair distribution
eta_r = 1; beta_r = 2;
pd_t_r = makedist('weibull','a',eta_r,'b',beta_r); 

% Simulation to calculate p(t) for each state, for 0<t<T_max
evaluation_horizon = 1e3;
t = linspace(.1,evaluation_horizon,20); % Evaluate time

r = [4,3,2,1];
y_th = 3000;
F_ref = zeros(1,length(t));

ns = 1e6; % Sample size

%% Create benchmark
sample_path = cell(1,ns); % Store all the sample path up to evaluation_horizon.
sample_jump_time = cell(1,ns); % Store all the jump times up to evaluation_horizon.

tic;
accumulated_reward = zeros(1,ns);
for j = 1:ns
    t_cur = 0; % Initial value of time
    state = 0; % Start from working
    
    % Initialize
    temp_sample_path = -1*ones(1,1000); % Vector for the sample path.
    temp_sample_jump_time = -1*ones(1,1000); % Vector for the sample path.
    
    index = 1; % Intial value for the index of the states.
    temp_sample_path(index) = 0; 
    temp_sample_jump_time(index) = 0; 
    
    while 1
        index = index + 1; % Next state.
        switch state
            case 0 % If current state is perfect state
                theta = pd_disp.random(); % Interarrival time
                t_next = t_cur + theta; % Calculate the event occurrence tiem
                if t_next > evaluation_horizon % Beyond evaluation horizon
                    t_cur = evaluation_horizon;
                    
                    temp_sample_path(index) = state;
                    temp_sample_jump_time(index) = t_cur;
                    accumulated_reward(j) = accumulated_reward(j) + r(state+1)*(evaluation_horizon-t_next+theta);
                    break;
                else
                    accumulated_reward(j) = accumulated_reward(j) + r(state+1)*theta;
                    t_cur = t_next;
                    % Decide what is the next state
                    is_large_damage = binornd(1,p_large);
                    is_fail_sb = binornd(1,pd_t_f_s.cdf(t_cur));                     
                    if is_large_damage == 1 % If large damage happens
                        state = 3; % Absorption state
                    else
                        if is_fail_sb == 1 % If safety barrier fails
                            state = 2; % Degradation state
                        else % If safety barrier is OK
                            state = 1; % Dummy state 0
                        end
                    end
                                        
                    temp_sample_path(index) = state;
                    temp_sample_jump_time(index) = t_cur;                    
                end
            case 1 % Dummy perfect state: it is equivalant to the perfect state
                theta = pd_disp.random(); % Interarrival time
                t_next = t_cur + theta; % Calculate the event occurrence tiem
                if t_next > evaluation_horizon % Beyond evaluation horizon
                    t_cur = evaluation_horizon;                   
                                        
                    temp_sample_path(index) = state;
                    temp_sample_jump_time(index) = t_cur;
                    accumulated_reward(j) = accumulated_reward(j) + r(state+1)*(evaluation_horizon-t_next+theta);
                    break;
                else
                    accumulated_reward(j) = accumulated_reward(j) + r(state+1)*theta;
                    t_cur = t_next;
                    % Decide what is the next state
                    is_large_damage = binornd(1,p_large);
                    is_fail_sb = binornd(1,pd_t_f_s.cdf(t_cur));                     
                    if is_large_damage == 1 % If large damage happens
                        state = 3; % Absorption state
                    else
                        if is_fail_sb == 1 % If safety barrier fails
                            state = 2; % Degradation state
                        else % If safety barrier is OK
                            state = 0; % Dummy state 0
                        end
                    end                    
                                        
                    temp_sample_path(index) = state;
                    temp_sample_jump_time(index) = t_cur;                    
                end
            case 2
                theta_d = pd_disp.random(); % The next degradation time
                theta_r = pd_t_r.random(); % The next recovery time
                if theta_d < theta_r % If next transition is degradation
                    theta = theta_d;
                    t_next = t_cur + theta; % Calculate the event occurrence tiem
                    if t_next > evaluation_horizon % Beyond evaluation horizon
                        t_cur = evaluation_horizon;                                               
                                            
                        temp_sample_path(index) = state;
                        temp_sample_jump_time(index) = t_cur;
                        accumulated_reward(j) = accumulated_reward(j) + r(state+1)*(evaluation_horizon-t_next+theta);
                        break;
                    else
                        accumulated_reward(j) = accumulated_reward(j) + r(state+1)*theta;
                        t_cur = t_next;
                        state = 3;                        
                                            
                        temp_sample_path(index) = state;
                        temp_sample_jump_time(index) = t_cur;                        
                    end
                else % If next state is recovery
                    theta = theta_r;
                    t_next = t_cur + theta; % Calculate the event occurrence tiem
                    if t_next > evaluation_horizon % Beyond evaluation horizon
                        t_cur = evaluation_horizon;                      
                                            
                        temp_sample_path(index) = state;
                        temp_sample_jump_time(index) = t_cur;
                        accumulated_reward(j) = accumulated_reward(j) + r(state+1)*(evaluation_horizon-t_next+theta);
                        break;
                    else
                        accumulated_reward(j) = accumulated_reward(j) + r(state+1)*theta;
                        t_cur = t_next;
                        state = 0;
                                            
                        temp_sample_path(index) = state;
                        temp_sample_jump_time(index) = t_cur;                        
                    end
                end
            case 3
                t_cur = evaluation_horizon;                             

                temp_sample_path(index) = state;
                temp_sample_jump_time(index) = t_cur;
                accumulated_reward(j) = accumulated_reward(j) + r(state+1)*(evaluation_horizon-temp_sample_jump_time(index-1));
                break;                    
            otherwise
                error('Undefined state!')
        end            
    end
    
    temp_sample_path = temp_sample_path(1:index);
    temp_sample_jump_time = temp_sample_jump_time(1:index);
    sample_path{j} = temp_sample_path;
    sample_jump_time{j} = temp_sample_jump_time;
end

%% Post processing.
for i = 1:length(t)
    tt = t(i);
    count = 0;
    for j = 1:ns
        temp_temp_sample_path = sample_path{j};
        temp_temp_sample_jump_time = sample_jump_time{j};
        
        index = temp_temp_sample_jump_time<tt;
        temp_temp_sample_path = temp_temp_sample_path(index);
        temp_temp_sample_jump_time = temp_temp_sample_jump_time(index);
        temp_sample_path = [temp_temp_sample_path,temp_temp_sample_path(end)];
        temp_sample_jump_time = [temp_temp_sample_jump_time,tt];
        
        time_interval = temp_sample_jump_time(2:end)-temp_sample_jump_time(1:end-1); % Get time interval.
        state_cur = temp_sample_path(1:end-1); % Get the current state.
        reward = sum(r(state_cur+1).*time_interval);
        count = count + (reward<=y_th);
    end
    F_ref(i) = count/ns;
end

elapsed_time_1 = toc;

save('result_ref_monte_carlo.mat');