% This is a test of the simulation algorithm for inhomogeneous semi-Markov process.
% Invere cdf based on NN is used to sample the holding time distribution.
% Benchmark: A four state system, where state transitions are caused by
% arrival of disruptive events, performance of safety barriers, and
% recovery time distributions.
%            State 0: Perfect
%            State 1: Dummy perfect state.
%            State 2: Performance degradation.
%            State 3: Total failure (absorption state).
% Last edit: 20191223 - Created by ZZ
clear; clc;
%% Create benchmark
% Parameter definition
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
T_max = 1e3;
t = linspace(.1,T_max,10); % Evaluate time
ns = 1e5; % Sample size
p_t_ref = zeros(4,length(t)); % p(t): 4*1, each column corresponds to one state

tic;

for i = 1:length(t)
    fprintf('%d/%d\n',i,length(t));
    count_0 = 0; count_1 = 0; count_2 = 0; count_3 = 0; % Counter for ending at each state
    for j = 1:ns
        t_cur = 0; % Initial value of time
        state = 0; % Start from working
        while 1
            switch state
                case 0 % If current state is perfect state
                    theta = pd_disp.random(); % Interarrival time
                    t_next = t_cur + theta; % Calculate the event occurrence tiem
                    if t_next > t(i) % Beyond evaluation horizon
                        t_cur = t(i);
                        count_0 = count_0 + 1;
                        break;
                    else
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
                    end
                case 1 % Dummy perfect state: it is equivalant to the perfect state
                    theta = pd_disp.random(); % Interarrival time
                    t_next = t_cur + theta; % Calculate the event occurrence tiem
                    if t_next > t(i) % Beyond evaluation horizon
                        t_cur = t(i);
                        count_1 = count_1 + 1;
                        break;
                    else
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
                    end
                case 2
                    theta_d = pd_disp.random(); % The next degradation time
                    theta_r = pd_t_r.random(); % The next recovery time
                    if theta_d < theta_r % If next transition is degradation
                        theta = theta_d;
                        t_next = t_cur + theta; % Calculate the event occurrence tiem
                        if t_next > t(i) % Beyond evaluation horizon
                            t_cur = t(i);
                            count_2 = count_2 + 1;
                            break;
                        else
                            t_cur = t_next;
                            state = 3;
                        end
                    else % If next state is recovery
                        theta = theta_r;
                        t_next = t_cur + theta; % Calculate the event occurrence tiem
                        if t_next > t(i) % Beyond evaluation horizon
                            t_cur = t(i);
                            count_2 = count_2 + 1;
                            break;
                        else
                            t_cur = t_next;
                            state = 0;
                        end
                    end
                case 3
                    t_cur = t(i);
                    count_3 = count_3 + 1;
                    break;                    
                otherwise
                    error('Undefined state!')
            end            
        end       
    end
    p_t_ref(1,i) = count_0/ns;
    p_t_ref(2,i) = count_1/ns;
    p_t_ref(3,i) = count_2/ns;
    p_t_ref(4,i) = count_3/ns;
end

elapsed_time_1 = toc;

% Display results
h_p_1 = figure();
plot(t,p_t_ref(1,:),'-o')
hold on
h_p_2 = figure();
plot(t,p_t_ref(2,:),'-o')
hold on
h_p_3 = figure();
plot(t,p_t_ref(3,:),'-o')
hold on
h_p_4 = figure();
plot(t,p_t_ref(4,:),'-o')
hold on

%% Simulate using holding time distribution and embedded chain. Sampling with rejection method.
load('inverse_cdf_f_20.mat')
load('inverse_cdf_f_02.mat')

% Define the embedded chain and the holding time dist.
pi_0 = [1,0,0,0]; % Initial distribution
p_01 = @(tau) (-0.887024E-1).*exp(1).^((1/1000).*tau).*...
    ((-0.1E1)+erf((1/100).*(5+tau)));
p_02 = @(tau) 0.9984E0+0.887024E-1.*exp(1).^((1/1000).*tau).*((-0.1E1)+...
    erf(0.5E-1+ 0.1E-1.*tau));
p_03 = @(tau) .16e-2;
p_10 = p_01;
p_12 = p_02;
p_13 = p_03;
p_20 = @(tau) .999114;
p_23 = @(tau) .000885727;
handle_zero = @(tau) 0;
handle_one = @(tau) 1;
P = {handle_zero, p_01, p_02, p_03;...
    p_10, handle_zero, p_12, p_13;...
    p_20, handle_zero, handle_zero, p_23;...
    handle_zero, handle_zero, handle_zero, handle_one}; % Transition probabability matrix

% Inverse functions of the CDF.
inv_01 = @(p,tau) -5 - tau + 100*erfinv(p + (1-p).*erf(.05+.01*tau));
inv_23 = @(p,tau) (-1 + 2000*erfinv(p + erf(1/2000) - p.*erf(1/2000)))/2000;
inv_03 = @(p,tau) -1e3*log(1-p);

inv_02 = @(p,tau) inv_cdf_f_02(p,tau,net_f_02_upper,net_f_02_tail,p_th_f_02);
inv_20 = @(p,tau) inv_cdf_f_20(p,net_f_20_upper,net_f_20_tail,p_th);

inv_absorb = @(p,tau) 1e8; % A very large number

inv_F_t = {handle_zero, inv_01, inv_02, inv_03;...
    inv_01, handle_zero, inv_02, inv_03;...
    inv_20, handle_zero, handle_zero, inv_23;...
    inv_absorb, inv_absorb, inv_absorb,inv_absorb};

% Calculate the same probability distribution
fprintf('\n PDF approach\n')
p_t_sim = zeros(4,length(t)); % p(t): 4*1, each column corresponds to one state

tic;

for i = 1:length(t)
    fprintf('%d/%d\n',i,length(t));
    count_0 = 0; count_1 = 0; count_2 = 0; count_3 = 0; % Counter
    for j = 1:ns
        [temp_t,y] = simulate_semi_markov_NNinv(inv_F_t,P,pi_0,T_max);
        switch y(end)
            case 1
                count_0 = count_0 + 1;
            case 2
                count_1 = count_1 + 1;
            case 3
                count_2 = count_2 + 1;
            case 4
                count_3 = count_3 + 1;
        end
    end
    p_t_sim(1,i) = count_0/ns;
    p_t_sim(2,i) = count_1/ns;
    p_t_sim(3,i) = count_2/ns;
    p_t_sim(4,i) = count_3/ns;
end

elapsed_time_2 = toc;

% Display results
figure(h_p_1);
plot(t,p_t_sim(1,:),'-dr')
figure(h_p_2);
plot(t,p_t_sim(2,:),'-dr')
figure(h_p_3);
plot(t,p_t_sim(3,:),'-dr')
figure(h_p_4);
plot(t,p_t_sim(4,:),'-dr')

fprintf('The elapsed time for 1 is %f seconds\n',elapsed_time_1);
fprintf('The elapsed time for 2 is %f seconds\n',elapsed_time_2);