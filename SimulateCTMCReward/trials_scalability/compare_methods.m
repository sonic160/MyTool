% This is to compare different methods on the multi-processor systems from
% Cloth (2006).
% Considered methods:
% - Extended Markov approach
% - Original numerical integration by Tijms
% - Monte Carlo

function [result_numit_ext,et_numit_ext,...
    result_numit,et_numit,...
    result_MC,et_MC,...
    result_vector,et_vector,nzz] = compare_methods(m,Delta,n_s)
    para.lambda = 6.89e-2; % Failure rate for the processor.
    para.gamma = 2.241e-1; % Failure rate for the memory.
    para.delta = 2.024e-1; % Failure rate for the switch.
    para.nu = 2; % Local repair rate for the processor.
    para.eta = 1; % Local repair rate for the memory.
    para.epsilon = .5; % Local repair rate for the switch.
    para.mu = .2; % Global repair rate.
    para.c = .9; % Coverage probability.
    % Initial state.
    x_0.i = m;
    x_0.j = m;
    x_0.k = 1;
    % Get the Q_matrix, r and pi_0.
    [~, Q_full, r, pi_0, nzz] = Q_matrix(m, para, x_0);
    T_max = 10; % Evaluation horizion.
    x_max = 10; % Calculate P(AR<x_max).
    % Scale the r and x_max to make them integral.
    r = floor(r*10);
    x_max = x_max*10;

    %% Numerical integration by extended Markov method.
    [result_numit_ext,et_numit_ext] = method_ext_markov(m,Delta,Q_full,pi_0,r,T_max,x_max);

    %% Numerical integration by original Tijm's approach.
    [result_numit,et_numit] = method_tijms(m,Delta,Q_full,pi_0,r,T_max,x_max);

    %% Benchmark MC.
    [result_MC,et_MC] = method_mc(n_s,Q_full,pi_0,r,T_max,x_max);

    %% Benchmark vectorized MC.
    [result_vector, et_vector] = ...
        method_vec_mc(n_s,Q_full,pi_0,r,T_max,x_max);
end



