% Script to implement the extended Markov approach.

function [result_numit_ext,et_numit_ext] = method_ext_markov(m,Delta,Q_full,pi_0,r,T_max,x_max)
    addpath('../');
    addpath('../compare_different_codings/')
    para_1.S = 1:(m-2)^2+(m-3)^2;
    para_1.pi = pi_0;
    para_1.Q = Q_full;
    para_1.r = r;
    tic;
    result_numit_ext = accumulated_reward_numit_ext_markov(T_max,x_max,Delta,para_1);
    et_numit_ext = toc;
end