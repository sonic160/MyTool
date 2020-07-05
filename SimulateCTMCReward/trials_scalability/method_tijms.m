% Script to implement Tijms' method.

function [result_numit,et_numit] = method_tijms(m,Delta,Q_full,pi_0,r,T_max,x_max)
    addpath('../');
    para_1.S = 1:(m-2)^2+(m-3)^2;
    para_1.pi = pi_0;
    para_1.Q = Q_full;
    para_1.r = r;
    tic;
    result_numit = accumulated_reward_numit_tijm(T_max,x_max,Delta,para_1);
    et_numit = toc;
end