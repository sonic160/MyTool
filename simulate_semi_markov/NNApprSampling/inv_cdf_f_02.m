% This is to generate a random sample from the cdf of f_02, using the
% NN-based approximation.

function t = inv_cdf_f_02(p,tau,net_upper,net_tail,p_th)
    if p<p_th
        t = net_upper([p;tau]);
    else
        t = net_tail([p;tau]);
    end
end