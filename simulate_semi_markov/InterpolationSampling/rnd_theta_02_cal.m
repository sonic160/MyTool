function theta = rnd_theta_02_cal(tau,inv_cdf_f_02)
    p = rand(size(tau));
    theta = inv_cdf_f_02([p,tau]);
end

