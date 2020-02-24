function theta = rnd_theta_20_cal(tau,inv_cdf_f_20)
    p = rand(size(tau));
    theta = inv_cdf_f_20(p);
end

