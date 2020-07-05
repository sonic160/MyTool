function theta = rnd_theta_10_cal(tau,inv_cdf_f_20,evaluation_horizon)
    p = rand(size(tau));
    theta = inv_cdf_f_20(p);
    theta(isnan(theta)) = evaluation_horizon;
end

