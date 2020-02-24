function theta = rnd_theta_03(tau)
    p = rand(size(tau));
    theta = -1e3*log(1-p);
end

