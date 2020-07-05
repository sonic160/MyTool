function theta = rnd_theta_02(tau)
    p = rand(size(tau));
    theta = -1e3*log(1-p);
end

