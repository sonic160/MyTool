function theta = rnd_theta_00(tau)
    p = rand(size(tau));
    theta = -5 - tau + 100*erfinv(p + (1-p).*erf(.05+.01*tau));
end

