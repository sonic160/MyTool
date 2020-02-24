function theta = rnd_theta_23(tau)
    p = rand(size(tau));
    theta = (-1 + 2000*erfinv(p + erf(1/2000) - p.*erf(1/2000)))/2000;
end

