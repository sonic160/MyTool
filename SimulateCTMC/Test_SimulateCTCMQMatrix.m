% This is to test how to simulate a CTMC using Q matrix
% Test case: labmda = .1, mu = .2, test steady state availability
clear; clc;
S = [1 2]; % state space
T_max = 40; % maximum time
pi = [1 0]; % initial distribution
lambda = .1;
mu = 1;
Q = [-lambda,lambda;mu,-mu];
[a,P] = CalculatePEmbedded(Q);
NS = 1e5;
count_success = 0;
for i = 1:NS
    disp([num2str(i) '/' num2str(NS)])
    [t,y] = SimulateCTMCEmbedded(a,P,pi,T_max);
    n = length(t);
    for j = n:-1:1
        if t(j) < T_max
            if y(j) == 1
                count_success = count_success + 1;
            end
            break;
        end
    end
end
A_sim = count_success/NS
A_real = mu/(mu+lambda)*(1+exp(-1*(mu+lambda)*T_max))