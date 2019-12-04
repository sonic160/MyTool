% This is to test how to simulate a CTMC using Q matrix
% Test case: labmda = .1, mu = .2, test availability
clear; clc;
S = [1 2]; % state space
pi = [1 0]; % initial distribution
lambda = .1;
mu = 1;
Q = [-lambda,lambda;mu,-mu];
NS = 1e5;
time = 0:10;
p_state = SolveMarkovLaplace(Q,pi);
A_eval = p_state(1);
A_real = mu/(mu+lambda)+lambda/(mu+lambda)*exp(-1*(mu+lambda).*time);
figure
fplot(A_eval,[min(time),max(time)],'r-');
hold on;
plot(time,A_real,'k-d')