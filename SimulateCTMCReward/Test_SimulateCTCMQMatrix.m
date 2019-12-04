% This is to test how to simulate a CTMC using Q matrix
% Test case: labmda = .1, mu = .2, test availability
clear; clc;
S = [1 2]; % state space
T_max = 40; % maximum time
pi = [1 0]; % initial distribution
lambda = .1;
mu = 1;
Q = [-lambda,lambda;mu,-mu];
NS = 1e5;
time = 0:10;
A_sim = zeros(size(time));
for j = time+1
    count_success = 0;
    fprintf('%d / %d interations\n',j,11);
    [t,y] = simulateCTMCQMatrix(Q,pi,time(j),NS);
    for i = 1:NS        
        y_jump = cell2num(y(i));
        if y_jump(end) == 1
            count_success = count_success + 1;
        end
    end
    A_sim(j) = count_success/NS;
end
A_real = mu/(mu+lambda)+lambda/(mu+lambda)*exp(-1*(mu+lambda).*time);
figure
plot(time,A_real,'k-d',...
    time,A_sim,'r-o')