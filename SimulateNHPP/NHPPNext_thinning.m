% Simulates the first arrival time of an NHPP after t = tau, using thinning
% method
% Input: handle-lambda - function handle to the intensity function
%        t_0: starting time
%        T_max: maximum time in the thinning method
%        t_next - Arrical time of next event of the NHPP
% History: 20170204 - Created by ZZ
function t_next = NHPPNext_thinning(handle_lambda,t_0,T_max)
if nargin == 2
    T_max = 4/handle_lambda(t_0); % default value
end
t = t_0; % current time
T_max = t_0 + T_max; % shift the origin of T_max
lambda = handle_lambda(T_max); % intensity at T_max
while 1
    % Generate arrival time of the next event of HPP with lambda
    t = t + exprnd(1/lambda);
    if t < T_max % if the next arrival is before T_max, do the accept/reject
        % Reject or accept
        v = rand;
        if v < handle_lambda(t)/lambda % accept 
            t_next = t;
            break;        
        else
            continue;
        end
    else % if the next arrival is beyond T_max, increase T_max and start again with t = T_max
        t = T_max;
        T_max = T_max + 1/handle_lambda(t_0);
        lambda = handle_lambda(T_max); % intensity at T_max
    end
end