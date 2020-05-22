% This is to implement the fast algorithm in Janssen and Manca (2001).
% Input parameters: evaluation_horizon - time limit for the evalutation
%                   Delta - step size for the approximation
%                   para - parameter structure: 
%                          para.pi_0 - Initial probability distribution
%                          para.cal_d_Q - Handle to calculate the d_Q
%                          para.cal_D - Handle to calculate the first item
%                                       in the iterative equation.
% Output parameter: state_probability - The state probability at
%                                       evaluation_horizion
% Version history: 12/05/2020: Created by ZZ
function state_probability = cal_state_p_numerical_int_use_Q(evaluation_horizon,Delta,para)
%% Initialize parameters
pi_0 = para.pi_0; % Initial distribution.
cal_Q = para.cal_Q; % Handle to calculate d_Q.
cal_D = para.cal_D; % Handle to calculate D: The first item in the iterative equation.

%% Algorithm begin: Initialization
n_t = floor(evaluation_horizon/Delta); % Number of steps in t.
Phi = cell(1,n_t+1); % Cell vector to store the result of each iteration.  The first step is zero.
Phi{n_t+1} = eye(length(pi_0));
%% Numerical integration starts.
for t = (n_t-1):-1:0
    % The first item in the iterative equation.
    Phi_cur = cal_D((n_t-t)*Delta, t*Delta);   
    % Making the consecutive summation.      
    for k = (t+1):n_t
        Phi_cur = Phi_cur + ...
            cal_Q((k-t-1)*Delta,(k-t)*Delta,t*Delta)*Phi{k+1};
    end
    Phi{t+1} = Phi_cur;
end

%% Calculate the state probability vector.
state_probability = pi_0*Phi{1};
