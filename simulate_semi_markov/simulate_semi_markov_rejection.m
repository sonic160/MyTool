% This is to simulate a trajectory of a homogeneous semi-Markov process.
% Rejection method is used for generating random samples
% Input: f_t - A handle to a function that generate values of f_{i,j}(t).
%        f_{i,j}(t) is the pdf of the corresponding holding time
%        distributions
%        max_value - Maximum of the pdf in its support. A n*n matrix.
%        support - A vector [lower,upper] that defines the support of the
%        pdf. A n*n cell matrix.
%        P - Transition probability matrix of the embedded Markov chain.
%        pi - Initial distribution.
%        T_max - Evaluation horizon.
% Output: t - times of the jumps.
%         y - the state after the jumps.
% History:  20191221 - Created by ZZ

function [t,y] = simulate_semi_markov_rejection(f_t,max_value,support,P,pi,T_max)
t(1) = 0; % start times at 0
y(1) = rando(pi); % generate the initial y

i = 1; % Index of the current state
while 1
    % Generate next state
    y_next = rando(P(y(i),:)); % Simulate next state using the transitin matrix of the embedded chain
    % Generate holding time
    pdf_holding_time = f_t{y(i),y_next}; % handle to the holding time distribution pdf
    upper_bound = max_value(y(i),y_next); % Get the upper bound of the pdf.
    range = support{y(i),y_next}; % Get the range of the pdf.
    t_next = sample_from_pdf(pdf_holding_time,upper_bound,1,range); % Generate random number from pdf.
    t_cur = t(i) + t_next; % Update time
    if t_cur < T_max % If this transition happens in the evaluation horizon
        % Update the parameters
        t(i+1) = t_cur;
        y(i+1) = y_next;
        i = i+1;
    else
        break;
    end
end

%  rando.m generates a random variable in 1, 2, ..., n given a distribution 
%  vector. 
function index = rando(p)
u = rand;
i = 1;
s = p(1);
while ((u > s) && (i < length(p)))
    i=i+1;
    s=s+p(i);
end
index=i;

% This is to generate a random number given a pdf.
function X = sample_from_pdf(f,M,N,b)
% SAMPLEDIST  Sample from an arbitrary distribution
%     sampleDist(f,M,N,b) retruns an array of size X of random values 
%     sampled from the distribution defined by the probability density 
%     function refered to by handle f, over the range b = [min, max].  
%     M is the threshold value for the proposal distribution, such that 
%     f(x) < M for all x in b.
%  
%     sampleDist(...,true) also generates a histogram of the results
%     with an overlay of the true pdf.
%  
%     Examples: 
%     %Sample from a step function over [0,1]:
%     X = sampleDist(@(x)1.3*(x>=0&x<0.7)+0.3*(x>=0.7&x<=1),...
%                    1.3,1e6,[0,1],true);
%     %Sample from a normal distribution over [-5,5]:
%     X = sampleDist(@(x) 1/sqrt(2*pi) *exp(-x.^2/2),...
%                    1/sqrt(2*pi),1e6,[-5,5],true);
%
% Dmitry Savransky (dsavrans@princeton.edu)
% May 11, 2010
n = 0;
X = zeros(N,1);
counter = 0;
while n < N && counter < 1000
    x = b(1) + rand(2*N,1)*diff(b);
    uM = M*rand(2*N,1);
    x = x(uM < f(x));
    X(n+1:min([n+length(x),N])) = x(1:min([length(x),N - n]));
    n = n + length(x);
    counter = counter+1;
end
if counter == 1000
    error('No sample generated.')
end