% This is to derive the expressions concerning the trainsition
% probabilities and test the subfunction that generates transition
% probabilities.

clear; clc;
%% Derive the expressions
syms x tau lambda_D eta_s beta_s eta_r beta_r p_large

% Disruption occurence
pdf_disp = lambda_D*exp(-1*lambda_D*x);
cdf_disp = 1-exp(-1*lambda_D*x);
% Failure probability of LPM
R_s = exp(-1*((x+tau)/eta_s)^beta_s);
% Recovery time cdf
cdf_t_R = 1-exp(-1*(x/eta_r)^beta_r);
pdf_t_R = diff(cdf_t_R,x);

% 0->0
d_Q_01 = pdf_disp*R_s*(1-p_large)
int(d_Q_01,x)
% 0->1
d_Q_02 = pdf_disp*(1-R_s)*(1-p_large)
int(d_Q_02,x)
% 0->3
d_Q_03 = pdf_disp*p_large
int(d_Q_03,x)
% 2->3
d_Q_23 = pdf_disp*(1-cdf_t_R)
int(d_Q_23,x)
% 2->0
d_Q_20 = pdf_t_R*(1-cdf_disp)
int(d_Q_20,x)

D_00 = 1-cdf_disp
D_11 = 1-cdf_disp
D_22 = (1-cdf_disp)*(1-cdf_t_R)