function [LDDs, slope, LD50_emp] = haber_slope_fit(surv, concs, yF)
%% Function to fit and report LD50 values when GUTS cannot fit the survival data

yInit = surv(2:2, 2:end); % initial survival number
resp = surv(end, 2:end); % Assuming last row is day 10 for LDDx,10

% Exclude control from fitting
xDRC = concs(2,3:end); yInitDRC = yInit(2:end); yEndDRC = resp(2:end); 

% Do the DRC fitting:
[LD50_emp, slope] = LDD_binomial(xDRC, yInitDRC, yEndDRC);
% Find 10% effect level at day 10
% yF = 0.9; % 1 - effect
LDD10_10d = ((1-yF)/yF)^(1/slope)*LD50_emp;
% Extrapolate out to 27 and 182d LDD10 using eq. in Annex G S7.1.2:
LDD10_27d = LDD10_10d/(2.7)^2; % 2.7 = 27/10
LDD10_182d = LDD10_10d/(18.2)^2; % 18.2 = 182/10

LDDs = [LDD10_27d, LDD10_182d];