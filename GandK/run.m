%% Performing BSL, uBSL and ABC using the skew t score of the simulated data (with skew t MLEs estimated from observed data).

warning('off','all'); % suppress warnings in Matlab

%% Usual synthetic likelihood (MCMC BSL)
load('data_gandk.mat');
load('Tskew_to_GandK_MLEs.mat'); %skew t MLEs based on observed data
cov_rw = [0.0001,0.0001,-0.00015,-9e-05;0.0001,0.0004,0.0003,-0.00048;-0.00015,0.0003,0.0025,0;-9e-05,-0.00048,0,0.0009]; % random walk covariance matrix
n = 10; % change the value of n as desired
%parpool(16); % for parallel computing
theta_BSL = bayes_sl_gandk(y,500000,n,cov_rw,a,b,mut,sigmas); % NOTE: this will take many hours to run (reduce number of iterations if desired and/or use parallel computing)
%delete(gcp); % also for parallel computing

%% Unbiased synthetic likelihood (MCMC uBSL)
load('data_gandk.mat');
load('Tskew_to_GandK_MLEs.mat'); %skew t MLEs based on observed data
cov_rw = [0.0001,0.0001,-0.00015,-9e-05;0.0001,0.0004,0.0003,-0.00048;-0.00015,0.0003,0.0025,0;-9e-05,-0.00048,0,0.0009]; % random walk covariance matrix
n = 10; % change the value of n as desired
%parpool(16); % for parallel computing
theta_uBSL = bayes_sl_gandk_go(y,500000,n,cov_rw,a,b,mut,sigmas); % NOTE: this will take many hours to run (reduce number of iterations if desired and/or use parallel computing)
%delete(gcp); % also for parallel computing


%% ABC analysis
load('data_gandk.mat');
load('Tskew_to_GandK_MLEs.mat'); %skew t MLEs based on observed data;
load('cov_abc_gandk.mat'); % pre-determined Matrix for Mahalanobis distance for ABC discrepancy function
cov_rw = [0.0001,0.0001,-0.00015,-9e-05;0.0001,0.0004,0.0003,-0.00048;-0.00015,0.0003,0.0025,0;-9e-05,-0.00048,0,0.0009]; % random walk covariance matrix
theta_curr = [3 1 2 0.5];
tol = 0.35; %change the value of tol as desired
[theta_ABC, dists_ABC, summ_stats_ABC] = abc_gandk(y,2e6,theta_curr,tol,cov_abc,cov_rw,a,b,mut,sigmas); % NOTE: this will take many hours to run (reduce number of iterations if desired)
