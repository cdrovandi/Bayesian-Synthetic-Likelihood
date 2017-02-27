%% Performing BSL, uBSL and ABC using the same summary statistics as in Wood (2010)

warning('off','all'); % suppress warnings in Matlab

%% Usual synthetic likelihood (MCMC BSL)
load('data_ricker.mat');
std_rw = [0.15 0.5 0.15];
corr_rw = [1 -0.7 -0.6; -0.7 1 0.4; -0.6 0.4 1];
cov_rw = corr2cov(std_rw,corr_rw); % random walk covariance matrix
n = 50; % change the value of n as desired
%parpool(16); % for parallel computing
theta_BSL = bayes_sl_ricker_wood(y,1,500000,n,cov_rw);  % NOTE: this will take many hours to run (reduce number of iterations if desired and/or use parallel computing)
%delete(gcp); % also for parallel computing

%% Unbiased synthetic likelihood (MCMC uBSL)
load('data_ricker.mat');
std_rw = [0.15 0.5 0.15];
corr_rw = [1 -0.7 -0.6; -0.7 1 0.4; -0.6 0.4 1];
cov_rw = corr2cov(std_rw,corr_rw); % random walk covariance matrix
n = 50; % change the value of n as desired
%parpool(16); % for parallel computing
theta_uBSL = bayes_sl_ricker_wood_go(y,1,500000,n,cov_rw); % NOTE: this will take many hours to run (reduce number of iterations if desired and/or use parallel computing)
%delete(gcp); % also for parallel computing


%% ABC analysis
load('data_ricker.mat');
load('cov_abc_ricker_wood.mat'); % pre-determined Matrix for Mahalanobis distance for ABC discrepancy function
std_rw = [0.15 0.5 0.15];
corr_rw = [1 -0.7 -0.6; -0.7 1 0.4; -0.6 0.4 1];
cov_rw = corr2cov(std_rw,corr_rw); % random walk covariance matrix
theta_curr = [3.8 10 0.3];
tol = 5; %change the value of tol as desired
[theta_ABC, dists_ABC, summ_stats_ABC] = abc_ricker_wood(y,1,25000000,theta_curr,tol,cov_abc,cov_rw); % NOTE: this will take many hours to run (reduce number of iterations if desired)

