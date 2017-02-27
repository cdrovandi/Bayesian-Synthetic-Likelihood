%% Performing BSL, uBSL and ABC

% The function 'mt19937ar.c' is required to run the below code. This code can be downloaded from http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/CODES/mt19937ar.c .

%To compile 'simulate_mex.c' on a Windows machine, use the line
mex simulate_mex.c mt19937ar.c

warning('off','all'); % suppress warnings in Matlab

%% Usual synthetic likelihood (MCMC BSL)
load('all_locations_simulated.mat');
obstimes = (1/12):(1/12):12;
cov_rw = [1e-4 0; 0 5e-8];
n = 5000; % change the value of n as desired
%parpool(16); % for parallel computing
theta_BSL = bayes_sl_cell(Yinit,S,obstimes,50000,n,cov_rw);  % NOTE: this will take many hours to run (reduce number of iterations if desired and/or use parallel computing)
%delete(gcp); % also for parallel computing

%% Unbiased synthetic likelihood (MCMC uBSL)
load('all_locations_simulated.mat');
obstimes = (1/12):(1/12):12;
cov_rw = [1e-4 0; 0 5e-8];
n = 5000; % change the value of n as desired
%parpool(16); % for parallel computing
theta_uBSL = bayes_sl_cell_go(Yinit,S,obstimes,50000,n,cov_rw);  % NOTE: this will take many hours to run (reduce number of iterations if desired and/or use parallel computing)
%delete(gcp); % also for parallel computing


%% ABC analysis
load('all_locations_simulated.mat');
load('cov_abc_cell.mat'); % pre-determined Matrix for Mahalanobis distance for ABC discrepancy function
obstimes = (1/12):(1/12):12;
cov_rw = [2.5e-4 0; 0 1e-7]; % random walk covariance matrix
theta_curr = [0.3 0.001];
tol = 1100; % change the value of tol as desired
%parpool(16); % for parallel computing
[theta_ABC, dists_ABC, summ_stats_ABC] = abc_cell_nkern(Yinit,S,obstimes,2e6,theta_curr,tol,cov_abc,cov_rw); % NOTE: this will take many hours to run (reduce number of iterations if desired)
%delete(gcp); % also for parallel computing
