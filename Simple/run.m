%% Performing BSL, uBSL and ABC

warning('off','all'); % suppress warnings in Matlab

%% Usual synthetic likelihood (MCMC BSL)
load('data_poisson.mat');
n = 5; % change the value of n as desired
theta_BSL = bayes_sl_simple(y,100000,n); % NOTE: this will take many hours to run (reduce number of iterations if desired)

%% Unbiased synthetic likelihood (MCMC uBSL)
load('data_poisson.mat');
n = 5; % change the value of n as desired
theta_uBSL = bayes_sl_simple_go(y,100000,n); % NOTE: this will take many hours to run (reduce number of iterations if desired)

%% ABC analysis
load('data_poisson.mat');
theta_curr = [3.8 10 0.3];
tol = 0.001; %change the value of tol as desired
theta_ABC = abc_simple(y,500000,30,tol); % NOTE: this will take many hours to run (reduce number of iterations if desired)

