function [theta] = bayes_sl_gandk(y,M,n,cov_rw,a,b,mut,sigmas)
% bayes_sl_gandk performs MCMC BSL on the g-and-k example
%
% INPUT:
% y - the observed data
% M - the number of iterations of BSL
% n - the number of simulated data sets SL estimation
% cov_rw - the covariance matrix of the random walk
% a - MLE estimate for skew t parameter 'a' based on observed data
% b - MLE estimate for skew t parameter 'b' based on observed data
% mut - MLE estimate for skew t parameter 'mu' based on observed data
% sigmas - MLE estimate for skew t parameter 'sigma' based on observed data
%
% OUTPUT:
% theta - MCMC samples from the BSL target


theta_curr = [3 1 2 0.5];
T = length(y);

theta = zeros(M,4);
score = zeros(n,4);

% simulating n data sets
%parfor k = 1:n % for parallel computing
for k = 1:n
    x = simulate_gandk(T,theta_curr);
    score(k,:) = Scores(a,b,mut,sigmas,x); %skew t score for simulated data (with skew t maximum likelihood estimates from observed data)
end

the_mean = mean(score);
the_cov = cov(score);

% estimating the SL for current value
loglike_ind_curr = -2*log(2*pi) - 0.5*log(det(the_cov)) - 0.5*the_mean*inv(the_cov)*the_mean';

for i = 1:M
    %i  % print out iteration number if desired
    theta_prop = mvnrnd(theta_curr,cov_rw);
	
	%simulating n data sets using the proposed parameters
    %parfor k = 1:n % for parallel computing
    for k = 1:n
        x = simulate_gandk(T,theta_prop);
        score(k,:) = Scores(a,b,mut,sigmas,x); %skew t score for simulated data (with skew t maximum likelihood estimates from observed data)
    end
	
    the_mean = mean(score);
    the_cov = cov(score);
	
	% estimating the SL for proposed value
    loglike_ind_prop = -2*log(2*pi) - 0.5*log(det(the_cov)) - 0.5*the_mean*inv(the_cov)*the_mean';
	
	% MH accept-reject step
    if (exp(loglike_ind_prop - loglike_ind_curr) > rand)
        theta_curr = theta_prop;
        loglike_ind_curr = loglike_ind_prop;
    end
    theta(i,:) = theta_curr;
	
end

end

