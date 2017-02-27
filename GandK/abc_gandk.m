function [theta, dists, summ_stats] = abc_gandk(y,M,init,tol,cov_abc,cov_rw,a,b,mut,sigmas)
% abc_gandk performs approximate Bayesian computation (ABC) on the g-and-k example.
%
% INPUT:
% y - the observed data
% M - the number of iterations of ABC
% init - the initial values of the parameters
% tol - the ABC tolerance
% cov_abc - estimate of the covariance matrix of the summary statistics used in Mahalanobis discrepancy function
% cov_rw - the covariance matrix of the random walk
% a - MLE estimate for skew t parameter 'a' based on observed data
% b - MLE estimate for skew t parameter 'b' based on observed data
% mut - MLE estimate for skew t parameter 'mu' based on observed data
% sigmas - MLE estimate for skew t parameter 'sigma' based on observed data
%
% OUTPUT:
% theta - MCMC samples from the ABC target
% dists - accepted mahalanobis distances
% summ_stats - accepted summary statistics

theta_curr = init;
T = length(y);

W_abc = inv(cov_abc);

theta = zeros(M,4);

dists = zeros(M,1);
summ_stats = zeros(M,4);

x = simulate_gandk(T,theta_curr); % simulating one data set from the model
ssx = Scores(a,b,mut,sigmas,x); %calculating the summary statistics (score) for the current data set

dist_curr = ssx*W_abc*ssx'; %Mahalanobis distance

loglike_curr = -0.5/tol*dist_curr^2; %Gaussian kernel weighting function
summ_stat_curr = ssx';

        
for i = 1:M
    %i  % print out iteration number if desired
    theta_prop = mvnrnd(theta_curr,cov_rw);
    if (any(theta_prop<0) || any(theta_prop>10))
        theta(i,:) = theta_curr;
        summ_stats(i,:) = summ_stat_curr;
        dists(i) = dist_curr;
        continue;
    end
    
    x = simulate_gandk(T,theta_prop); % simulating one data set from the model based on proposed parameters
    ssx = Scores(a,b,mut,sigmas,x); % calculating the summary statistics (score) for the proposed data set
        
    dist_prop = ssx*W_abc*ssx'; %Mahalanobis distance
    
    loglike_prop = -0.5/tol*dist_prop^2; %Gaussian kernel weighting function
    
    % MCMC ABC accept/reject step
    if (rand<exp(loglike_prop - loglike_curr))
        theta_curr = theta_prop;
        dist_curr = dist_prop;
        summ_stat_curr = ssx';
        loglike_curr = loglike_prop;
    end
    theta(i,:) = theta_curr;
    summ_stats(i,:) = summ_stat_curr;
    dists(i) = dist_curr;
    
end

end
