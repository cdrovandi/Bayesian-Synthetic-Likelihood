function [theta, dists, summ_stats] = abc_ricker_wood(y,N,M,init,tol,cov_abc,cov_rw)
% abc_ricker_wood performs approximate Bayesian computation (ABC) on the Ricker example, using the same summary statistics as Wood (2010).
%
% INPUT:
% y - the observed data
% N - the starting population size (this will be 1 for our application)
% M - the number of iterations of ABC
% init - the initial values of the parameters
% tol - the ABC tolerance
% cov_abc - estimate of the covariance matrix of the summary statistics used in Mahalanobis discrepancy function
% cov_rw - the covariance matrix of the random walk
%
% OUTPUT:
% theta - MCMC samples from the ABC target
% dists - accepted mahalanobis distances
% summ_stats - accepted summary statistics

theta_curr = init;
ssy = ricker_summstats(y,y); %calculating the summary statistics of the observed data, using Wood (2010) summary statistics
ns = length(ssy);
T=length(y);

W_abc = inv(cov_abc);

theta = zeros(M,3);

dists = zeros(M,1);
summ_stats = zeros(M,ns);

x = simulate_ricker(theta_curr,N,T); % simulating one data set from the model
ssx = ricker_summstats(x,y); %calculating the summary statistics for the current data set
 
dist_curr = (ssy-ssx)'*W_abc*(ssy-ssx); %Mahalanobis distance

loglike_curr = -0.5/tol*dist_curr^2; %Gaussian kernel weighting function
summ_stat_curr = ssx';
       
for i = 1:M
    %i  % print out iteration number if desired
    theta_prop = mvnrnd(theta_curr,cov_rw);
    if (theta_prop(3)<0) %can't have negative sigma_e
        theta(i,:) = theta_curr;
        summ_stats(i,:) = summ_stat_curr;
        dists(i) = dist_curr;
        continue;
    end
    
    x = simulate_ricker(theta_prop,N,T); % simulating one data set from the model based on proposed parameters
    ssx = ricker_summstats(x,y); % calculating the summary statistics for the proposed data set
    
    dist_prop = (ssy-ssx)'*W_abc*(ssy-ssx); %Mahalanobis distance
    
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
