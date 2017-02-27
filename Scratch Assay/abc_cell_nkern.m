function [theta, dists_firsts, summ_stats_firsts] = abc_cell_nkern(Yinit,Y,obstimes,M,init,tol,cov_abc,cov_rw)
% abc_cell_nkern performs approximate Bayesian computation (ABC) and takes
% advantage of 16 cores by averaging the kernel weighting function values
% from 16 independent model simulations per proposed parameter value.
% Note the number of model simulations can be changed as desired, and
% paralell computing can also be used by using 'parfor' instead of 'for'.
%
% INPUT:
% Yinit - the initial binary matrix for cell presence
% Y - the remaining matrices for cell presence
% obstimes - the time intervals at which the cell locations were observed (in hours)
% M - the number of iterations of ABC
% init - the initial values of the parameters
% tol - the ABC tolerance
% cov_abc - estimate of the covariance matrix of the summary statistics used in Mahalanobis discrepancy function
% cov_rw - the covariance matrix of the random walk
%
% OUTPUT:
% theta - MCMC samples from the ABC target
% dists - accepted mahalanobis distances (the first of 16 draws) for regression adjustment.
% summ_stats - accepted summary statistics (the first of 16 draws) for regression adjustment.


theta_curr = init;
tau = 1/24; %Duration of each time step
ssy = cell_summstats(Y,Yinit)'; %Summary statistics of observed data
num_obs = size(Y,3);

W_abc = inv(cov_abc);

sampling_rate = max(obstimes)/num_obs; %Time between consecutive observations
sim_iters = sampling_rate/tau;

[rows, cols] = find(Yinit); %Finds non-zero indices (positions of cells)
rows = int32(rows-1);
cols = int32(cols-1);
X = int32(Yinit); %signed 32 bit integers

theta = zeros(M,2);
dists = zeros(16,1);
loglikes = zeros(16,1);
summ_stats_mult = zeros(16,145);

summ_stats_firsts = zeros(M,145);
dists_firsts = zeros(M,1);

% simulating 16 data sets
%parfor j = 1:16 % for parallel computing
for j = 1:16
    S = simulate_cell(theta_curr,X,rows,cols,num_obs,sim_iters);
    ssx = cell_summstats(S,X)';
    dists(j) = (ssx-ssy)*W_abc*(ssx-ssy)';  %Mahalanobis distance
    summ_stats_mult(j,:) = ssx;
end
dist_curr = dists(1);
summ_stats_curr = summ_stats_mult(1,:);

for j = 1:16
    loglikes(j) = -0.5/tol*dists(j)^2; %Gaussian kernel weighting function
end
loglike_curr = -log(16) + logsumexp(loglikes); %Average kernel weighting function value across 16 runs

for i = 1:M
    i  % print out iteration number if desired
    theta_prop = mvnrnd(theta_curr,cov_rw);
    
    if (any(theta_prop<0) || any(theta_prop>1)) %Rejecting negative and >1 proportions and going to next iteration
        theta(i,:) = theta_curr;
        dists_firsts(i) = dist_curr;
        summ_stats_firsts(i,:) = summ_stats_curr;
        continue;
    end
    
    % simulating 16 data sets using the proposed parameters
    %parfor j = 1:16 % for parallel computing
    for j = 1:16
        S = simulate_cell(theta_prop,X,rows,cols,num_obs,sim_iters);
        ssx = cell_summstats(S,X)';
        dists(j) = (ssx-ssy)*W_abc*(ssx-ssy)';
        summ_stats_mult(j,:) = ssx;
    end
    dist_prop = dists(1);
    summ_stats_prop = summ_stats_mult(1,:);
    
    for j = 1:16
        loglikes(j) = -0.5/tol*dists(j)^2; 
    end
    loglike_prop = -log(16) + logsumexp(loglikes);

    % MH accept-reject step
    if (rand < exp(loglike_prop - loglike_curr))
        theta_curr = theta_prop;
		dist_curr = dist_prop;
        loglike_curr = loglike_prop;
        summ_stats_curr = summ_stats_prop;
    end
    
    theta(i,:) = theta_curr;
	dists_firsts(i) = dist_curr;
    summ_stats_firsts(i,:) = summ_stats_curr;
    
end

end







