function [theta] = bayes_sl_cell_go(Yinit,Y,obstimes,M,n,cov_rw)
% bayes_sl_cell_go performs MCMC uBSL for the cell biology example
%
% INPUT:
% Yinit - the initial binary matrix for cell presence
% Y - the remaining matrices for cell presence
% obstimes - the time intervals at which the cell locations were observed (in hours)
% M - the number of iterations of uBSL
% n - the number of simulated data sets uSL estimation
% cov_rw - the covariance matrix of the random walk
%
% OUTPUT:
% theta - MCMC samples from the uBSL target

theta_curr = [0.3 0.001]; %Initial guesses for parameters
tau = 1/24; %Duration of each time step
ssy = cell_summstats(Y,Yinit)'; %Summary statistics of observed data
num_obs = size(Y,3);

sampling_rate = max(obstimes)/num_obs; %Time between consecutive observations
sim_iters = sampling_rate/tau;

[rows, cols] = find(Yinit); %Finds non-zero indices (positions of cells)
rows = int32(rows-1);
cols = int32(cols-1);
X = int32(Yinit); %signed 32 bit integers

theta = zeros(M,2);

ns = length(ssy);
ssx = zeros(n,ns);

% simulating n data sets
%parfor k = 1:n % for parallel computing
for k = 1:n
    S = simulate_cell(theta_curr,X,rows,cols,num_obs,sim_iters);
    ssx(k,:) = cell_summstats(S,X)';
end

the_mean = mean(ssx);
the_cov = cov(ssx);

% estimating the unbiased SL (under multivariate normality) for current value using results from Ghurye & Olkin (1969)
loglike_ind_curr = sl_log_like_ghuryeolkin(ssy,the_mean,the_cov,n);

for i = 1:M
    i  % print out iteration number if desired
    theta_prop = mvnrnd(theta_curr,cov_rw);
    if (any(theta_prop<0) || any(theta_prop>1)) %Rejecting negative and >1 proportions and going to next iteration
        theta(i,:) = theta_curr;
        continue;
    end
    
    %simulating n data sets using the proposed parameters
    %parfor k = 1:n % for parallel computing
    for k = 1:n
        S = simulate_cell(theta_prop,X,rows,cols,num_obs,sim_iters);
        ssx(k,:) = cell_summstats(S,X)';
    end
    
    the_mean = mean(ssx);
    the_cov = cov(ssx);
    
    % estimating the unbiased SL (under multivariate normality) for proposed value using results from Ghurye & Olkin (1969)
    loglike_ind_prop = sl_log_like_ghuryeolkin(ssy,the_mean,the_cov,n);
    
    % MH accept-reject step
    if (exp(loglike_ind_prop - loglike_ind_curr) > rand)
        theta_curr = theta_prop;
        loglike_ind_curr = loglike_ind_prop;
    end
    theta(i,:) = theta_curr;
    
end

end

