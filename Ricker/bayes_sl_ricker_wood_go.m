function theta = bayes_sl_ricker_wood_go(y,N,M,n,cov_rw)
% bayes_sl_ricker_wood_go performs MCMC uBSL on the Ricker example, using the same summary statistics as Wood (2010).
%
% INPUT:
% y - the observed data
% N - the starting population size (this will be 1 for our application)
% M - the number of iterations of BSL
% n - the number of simulated data sets uSL estimation
% cov_rw - the covariance matrix of the random walk
%
% OUTPUT:
% theta - approximate samples from the uBSL target


theta_curr = [3.8 10 0.3];
ssy = ricker_summstats(y,y)';
ns = length(ssy);
T=length(y);

theta = zeros(M,3);
ssx = zeros(n,ns);

% simulating n data sets
%parfor k = 1:n % for parallel computing
for k = 1:n
   x = simulate_ricker(theta_curr,N,T);
   ssx(k,:) = ricker_summstats(x,y)';
end

the_mean = mean(ssx);
the_cov = cov(ssx);

% estimating the unbiased SL (under multivariate normality) for current value using results from Ghurye & Olkin (1969)
loglike_ind_curr = sl_log_like_ghuryeolkin(ssy,the_mean,the_cov,n);
        
for i = 1:M
    %i  % print out iteration number if desired
    theta_prop = mvnrnd(theta_curr,cov_rw);
    if (theta_prop(3)<0) %sigma_e can't be negative
        theta(i,:) = theta_curr;
        continue;
    end
    
    % simulating n data sets using the proposed parameter value
	%parfor k = 1:n % for parallel computing
    for k = 1:n
        x = simulate_ricker(theta_prop,N,T);
        ssx(k,:) = ricker_summstats(x,y)';
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
