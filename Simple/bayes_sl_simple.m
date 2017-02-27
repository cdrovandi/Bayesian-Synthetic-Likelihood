function [theta] = bayes_sl_simple(y,M,n)
% bayes_sl_simple performs MCMC BSL for the toy example
%
% INPUT:
% y - the observed data
% M - the number of iterations of BSL
% n - the number of simulated data sets SL estimation
%
% OUTPUT:
% theta - MCMC samples from the BSL target


%For Gamma distribution
alpha = 0.001;
beta = 0.001;

theta_curr = 30;
ssy = mean(y);
T = length(y);
sigma = sqrt(sum(y)+alpha)/(T+beta); %Best to use exact sigma!

theta = zeros(M,1);

% simulating n data sets
x = poissrnd(theta_curr,T,n);
ssx = mean(x);

the_mean = mean(ssx);
the_var = var(ssx);

% estimating the SL for current value
loglike_ind_curr = -normlike([the_mean sqrt(the_var)],ssy);
logprior_curr = (alpha-1)*log(theta_curr)-beta*theta_curr;

for i = 1:M
    %i  % print out iteration number if desired
    theta_prop = normrnd(theta_curr,sigma);
    while theta_prop<0 %cannot have negative parameter values
        theta_prop = normrnd(theta_curr,sigma);
    end
    
	% simulating n data sets using the proposed parameters
    x = poissrnd(theta_prop,T,n);
    ssx = mean(x);
    
    the_mean = mean(ssx);
    the_var = var(ssx);
    
    % estimating the SL for the proposed parameter values
    loglike_ind_prop = -normlike([the_mean sqrt(the_var)],ssy);
    logprior_prop = (alpha-1)*log(theta_prop)-beta*theta_prop;
    
    if (exp(loglike_ind_prop + logprior_prop - loglike_ind_curr - logprior_curr) > rand)
        theta_curr = theta_prop;
        loglike_ind_curr = loglike_ind_prop;
        logprior_curr = logprior_prop;
    end
    theta(i,:) = theta_curr;
    
end

end

