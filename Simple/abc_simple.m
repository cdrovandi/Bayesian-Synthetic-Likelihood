function theta = abc_simple(y,M,init,tol)
% abc_simple performs approximate Bayesian computation (ABC) for toy
% example
%
% INPUT:
% y - the observed data
% M - the number of iterations of ABC
% init - the initial values of the parameters
% tol - the ABC tolerance
%
% OUTPUT:
% theta - MCMC samples from the ABC target

%For Gamma distribution
alpha = 0.001;
beta = 0.001;

theta_curr = init;
ssy = mean(y);
T = length(y); %Size of initial data set, y
sigma = sqrt(sum(y)+alpha)/(T+beta); %Best to use exact sigma!

theta = zeros(M,1);

ssx = mean(poissrnd(theta_curr,T,1));
dist_curr = (ssy-ssx)^2;

loglike_curr = -0.5/tol*dist_curr^2;
logprior_curr = (alpha-1)*log(theta_curr)-beta*theta_curr;
       
for i = 1:M
    %i  % print out iteration number if desired
    theta_prop = normrnd(theta_curr,sigma);
    while theta_prop<0 %cannot have negative parameter values
        theta_prop = normrnd(theta_curr,sigma);
    end
        
    ssx = mean(poissrnd(theta_prop,T,1));
    dist_prop = (ssy-ssx)^2;
    
    loglike_prop = -0.5/tol*dist_prop^2;
    logprior_prop = (alpha-1)*log(theta_prop)-beta*theta_prop;
    
    % MCMC ABC accept/reject step
    if (rand<exp(loglike_prop + logprior_prop - loglike_curr - logprior_curr))
        theta_curr = theta_prop;
        loglike_curr = loglike_prop;
        logprior_curr = logprior_prop;
    end
    theta(i,:) = theta_curr;
    
end

end
