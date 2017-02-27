function y = simulate_ricker(theta,N,T)
% simulate_ricker simulates one data set from the model
%
% INPUT:
% N - the starting population (equal to 1 in our application)
% T - the length of the data set
%
% OUTPUT:
% y - the simulated data set

y = zeros(T,1);
r = exp(theta(1)); % the parameter we sample over is log(r)
phi = theta(2);
sigmae = theta(3);

for t = 1:T
    N = r*N*exp(-N+sigmae*randn); % population size
    y(t) = poissrnd(phi*N); % the actual observed random variable
end

end
    
    