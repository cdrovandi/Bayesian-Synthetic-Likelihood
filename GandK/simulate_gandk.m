function q = simulate_gandk(n,parms)
% simulate_gandk simulates a data set of length n from the g-and-k distribution
%
% INPUT:
% n - number of observations
% parms - vector [a b g k] of g-and-k distribution
%
% OUTPUT:
% q - n draws from g-and-k distribution with parameters parms.

zu = randn(n,1);
q = fun_gandk(parms,zu); %the g-and-k quantile function
end