function f = fun_gandk(parms,zu)
% fun_gandk is the g-and-k quantile function
%
% INPUT:
% parms - vector [a b g k] of g-and-k distribution
% zu - quantiles of the standard normal distribution (vector)
%
% OUTPUT:
% f - the calculated quantiles

a = parms(1); b = parms(2); c = 0.8; g = parms(3); k = parms(4);
f = a + b*(1 + c*(1-exp(-g*zu))./(1 + exp(-g*zu))).*(1 + zu.^2).^k.*zu;

end