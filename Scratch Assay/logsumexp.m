function f = logsumexp(x)
% logsumexp performs a stable calculation of log(sum(exp(x)))

the_max  = max(x);
x = x - the_max;
f = the_max + log(sum(exp(x)));

end
