function f = sl_log_like_ghuryeolkin(s,muhat,sighat,m)
% sl_log_like_ghuryeolkin gets the unbiased estimator for the synthetic likelihood (when the summary statistics are MVN)
%
% INPUT:
% s - summary statistics of the observed data
% muhat - mean of the simulated summary statistics
% sighat - covariance of the simulated summary statistics
% m - the number of simulated data sets used to estimate uSL
%
% OUTPUT:
% f - the unbiased estimate of the SL (under multivariate normality of the summary statistics)

p = length(s);
temp = (m-1)*sighat- (s-muhat)'*(s-muhat)/(1-1/m);
[~,a] = chol(temp);

if (a == 0) % then positive definite
    result = -p/2*log(2*pi)+wcon(p,m-2)-wcon(p,m-1)-p/2*log(1-1/m);
    result = result-(m-p-2)/2*(log(m-1)+logdet(sighat));
    f = result + (m-p-3)/2*logdet(temp);
else
    f = -Inf;
end


function f = wcon(k,nu)
%c(k,nu) from Ghurye & Olkin (1969)
f = -k*nu/2*log(2)-k*(k-1)/4*log(pi)-sum(gammaln(0.5*(nu-(1:k)+1)));


function y = logdet(A)
% calculating the log of the determinant
U = chol(A);
y = 2*sum(log(diag(U)));


