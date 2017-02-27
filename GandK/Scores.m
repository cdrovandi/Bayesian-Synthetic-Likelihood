function [score_results] = Scores(a,b,mut,sigmas,x)
% Scores computes the score for the simulated data, with the skew t
% parameters set to the MLEs for the observed data.
%
% INPUT:
% a - MLE estimate for skew t parameter 'a' based on observed data
% b - MLE estimate for skew t parameter 'b' based on observed data
% mut - MLE estimate for skew t parameter 'mu' based on observed data
% sigmas - MLE estimate for skew t parameter 'sigma' based on observed data
% x - simulated data
%
% OUTPUT:
% score_results - the vector of scores


m = length(x);

%dlogL_da = -m*log(2)-m*psi(a)+m*psi(a+b)-m/2/(a+b)+sum(d1_da+d2_da);
dlogL_da = -m*log(2)-m*psi(a)+m*psi(a+b)-m/2/(a+b)+sum(...
    log(1 - (mut - x)./(sigmas*(a + b + (mut - x).^2./sigmas^2).^(1/2))) - ((mut - x).*(a + 1/2))./(2*sigmas*((mut - x)./(sigmas*(a + b + (mut - x).^2./sigmas^2).^(1/2)) - 1).*(a + b + (mut - x).^2/sigmas^2).^(3/2))...
    -((mut - x).*(b + 1/2))./(2*sigmas*((mut - x)./(sigmas*(a + b + (mut - x).^2./sigmas^2).^(1/2)) + 1).*(a + b + (mut - x).^2./sigmas^2).^(3/2)));

%dlogL_db = -m*log(2)-m*psi(b)+m*psi(a+b)-m/2/(a+b)+sum(d1_db+d2_db);
dlogL_db = -m*log(2)-m*psi(b)+m*psi(a+b)-m/2/(a+b)+sum(...
    -((mut - x).*(a + 1/2))./(2*sigmas*((mut - x)./(sigmas*(a + b + (mut - x).^2./sigmas^2).^(1/2)) - 1).*(a + b + (mut - x).^2./sigmas^2).^(3/2))...
    +log((mut - x)./(sigmas*(a + b + (mut - x).^2./sigmas^2).^(1/2)) + 1) - ((mut - x).*(b + 1/2))./(2*sigmas*((mut - x)./(sigmas*(a + b + (mut - x).^2./sigmas^2).^(1/2)) + 1).*(a + b + (mut - x).^2./sigmas^2).^(3/2)));

%dlogL_dmu = sum(d1_dmu+d2_dmu);
dlogL_dmu = sum(((a + 1/2)*(1./(sigmas*(a + b + (mut - x).^2./sigmas^2).^(1/2)) - ((2*mut - 2.*x).*(mut - x))./(2*sigmas^3*(a + b + (mut - x).^2./sigmas^2).^(3/2))))./((mut - x)./(sigmas*(a + b + (mut - x).^2./sigmas^2).^(1/2)) - 1)...
    +((b + 1/2)*(1./(sigmas*(a + b + (mut - x).^2/sigmas^2).^(1/2)) - ((2*mut - 2.*x).*(mut - x))./(2*sigmas^3*(a + b + (mut - x).^2/sigmas^2).^(3/2))))./((mut - x)./(sigmas*(a + b + (mut - x).^2./sigmas^2).^(1/2)) + 1));

%dlogL_dsigma = -m/sigmas+sum(d1_dsigma+d2_dsigma);
dlogL_dsigma = -m/sigmas+sum(-((a + 1/2)*((mut - x)./(sigmas^2*(a + b + (mut - x).^2./sigmas^2).^(1/2)) - (mut - x).^3./(sigmas^4*(a + b + (mut - x).^2./sigmas^2).^(3/2))))./((mut - x)./(sigmas*(a + b + (mut - x).^2./sigmas^2).^(1/2)) - 1)...
    -((b + 1/2)*((mut - x)./(sigmas^2*(a + b + (mut - x).^2./sigmas^2).^(1/2)) - (mut - x).^3./(sigmas^4*(a + b + (mut - x).^2./sigmas^2).^(3/2))))./((mut - x)./(sigmas*(a + b + (mut - x).^2./sigmas^2).^(1/2)) + 1));

score_results = [dlogL_da dlogL_db dlogL_dmu dlogL_dsigma];

end