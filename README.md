# Bayesian-Synthetic-Likelihood

Approximate Bayesian computation (ABC, see Sisson and Fan (2011) for example) is a simulation-based method to approximate the posterior distribution in Bayesian inference where the likelihood function for a statistical model is difficult to compute in some way.  Likelihood evaluations are avoided by generating simulations from the model for many proposed parameter values and keeping those that generate simulated data 'close' to the observed data in some sense.  The distance measure is usually based on comparing some summarisation of the data.  However, in order to lose not too much information, it might be necessary to use a high-dimensional summary statistic.   Standard ABC methods are known not to scale efficiently with the dimension of the summary statistic.

Synthetic likelihood (Wood, 2010) is an alternative simulation-based approach that assumes that the model summary statistic follows a multivariate normal distribution with mean and covariance matrix that can depend on the parameter.  In the paper published in the Journal of Computational and Graphical Statistics (Price et al, 2016), we consider the synthetic likelihood within a Bayesian framework, which we refer to as BSL.  In short, we find the following:

- The BSL approximation does not seem to depend much on its tuning parameter, the number of simulations, n, to perform to estimate the mean and covariance matrix per proposed parameter value.
- We use an unbiased estimator of the multivariate normal density from the literature (Ghurye and Olkin, 1969) to construct a BSL algorithm that is theoretically unaffected by n when the model summary statistic does indeed have a multivariate normal distribution.
- BSL becomes increasingly more computationally efficient than ABC with an increase in the dimension of the summmary statistic.
- BSL can be a useful approach in interesting applications despite the multivariate normal assumption.  However, strong departures away from normality can lead to approximations of very poor quality.

In the paper we consider four examples.  Computer code is provided here for all examples, with one example per folder.  A README file is provided in each folder containing instructions pertaining to each of the examples.

## Key References
- Ghurye, S. G. and Olkin, I. (1969). Unbiased estimation of some multivariate probability densities and related functions. *The Annals of Mathematical Statistics*, 40(4):1261-1271.
- Price, L. F., Drovandi, C. C., Lee, A., and Nott, D. J. (2016). Bayesian synthetic likelihood. To appear in *Journal of Computational and Graphical Statistics*. https://eprints.qut.edu.au/92795/
- Sisson, S. A. and Fan, Y. (2011). *MCMC handbook*, chapter Likelihood-free Markov chain Monte Carlo, pages 313-335. Chapman & Hall.
- Wood, S. N. (2010). Statistical inference for noisy nonlinear ecological dynamic systems. *Nature*, 466:1102-1107.
