The code contained in this folder implements BSL and uBSL for the Poisson example from Section 4.1 of the main paper.
For completeness, we have also included our ABC implementation.

For more details on model, refer to the main paper.

The 'run.m' file contains everything required to run BSL, uBSL and ABC for this example.

The files contained here are:

sl_log_like_ghuryeolkin	-	gets the unbiased estimator for the synthetic likelihood (when the summary statistics are MVN)
bayes_sl_simple			-	performs MCMC BSL
bayes_sl_simple_go		-	performs MCMC uBSL
abc_simple				-	performs ABC
