The code contained in this folder implements BSL and uBSL for the g-and-k example from the supplementary materials.
For completeness, we have also included our ABC implementation.

For more details on model, refer to the Appendix E of the supplementary materials.

The 'run.m' file contains everything required to run BSL, uBSL and ABC for this example.

The files contained here are:

simulate_gandk			-	simulates from the g-and-k model
fun_gandk				-	contains the g-and-k quantile function for simulating from the model
Scores					-	computes the score (the summary statistic for this example)
sl_log_like_ghuryeolkin	-	gets the unbiased estimator for the synthetic likelihood (when the summary statistics are MVN)
bayes_sl_gandk			-	performs MCMC BSL
bayes_sl_gandk_go		-	performs MCMC uBSL
abc_gandk				-	performs ABC

The BSL and uBSL implementations can easily be adjusted to take advantage of parallel
computing architectures by uncommenting the relevant lines of code. Under 'run.m',
uncomment the line 'parpool(16)' and change 16 to the desired number of cores.
Also uncomment delete(gcp). In 'bayes_sl_gandk' and 'bayes_sl_gandk_go', replace
the for loops with parallel for loops by uncommenting the lines which begin with
'parfor' and commenting out the next lines which are for standard for loops.