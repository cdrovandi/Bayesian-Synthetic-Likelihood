The code contained in this folder implements BSL and uBSL for the Ricker model example from Section 4.2 of the main paper.
For completeness, we have also included our ABC implementation.
The Ricker model was presented in Wood (2010), and the model simulation and summary statistics have been replicated here.

For more details on model, refer to the main paper or to Wood (2010).

The 'run.m' file contains everything required to run BSL, uBSL and ABC for this example.

The files contained here are:

simulate_ricker			-	simulates from the Ricker model
ricker_summstats		-	computes the same summary statistics as in Wood (2010)
sl_log_like_ghuryeolkin	-	gets the unbiased estimator for the synthetic likelihood (when the summary statistics are MVN)
bayes_sl_ricker_wood	-	performs MCMC BSL
bayes_sl_ricker_wood_go	-	performs MCMC uBSL
abc_ricker_wood			-	performs ABC

The BSL and uBSL implementations can easily be adjusted to take advantage of parallel
computing architectures by uncommenting the relevant lines of code. Under 'run.m',
uncomment the line 'parpool(16)' and change 16 to the desired number of cores.
Also uncomment delete(gcp). In 'bayes_sl_ricker_wood' and 'bayes_sl_ricker_wood_go',
replace the for loops with parallel for loops by uncommenting the lines which begin with
'parfor' and commenting out the next lines which are for standard for loops.