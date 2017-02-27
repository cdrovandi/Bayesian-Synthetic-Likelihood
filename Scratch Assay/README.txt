The code contained in this folder implements BSL and uBSL for the cell biology model example from Section 4.3 of the main paper.
For completeness, we have also included our ABC implementation.

For more details on model, refer to the main paper or to Johnston et al. (2014).

The 'run.m' file contains the code required to run BSL, uBSL and ABC for this example.

The function 'simulate_mex.c' uses the function 'mt19937ar.c' to get random numbers. This code can be downloaded from http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/CODES/mt19937ar.c .

To compile 'simulate_mex.c' on a Windows machine, use the line
mex simulate_mex.c mt19937ar.c

The files contained here are:

simulate_cell			-	performs a call to the C function 'simulate_mex'
simulate_mex			-	simulates binary matrices with cell locations
cell_summstats			-	computes the summary statistics (Hamming distances and number of cells at final time)
sl_log_like_ghuryeolkin	-	gets the unbiased estimator for the synthetic likelihood (when the summary statistics are MVN)
bayes_sl_cell			-	performs MCMC BSL
bayes_sl_cell_go		-	performs MCMC uBSL
abc_cell_nkern			-	performs ABC (taking advantage of 16 cores)
logsumexp				-	performs a stable calculation of log(sum(exp(x)))


These implementations can easily be adjusted to take advantage of parallel
computing architectures by uncommenting the relevant lines of code. Under 'run.m',
uncomment the line 'parpool(16)' and change 16 to the desired number of cores.
Also uncomment delete(gcp). In 'bayes_sl_cell', 'bayes_sl_cell_go' and
'abc_cell_nkern', replace the for loops with parallel for loops by uncommenting
the lines which begin with 'parfor' and commenting out the next lines which are
for standard for loops.