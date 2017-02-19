These files provide a friendly way of computing the P-value Prediction Interval for future P-values given an observed P-value and sample size for the hypothesis test. The test statistic follows either the standard normal distribution or a student's t-distribution with N-1 degrees of freedome, while the prior distribution on the test statistic's true mean follows a discretized normal distribution according to the parameters specified by the user. The user should edit the entries in `computerInterval.R` as necessary.

In order to compute these intervals, the values assigned to `p.obt` and `N` in the file `computeInterval.R` should be changed to the observed p-value and sample size, respectively. After these and other edits that include specifying the parameters for the prior distribution and type of test statistic, the entire script, `computeInterval.R`, should be run. The data.frame object `summaryInfo` stores the following two Prediction Intervals: 

1. Prediction interval according to Vsevolozhskaya et al. ("Discrete Bayesian") [1]
2. P-interval according to Cumming [2]; Lazzeroni et al. ("P-interval") [3]

Additionally, this information is output to "summaryInfo.csv".

The file "helperScript.R" provides the helper functions for computing this information. 

References: 
1. Vsevolozhskaya OA, Ruiz G, Zaykin DV. Assessment of P-value variability in the current replicability crisis." arXiv preprint arXiv:1609.01664, 2016. 
2. Cumming G. Replication and p intervals: p values predict the future only vaguely, but confidence intervals do much better. Perspect Psychol Sci. 2008;3(4):286–300.
3. Lazzeroni L, Lu Y, Belitskaya-Levy I. P-values in genomics: apparent precision masks high uncertainty. Molecular Psychiatry. 2014;19(12):1336–1340.
