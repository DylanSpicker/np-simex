# Nonparametric SIMEX (NP-SIMEX)
This repository contains an implementation of the nonparametric simulation extrapolation (NP-SIMEX) methodology, as described in [Spicker, Wallace, Yi (2022)](https://arxiv.org/abs/2111.02863). The implementation is contained entirely in `np_simex.R` with the primary function of interest being `np_simex()`. 

This function takes:
* `replicates` - This contains either the error prone variable (as a numeric vector), or else a list of replicated versions of the error-prone vector. Note that if the `replicates` are not a list, then `U` must be provided.
* `estimator` - This is a callable option which takes as input a named argument (`X`) for the error-prone proxies, and then any additional arguments (passable to `np_simex()`), and returns the estimated value of interest. 
* `U` (defaults to `NULL`) - This is an optional parameter containing the set of observed errors. If a validation sample is to be used, this should be computed first and passed to the argument. If `replicates` is a list, then this argument will be ignored.
* `val_Xs` (defaults to `NULL`) - This is an optional parameter, used only when `het = TRUE` and when a validation sample is used. This contains the vector of error-prone variables from the validation sample. Note that this must correspond to `U` directly. 
* `M` (defaults to `10`) - The grid size for lambda. Will consider values of 1...M.
* `B` (defaults to `50`) - The number of re-sampled iterations to average over for each value of lambda.
* `parallel` (defaults to `TRUE`) - Whether the re-sampling should be parallelized or not. Implementations exist for both the `foreach` package and `parallel::mclapply`.
* `numCores` (defaults to `parallel::detectCores()/2`) - The number of cores to be used, if parallelization occurs.
* `est.variance` (defaults to `"none"`) - The method for estimating the variance. If it is provided as `"jackknife"` then the modified Jackknife procedure is used, otherwise no asymptotic variances are estimated.
* `parPackage` (defaults to `"foreach"`) - Which method for parallelizing is used. If this is anything other than `"foreach"`, and `parallel = TRUE`, then `parallel::mclapply` will be used. 
* `smoothed` (defaults to `FALSE`) - Should smoothed density estimators be used. If so samples are drawn from the KDE estimate of the distribution of `U`, rather than from the empirical distribution.
* `het` (defaults to `FALSE`) - Are the errors heterogenous, in that `U` and `X` are dependent. If so conditional KDEs are used in place of the empirical error distribution.
