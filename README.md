# mutualism-diversification

## HiSSE_ML

All files correspond to running HiSSE analyses in R using Maximum Likelihood methods.

``HiSSE_ML.R`` 7 types of HiSSE models

``corHMM_and_simmap.R`` Uses output of HiSSE models as input for corHMM and subsequently takes the corHMM model to run a simmap. Also includes code to calulcate the time spent in each state, to be used to calculate weighted diversification rates.

``corHMM_figure.R`` plots for weighted diversification rates from simmaps

``diversification_rates_HiSSE_ML.R`` Converts turnover and extinction fraction rates from HiSSE output into speciation and extinction using equations from Caetano et al. 2018 (https://doi.org/10.1111/evo.13602), to be used to calculate weighted diversification rates.


## HiSSE_RevBayes

All files correspond to running HiSSE analyses in RevBayes (and downstream analyses for RevBayes output conducted in R).

``HiSSE_RB_constrained.Rev`` code to run HiSSE in RevBayes with extinction constrained to 0.8(speciation)

``HiSSE_RB_estimated.Rev`` code to run HiSSE in RevBayes with extinction estimated

``RevBayes_figures.R``plots for weighted diversification rates from RevBayes

``RevBayes_weighted_rates.R``used to calculate weighted diversification rates from RevBayes

## Other Files

``sister_contrasts.R`` runs sister contrasts
