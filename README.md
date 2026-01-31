# Tail risk in the tail: Estimating high quantiles when a related variable is extreme

This repository contains the code for replication of the results in the paper "Tail risk in the tail: Estimating high quantiles when a related variable is extreme".


## Files
- The `code` directory: contains the code and data to reproduce the simulation studies in Sections 3.3 and application studied in Section 4.
  - `functions.R`: The functions used to simulation and real data analysis.
  - `code_simulation.R`: The main code to generate simulation results in Section 3.3.
  - `code_application.R`: The main code to conduct dynamic forecasting in Section 4.
- The `manuscript` directory: contains the source files for the manuscript.
- The `data` directory: contains original data (Snp500_dataset.csv) used in Section 4.


## Quick start — which functions to call

Recommended order to run the pipeline (fast path):

1. Pre-filter returns (GARCH) to obtain standardized residuals:
   - filtering(...)

2. Estimate system tail index and VaR:
   - sys_tail <- tail_estimate(...)
   - VaR_sys <- Qvar_esti(...)

3. Estimate CoVaR:
   - Proposed method: CoVaR_estimate(...)
   - Fully parametric method: CoVaR_FP(...)
   - Nolde & Zhang EVT: CoVaR_NZ(...)

4. Backtesting and comparison:
   - average_score(...)
   - unconditional_test(...)
   - comparative.test(...)

Example:
```r
# 1) Pre-filter returns -> obtain standardized residuals:
res_info = filtering(dat = my_returns, forecast = FALSE)
residuals = res_info  # standardized residuals

# 2) Estimate tail index and VaR for market (system):
gamma_sys = tail_estimate(dat = residuals[,2], k = round(0.1 * nrow(residuals)))
VaR_sys = Qvar_esti(dat = residuals[,2], gamma = gamma_sys, p = 0.05, k = round(0.01 * nrow(residuals)))

# 3) Use CoVaR_estimate for proposed method:
covar_est = CoVaR_estimate(Data = residuals, group = c("t","alog"), m = c(200,200),k = c(300,100), p = c(0.01, 0.05))
```


## Quick reference - functions.R

This file documents the main user-facing functions and suggested usage for the code in `code/functions.R`.

It contains the R comment block recommended to help users distinguish high-level API functions from helpers.

### User-facing / Main functions
```r
# These are the high-level functions most users will call directly:
#
#  - CoVaR_estimate(Data, par_hat=NULL, group=NULL, m=NULL, k, p)
#      Main proposed method to estimate CoVaR. If par_hat is NULL,
#      M_estimate() will be run for each distribution in `group`.
#      Inputs:
#        Data: two-column data.frame/matrix (institution, market)
#        group: character vector of parametric TD families (e.g. "t", "log", "alog", "hr")
#        m: sample fractions corresponding to `group` (required if par_hat is NULL)
#        k: vector of k values for tail estimation and VaR estimation (e.g. c(k_tail, k_var))
#        p: two-element risk vector (VaR level, CoVaR level)
#      Returns: numeric vector of CoVaR estimates (one per `group` or per par_hat)
#
#  - CoVaR_FP(Data, fit_par = NULL, VaR = NULL, level)
#      Parametric Skew-t based CoVaR (Girardi 2013 style). Fits multivariate skew-t
#      (mst.mle) then finds CoVaR via CoVaR_Skew objective.
#
#  - CoVaR_NZ(Data, VaR = NULL, k = NULL, level)
#      Nolde & Zhang EVT-based CoVaR: uses New_cov_update / New_update to fit EVT parameters
#      and then solves for CoVaR.
#
#  - filtering(dat, forecast = TRUE, n_out = 1, model = "sGARCH")
#      Pre-filtering workflow (GARCH fit + standardized residuals). Use this before tail
#      estimation / VaR forecasting when working with time series returns.
#
#  - average_score(forecast_CoVaR, forecast_VaR, test, level)
#  - unconditional_test(forecast_CoVaR, forecast_VaR, test, level)
#  - comparative.test(VaR_fr, CoVaR_fr, CoVaR_fr_r, test_ins, test_sys, p)
#      Backtesting / scoring utilities for evaluating forecasts and comparing methods.
#
#  - New_cov_update(sim_y, method="combine", prob = 0.85)
#  - New_update(sim_y, method="combine", prob = 0.85)
#      EVT estimation routines (New_cov_update is the most general).
#
#  - M_estimate(X, family, m, start = NULL)
#      M-estimation for tail-dependence parameters (used internally by CoVaR_estimate
#      but can be called standalone to get estimates).
#
#  - tail_estimate(dat, k), Qvar_esti(dat, gamma, p, k), eta_estimate(par_hat, p, family)
#      Tail index and VaR helpers: compute Hill tail-index, nonparametric VaR estimator,
#      and the eta_star solving routine used in the CoVaR formula.
#
#  - real_compute / data_generate
#      Utilities used for simulation studies and computing "true" values.
```


### Helper / internal functions (usually not called directly)
```r
#  - tail_dependence, generate_gfun, nonpar_tail, optim_fun
#      Core math used by M_estimate.
#
#  - neglogfn, neglogfn1, q_x, q_x_cov, func, func_cov, K, Thrm4_1
#      Numerical likelihoods / integrals used in EVT estimation (New_update / New_cov_update).
#
#  - Hill_small, f_theta, f_theta_plot
#      Small helpers for Hill estimator and plotting spectral measures.
#
#  - mst.fit, mst.mle, mst.dev, mst.dev.grad, solvePD, log.pt, num.deriv2, etc.
#      Skew-t fitting and internal likelihood machinery (copied/derived from 'sn' internals).
```

### Notes
- Many of the functions accept k/m values as either fractions (<1) or absolute counts.Passing fractional values retains backward compatibility; the functions convert to counts.

