<!-- README.md is generated from README.Rmd. Please edit that file -->
synthate package: Synthetic estimators for the average treatment effect
=======================================================================

This package comprises a suite of tools to simulate various data-generating processes. Its purpose is to check the operating characteristics of various causal estimators and to optimally combine them.

Installation
------------

The `synthate` package can currently be obtained from GitLab: <https://code.rand.org/agniel-projects/personal-packages/synthate> or from GitHub at <https://github.com/denisagniel/synthate>. It may be downloaded by running

``` r
devtools::install_github('denisagniel/synthate')
```

Overview
--------

There are three main features of `synthate`:

-   Simulating data
-   Causal estimation
-   Combining estimators

### Simulating data

There are currently seven DGPs derived from six different papers to choose from:

| Short-hand | Paper                                            |
|------------|--------------------------------------------------|
| `ks`       | Kang and Schafer (Kang and Schafer (2007))       |
| `ld`       | Lunceford and Davidian (Lunceford et al. (2004)) |
| `iw`       | Waernbaum (Waernbaum (2012))                     |
| `fi`       | Fan and Imai (Fan (2016))                        |
| `pa`       | Austin (Austin (2009))                           |
| `ls`       | Leacy and Stuart (Leacy and Stuart (2014))       |
| `ik`       | Iacus and King (Iacus, King, and Porro (2009))   |

References
----------

Austin, Peter C. 2009. “Some methods of propensity-score matching had superior performance to others: results of an empirical investigation and monte carlo simulations.” *Biometrical Journal* 51 (1): 171–84. doi:[10.1002/bimj.200810488](https://doi.org/10.1002/bimj.200810488).

Fan, Jianqing. 2016. “Improving Covariate Balancing Propensity Score : A Doubly Robust and Efficient Approach ∗,” 1–47.

Iacus, Stefano M., Gary King, and Giuseppe Porro. 2009. “&lt;b&gt;cem&lt;/b&gt; : Software for Coarsened Exact Matching.” *Journal of Statistical Software* 30 (9). doi:[10.18637/jss.v030.i09](https://doi.org/10.18637/jss.v030.i09).

Kang, Joseph D. Y., and Joseph L. Schafer. 2007. “Demystifying Double Robustness: A Comparison of Alternative Strategies for Estimating a Population Mean from Incomplete Data.” *Statistical Science* 22 (4): 523–39. doi:[10.1214/07-STS227](https://doi.org/10.1214/07-STS227).

Leacy, Finbarr P., and Elizabeth A. Stuart. 2014. “On the joint use of propensity and prognostic scores in estimation of the ATT: A simulation study.” *Statistics in Medicine* 33 (20): 161–69. doi:[10.3851/IMP2701.Changes](https://doi.org/10.3851/IMP2701.Changes).

Lunceford, Jared K, Jared K Lunceford, Marie Davidian, and Marie Davidian. 2004. “Stratification and weighting via the propensity score in estimation of causal treatment e ects: a comparative study.” *Statistics in Medicine* 2960 (19): 2937–60. doi:[10.1002/sim.1903](https://doi.org/10.1002/sim.1903).

Waernbaum, Ingeborg. 2012. “Model misspecification and robustness in causal inference: Comparing matching with doubly robust estimation.” *Statistics in Medicine* 31 (15): 1572–81. doi:[10.1002/sim.4496](https://doi.org/10.1002/sim.4496).
