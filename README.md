# Bayesian Estimation of Change-Points in the Slope of Multivariate Time-Series

Assume that a temporal process is composed of contiguous segments with differing slopes and replicated noise-corrupted time series measurements are observed. The unknown mean of the data generating process is modelled as a piecewise linear function of time with an unknown number of change-points. The package infers the joint posterior distribution of the number and position of change-points as well as the unknown mean parameters per time-series by MCMC sampling. A-priori, the proposed model uses an overfitting number of mean parameters but, conditionally on a set of change-points, only a subset of them influences the likelihood. An exponentially decreasing prior distribution on the number of change-points gives rise to a posterior distribution concentrating on sparse representations of the underlying sequence, but also available is the Poisson distribution. See  [Papastamoulis et al (2019)](https://arxiv.org/abs/1709.06111) for a detailed presentation of the method.


## The software is available as an [R package on CRAN](https://CRAN.R-project.org/package=beast)


## Replication script is available [here](https://github.com/mqbssppe/growthPhaseMCMC/blob/master/replication_scripts/)

