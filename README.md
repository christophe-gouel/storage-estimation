## Dependencies

* MATLAB R2013b
* RECS v0.7
* R, with the package knitr

## Files organization

## First guess of the parameters for the estimation

We propose a simple heuristic, relying on a method of moments, to determine a
first guess for the parameters to be estimated.

Considering the inverse demand function $P=a+b D$, we take its first two
moments:

$$\mathrm{E}(P)=a+b\mathrm{E}(D),$$

$$\sigma_P^2=b^2\sigma_D^2.$$

The mean and variance of price, $\mathrm{E}(P)$ and $\sigma_P^2$, can be
estimated from the sample. Mean demand, $\mathrm{E}(D)$ is equal to mean supply,
as in expectations stock accumulations are compensated by stock releases. Given
that the production shocks have a zero mean, we have $\mathrm{E}(D)=0$. The
variance of demand is, because of the smoothing effect of storage, inferior to
the variance of production shocks, so inferior to 1, in our calibration.

Storage costs $k$ can be assumed to be a share of mean price.

Following these considerations, providing a first guess for the parameters
amounts to guess two dimensionless values: the proportion of reduction in the
variance of demand with respect to the variance of production and storage costs
expressed as percentage of mean price:

$$\left\{a=\mathrm{E}(P),b=-\frac{\sigma_P}{\sigma_D},k=sh_{k/P}\mathrm{E}(P)\right\}.$$

From previous work estimating the storage model, we know that to reproduce the
observed price autocorrelation, one needs that storage plays a large role in
dynamics. So we assume $\sigma_D=0.3$ and $sh_{k/P}=0.02$.
