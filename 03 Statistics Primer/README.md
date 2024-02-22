# Statistics Primer

### Random Numbers

- `draw.py` - draw a specified number of IID samples from a distribution
- `erfinv.py` - inverse error function approximation according to [Giles (2012)](https://people.maths.ox.ac.uk/gilesm/files/gems_erfinv.pdf)
- `bin.py` - bin IID random samples to estimate PDF
- `cdf.py` - rank IID random samples to estimate CDF and PDF
- `moments.py` - estimate moments and L-moments of sampled data
- `bootstrap.py` - re-draw samples (sampled from IID) randomly using bootstrap

### Multivariate Distributions

- `stocks/`
	- `*.csv` - market stock data from [NASDAQ](https://www.nasdaq.com/market-activity/quotes/historical)
	- `stocks.dat` - 13 tech stocks in `NumPy`-ready format
	- `parse.py` - parse `CSV` files and merge common data into `NumPy` array
	- `plot.py` - plot stock prices
- `scatter.py` - scatter plot of bivariate PDF, with a few variations
- `covariance.py` - compute the covariance of supplied data, and draw Gaussian sample with the same...
- `pca.py` - compute SVD of supplied data, and keep only a few components

### Random Processes

- `walk.py` - stochastic processes (aka random walk) demo
- `mcmc-?.py` - MCMC samplers using Metropolis-Hastings algorithm
- `field-?.py` - Gaussian random field generators using FFTs

