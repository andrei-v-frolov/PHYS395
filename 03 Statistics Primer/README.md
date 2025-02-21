# Statistics Primer

### Random Numbers

- `rand.py` - random number generator demo
- `draw.py` - draw a specified number of IID samples from a distribution
- `erfinv.py` - inverse error function approximation by [Giles (2012)](https://people.maths.ox.ac.uk/gilesm/files/gems_erfinv.pdf)
- `bin.py` - bin IID random samples to estimate PDF
- `cdf.py` - sort IID random samples to estimate CDF

### Multivariate Distributions

- `stocks/`
	- `*.csv` - market stock data from [NASDAQ](https://www.nasdaq.com/market-activity/quotes/historical)
	- `stocks.dat` - 13 tech stocks in `NumPy`-ready format
	- `parse.py` - parse `CSV` files and merge common data into `NumPy` array
	- `plot.py` - plot stock prices
