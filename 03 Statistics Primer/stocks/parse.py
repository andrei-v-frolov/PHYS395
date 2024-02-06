#!/usr/bin/env python
# -*- coding: utf8 -*-
#######################################################################
# parse CSV stock prices and merge common data into NumPy array
# https://www.nasdaq.com/market-activity/quotes/historical
#######################################################################

from functools import reduce
import csv, datetime, locale
from dateutil import parser

# currency formatting for US locale
locale.setlocale(locale.LC_ALL, 'en_US')
currency = locale.localeconv()['currency_symbol']

# parse CSV stock data, returning a dictionary in {epoch: price} format
def stock(filename):
	data = {}
	
	with open(filename,'rt') as csvfile:
		for row in csv.reader(csvfile):
			if row[0] == 'Date': continue
			
			date = parser.parse(row[0])
			epoch = int(date.timestamp()/(24*60*60))
			price = locale.atof(row[1].replace(currency,''))
			
			data[epoch] = price
	
	return data

#######################################################################

# NASDAQ stock symbols we grabbed data for...
stocks = ["AAPL", "AMD", "AMZN", "CSCO", "GOOGL", "INTC",
"META", "MSFT", "NFLX", "NVDA", "QCOM", "SBUX", "TSLA"]

data = []

# load up stock prices data from Excel CSV files
for s in stocks:
	data.append(stock(s+'.csv'))

# set of dates when we have prices for all companies
common = reduce(lambda x, y: x & y, (set(stock.keys()) for stock in data))

#######################################################################

for epoch in sorted(common):
	print(*[stock[epoch] for stock in data])
