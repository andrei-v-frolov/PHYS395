#!/usr/bin/env python
# -*- coding: utf8 -*-

import csv, datetime
from dateutil import parser

import numpy as np
import scipy.linalg as linalg

def stock(filename):
	"""Parse CSV stock data, return a dictionary in {date: price} format"""
	
	data = {}
	
	with open(filename,'rb') as csvfile:
		for row in csv.reader(csvfile):
			if row[0] == 'Date': continue
			
			date = row[0]
			price = float(row[1])
			
			data[date] = price
	
	return data

# load up stock prices data from Excel CSV files
apple = stock('APPL.csv')
google = stock('GOOG.csv')
microsoft = stock('MSFT.csv')

# set of dates when we have prices for all three companies
common = set(apple.keys()) & set(google.keys()) & set(microsoft.keys())

data = []

for date in sorted(common):
	# date as seconds elapsed after epoch marker
	secs = int(parser.parse(date).strftime('%s'))
	#print secs, apple[date], google[date], microsoft[date]
	data.append([apple[date], google[date], microsoft[date]])

# compute singular value decomposition of the stock data
U, S, V = linalg.svd(data, compute_uv=True, overwrite_a=False)

n = U.shape[0]; m = V.shape[0]; S[-1] = 0.0

# recombine the matrix with some of the eigenvalues zeroed
A = np.dot(U, np.dot(linalg.diagsvd(S,n,m), V))

for row in A:
	print row[0], row[1], row[2]
