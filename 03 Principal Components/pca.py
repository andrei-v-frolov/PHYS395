#!/usr/bin/env python
# -*- coding: utf8 -*-

import csv, datetime
from dateutil import parser

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

# set of dates when we have proces for all three companies
common = set(apple.keys()) & set(google.keys()) & set(microsoft.keys())

for date in common:
	# date as seconds elapsed after epoch marker
	secs = int(parser.parse(date).strftime('%s'))
	print secs, apple[date], google[date], microsoft[date]
