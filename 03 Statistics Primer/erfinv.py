# inverse error function approximation according to Giles (2012)

from numpy import log, sqrt

# evaluate polynomial in Horner's form
def horner(coeff,x):
	p = 0.0
	for c in coeff:
		p = c + p*x
	return p

# erfinv approximation coefficients
C_ERF1 = [2.81022636e-08, 3.43273939e-07, -3.5233877e-06, -4.39150654e-06, 0.00021858087, -0.00125372503, -0.00417768164, 0.246640727, 1.50140941]
C_ERF2 = [-0.000200214257, 0.000100950558, 0.00134934322, -0.00367342844, 0.00573950773, -0.0076224613, 0.00943887047, 1.00167406, 2.83297682]

# inverse error function approximation, good to single precision
def erfinv(x):
	w = log((1.0-x)*(1.0+x))
	p = horner(C_ERF1,w-2.5) if w < 5.0 else horner(C_ERF2,sqrt(w)-3.0)
	return p*x
