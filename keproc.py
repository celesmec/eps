# -*- coding: utf-8 -*-
##################################################
# Copyright (C) 2013-2018 by Alexander Perminov
# perminov12@yandex.ru
#
##################################################
# THIS IS KEPLERIAN PROCESSOR
##################################################
def besselJ(n, x, order):
	"""
	Bessel functions.
	"""
	from pyranha.math import factorial
	from fractions import Fraction as F
	temp = 0
	for m in range(0, int((order-n)/2) + 1):
		temp = temp + F((-1)**m, factorial(m)*factorial(m+n)) * (x/2)**(2*m+n)
	return temp
##################################################
def binomial_exp(x, y, r, order):
	"""
	Binomial exponentiations.
	"""
	from pyranha.math import binomial
	temp = 0
	for k in range(0, order + 1):
		temp = temp + binomial(r, k) * x**(r-k) * y**k
	return temp
##################################################
def LC(n, k, alpha, imax = 200):
	"""
	Laplace coefficients.
	"""
	temp = 1
	if k == 0:
		for i in range(1, imax+1):
			C = 1
			for j in range(1, i+1):
				C = C * ((n+2*(j-1))/(2.*j))
			temp = temp + C**2*alpha**(2*i)
		temp = 2*temp
	if k != 0:
		for i in range(1, imax+1):
			C = 1
			for j in range(1, i+1):
				C = C * ((n+2*(j-1))/(2.*j))*((n+2*(j-1)+2*k)/(2.*j+2.*k))
			temp = temp + C*alpha**(2*i)
		C = 1
		for j in range(1, k+1):
			C = C * ((n+2*(j-1))/(2.*j))
		temp = 2*C*alpha**k*temp
	return temp
##################################################
def legendrePn(n, x):
	"""
	Legendre polynomials.
	"""
	from pyranha.math import factorial
	from fractions import Fraction as F
	temp = 0
	for k in range(0, int(n/2) + 1):
		temp = temp + (-1)**k * F(factorial(2*(n - k)), 2**n * factorial(k) * factorial(n - k) * factorial(n - 2*k)) * x**(n - 2*k)
	return temp
##################################################
def r_a(e, M, order):
	"""
	Celestial mechanics classical expansion of r/a.
	"""
	from pyranha.math import cos
	from fractions import Fraction as F
	temp = 1 + F(1, 2) * e**2
	for k in range(1, order + 1):
	 	temp = temp - e * F(1, k) * (besselJ(k-1, k*e, order-1) - besselJ(k+1, k*e, order-1)) * cos(k*M)
	return temp
##################################################
def a_r(e, M, order):
	"""
	Celestial mechanics classical expansion of a/r.
	"""
	from pyranha.math import cos
	from fractions import Fraction as F
	temp = 1
	for k in range(1, order + 1):
		temp = temp + F(2, 1) * besselJ(k, k*e, order) * cos(k*M)
	return temp
##################################################
def cos_E(e, M, order):
	"""
	Celestial mechanics classical expansion of cosE.
	"""
	from pyranha.math import cos
	from fractions import Fraction as F
	temp = F(-1, 2) * e
	for k in range(1, order + 2):
		temp = temp + F(1, k) * (besselJ(k-1, k*e, order) - besselJ(k+1, k*e, order)) * cos(k*M)
	return temp
##################################################
def sin_E(e, M, order):
	"""
	Celestial mechanics classical expansion of sinE.
	"""
	from pyranha.math import sin
	from fractions import Fraction as F
	temp = 0
	for k in range(1, order + 2):
		temp = temp + F(1, k) * (besselJ(k-1, k*e, order) + besselJ(k+1, k*e, order)) * sin(k*M)
	return temp
##################################################
def cos_f(e, M, order):
	"""
	Celestial mechanics classical expansion of cosf.
	"""
	from pyranha.math import cos
	temp = 0
	for k in range(1, order + 2):
		temp = temp + (besselJ(k-1, k*e, order) + besselJ(k+1, k*e, order)) * cos(k*M)
	temp = temp * (1 - e**2) - e
	return temp.transform(lambda t: (t[0].filter(lambda u: u[1].degree(['e']) <= order), t[1]))
##################################################
def sin_f(e, M, order):
	"""
	Celestial mechanics classical expansion of sinf.
	"""
	from pyranha.math import sin
	from fractions import Fraction as F
	temp = 0
	for k in range(1, order + 2):
		temp = temp + (besselJ(k-1, k*e, order) - besselJ(k+1, k*e, order)) * sin(k*M)
	temp = temp * binomial_exp(1, -e**2, F(1,2), int(order/2))
	return temp.transform(lambda t: (t[0].filter(lambda u: u[1].degree(['e']) <= order), t[1]))
##################################################
def EE(e, M, order):
	"""
	Celestial mechanics classical expansion of E.
	"""
	from pyranha.math import sin
	from fractions import Fraction as F
	temp = M
	for k in range(1, order + 1):
		temp = temp + F(2, k) * besselJ(k, k*e, order) * sin(k*M)
	return temp
##################################################

