# -*- coding: utf-8 -*-
##################################################
# Copyright (C) 2013-2018 by Alexander Perminov
# perminov12@yandex.ru
#
##################################################
# THIS IS ALGORITHM FOR THE HAMILTONIAN EXPANSION
##################################################
import os, time
from eps.config import *
##################################################
# THIS IS A SET OF BASE EXPANSION ALGORITHMS
##################################################
def kepler_trig_expand(s):
	"""
	Transform trigonometric part of series s [sin(nM) and cos(nM)]
	to sum of items like cos(M)^k1 * sin(M)^k2, where n=k1+k2.
	- s is series of e and M,
	- n is number of terms.
	"""
	from math import ceil
	from pyranha.math import binomial, cos, sin, degree
	temp, e, M, ecosM, esinM = 0, epst('e'), epst('M'), epst('ecosM'), epst('esinM')
	n = degree(s).numerator
	s_list = s.list
	for i in range(len(s_list)):
		for j in range(n+1):
			trig_cos, trig_sin = 0, 0
			# calculate cos(nM):
			if s_list[i][1] == cos(j*M):
				for k in range(int(ceil(j/2))+1):
					trig_cos = trig_cos + (-1)**k * binomial(j, 2*k) * ecosM**(j-2*k) * esinM**(2*k) * e**-j
				temp = temp + s_list[i][0] * trig_cos
			# calculate sin(nM):
			if s_list[i][1] == sin(j*M):
				for k in range(int(ceil((j-1)/2))+1):
					trig_sin = trig_sin + (-1)**k * binomial(j, 2*k+1) * ecosM**(j-2*k-1) * esinM**(2*k+1) * e**-j
				temp = temp + s_list[i][0] * trig_sin
	if type(temp) != type(1): temp = temp.trim()
	return temp
##################################################
def xyz(p, n):
	"""
	Calculation expansions for vector of x, y, z-coordinates.
	- n is maximum degree of variables.
	"""
	from pyranha.math import cos, sin
	from eps.keproc import binomial_exp, cos_E, sin_E
	deg_max, e_max = n/2 + 1, n
	MM, LL, xx, yy, uu, vv, qq = [epst(i+str(p)) for i in ('m', 'L', 'x', 'y', 'u', 'v', 'q')]
	K0, ee, MA, ecosM, esinM = epst('K0'), epst('e'), epst('M'), epst('ecosM'), epst('esinM')
	S1 = epst('s'+str(p-1)) if p != 1 else F(1, 1)
	S2 = epst('s'+str(p))
	##################################################
	temp1 = binomial_exp(epst(1), -epst('eps'), -1, deg_max)
	temp2 = binomial_exp(epst(1), -epst('eps'), F(1,2), deg_max)
	temp3 = F(1,4)*LL**-1 * temp1.subs('eps', F(1,2)*LL**-1 *    (xx**2 + yy**2))
	temp5 = LL**F(-1,2)   * temp2.subs('eps', F(1,4)*LL**-1 *    (xx**2 + yy**2))
	temp4 = LL**F(-1,2)   * temp1.subs('eps', F(1,2)*LL**-1 *    (xx**2 + yy**2)) * \
	                        temp2.subs('eps', F(1,4)*LL**-1 * (2*(xx**2 + yy**2) + (uu**2 + vv**2)))
	cos2I  =  1 - (uu**2 + vv**2) * temp3
	sin2Icos2Om = (uu**2 - vv**2) * temp3
	sin2Isin2Om =    -2*(uu * vv) * temp3
	sinIsinOm   = -vv * temp4
	sinIcosOm   =  uu * temp4
	ecospi      =  xx * temp5
	esinpi      = -yy * temp5
	##################################################
	temp = binomial_exp(epst(1), -ee**2, F(1,2), e_max)
	sqrt1 = F(1,2)*(1 + temp)
	sqrt2 = F(1,2)*(1 - temp)
	cosEM1 = cos_E(ee, MA, e_max) * cos(MA) + sin_E(ee, MA, e_max) * sin(MA)
	cosEM2 = cos_E(ee, MA, e_max) * cos(MA) - sin_E(ee, MA, e_max) * sin(MA)
	sinEM1 = sin_E(ee, MA, e_max) * cos(MA) - cos_E(ee, MA, e_max) * sin(MA)
	sinEM2 = sin_E(ee, MA, e_max) * cos(MA) + cos_E(ee, MA, e_max) * sin(MA)
	cosEM1 = (sqrt1 * cosEM1).truncate_degree(e_max, ['e'])
	cosEM2 = (sqrt2 * cosEM2).truncate_degree(e_max, ['e'])
	sinEM1 = (sqrt1 * sinEM1).truncate_degree(e_max, ['e'])
	sinEM2 = (sqrt2 * sinEM2).truncate_degree(e_max, ['e'])
	cosEM1 = kepler_trig_expand(cosEM1)
	cosEM2 = kepler_trig_expand(cosEM2)
	sinEM1 = kepler_trig_expand(sinEM1)
	sinEM2 = kepler_trig_expand(sinEM2)
	X_a = -ecosM + cosEM1 + cosEM2
	Y_a =  esinM + sinEM1 - sinEM2
	X_a = X_a.subs('e', ee**F(1,2))
	Y_a = Y_a.subs('e', ee**F(1,2))
	subs1 = ecospi * cos(qq) + esinpi * sin(qq)
	subs2 = ecospi * sin(qq) - esinpi * cos(qq)
	subs3 = LL**-1 * (xx**2 + yy**2) * (1 - F(1,4)*LL**-1 * (xx**2 + yy**2))
	##################################################
	pt.set_auto_truncate_degree(F(n), ['x'+str(p), 'y'+str(p), 'u'+str(p), 'v'+str(p)])
	X_a = X_a.subs('ecosM', subs1)
	X_a = X_a.subs('esinM', subs2)
	X_a = X_a.subs('e', subs3)
	Y_a = Y_a.subs('ecosM', subs1)
	Y_a = Y_a.subs('esinM', subs2)
	Y_a = Y_a.subs('e', subs3)
	aa = (K0*MM**2*S1)**-1*S2*LL**2
	xa = aa * (X_a * (cos2I * cos(qq) + sin2Icos2Om * cos(qq) + sin2Isin2Om * sin(qq)) -\
			   Y_a * (cos2I * sin(qq) + sin2Icos2Om * sin(qq) - sin2Isin2Om * cos(qq)))
	ya = aa * (X_a * (cos2I * sin(qq) - sin2Icos2Om * sin(qq) + sin2Isin2Om * cos(qq)) +\
			   Y_a * (cos2I * cos(qq) - sin2Icos2Om * cos(qq) - sin2Isin2Om * sin(qq)))
	za = aa * (X_a * (sinIcosOm * sin(qq) - sinIsinOm * cos(qq)) +\
			   Y_a * (sinIcosOm * cos(qq) + sinIsinOm * sin(qq)))
	pt.unset_auto_truncate_degree()
	return [xa.trim(), ya.trim(), za.trim()]
##################################################
def one_r(p, n, deg):
	"""
	Calculation expansions for 1/r and them powers.
	- n is maximum degree of variables.
	"""
	from pyranha.math import cos, sin
	from eps.keproc import binomial_exp, a_r
	deg_max, e_max = n/2 + 1, n
	MM, LL, xx, yy, uu, vv, qq = [epst(i+str(p)) for i in ('m', 'L', 'x', 'y', 'u', 'v', 'q')]
	K0, ee, MA = epst('K0'), epst('e'), epst('M')
	S1 = epst('s'+str(p-1)) if p != 1 else F(1, 1)
	S2 = epst('s'+str(p))
	##################################################
	temp2  = binomial_exp(epst(1), -epst('eps'), F(1,2), deg_max)
	temp5  = LL**F(-1,2) * temp2.subs('eps', F(1,4)*LL**-1 * (xx**2 + yy**2))
	ecospi =  xx * temp5
	esinpi = -yy * temp5
	##################################################
	ar = a_r(ee, MA, e_max)
	ar = kepler_trig_expand(ar)
	ar = ar.subs('e', ee**F(1,2))
	subs1 = ecospi * cos(qq) + esinpi * sin(qq)
	subs2 = ecospi * sin(qq) - esinpi * cos(qq)
	subs3 = LL**-1 * (xx**2 + yy**2) * (1 - F(1,4)*LL**-1 * (xx**2 + yy**2))
	##################################################
	pt.set_auto_truncate_degree(F(n), ['x'+str(p), 'y'+str(p)])
	ar = ar.subs('ecosM', subs1)
	ar = ar.subs('esinM', subs2)
	ar = ar.subs('e', subs3)
	ar = ar.trim()
	aa = (K0*MM**2*S1)**-1*S2*LL**2
	ar_deg = (ar*aa**-1)**deg
	pt.unset_auto_truncate_degree()
	return ar_deg.trim()
##################################################
def rrr(q, n, deg):
	"""
	Calculation expansions for r and them powers.
	- n is maximum degree of variables.
	"""
	from pyranha.math import cos, sin
	from eps.keproc import binomial_exp, r_a
	deg_max, e_max = n/2 + 1, n
	MM, LL, xx, yy, uu, vv, qq = [epst(i+str(q)) for i in ('m', 'L', 'x', 'y', 'u', 'v', 'q')]
	K0, ee, MA = epst('K0'), epst('e'), epst('M')
	S1 = epst('s'+str(q-1)) if q != 1 else F(1, 1)
	S2 = epst('s'+str(q))
	##################################################
	temp2  = binomial_exp(epst(1), -epst('eps'), F(1,2), deg_max)
	temp5  = LL**F(-1,2) * temp2.subs('eps', F(1,4)*LL**-1 * (xx**2 + yy**2))
	ecospi =  xx * temp5
	esinpi = -yy * temp5
	##################################################
	ra = r_a(ee, MA, e_max)
	ra = kepler_trig_expand(ra)
	ra = ra.subs('e', ee**F(1,2))
	subs1 = ecospi * cos(qq) + esinpi * sin(qq)
	subs2 = ecospi * sin(qq) - esinpi * cos(qq)
	subs3 = LL**-1 * (xx**2 + yy**2) * (1 - F(1,4)*LL**-1 * (xx**2 + yy**2))
	##################################################
	pt.set_auto_truncate_degree(F(n), ['x'+str(q), 'y'+str(q)])
	ra = ra.subs('ecosM', subs1)
	ra = ra.subs('esinM', subs2)
	ra = ra.subs('e', subs3)
	ra = ra.trim()
	aa = (K0*MM**2*S1)**-1*S2*LL**2
	ra_deg = (ra*aa)**deg
	pt.unset_auto_truncate_degree()
	return ra_deg.trim()
##################################################
def cosine(p, q, n, deg):
	"""
	Calculation expansions for cosine between two radius vectors.
	- n is maximum degree of variables.
	"""
	one_ra = one_r(p, n, deg)
	one_rb = one_r(q, n, deg)
	xa, ya, za = xyz(p, n)
	xb, yb, zb = xyz(q, n)
	pt.set_auto_truncate_degree(F(n), ['x'+str(p), 'y'+str(p), 'u'+str(p), 'v'+str(p), 'x'+str(q), 'y'+str(q), 'u'+str(q), 'v'+str(q)])
	cosine = (xa*xb + ya*yb + za*zb)**deg * one_ra*one_rb
	pt.unset_auto_truncate_degree()
	return cosine
##################################################
def one_delta(p, q, n, c_max, deg, isAver, isCalc, isPrint):
	"""
	Construction of the expansion 1/DELTA through cosines between radius vectors.
	- n is maximum degree of variables;
	- c_max is maximum power of cosines;
	- deg is degree of 1/DELTA.
	"""
	from pyranha.math import factorial
	start = time.time()
	ra = rrr(p, n, 1)
	br = one_r(q, n, 1)
	if isAver and deg == 1:
		for item in br.list:
			if item[1] == epst(1):
				delta = item[0]*item[1]
	else: delta = br
	cos_deg = [1]
	for i in range(1, c_max+1):
		if isCalc:
			cos_deg.append(cosine(p, q, n, i))
		else:
			cos_temp = epst()
			lf(cos_temp, PREFIX + 'COS/'+str(n)+'/'+str(i)+'/cosine_'+str(q)+str(p)+'.epst.boostp.bz2', df.boost_portable, cf.bzip2)
			cos_deg.append(cos_temp)
		if isPrint:
			if isCalc: print 'cosine with power of', i, 'is calculated for %.2f'%(time.time() - start), 's'
			else: print 'cosine with power of', i, 'is read for %.2f'%(time.time() - start), 's'
	##################################################
	pt.set_auto_truncate_degree(F(n), ['x'+str(p), 'y'+str(p), 'u'+str(p), 'v'+str(p), 'x'+str(q), 'y'+str(q), 'u'+str(q), 'v'+str(q)])
	rr = ra * br
	for c in range(1, c_max+1):
		legPn = 0
		for k in range(0, int(c/2) + 1):
			legPn = legPn + (-1)**k * F(factorial(2*(c - k)), 2**c * factorial(k) * factorial(c - k) * factorial(c - 2*k)) * cos_deg[c - 2*k]
		rr_n = br * rr**c
		if isAver and deg == 1:
			temp = 0
			for item in rr_n.list:
				for jtem in legPn.list:
					trig = item[1]*jtem[1]
					for ktem in trig.list:
						if ktem[1] == 1:
							temp = temp + item[0]*jtem[0]*ktem[0]
		else: temp = rr_n * legPn
		delta = delta + temp
		if isPrint: print 'Legendre polynomial with power of', c, 'is added for %.2f'%(time.time() - start), 's'
	if deg > 1: delta = delta**deg
	pt.unset_auto_truncate_degree()
	return delta
##################################################
# CONSTRUCTION OF THE UNDISTURBED HAMILTONIAN
##################################################
def h_0(num_pl = 4):
	"""
	Calculate undisturbed part of the Hamiltonian.
	"""
	start = time.time()
	print '> calc: undisturbed hamiltonian'
	L = [epst('L'+str(i+1)) for i in range(num_pl)]
	M = [epst('m'+str(i+1)) for i in range(num_pl)]
	S = [epst('s'+str(i+1)) for i in range(num_pl)]
	h0 = -F(1,2) * epst('K0')**2 * (M[0]**3 * S[0]**-1) * L[0]**-2
	for i in range(1, num_pl): h0 += -F(1,2) * epst('K0')**2 * (M[i]**3 * S[i-1] * S[i]**-1) * L[i]**-2
	file_out = 'H0.epst.boostp.bz2'
	print h0
	if not os.path.exists(PREFIX + 'HAM/'): os.makedirs(PREFIX + 'HAM/')
	sf(h0, PREFIX + 'HAM/'+file_out, df.boost_portable, cf.bzip2)
	print '> done! FILE:', file_out, 'on', '%.2f'%(time.time() - start), 's'
##################################################
# CONSTRUCTION OF THE MAJOR PART EXPANSION
##################################################
def main_1(p, q, n, c_max, isAver = True, isCalc = True, isPrint = True):
	"""
	Construction of the expansion 1/DELTA through cosines between radius vectors.
	- n is maximum degree of variables;
	- c_max is maximum power of a cosine.
	"""
	start = time.time()
	(p, q) = (q, p) if p > q else (p, q)
	print '> calc: main part (i='+str(q)+', j='+str(p)+', n='+str(n)+', c_max='+str(c_max)+') for 0th degree'
	delta = -epst('m'+str(p))*epst('m'+str(q)) * one_delta(p, q, n, c_max, 1, isAver, isCalc, isPrint)
	name = 'A' if isAver else 'H'
	file_out = name+'1_'+str(q)+str(p)+'_n'+str(n).rjust(2, '0')+'_c'+str(c_max).rjust(2, '0')+'.epst.boostp.bz2'
	if not os.path.exists(PREFIX + 'HAM/'): os.makedirs(PREFIX + 'HAM/')
	sf(delta.trim(), PREFIX + 'HAM/'+file_out, df.boost_portable, cf.bzip2)
	print '> done! FILE:', file_out, 'on', '%.2f'%(time.time() - start), 's'
##################################################
def main_2(p, q, n, c_max, isAver = True, isCalc = True, isPrint = True):
	"""
	Construction of the expansion A/DELTA^3 through cosines between radius vectors.
	- n is maximum degree of variables;
	- c_max is maximum power of a cosine.
	"""
	start = time.time()
	(p, q) = (q, p) if p > q else (p, q)
	print '> calc: main part (i='+str(q)+', j='+str(p)+', n='+str(n)+', c_max='+str(c_max)+') for 1st degree'
	pr = q - p
	mm = [epst('m'+str(p+i))*epst('s'+str(p+i))**-1 for i in range(pr)]
	r1 = [rrr(p+i, n, 1) for i in range(pr+1)]
	delta3 = one_delta(p, q, n, c_max, 3, False, isCalc, isPrint)
	cos_list = {str(q)+str(p):   cosine(q, p, n, 1),
				str(q)+str(q-1): cosine(q, q-1, n, 1), str(q)+str(q-2): cosine(q, q-2, n, 1),
				str(p+1)+str(p): cosine(p+1, p, n, 1), str(p+2)+str(p): cosine(p+2, p, n, 1)}
	pt.set_auto_truncate_degree(F(n), ['x'+str(p+i) for i in range(pr+1)] + ['y'+str(p+i) for i in range(pr+1)] +
									  ['u'+str(p+i) for i in range(pr+1)] + ['v'+str(p+i) for i in range(pr+1)])
	if isPrint: print '> calc: combinations of cosines'
	temp_cos = mm[0] * (r1[0]*r1[pr] * cos_list[str(q)+str(p)] - r1[0]*r1[0])
	if pr == 2:
		temp_cos += mm[1] * (r1[1]*r1[2] * cos_list[str(q)+str(q-1)] - r1[0]*r1[1] * cos_list[str(p+1)+str(p)])
	if pr == 3:
		temp_cos += mm[1] * (r1[1]*r1[3] * cos_list[str(q)+str(q-2)] - r1[0]*r1[1] * cos_list[str(p+1)+str(p)])
		temp_cos += mm[2] * (r1[2]*r1[3] * cos_list[str(q)+str(q-1)] - r1[0]*r1[2] * cos_list[str(p+2)+str(p)])
	if isPrint: print '> calc: multiplication of cosines and delta'
	if isAver:
		temp = 0
		for item in delta3.list:
			if isPrint: print 'delta #3:', item[1]
			for jtem in temp_cos.list:
				trig = item[1]*jtem[1]
				for ktem in trig.list:
					if ktem[1] == 1:
						temp = temp + item[0]*jtem[0]*ktem[0]
	else: temp = temp_cos * delta3
	temp = epst('m'+str(p))*epst('m'+str(q)) * temp
	name = 'A' if isAver else 'H'
	file_out = name+'2_'+str(q)+str(p)+'_n'+str(n).rjust(2, '0')+'_c'+str(c_max).rjust(2, '0')+'.epst.boostp.bz2'
	if not os.path.exists(PREFIX + 'HAM/'): os.makedirs(PREFIX + 'HAM/')
	sf(temp, PREFIX + 'HAM/'+file_out, df.boost_portable, cf.bzip2)
	print '> done! FILE:', file_out, 'on', '%.2f'%(time.time() - start), 's'
	pt.unset_auto_truncate_degree()
##################################################
def main_3(p, q, n, c_max1, c_max2, isAver = True, isCalc = True, isPrint = True):
	"""
	Construction of the expansion 1.5A^2/DELTA^5 - 0.5B/DELTA^3 through cosines between radius vectors.
	- n is maximum degree of variables;
	- c_max is maximum power of a cosine.
	"""
	start = time.time()
	(p, q) = (q, p) if p > q else (p, q)
	print '> calc: main part (i='+str(q)+', j='+str(p)+', n='+str(n)+', c_max#1='+str(c_max1)+', c_max#2='+str(c_max2)+') for 2nd degree'
	pr = q - p
	mm = [epst('m'+str(p+i))*epst('s'+str(p+i))**-1 for i in range(pr)]
	r1 = [rrr(p+i, n, 1) for i in range(pr+1)]
	delta3 = one_delta(p, q, n, c_max1, 3, False, isCalc, isPrint)
	delta5 = one_delta(p, q, n, c_max2, 5, False, isCalc, isPrint)
	cos_list = {str(q)+str(p):   cosine(q, p, n, 1),   str(q-1)+str(p):   cosine(q-1, p, n, 1),
				str(q-2)+str(p): cosine(q-2, p, n, 1), str(q-1)+str(p+1): cosine(q-1, p+1, n, 1),
				str(q)+str(q-1): cosine(q, q-1, n, 1), str(q)+str(q-2):   cosine(q, q-2, n, 1),
				str(p+1)+str(p): cosine(p+1, p, n, 1), str(p+2)+str(p):   cosine(p+2, p, n, 1)}
	pt.set_auto_truncate_degree(F(n), ['x'+str(p+i) for i in range(pr+1)] + ['y'+str(p+i) for i in range(pr+1)] +
									  ['u'+str(p+i) for i in range(pr+1)] + ['v'+str(p+i) for i in range(pr+1)])
	if isPrint: print '> calc: combinations of cosines'
	temp_cos1 = mm[0] * (r1[0]*r1[pr] * cos_list[str(q)+str(p)] - r1[0]*r1[0])
	temp_cos2 = mm[0]**2 * r1[0]*r1[0]
	if pr == 2:
		temp_cos1 += mm[1] * (r1[1]*r1[2] * cos_list[str(q)+str(q-1)] - r1[0]*r1[1] * cos_list[str(p+1)+str(p)])
		temp_cos2 += 2*mm[0]*mm[1] * r1[0]*r1[1] * cos_list[str(q-1)+str(p)] + mm[1]**2 * r1[1]*r1[1]
	if pr == 3:
		temp_cos1 += mm[1] * (r1[1]*r1[3] * cos_list[str(q)+str(q-2)] - r1[0]*r1[1] * cos_list[str(p+1)+str(p)])
		temp_cos1 += mm[2] * (r1[2]*r1[3] * cos_list[str(q)+str(q-1)] - r1[0]*r1[2] * cos_list[str(p+2)+str(p)])
		temp_cos2 += mm[1]**2 * r1[1]*r1[1] + mm[2]**2 * r1[2]*r1[2]
		temp_cos2 += 2*mm[0]*mm[1] * r1[0]*r1[1] * cos_list[str(q-2)+str(p)] + 2*mm[0]*mm[2] * r1[0]*r1[2] * cos_list[str(q-1)+str(p)] + \
					 2*mm[1]*mm[2] * r1[1]*r1[2] * cos_list[str(q-1)+str(p+1)]
	temp_cos1 = temp_cos1**2
	if isPrint: print '> calc: multiplication of cosines and delta'
	if isAver:
		temp1 = 0
		for item in delta5.list:
			if isPrint: print 'delta #5:', item[1]
			for jtem in temp_cos1.list:
				trig = item[1]*jtem[1]
				for ktem in trig.list:
					if ktem[1] == 1:
						temp1 = temp1 + item[0]*jtem[0]*ktem[0]
		temp2 = 0
		for item in delta3.list:
			if isPrint: print 'delta #3:', item[1]
			for jtem in temp_cos2.list:
				trig = item[1]*jtem[1]
				for ktem in trig.list:
					if ktem[1] == 1:
						temp2 = temp2 + item[0]*jtem[0]*ktem[0]
	else:
		temp1 = temp_cos1 * delta5
		temp2 = temp_cos2 * delta3
	temp1 = -F(3,2) * epst('m'+str(p))*epst('m'+str(q)) * temp1
	temp2 =  F(1,2) * epst('m'+str(p))*epst('m'+str(q)) * temp2
	temp = temp1 + temp2
	name = 'A' if isAver else 'H'
	c_max = c_max1 if c_max1 >= c_max2 else c_max2
	file_out = name+'3_'+str(q)+str(p)+'_n'+str(n).rjust(2, '0')+'_c'+str(c_max).rjust(2, '0')+'.epst.boostp.bz2'
	if not os.path.exists(PREFIX + 'HAM/'): os.makedirs(PREFIX + 'HAM/')
	sf(temp, PREFIX + 'HAM/'+file_out, df.boost_portable, cf.bzip2)
	print '> done! FILE:', file_out, 'on', '%.2f'%(time.time() - start), 's'
	pt.unset_auto_truncate_degree()
##################################################
# CONSTRUCTION OF THE SECOND PART EXPANSION
##################################################
def second_1(n, isAver = True):
	"""
	Construction of the second part expansion C/r^3.
	- n is maximum degree of variables.
	"""
	start = time.time()
	print '> calc: second part (n='+str(n)+') for 0th degree'
	mm = [epst('m'+str(1+i))*epst('s'+str(1+i))**-1 for i in range(3)]
	r1 = [rrr(1+i, n, 1) for i in range(4)]
	one_r3 = [one_r(2+i, n, 3) for i in range(3)]
	cos_list = {'P': cosine(2, 1, n, 1), 'Q': cosine(3, 2, n, 1), 'R': cosine(4, 3, n, 1),
				'S': cosine(3, 1, n, 1), 'T': cosine(4, 2, n, 1), 'U': cosine(4, 1, n, 1)}
	pt.set_auto_truncate_degree(F(n), ['x'+str(1+i) for i in range(4)] + ['y'+str(1+i) for i in range(4)] +
									  ['u'+str(1+i) for i in range(4)] + ['v'+str(1+i) for i in range(4)])
	temp = epst('m2') * (mm[0]*r1[0]*r1[1] * cos_list['P']) * one_r3[0]
	temp += epst('m3') * (mm[0]*r1[0]*r1[2] * cos_list['S'] + mm[1]*r1[1]*r1[2] * cos_list['Q']) * one_r3[1]
	temp += epst('m4') * (mm[0]*r1[0]*r1[3] * cos_list['U'] + mm[1]*r1[1]*r1[3] * cos_list['T'] + mm[2]*r1[2]*r1[3] *  cos_list['R']) * one_r3[2]
	if isAver:
		aver = epst(0)
		for item in temp.list:
			if item[1] == epst(1): aver = item[0]*item[1]
		temp = aver
	name = 'A' if isAver else 'H'
	file_out = name+'1_44_n'+str(n).rjust(2, '0')+'.epst.boostp.bz2'
	if not os.path.exists(PREFIX + 'HAM/'): os.makedirs(PREFIX + 'HAM/')
	sf(temp, PREFIX + 'HAM/'+file_out, df.boost_portable, cf.bzip2)
	print '> done! FILE:', file_out, 'on', '%.2f'%(time.time() - start), 's'
	pt.unset_auto_truncate_degree()
##################################################
def second_2(n, isAver = True):
	"""
	Construction of the second part expansions -1.5C^2/r^5 + 0.5D/r^3.
	- n is maximum degree of variables.
	"""
	start = time.time()
	print '> calc: second part (n='+str(n)+') for 1st degree'
	mm = [epst('m'+str(1+i))*epst('s'+str(1+i))**-1 for i in range(3)]
	r1 = [rrr(1+i, n, 1) for i in range(4)]
	r2 = [rrr(1+i, n, 2) for i in range(3)]
	one_r3 = [one_r(2+i, n, 3) for i in range(3)]
	one_r5 = [one_r(2+i, n, 5) for i in range(3)]
	cos_list = {'P': cosine(2, 1, n, 1), 'Q': cosine(3, 2, n, 1), 'R': cosine(4, 3, n, 1),
				'S': cosine(3, 1, n, 1), 'T': cosine(4, 2, n, 1), 'U': cosine(4, 1, n, 1)}
	pt.set_auto_truncate_degree(F(n), ['x'+str(1+i) for i in range(4)] + ['y'+str(1+i) for i in range(4)] +
									  ['u'+str(1+i) for i in range(4)] + ['v'+str(1+i) for i in range(4)])
	temp1 = -F(3,2) * (mm[0]*r1[0]*r1[1] * cos_list['P'])**2 * one_r5[0]
	temp2 = mm[0]**2 * r2[0]
	temp3 = epst('m2') * (temp1 + F(1,2) * temp2 * one_r3[0])
	temp1 = -F(3,2) * (mm[0]*r1[0]*r1[2] * cos_list['S'] + mm[1]*r1[1]*r1[2] * cos_list['Q'])**2 * one_r5[1]
	temp2 += mm[1]**2 * r2[1] + 2*mm[0]*mm[1] * r1[0]*r1[1] * cos_list['P']
	temp3 += epst('m3') * (temp1 + F(1,2) * temp2 * one_r3[1])
	temp1 = -F(3,2) * (mm[0]*r1[0]*r1[3] * cos_list['U'] + mm[1]*r1[1]*r1[3] * cos_list['T'] + mm[2]*r1[2]*r1[3] * cos_list['R'])**2 * one_r5[2]
	temp2 += mm[2]**2 * r2[2] + 2*mm[1]*mm[2] * r1[1]*r1[2] * cos_list['Q'] + 2*mm[0]*mm[2] * r1[0]*r1[2] * cos_list['S']
	temp3 += epst('m4') * (temp1 + F(1,2) * temp2 * one_r3[2])
	if isAver:
		aver = epst(0)
		for item in temp3.list:
			if item[1] == epst(1): aver = item[0]*item[1]
		temp3 = aver
	name = 'A' if isAver else 'H'
	file_out = name+'2_44_n'+str(n).rjust(2, '0')+'.epst.boostp.bz2'
	if not os.path.exists(PREFIX + 'HAM/'): os.makedirs(PREFIX + 'HAM/')
	sf(temp3, PREFIX + 'HAM/'+file_out, df.boost_portable, cf.bzip2)
	print '> done! FILE:', file_out, 'on', '%.2f'%(time.time() - start), 's'
	pt.unset_auto_truncate_degree()
##################################################
def second_3(n, isAver = True):
	"""
	Construction of the second part expansions 2.5C^3/r^7 - 1.5CD/r^5.
	- n is maximum degree of variables.
	"""
	start = time.time()
	print '> calc: second part (n='+str(n)+') for 2nd degree'
	mm = [epst('m'+str(1+i))*epst('s'+str(1+i))**-1 for i in range(3)]
	r1 = [rrr(1+i, n, 1) for i in range(4)]
	r2 = [rrr(1+i, n, 2) for i in range(3)]
	one_r5 = [one_r(2+i, n, 5) for i in range(3)]
	one_r7 = [one_r(2+i, n, 7) for i in range(3)]
	cos_list = {'P': cosine(2, 1, n, 1), 'Q': cosine(3, 2, n, 1), 'R': cosine(4, 3, n, 1),
				'S': cosine(3, 1, n, 1), 'T': cosine(4, 2, n, 1), 'U': cosine(4, 1, n, 1)}
	pt.set_auto_truncate_degree(F(n), ['x'+str(1+i) for i in range(4)] + ['y'+str(1+i) for i in range(4)] +
									  ['u'+str(1+i) for i in range(4)] + ['v'+str(1+i) for i in range(4)])
	temp1 = mm[0]*r1[0]*r1[1] * cos_list['P']
	temp2 = mm[0]**2 * r2[0]
	temp3 = epst('m2') * (F(5,2) * temp1**3 * one_r7[0] - F(3,2) * temp1 * temp2 * one_r5[0])
	temp1 = mm[0]*r1[0]*r1[2] * cos_list['S'] + mm[1]*r1[1]*r1[2] * cos_list['Q']
	temp2 += mm[1]**2 * r2[1] + 2*mm[0]*mm[1] * r1[0]*r1[1] * cos_list['P']
	temp3 += epst('m3') * (F(5,2) * temp1**3 * one_r7[1] - F(3,2) * temp1 * temp2 * one_r5[1])
	temp1 = mm[0]*r1[0]*r1[3] * cos_list['U'] + mm[1]*r1[1]*r1[3] * cos_list['T'] + mm[2]*r1[2]*r1[3] * cos_list['R']
	temp2 += mm[2]**2 * r2[2] + 2*mm[1]*mm[2] * r1[1]*r1[2] * cos_list['Q'] + 2*mm[0]*mm[2] * r1[0]*r1[2] * cos_list['S']
	temp3 += epst('m4') * (F(5,2) * temp1**3 * one_r7[2] - F(3,2) * temp1 * temp2 * one_r5[2])
	if isAver:
		aver = epst(0)
		for item in temp3.list:
			if item[1] == epst(1): aver = item[0]*item[1]
		temp3 = aver
	name = 'A' if isAver else 'H'
	file_out = name+'3_44_n'+str(n).rjust(2, '0')+'.epst.boostp.bz2'
	if not os.path.exists(PREFIX + 'HAM/'): os.makedirs(PREFIX + 'HAM/')
	sf(temp3, PREFIX + 'HAM/'+file_out, df.boost_portable, cf.bzip2)
	print '> done! FILE:', file_out, 'on', '%.2f'%(time.time() - start), 's'
	pt.unset_auto_truncate_degree()
##################################################
# TOOL FUNCTIONS
##################################################
def check_base(elmass, num = 0, p = 1, q = 2, c_max = 25, deg = 1):
	"""
	Checking accuracy of base series.
	"""
	from math import fabs
	from mpmath import mpf
	from pyranha.math import evaluate
	from eps.tools import kepler2cart, kepler2poincare, jacobiAMass
	elements, gmass, mu = elmass
	# -------- masses --------
	mp = [1./mu]+[gmass[i]/gmass[0]/mu for i in range(len(gmass))][1:]
	sp = [1. + mu*sum(mp[1:i+1]) for i in range(len(gmass))]
	gmj, mpj, km = jacobiAMass(gmass, mu)
	# -------- coordinates --------
	rlen1 = range(len(elements))
	x, y, z = [0. for i in rlen1], [0. for i in rlen1], [0. for i in rlen1]
	for i in rlen1:
		x[i], y[i], z[i] = kepler2cart(elements[i], km[i])[:3]
	r  = [(x[i]**2 + y[i]**2 + z[i]**2)**0.5  for i in rlen1]
	PE = [kepler2poincare(elements[i], km[i]) for i in rlen1]
	# pairs 21, 32, 43, 31, 42, 41:
	rr = {'21': r[0]/r[1], '32': r[1]/r[2], '43': r[2]/r[3], '31': r[0]/r[2], '42': r[1]/r[3], '41': r[0]/r[3]}
	cs = {'21': (x[0]*x[1] + y[0]*y[1] + z[0]*z[1])/r[0]/r[1], '32': (x[1]*x[2] + y[1]*y[2] + z[1]*z[2])/r[1]/r[2],
		  '43': (x[2]*x[3] + y[2]*y[3] + z[2]*z[3])/r[2]/r[3], '31': (x[0]*x[2] + y[0]*y[2] + z[0]*z[2])/r[0]/r[2],
		  '42': (x[1]*x[3] + y[1]*y[3] + z[1]*z[3])/r[1]/r[3], '41': (x[0]*x[3] + y[0]*y[3] + z[0]*z[3])/r[0]/r[3]}
	dd = {'21': ((x[1]-x[0])**2 + (y[1]-y[0])**2 + (z[1]-z[0])**2)**-0.5, '32': ((x[2]-x[1])**2 + (y[2]-y[1])**2 + (z[2]-z[1])**2)**-0.5,
		  '43': ((x[3]-x[2])**2 + (y[3]-y[2])**2 + (z[3]-z[2])**2)**-0.5, '31': ((x[2]-x[0])**2 + (y[2]-y[0])**2 + (z[2]-z[0])**2)**-0.5,
		  '42': ((x[3]-x[1])**2 + (y[3]-y[1])**2 + (z[3]-z[1])**2)**-0.5, '41': ((x[3]-x[0])**2 + (y[3]-y[0])**2 + (z[3]-z[0])**2)**-0.5}
	# -------- dictionaries --------
	evaldict = {'K0': mpf(gmass[0]), 'm1': mpf(mp[1]), 'm2': mpf(mp[2]), 'm3': mpf(mp[3]), 'm4': mpf(mp[4]), 's1': mpf(sp[1]), 's2': mpf(sp[2]), 's3': mpf(sp[3]), 's4': mpf(sp[4]),\
				'K1': mpf(gmj[0]**0.5), 'M1': mpf(mpj[0]), 'L1': mpf(PE[0][0]), 'x1': mpf(PE[0][1]), 'y1': mpf(PE[0][2]), 'u1': mpf(PE[0][3]), 'v1': mpf(PE[0][4]), 'q1': mpf(PE[0][5]),\
				'K2': mpf(gmj[1]**0.5), 'M2': mpf(mpj[1]), 'L2': mpf(PE[1][0]), 'x2': mpf(PE[1][1]), 'y2': mpf(PE[1][2]), 'u2': mpf(PE[1][3]), 'v2': mpf(PE[1][4]), 'q2': mpf(PE[1][5]),\
				'K3': mpf(gmj[2]**0.5), 'M3': mpf(mpj[2]), 'L3': mpf(PE[2][0]), 'x3': mpf(PE[2][1]), 'y3': mpf(PE[2][2]), 'u3': mpf(PE[2][3]), 'v3': mpf(PE[2][4]), 'q3': mpf(PE[2][5]),\
				'K4': mpf(gmj[3]**0.5), 'M4': mpf(mpj[3]), 'L4': mpf(PE[3][0]), 'x4': mpf(PE[3][1]), 'y4': mpf(PE[3][2]), 'u4': mpf(PE[3][3]), 'v4': mpf(PE[3][4]), 'q4': mpf(PE[3][5])}
	# -------- evaluation x, y, z --------
	if num in [0, 1]:
		for n in range(1, 17):
			start = time.time()
			coord_s = xyz(p, n)
			coord_e = [x[p-1], y[p-1], z[p-1]]
			accuracy = [0, 0, 0]
			sum_items = [0, 0, 0]
			finish = time.time()-start
			for i in range(len(coord_s)):
				for item in coord_s[i].list:
					for jtem in item[0].list:
						sum_items[i] += len(jtem[0].list)
				accuracy[i] = fabs((coord_e[i] - evaluate(coord_s[i], evaldict))/coord_e[i])
			print 'n =', str(n).rjust(2, ' '),  '|x  ', ('%.2f'%finish).rjust(5, ' '), str(sum_items[0]).rjust(5, ' '), '%.1e'%accuracy[0], ' ',\
												'|y  ', ('%.2f'%finish).rjust(5, ' '), str(sum_items[1]).rjust(5, ' '), '%.1e'%accuracy[1], ' ',\
												'|z  ', ('%.2f'%finish).rjust(5, ' '), str(sum_items[2]).rjust(5, ' '), '%.1e'%accuracy[2]
	# -------- evaluation r, 1/r', r/r' --------
	if num in [0, 2]:
		for n in range(1, 17):
			temp, finish, accuracy = [], [], [0, 0, 0]
			start = time.time()
			temp.append(rrr(p, n, deg))
			finish.append(time.time()-start)
			start = time.time()
			temp.append(one_r(q, n, deg))
			finish.append(time.time()-start)
			pt.set_auto_truncate_degree(F(n), ['x'+str(p), 'y'+str(p), 'u'+str(p), 'v'+str(p), 'x'+str(q), 'y'+str(q), 'u'+str(q), 'v'+str(q)])
			start = time.time()
			temp.append(rrr(p, n, deg) * one_r(q, n, deg))
			finish.append(time.time()-start)
			pt.unset_auto_truncate_degree()
			sum_items = [0, 0, 0]
			for i in range(3):
				for item in temp[i].list:
					for jtem in item[0].list:
						sum_items[i] += len(jtem[0].list)
			accuracy[0] = fabs((r[p-1]**deg - evaluate(temp[0], evaldict))/r[p-1]**deg)
			accuracy[1] = fabs((r[q-1]**-deg - evaluate(temp[1], evaldict))/r[q-1]**-deg)
			accuracy[2] = fabs(((rr[str(q)+str(p)])**deg - evaluate(temp[2], evaldict))/rr[str(q)+str(p)]**deg)
			print 'n =', str(n).rjust(2, ' '),  '|r  ', ('%.2f'%finish[0]).rjust(5, ' '), str(sum_items[0]).rjust(5, ' '), '%.1e'%accuracy[0], ' ',\
												'|1/r', ('%.2f'%finish[1]).rjust(5, ' '), str(sum_items[1]).rjust(5, ' '), '%.1e'%accuracy[1], ' ',\
												'|r/r', ('%.2f'%finish[2]).rjust(5, ' '), str(sum_items[2]).rjust(5, ' '), '%.1e'%accuracy[2]
	# -------- evaluation cos --------
	if num in [0, 3]:
		for n in range(1, 17):
			start = time.time()
			temp = cosine(p, q, n, deg)
			finish = time.time()-start
			sum_items = 0
			for item in temp.list:
				for jtem in item[0].list:
					sum_items += len(jtem[0].list)
			accuracy = fabs((cs[str(q)+str(p)]**deg - evaluate(temp, evaldict))/cs[str(q)+str(p)]**deg)
			print 'n =', str(n).rjust(2, ' '), '|cos', '%.2f'%finish, str(sum_items).rjust(8, ' '), '%.1e'%accuracy
	# -------- evaluation 1/d --------
	if num in [0, 4]:
		for i in range(1, 17):
			start = time.time()
			temp = one_delta(p, q, i, c_max, deg, isPrint = False)
			finish = time.time()-start
			sum_items = 0
			for item in temp.list:
				for jtem in item[0].list:
					sum_items += len(jtem[0].list)
			evaluated = evaluate(temp, evaldict)
			accuracy = fabs((dd[str(q)+str(p)]**deg - evaluated)/dd[str(q)+str(p)]**deg)
			print 'n =', str(n).rjust(2, ' '), '|1/d', '%.2f'%finish, str(sum_items).rjust(8, ' '), '%.16e'%dd[str(q)+str(p)]**deg, '%.16e'%evulated, '%.1e'%accuracy
##################################################
def check_items(elmass, series2, isPrint):
	"""
	Checking accuracy of items of the hamiltonian.
	"""
	from math import fabs
	from mpmath import mpf
	from pyranha.math import evaluate
	from eps.tools import kepler2cart, kepler2poincare, poincare2kepler, jacobiAMass
	elements, gmass, mu = elmass
	# -------- masses --------
	mp = [1./mu]+[gmass[i]/gmass[0]/mu for i in range(len(gmass))][1:]
	sp = [1. + mu*sum(mp[1:i+1]) for i in range(len(gmass))]
	gmj, mpj, km = jacobiAMass(gmass, mu)
	# -------- coordinates --------
	rlen1 = range(len(elements))
	kepler = [poincare2kepler(elements[i], km[i]) for i in rlen1]
	x, y, z = [0. for i in rlen1], [0. for i in rlen1], [0. for i in rlen1]
	for i in rlen1:
		x[i], y[i], z[i] = kepler2cart(kepler[i], km[i])[:3]
	r  = [(x[i]**2 + y[i]**2 + z[i]**2)**0.5  for i in rlen1]
	PE = [kepler2poincare(kepler[i], km[i]) for i in rlen1]
	# pairs 21, 32, 43, 31, 42, 41:
	rr = {'21': r[0]/r[1], '32': r[1]/r[2], '43': r[2]/r[3], '31': r[0]/r[2], '42': r[1]/r[3], '41': r[0]/r[3]}
	sc = {'21': (x[0]*x[1] + y[0]*y[1] + z[0]*z[1]), '32': (x[1]*x[2] + y[1]*y[2] + z[1]*z[2]),
		  '43': (x[2]*x[3] + y[2]*y[3] + z[2]*z[3]), '31': (x[0]*x[2] + y[0]*y[2] + z[0]*z[2]),
		  '42': (x[1]*x[3] + y[1]*y[3] + z[1]*z[3]), '41': (x[0]*x[3] + y[0]*y[3] + z[0]*z[3])}
	cs = {'21': (x[0]*x[1] + y[0]*y[1] + z[0]*z[1])/r[0]/r[1], '32': (x[1]*x[2] + y[1]*y[2] + z[1]*z[2])/r[1]/r[2],
		  '43': (x[2]*x[3] + y[2]*y[3] + z[2]*z[3])/r[2]/r[3], '31': (x[0]*x[2] + y[0]*y[2] + z[0]*z[2])/r[0]/r[2],
		  '42': (x[1]*x[3] + y[1]*y[3] + z[1]*z[3])/r[1]/r[3], '41': (x[0]*x[3] + y[0]*y[3] + z[0]*z[3])/r[0]/r[3]}
	dd = {'21': ((x[1]-x[0])**2 + (y[1]-y[0])**2 + (z[1]-z[0])**2)**-0.5, '32': ((x[2]-x[1])**2 + (y[2]-y[1])**2 + (z[2]-z[1])**2)**-0.5,
		  '43': ((x[3]-x[2])**2 + (y[3]-y[2])**2 + (z[3]-z[2])**2)**-0.5, '31': ((x[2]-x[0])**2 + (y[2]-y[0])**2 + (z[2]-z[0])**2)**-0.5,
		  '42': ((x[3]-x[1])**2 + (y[3]-y[1])**2 + (z[3]-z[1])**2)**-0.5, '41': ((x[3]-x[0])**2 + (y[3]-y[0])**2 + (z[3]-z[0])**2)**-0.5}
	# -------- dictionaries --------
	evaldict = {'K0': mpf(gmass[0]), 'm1': mpf(mp[1]), 'm2': mpf(mp[2]), 'm3': mpf(mp[3]), 'm4': mpf(mp[4]), 's1': mpf(sp[1]), 's2': mpf(sp[2]), 's3': mpf(sp[3]), 's4': mpf(sp[4]),\
				'K1': mpf(gmj[0]**0.5), 'M1': mpf(mpj[0]), 'L1': mpf(PE[0][0]), 'x1': mpf(PE[0][1]), 'y1': mpf(PE[0][2]), 'u1': mpf(PE[0][3]), 'v1': mpf(PE[0][4]), 'q1': mpf(PE[0][5]),\
				'K2': mpf(gmj[1]**0.5), 'M2': mpf(mpj[1]), 'L2': mpf(PE[1][0]), 'x2': mpf(PE[1][1]), 'y2': mpf(PE[1][2]), 'u2': mpf(PE[1][3]), 'v2': mpf(PE[1][4]), 'q2': mpf(PE[1][5]),\
				'K3': mpf(gmj[2]**0.5), 'M3': mpf(mpj[2]), 'L3': mpf(PE[2][0]), 'x3': mpf(PE[2][1]), 'y3': mpf(PE[2][2]), 'u3': mpf(PE[2][3]), 'v3': mpf(PE[2][4]), 'q3': mpf(PE[2][5]),\
				'K4': mpf(gmj[3]**0.5), 'M4': mpf(mpj[3]), 'L4': mpf(PE[3][0]), 'x4': mpf(PE[3][1]), 'y4': mpf(PE[3][2]), 'u4': mpf(PE[3][3]), 'v4': mpf(PE[3][4]), 'q4': mpf(PE[3][5])}
	# -------- evaluation --------
	temp1_all, temp2_all, length1_all = 0, 0, 0
	tmp1, tmp2 = [0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0]
	for file2, series1, length1 in series2:
		if file2[0] != 'A':
			temp1, temp2 = 0, 0
			length1_all += length1
			if isPrint: print file2[:-16]+':',
			if file2[:2] == 'H0':
				temp1 = evaluate(series1, evaldict)
				temp2 = -F(1,2) * (mpj[0]**3 * gmj[0]**2 * PE[0][0]**-2 + mpj[1]**3 * gmj[1]**2 * PE[1][0]**-2 + mpj[2]**3 * gmj[2]**2 * PE[2][0]**-2 + mpj[3]**3 * gmj[3]**2 * PE[3][0]**-2)
				temp1_h0, temp2_h0 = temp1, temp2
			if '_21_' in file2:
				temp1 = evaluate(series1, evaldict)
				if 'H1' in file2:
					temp2 = -mp[1]*mp[2]*dd['21']
					temp1_all += temp1; temp2_all += temp2
					tmp1[0] += temp1; tmp2[0] += temp2
				if 'H2' in file2:
					temp2 = mp[1]*mp[2]*(mm[1]*sc['21'] - mm[1]*r[0]**2) * dd['21']**3
					temp1_all += mu*temp1; temp2_all += mu*temp2
					tmp1[0] += mu*temp1; tmp2[0] += mu*temp2
				if 'H3' in file2:
					temp2 = -1.5*mp[1]*mp[2]*(mm[1]*sc['21'] - mm[1]*r[0]**2)**2 * dd['21']**5 + 0.5*mp[1]*mp[2]*(mm[1]*r[0])**2 * dd['21']**3
					temp1_all += mu**2*temp1; temp2_all += mu**2*temp2
					tmp1[0] += mu**2*temp1; tmp2[0] += mu**2*temp2
			if '_32_' in file2:
				temp1 = evaluate(series1, evaldict)
				if 'H1' in file2:
					temp2 = -mp[2]*mp[3]*dd['32']
					temp1_all += temp1; temp2_all += temp2
					tmp1[1] += temp1; tmp2[1] += temp2
				if 'H2' in file2:
					temp2 = mp[2]*mp[3]*(mm[2]*sc['32'] - mm[2]*r[1]**2) * dd['32']**3
					temp1_all += mu*temp1; temp2_all += mu*temp2
					tmp1[1] += mu*temp1; tmp2[1] += mu*temp2
				if 'H3' in file2:
					temp2 = -1.5*mp[2]*mp[3]*(mm[2]*sc['32'] - mm[2]*r[1]**2)**2 * dd['32']**5 + 0.5*mp[2]*mp[3]*(mm[2]*r[1])**2 * dd['32']**3
					temp1_all += mu**2*temp1; temp2_all += mu**2*temp2
					tmp1[1] += mu**2*temp1; tmp2[1] += mu**2*temp2
			if '_43_' in file2:
				temp1 = evaluate(series1, evaldict)
				if 'H1' in file2:
					temp2 = -mp[3]*mp[4]*dd['43']
					temp1_all += temp1; temp2_all += temp2
					tmp1[2] += temp1; tmp2[2] += temp2
				if 'H2' in file2:
					temp2 = mp[3]*mp[4]*(mm[3]*sc['43'] - mm[3]*r[2]**2) * dd['43']**3
					temp1_all += mu*temp1; temp2_all += mu*temp2
					tmp1[2] += mu*temp1; tmp2[2] += mu*temp2
				if 'H3' in file2:
					temp2 = -1.5*mp[3]*mp[4]*(mm[3]*sc['43'] - mm[3]*r[2]**2)**2 * dd['43']**5 + 0.5*mp[3]*mp[4]*(mm[3]*r[2])**2 * dd['43']**3
					temp1_all += mu**2*temp1; temp2_all += mu**2*temp2
					tmp1[2] += mu**2*temp1; tmp2[2] += mu**2*temp2
			if '_31_' in file2:
				temp1 = evaluate(series1, evaldict)
				if 'H1' in file2:
					temp2 = -mp[1]*mp[3]*dd['31']
					temp1_all += temp1; temp2_all += temp2
					tmp1[3] += temp1; tmp2[3] += temp2
				if 'H2' in file2:
					temp2 = mp[1]*mp[3]*(mm[1]*sc['31'] - mm[1]*r[0]**2 + mm[2]*sc['32'] - mm[2]*sc['21']) * dd['31']**3
					temp1_all += mu*temp1; temp2_all += mu*temp2
					tmp1[3] += mu*temp1; tmp2[3] += mu*temp2
				if 'H3' in file2:
					temp2 = -1.5*mp[1]*mp[3]*(mm[1]*sc['31'] - mm[1]*r[0]**2 + mm[2]*sc['32'] - mm[2]*sc['21'])**2 * dd['31']**5 +\
							 0.5*mp[1]*mp[3]*((mm[1]*r[0])**2 + (mm[2]*r[1])**2 + 2*mm[1]*mm[2]*sc['21']) * dd['31']**3
					temp1_all += mu**2*temp1; temp2_all += mu**2*temp2
					tmp1[3] += mu**2*temp1; tmp2[3] += mu**2*temp2
			if '_42_' in file2:
				temp1 = evaluate(series1, evaldict)
				if 'H1' in file2:
					temp2 = -mp[2]*mp[4]*dd['42']
					temp1_all += temp1; temp2_all += temp2
					tmp1[4] += temp1; tmp2[4] += temp2
				if 'H2' in file2:
					temp2 = mp[2]*mp[4]*(mm[2]*sc['42'] - mm[2]*r[1]**2 + mm[3]*sc['43'] - mm[3]*sc['32']) * dd['42']**3
					temp1_all += mu*temp1; temp2_all += mu*temp2
					tmp1[4] += mu*temp1; tmp2[4] += mu*temp2
				if 'H3' in file2:
					temp2 = -1.5*mp[2]*mp[4]*(mm[2]*sc['42'] - mm[2]*r[1]**2 + mm[3]*sc['43'] - mm[3]*sc['32'])**2 * dd['42']**5 +\
							 0.5*mp[2]*mp[4]*((mm[2]*r[1])**2 + (mm[3]*r[2])**2 + 2*mm[2]*mm[3]*sc['32']) * dd['42']**3
					temp1_all += mu**2*temp1; temp2_all += mu**2*temp2
					tmp1[4] += mu**2*temp1; tmp2[4] += mu**2*temp2
			if '_41_' in file2:
				temp1 = evaluate(series1, evaldict)
				if 'H1' in file2:
					temp2 = -mp[1]*mp[4]*dd['41']
					temp1_all += temp1; temp2_all += temp2
					tmp1[5] += temp1; tmp2[5] += temp2
				if 'H2' in file2:
					temp2 = mp[1]*mp[4]*(mm[1]*sc['41'] - mm[1]*r[0]**2 + mm[2]*sc['42'] - mm[2]*sc['21'] + mm[3]*sc['43'] - mm[3]*sc['31']) * dd['41']**3
					temp1_all += mu*temp1; temp2_all += mu*temp2
					tmp1[5] += mu*temp1; tmp2[5] += mu*temp2
				if 'H3' in file2:
					temp2 = -1.5*mp[1]*mp[4]*(mm[1]*sc['41'] - mm[1]*r[0]**2 + mm[2]*sc['42'] - mm[2]*sc['21'] + mm[3]*sc['43'] - mm[3]*sc['31'])**2 * dd['41']**5 +\
							 0.5*mp[1]*mp[4]*((mm[1]*r[0])**2 + (mm[2]*r[1])**2 + (mm[3]*r[2])**2 + 2*mm[1]*mm[2]*sc['21'] + 2*mm[1]*mm[3]*sc['31'] + 2*mm[2]*mm[3]*sc['32']) * dd['41']**3
					temp1_all += mu**2*temp1; temp2_all += mu**2*temp2
					tmp1[5] += mu**2*temp1; tmp2[5] += mu**2*temp2
			if 'H1_44_' in file2:
				temp1 = evaluate(series1, evaldict)
				temp2 = mp[2]*(mm[1]*sc['21'])*r[1]**-3 + mp[3]*(mm[1]*sc['31'] + mm[2]*sc['32']) * r[2]**-3 + mp[4]*(mm[1]*sc['41'] + mm[2]*sc['42'] + mm[3]*sc['43']) * r[3]**-3
				temp1_all += temp1; temp2_all += temp2
				tmp1[6] += temp1; tmp2[6] += temp2
			if 'H2_44_' in file2:
				temp1 = evaluate(series1, evaldict)
				temp2 = mp[2]*(-1.5*(mm[1]*sc['21'])**2*r[1]**-5 + 0.5*(mm[1]*r[0])**2 * r[1]**-3) + \
						mp[3]*(-1.5*(mm[1]*sc['31'] + mm[2]*sc['32'])**2 * r[2]**-5 + 0.5*((mm[1]*r[0])**2 + (mm[2]*r[1])**2 + 2*mm[1]*mm[2]*sc['21']) * r[2]**-3) + \
						mp[4]*(-1.5*(mm[1]*sc['41'] + mm[2]*sc['42'] + mm[3]*sc['43'])**2 * r[3]**-5 + \
								0.5*((mm[1]*r[0])**2 + (mm[2]*r[1])**2 + (mm[3]*r[2])**2 + 2*mm[1]*mm[2]*sc['21'] + 2*mm[2]*mm[3]*sc['32'] + 2*mm[1]*mm[3]*sc['31']) * r[3]**-3)
				temp1_all += mu*temp1; temp2_all += mu*temp2
				tmp1[6] += mu*temp1; tmp2[6] += mu*temp2
			if 'H3_44_' in file2:
				temp1 = evaluate(series1, evaldict)
				temp2 = mp[2]*(2.5*(mm[1]*sc['21'])**3 * r[1]**-7 - 1.5*(mm[1]*sc['21'])*(mm[1]*r[0])**2 * r[1]**-5) + \
						mp[3]*(2.5*(mm[1]*sc['31'] + mm[2]*sc['32'])**3 * r[2]**-7 - 1.5*(mm[1]*sc['31'] + mm[2]*sc['32'])*((mm[1]*r[0])**2 + (mm[2]*r[1])**2 + 2*mm[1]*mm[2]*sc['21']) * r[2]**-5) + \
						mp[4]*(2.5*(mm[1]*sc['41'] + mm[2]*sc['42'] + mm[3]*sc['43'])**3 * r[3]**-7 - \
							   1.5*(mm[1]*sc['41'] + mm[2]*sc['42'] + mm[3]*sc['43'])* \
								  ((mm[1]*r[0])**2 + (mm[2]*r[1])**2 + (mm[3]*r[2])**2 + 2*mm[1]*mm[2]*sc['21'] + 2*mm[2]*mm[3]*sc['32'] + 2*mm[1]*mm[3]*sc['31']) * r[3]**-5)
				temp1_all += mu**2*temp1; temp2_all += mu**2*temp2
				tmp1[6] += mu**2*temp1; tmp2[6] += mu**2*temp2
			if isPrint:
				if file2[:2] != 'H0': print '%.12e'%temp1, '%.12e'%temp2, '%.1e'%fabs((temp1-temp2)/temp2), length1
				else:  				  print '%.12e'%temp1, '%.12e'%temp2, '%.1e'%0, length1
	print '\n'
	for i in range(7):
		if isPrint: print i+1, '%.12e'%tmp1[i], '%.12e'%tmp2[i], '%.1e'%fabs((tmp1[i]-tmp2[i])/tmp2[i])
	h0_1 = temp1_h0 + gmass[0]*mu*temp1_all
	h0_2 = temp2_h0 + gmass[0]*mu*temp2_all
	prec = fabs((h0_2-h0_1)/h0_2)
	if isPrint: print '%.12e'%temp1_all, '%.12e'%temp2_all, '%.1e'%fabs((temp1_all-temp2_all)/temp2_all), length1_all
	print '%.12e'%h0_1, '%.12e'%h0_2, '%.1e'%prec, length1_all+4
	return h0_1, h0_2, prec
##################################################
def loop(elmass, N = 1, isPrint = True):
	"""
	Checking accuracy of the hamiltonian in the loop.
	"""
	elements, gmass, mu = elmass
	oe, series2, RESULT = [], [], []
	if N > 1:
		ln = len(elements)
		for i in range(ln):		# loop for planets
			temp = []
			for M in [2*j*pi/N for j in range(0, N)]:		# loop for mean anomalies
				temp.append([elements[i][0], elements[i][1], elements[i][2], elements[i][3], elements[i][4], M])
			oe.append(temp)
		KEYS = [str(i)+str(j)+str(k)+str(m) for i in range(N) for j in range(N) for k in range(N) for m in range(N)]
		ELEMS = {key: [oe[i][int(key[i])] for i in range(ln)] for key in KEYS}
	else:
		ELEMS = {'elem': elements}
	for file2 in os.listdir(PREFIX + 'HAM/'):
		series1 = epst()
		lf(series1, PREFIX + 'HAM/'+file2, df.boost_portable, cf.bzip2)
		length1 = sum([len(item.list) for item in [jtem[0].list[0][0] for jtem in series1.list]])
		series2.append([file2, series1, length1])
	for KEY in ELEMS.keys():
		print '\n'+KEY+':'
		result = check_items([ELEMS[KEY], gmass, mu], series2, isPrint)
		RESULT.append([KEY, result])
		if result[2] >= 1e-10: print KEY, result[0], '%.1e'%result[2]
	return RESULT
##################################################
def cosine_save(p, q, n, deg):
	"""
	Saving of cosine's expansion.
	"""
	(p, q) = (q, p) if p > q else (p, q)
	PATH_COS = PREFIX + 'COS/'+str(n)+'/'+str(deg)+'/'
	if not os.path.exists(PATH_COS): os.makedirs(PATH_COS)
	sf(cosine(p, q, n, deg), PATH_COS + 'cosine_'+str(q)+str(p)+'.epst.boostp.bz2', df.boost_portable, cf.bzip2)
##################################################

