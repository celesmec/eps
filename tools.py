# -*- coding: utf-8 -*-
##################################################
# Copyright (C) 2013-2018 by Alexander Perminov
# perminov12@yandex.ru
#
##################################################
	"""
	Converting JD to day, month, year.
	"""
	JDZ = int(JD)
	JDQ = JD - int(JD)
	# --------------------------------
	a = JDZ + 32044
	b = (4*a + 3)/146097
	c = a - 146097*b/4
	d = (4*c + 3)/1461
	e = c - 1461*d/4
	m = (5*e + 2)/153
	# --------------------------------
	day = e - (153*m + 2)/5 + 1
	month = m + 3 - 12*(m/10)
	year = 100*b + d - 4800 + m/10
	hourq = JDQ*24 + 12
	hour = int(hourq)
	minq = (hourq - hour)*60
	min = int(minq)
	secq = (minq - min)*60
	sec = int(secq)
	# --------------------------------
	m30 = [4, 6, 9, 11]
	m31 = [1, 3, 5, 7, 8, 10, 12]
	if hour >= 24:
		hour = hour - 24
		if month in m30:
			if day == 30: day, month = 1, month + 1
			else: day = day + 1
		elif month in m31:
			if day == 31:
				if month == 12: day, month, year = 1, 1, year + 1
				else: day, month = 1, month + 1
			else: day = day + 1
		elif month == 2:
			year = int(year)
			if year%400 == 0:	isLeapYear = True
			elif year%100 == 0:	isLeapYear = False
			elif year%4 == 0:	isLeapYear = True
			else:				isLeapYear = False
			if isLeapYear:
				if day == 29: day, month = 1, month + 1
				else: day = day + 1
			else: day = day + 1
	return year, month, day, hour, min, secq
def kepler2poincare(elem_k, kM):
	"""
	Transform Kepler orbital elements to Poincare elements.
	- 'kM' is sqrt(kappa2) * Mass,
	- 'aa' is semi-major axis,
	- 'ee'  is eccentricity,
	- 'ii'  is inclination,
	- 'om'  is argument of pericenter,
	- 'OM'  is longitude of ascending node,
	- 'MM'  is mean anomaly.
	"""
	from math import sin, cos, sqrt, fmod, pi
	aa, ee, ii, om, OM, MM = elem_k
	e2 = sqrt(1. - ee**2)
	LL = kM*sqrt(aa)
	PP = sqrt(2*LL * (1. - e2))
	QQ = sqrt(2*LL * e2 * (1. - cos(ii)))
	xx =  PP * cos(om + OM)
	yy = -PP * sin(om + OM)
	uu =  QQ * cos(OM)
	vv = -QQ * sin(OM)
	ll = om + OM + MM
	return [LL, xx, yy, uu, vv, ll]
##################################################
def poincare2kepler(elem_p, kM):
	"""
	Transform Poincare elements to Kepler orbital elements.
	"""
	from math import acos, atan2, sqrt, fmod, pi
	LL, xx, yy, uu, vv, ll = elem_p
	aa = (LL/kM)**2
	ee = sqrt(1. - (1. - 0.5*(xx**2 + yy**2)/LL)**2)
	ii = acos(1. - 0.5*(uu**2 + vv**2)/LL/sqrt(1. - ee**2))
	OM = atan2(-vv, uu)
	UU = atan2(-yy, xx)
	om = UU - OM
	MM = ll - om - OM
	return [aa, ee, ii, om, OM, MM]
##################################################
def cart2kepler(cart, k2):
	"""
	Calculation values of Kepler orbital elements.
	- 'k2' is gravitational parameter,
	- 'x, y, z' are coordinates,
	- 'vx, vy, vz' are velocities.
	"""
	from math import acos, atan, atan2, sin, cos, tan, sqrt
	x, y, z, vx, vy, vz = cart
	r = sqrt(x**2 + y**2 + z**2)
	v = sqrt(vx**2 + vy**2 + vz**2)
	a = (2./r - v**2/k2)**-1
	if a < 0: a = -1*a
	if a > 0:
		s = x*vx + y*vy + z*vz
		sx = y*vz - z*vy
		sy = x*vz - z*vx
		sz = x*vy - y*vx
		p = (sx**2 + sy**2 + sz**2)/k2
		ra = 1. - r/a
		e = sqrt(s**2/(k2*a) + ra**2)
		if e > 1: e = 0.
		if z == 0.: i = 0.
		else: i = acos(sz/sqrt(k2*p))
		if i == 0:
			i = 1e-14
		OM = atan2(sx, sy)
		th = atan2(sqrt(p)*s, sqrt(k2)*(p - r))
		uu = atan2(z/sin(i), x*cos(OM) + y*sin(OM))
		om = uu - th
		EE = 2.*atan(sqrt((1. - e)/(1. + e))*tan(th/2.))
		MM = EE - e*sin(EE)
		return [a, e, i, om, OM, MM]
	else:
		print a, 'WARNING: hyperbolic or parabolic orbit!'
		return [a, 0, 0, 0, 0, 0]
##################################################
def kepler2cart(elem_k, k2):
	"""
	Calculation values of x, y, z.
	- 'k2' is gravitational parameter,
	- 'aa' is semi-major axis,
	- 'ee', 'ii', 'om', 'OM', 'MM' see above.
	"""
	from math import atan, tan, sin, cos, sqrt, fabs
	aa, ee, ii, om, OM, MM = elem_k
	EE, eps = MM, 1000
	while eps >= 1e-14:
		EM = MM + ee*sin(EE)
		eps = fabs(EM - EE)
		EE = EM
	th = 2.*atan(sqrt((1. + ee)/(1. - ee))*tan(EE/2.))
	uu = th + om
	pp = aa*(1. - ee**2)
	rr = pp/(1. + ee*cos(th))
	xx = rr * (cos(uu)*cos(OM) - sin(uu)*sin(OM)*cos(ii))
	yy = rr * (cos(uu)*sin(OM) + sin(uu)*cos(OM)*cos(ii))
	zz = rr * sin(uu)*sin(ii)
	vr = sqrt(k2/pp)*ee*sin(th)
	vn = sqrt(k2/pp)*(1. + ee*cos(th))
	vx = xx*vr/rr + (-sin(uu)*cos(OM) - cos(uu)*sin(OM)*cos(ii)) * vn
	vy = yy*vr/rr + (-sin(uu)*sin(OM) + cos(uu)*cos(OM)*cos(ii)) * vn
	vz = zz*vr/rr + cos(uu)*sin(ii) * vn
	return [xx, yy, zz, vx, vy, vz]
##################################################
def barySolar(cart, mass):
	"""
	Calculation barycentric Cartesian coordinates of the Sun.
	- 'mass' is list of planets masses (by Solar mass);
	- 'cart' is list of Cartesian coordinates (or velocities).
	"""
	N = len(mass)
	solc = [sum([-mass[k]*cart[k][j] for k in range(1, N)]) for j in range(3)]
	return solc
##################################################
def cart2jacobi(cart, mass):
	"""
	Transform barycentric Cartesian coordinates to Jacobi coordinates.
	- 'mass' is list of planets masses (by Solar mass);
	- 'cart' is list of Cartesian coordinates (or velocities).
	"""
	N = len(mass)
	sum_m = [sum([mass[i] for i in range(n+1)]) for n in range(N)]
	jc = []
	jc.append([sum([(mass[k]/sum_m[N-1])*cart[k][j] for k in range(1, N)], cart[0][j]/sum_m[N-1]) for j in range(3)])
	for i in range(1, N):
		jc.append([sum([-(mass[k]/sum_m[i-1])*cart[k][j] for k in range(1, i)], cart[i][j]-cart[0][j]/sum_m[i-1]) for j in range(3)])
	return jc
##################################################
def jacobi2cart(jc, mass):
	"""
	Transform Jacobi coordinates to barycentric Cartesian coordinates.
	- 'mass' is list of planets masses (by Solar mass);
	- 'jc' is list of Jacobi coordinates (or velocities).
	"""
	N = len(mass)
	sum_m = [sum([mass[i] for i in range(n+1)]) for n in range(N)]
	cart = []
	cart.append([sum([-(mass[k]/sum_m[k])*jc[k][j] for k in range(1, N)], jc[0][j]) for j in range(3)])
	for i in range(1, N):
		cart.append([sum([-(mass[k]/sum_m[k])*jc[k][j] for k in range(i+1, N)], jc[0][j]+(sum_m[i-1]/sum_m[i])*jc[i][j]) for j in range(3)])
	return cart
##################################################
def jacobiAMass(gmass, mu):
	"""
	Calculation mass-parameters in Jacobi coordinates.
	"""
	rlen = range(len(gmass))
	mp = [gmass[i]/gmass[0]/mu for i in rlen]
	sum_mp = [1.+mu*sum(mp[1:i+1]) for i in rlen]
	gmj = [gmass[0] if i == 0 else gmass[0]*sum_mp[i]/sum_mp[i-1] for i in rlen]
	mpj = [1. if i == 0 else gmass[0]*mp[i]/gmj[i] for i in rlen]
	km = [mpj[i]*gmj[i]**0.5 for i in rlen]
	return gmj[1:], mpj[1:], km[1:]
##################################################

