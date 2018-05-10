# -*- coding: utf-8 -*-
##################################################
# Copyright (C) 2013-2018 by Alexander Perminov
# perminov12@yandex.ru
#
##################################################
# THIS IS ALGORITHM FOR HORI-DEPRIT METHOD
##################################################
import os, time
from eps.config import *
##################################################
# CONSTRUCTION OF AVERAGED MOTION EQUATIONS
##################################################
def ham_1(isAver = True):
	"""
	"""
	PATH_PHI = PREFIX+'HAM/'
	NAME = ['PHI', 'HAM', 'TRG', 'INT']
	for file2 in os.listdir(PATH_PHI):
		order = int(file2[1])
		if order > 0:
			start = time.time()
			series = {NAME[i]: epst(0) for i in range(4)}
			lf(series['PHI'], PATH_PHI + file2, df.boost_portable, cf.bzip2)
			if file2[0] == 'H':
				for item in series['PHI'].list:
					if item[1] == epst(1):
						series['HAM'] = item[0]*item[1]
				series['TRG'] = series['PHI'] - series['HAM']
				series['INT'] = series['TRG'].t_integrate()
			length = {}
			for key in series.keys():
				key2 = 'HAM' if (file2[0] == 'A' and key == 'PHI') else key
				PATH_SAVE = PREFIX + 'HDM/'+key2+str(order)+'/H'+str(order)+'/'
				if not os.path.exists(PATH_SAVE): os.makedirs(PATH_SAVE)
				if isAver:
					if file2[0] == 'A':
						is_write = True if key == 'PHI' else False
					else:
						is_write = True if key != 'HAM' else False
				else:
					is_write = True if file2[0] == 'H' else False
				if is_write:
					sf(series[key], PATH_SAVE + key2[0].lower()+str(order)+'_'+file2[3:5]+'.epst.boostp.bz2', df.boost_portable, cf.bzip2)
				length.update({key: sum([len(item.list) for item in [jtem[0].list[0][0] for jtem in series[key].list]])})
			print file2[:-16].ljust(16, ' '),
			for i in range(4): print NAME[i]+':', str(length[NAME[i]]).rjust(8, ' ')+',',
			print str('%.1f'%(time.time()-start)).rjust(4, ' '), 's'
##################################################
def ham_s(max_deg, order = 2, num_item = 1, t_order_max = 15, bracket = '', prefix = '', save_bracket = True):
	"""
	"""
	from pyranha.math import truncate_degree
	if order == 1:
		PATH, index = PREFIX + 'HDM/HAM1/H1/', ''
	if order > 1:
		pbrackets = ['T1H1'] if order == 2 else ['T1H2', 'T2H1_H2', 'T2H1_T1A1', 'T2H1_T1H1', 'T1T1H1', 'T1T1A1']
		if num_item == 1: PATH, index = PREFIX + 'HDM/HAM'+str(order)+'/H'+str(order)+'/', '_'+str(order)
		if num_item > 1:
			PATH = PREFIX + 'HDM/HAM'+str(order)+'/'+pbrackets[num_item-2]+'_TOR'+str(t_order_max).rjust(2, '0')+'/'
			index = '_'+pbrackets[num_item-2].lower()+'_tor'+str(t_order_max).rjust(2, '0')
	PATH_SAVE = PREFIX + 'HDM/SUMS/HAM/HAM'+str(order)+prefix+'/'
	if not os.path.exists(PATH_SAVE): os.makedirs(PATH_SAVE)
	print PATH[len(PREFIX):], '->', PATH_SAVE[len(PREFIX):]
	list_dirs = os.listdir(PATH) if num_item != 1 else ['']
	series = epst(0)
	for subdir in list_dirs:
		series_temp = epst(0)
		for file2 in os.listdir(PATH + subdir+'/'):
			temp = epst(0)
			if bracket == '':
				lf(temp, PATH + subdir+'/'+file2, df.boost_portable, cf.bzip2)
				if max_deg > 0: temp = truncate_degree(temp, max_deg, pq_list[8:])
				if save_bracket: series_temp += temp
				else:			 series 	 += temp
				print '>', file2[:-16]
			if file2[4:9] == bracket:
				lf(temp, PATH + subdir+'/'+file2, df.boost_portable, cf.bzip2)
				if max_deg > 0: temp = truncate_degree(temp, max_deg, pq_list[8:])
				series_temp += temp
				print '>', file2[:-16]
		if bracket != '':
			sf(series_temp.trim(), PATH_SAVE + 'h'+str(order)+index+'_('+subdir+')_'+bracket+'.epst.boostp.bz2', df.boost_portable, cf.bzip2)
		if save_bracket and bracket == '':
			sf(series_temp.trim(), PATH_SAVE + 'h'+str(order)+index+'_('+subdir+').epst.boostp.bz2', df.boost_portable, cf.bzip2)
	if (not save_bracket) and (bracket == ''):
		sf(series.trim(), PATH_SAVE + 'h'+str(order)+index+'.epst.boostp.bz2', df.boost_portable, cf.bzip2)
##################################################
def eqn_s(max_deg, order = 2, num_item = 1, t_order_max = 15, prefix = ''):
	"""
	"""
	from pyranha.math import partial, truncate_degree
	if order == 1:
		index = ''
	if order > 1:
		pbrackets = ['t1h1'] if order == 2 else ['t1h2', 't2h1_h2', 't2h1_t1a1', 't2h1_t1h1', 't1t1h1', 't1t1a1']
		if num_item == 1: index = '_'+str(order)
		if num_item > 1: index = '_'+pbrackets[num_item-2]+'_tor'+str(t_order_max).rjust(2, '0')
	PATH = PREFIX + 'HDM/SUMS/HAM/HAM'+str(order)+prefix+'/'
	for file2 in os.listdir(PATH):
		if index in file2 and file2[0] == 'h':
			series = epst(0)
			lf(series, PATH + file2, df.boost_portable, cf.bzip2)
			for i in range(8, len(pq_list)):
				motion = -partial(series, pq_list[i]) if (pq_list[i][0] in 'qyv') else partial(series, pq_list[i])
				if max_deg > 0: motion = truncate_degree(motion, max_deg, pq_list[8:])
				PATH_SAVE = PREFIX + 'HDM/SUMS/EQN/EQN'+str(order)+prefix+'/'
				sf(motion.trim(), PATH_SAVE + 'e'+file2[1:-16]+'_'+qp_list[i]+'.epst.boostp.bz2', df.boost_portable, cf.bzip2)
##################################################
# CONSTRUCTION OF FUNCTIONS FOR THE CHANGE
##################################################
def int_N(N, num_item = 1, t_order_max = 15, isDouble = False):
	"""
	"""
	BRACKETS = ['T1H1', 'T1A1'] if N == 2 else ['T1H2', 'T1A2', 'T2H1_H2', 'T2H1_T1A1', 'T2H1_T1H1', 'T2A1_H2', 'T2A1_T1A1', 'T2A1_T1H1', 'T1T1H1', 'T1T1A1']
	PATH_PHI = PREFIX + 'HDM/PHI'+str(N)+'/'+BRACKETS[num_item-1]+'_TOR'+str(t_order_max).rjust(2, '0')
	NAME = ['PHI', 'HAM', 'TRG', 'INT']
	(epst2, type2) = (epsd, 'd') if isDouble else (epst, 't')
	for subdir in os.listdir(PATH_PHI):
		for file2 in os.listdir(PATH_PHI + subdir+'/'):
			start = time.time()
			order = int(file2[1])
			series = {NAME[i]: epst2(0) for i in range(4)}
			lf(series['PHI'], PATH_PHI + subdir+'/'+file2, df.boost_portable, cf.bzip2)
			for item in series['PHI'].list:
				if item[1] == epst2(1):
					series['HAM'] = item[0]*item[1]
			series['TRG'] = series['PHI'] - series['HAM']
			series['INT'] = series['TRG'].t_integrate()
			PATH_SAVE = PREFIX + 'HDM/INT2/'+PHI_DIRS+'/'+subdir+'/'
			if not os.path.exists(PATH_SAVE): os.makedirs(PATH_SAVE)
			sf(series['INT'], PATH_SAVE + 'i'+str(order)+file2[2:-16]+'.eps'+type2+'.boostp.bz2', df.boost_portable, cf.bzip2)
			print subdir+'/'+file2[:-16].ljust(12, ' '), str('%.1f'%(time.time()-start)).rjust(4, ' '), 's'
##################################################
def chn_N(N, num_item = 1, t_order_max = 15, isDouble = False):
	"""
	"""
	from pyranha.math import partial
	if N == 1:
		PATH_INT = [PREFIX + 'HDM/INT'+str(i+1)+'/H'+str(i+1)+'/' for i in range(3)]
	else:
		BRACKETS = ['T1H1', 'T1A1'] if N == 2 else ['T1H2', 'T1A2', 'T2H1_H2', 'T2H1_T1A1', 'T2H1_T1H1', 'T2A1_H2', 'T2A1_T1A1', 'T2A1_T1H1', 'T1T1H1', 'T1T1A1']
		PATH_INT = PREFIX + 'HDM/INT'+str(N)+'/'+BRACKETS[num_item-1]+'_TOR'+str(t_order_max).rjust(2, '0')
	(type2, epst2) = ('d', epsd) if isDouble else ('t', epst)
	for j in range(len(PATHS_INT)):
		list_dir = os.listdir(PATHS_INT[j]) if N > 1 else ['']
		for subdir in list_dir:
			for file2 in os.listdir(PATHS_INT[j] + subdir+'/'):
				start = time.time()
				order = int(file2[1])
				series = epst2(0)
				lf(series, PATHS_INT[j] + subdir+'/'+file2, df.boost_portable, cf.bzip2)
				for i in range(len(pq_list)):
					if i in range(4):
						s_prev = epst2('s'+str(i)) if i != 0 else 1
						part_diff_nu = -3 * epst2('K0')**2 * (epst2('m'+str(i+1))**3 * s_prev * epst2('s'+str(i+1))**-1) * epst2('L'+str(i+1))**-4
						epst2.register_custom_derivative(pq_list[i], lambda temp: temp.partial(pq_list[i]) + temp.partial(r'\nu_{'+qp_list[i]+r'}') * part_diff_nu)
						change = partial(series, pq_list[i])
						epst2.unregister_all_custom_derivatives()
					else:
						change = -partial(series, pq_list[i]) if (pq_list[i][0] in 'qyv') else partial(series, pq_list[i])
					if change != epst2(0):
						PATH_SAVE = PREFIX + 'HDM/CHN'+str(order)+'/'+INT_DIRS[j]+'/'+qp_list[i]+'/'+subdir+'/'
						if not os.path.exists(PATH_SAVE): os.makedirs(PATH_SAVE)
						sf(change.trim(), PATH_SAVE + 'c'+file2[1:-16]+'_'+pq_list[i]+'.eps'+type2+'.boostp.bz2', df.boost_portable, cf.bzip2)
				print PATHS_INT[j]+subdir+'/'+file2[:-16], str('%.1f'%(time.time()-start)).rjust(4, ' '), 's'
##################################################
# CONSTRUCTION OF POISSON BRACKETS
##################################################
def d_phi1():
	"""
	"""
	from pyranha.math import partial
	PATHS_PHI = [PREFIX + 'HDM/PHI'+str(i+1)+'/H'+str(i+1)+'/' for i in range(3)]
	for PATH_PHI in PATHS_PHI:
		if os.path.isdir(PATH_PHI):
			for file2 in os.listdir(PATH_PHI):
				start = time.time()
				order = int(file2[1])
				series = epst(0)
				lf(series, PATH_PHI + file2, df.boost_portable, cf.bzip2)
				for i in range(len(pq_list)):
					diff = partial(series, pq_list[i])
					if diff != epst(0):
						PATH_SAVE = PREFIX + 'HDM/DIF1/PHI'+str(order)+'/H'+str(order)+'/'+pq_list[i]+'/'
						if not os.path.exists(PATH_SAVE): os.makedirs(PATH_SAVE)
						sf(diff.trim(), PATH_SAVE + file2[:-16]+'_'+pq_list[i]+'.epst.boostp.bz2', df.boost_portable, cf.bzip2)
				print '/PHI'+str(order)+'/H'+str(order)+'/'+file2[:-16], str('%.1f'%(time.time()-start)).rjust(4, ' '), 's'
##################################################
def pb_t1hN(N = 1, m = 6, t_order_max = 20, isAver = True, isDouble = False, isPrint = True):
	"""
	"""
	from pyranha.math import cos, sin, degree, ldegree, t_degree, t_ldegree, t_order, t_lorder, truncate_degree
	PATH_CHN = PREFIX + 'HDM/CHN1/H1/'
	PATH_PHI = PREFIX + 'HDM/DIF1/PHI'+str(N)+'/H'+str(N)+'/'
	if N == 1: c =  0.5 if isDouble else F( 1,2)
	if N == 2: c = -0.5 if isDouble else F(-1,2)
	(one, max_m, type2, pt2, epst2) = (1., m, 'd', pd, epsd) if isDouble else (1, F(m), 't', pt, epst)
	for i in range(len(pq_list)):
		for file1 in os.listdir(PATH_CHN + qp_list[i]):
			diff1 = epst(0)
			lf(diff1, PATH_CHN + qp_list[i]+'/'+file1, df.boost_portable, cf.bzip2)
			diff1 = truncate_degree(diff1, F(m), pq_list[8:])
			new_diff1 = epst2(0)
			for item1 in diff1.list:
				if t_order(item1[1]) <= t_order_max:
					new_diff1 += one*(item1[0]*item1[1])
			if isPrint: print file1[:-16]+':', sum([len(item.list) for item in [jtem[0].list[0][0] for jtem in new_diff1.list]]), len(new_diff1.list)
			for file2 in os.listdir(PATH_PHI + qp_list[i]):
				start = time.time()
				diff2 = epst(0)
				lf(diff2, PATH_PHI + qp_list[i]+'/'+file2, df.boost_portable, cf.bzip2)
				diff2 = truncate_degree(diff2, F(m), pq_list[8:])
				new_diff2 = epst2(0)
				for item1 in diff2.list:
					if t_order(item1[1]) <= t_order_max:
						new_diff2 += one*(item1[0]*item1[1])
				if isPrint: print file2[:-16]+':', sum([len(item.list) for item in [jtem[0].list[0][0] for jtem in new_diff2.list]]), len(new_diff2.list)
				if isAver:
					count = 0
					new_result = epst2(0)
					pt2.set_auto_truncate_degree(max_m, pq_list[8:])
					for item in new_diff1.list:
						for jtem in new_diff2.list:
							trig = item[1]*jtem[1]
							for ktem in trig.list:
								if ktem[1] == 1:
									new_result = new_result + c*item[0]*jtem[0]*ktem[0]
						count += 1
						print count, item[1]
					pt2.unset_auto_truncate_degree()
				else:
					pt2.set_auto_truncate_degree(max_m, pq_list[8:])
					result = c*new_diff1*new_diff2
					pt2.unset_auto_truncate_degree()
				if not isAver:
					new_result = epst2(0)
					for item1 in result.list:
						if t_order(item1[1]) <= t_order_max:
							new_result += item1[0]*item1[1]
				if isPrint: print 'result:', sum([len(item.list) for item in [jtem[0].list[0][0] for jtem in new_result.list]]), len(new_result.list)
				if new_result != epst2(0):
					TYPE = 'HAM' if isAver else 'PHI'
					PATH_SAVE = PREFIX + 'HDM/'+TYPE+str(N+1)+'/T1H'+str(N)+'_TOR'+str(t_order_max).rjust(2, '0')+'/'+pq_list[i]+'_'+qp_list[i]+'/'
					if not os.path.exists(PATH_SAVE): os.makedirs(PATH_SAVE)
					sf(new_result.trim(), PATH_SAVE + '/'+TYPE[0].lower()+str(N+1)+'_'+file1[3:5]+'_'+file2[3:5]+'.eps'+type2+'.boostp.bz2', df.boost_portable, cf.bzip2)
				print 'bracket:', pq_list[i]+'_'+qp_list[i], file1[:-16], file2[:-16], str('%.1f'%(time.time()-start)).rjust(4, ' '), 's'
##################################################
def pb_t1a1(m = 6, t_order_max = 20, isDouble = False, isPrint = True):
	"""
	"""
	from pyranha.math import cos, sin, degree, ldegree, t_degree, t_ldegree, t_order, t_lorder, truncate_degree
	PATH_CHN = PREFIX + 'HDM/CHN1/H1/'
	PATH_EQN = PREFIX + 'HDM/EQN1/H1/'
	(c, one, max_m, type2, pt2, epst2) = (0.5, 1., m, 'd', pd, epsd) if isDouble else (F(1,2), 1, F(m), 't', pt, epst)
	for i in range(len(pq_list)):
		if pq_list[i][0] != 'L':
			for file1 in os.listdir(PATH_CHN + qp_list[i]):
				diff1 = epst(0)
				lf(diff1, PATH_CHN + qp_list[i]+'/'+file1, df.boost_portable, cf.bzip2)
				diff1 = truncate_degree(diff1, F(m), pq_list[8:])
				new_diff1 = epst2(0)
				for item1 in diff1.list:
					if t_order(item1[1]) <= t_order_max:
						new_diff1 += one*(item1[0]*item1[1])
				if isPrint: print file1[:-16]+':', sum([len(item.list) for item in [jtem[0].list[0][0] for jtem in new_diff1.list]]), len(new_diff1.list)
				for file2 in os.listdir(PATH_EQN + pq_list[i]):
					start = time.time()
					diff2 = epst(0)
					lf(diff2, PATH_EQN + pq_list[i]+'/'+file2, df.boost_portable, cf.bzip2)
					diff2 = one*diff2
					if isPrint: print file2[:-16]+':', sum([len(item.list) for item in [jtem[0].list[0][0] for jtem in diff2.list]]), len(diff2.list)
					pt2.set_auto_truncate_degree(max_m, pq_list[8:])
					result = c*new_diff1*diff2 if (pq_list[i][0] in 'qyv') else -c*new_diff1*diff2
					pt2.unset_auto_truncate_degree()
					if isPrint: print 'result:', sum([len(item.list) for item in [jtem[0].list[0][0] for jtem in result.list]]), len(result.list)
					if result != epst2(0):
						PATH_SAVE = PREFIX + 'HDM/PHI2/T1A1_TOR'+str(t_order_max).rjust(2, '0')+'/'+pq_list[i]+'_'+qp_list[i]+'/'
						if not os.path.exists(PATH_SAVE): os.makedirs(PATH_SAVE)
						sf(result.trim(), PATH_SAVE + '/p2'+'_'+file1[3:5]+'_'+file2[3:5]+'.eps'+type2+'.boostp.bz2', df.boost_portable, cf.bzip2)
					print 'bracket:', pq_list[i]+'_'+qp_list[i], file1[:-16], file2[:-16], str('%.1f'%(time.time()-start)).rjust(4, ' '), 's'
##################################################

