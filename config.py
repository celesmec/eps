# -*- coding: utf-8 -*-
##################################################
# Copyright (C) 2013-2018 by Alexander Perminov
# perminov12@yandex.ru
#
##################################################
# THAT ARE GLOBAL VARIABLES AND FUNCTIONS
##################################################
from math import pi
from fractions import Fraction as F
from eps import pt, epst, pd, epsd, cf, df, sf, lf

NUM_PL = 4
PREFIX = ''
if PREFIX == '': raw_input('Please, input path to work directory: ')

pq_list = [item+str(i+1) for item in 'Lqxyuv' for i in range(NUM_PL)]
qp_list = [item+str(i+1) for item in 'qLyxvu' for i in range(NUM_PL)]
names_ke = [item+str(i+1) for i in range(NUM_PL) for item in ['a', 'e', 'i', 'om', 'Om']]
names_pe = [item+str(i+1) for i in range(NUM_PL) for item in ['L', 'x', 'y', 'u', 'v']]	

