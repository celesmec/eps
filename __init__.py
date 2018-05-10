# -*- coding: utf-8 -*-
########################################################################
# Copyright (C) 2013-2018 by Alexander Perminov:)
# perminov12@yandex.ru
#
# This program is free software. It uses CAS Piranha
# (F. Biscani, https://github.com/bluescarni/piranha)
########################################################################
# THIS IS PROGRAM FOR MODELING OF PLANETARY SYSTEMS ORBITAL EVOLUTION
########################################################################
# piranha file compression and data format:
from pyranha import data_format as df, compression as cf

# piranha save/load functions:
from pyranha import save_file as sf, load_file as lf

# piranha series types:
from pyranha.types import divisor, divisor_series, monomial, polynomial, poisson_series, int16, rational, double

pt = polynomial[rational,monomial[rational]]()
pd = polynomial[double,monomial[rational]]()
epst = poisson_series[divisor_series[polynomial[rational,monomial[rational]],divisor[int16]]]()
epsd = poisson_series[divisor_series[polynomial[double,monomial[rational]],divisor[int16]]]()

del divisor, divisor_series, monomial, polynomial, poisson_series, int16, rational, double

# package modules:
__name__ = 'eps'
__all__ = ['ham', 'hdm', 'keproc', 'tools', 'config', 'pt', 'epst', 'cf', 'df', 'lf', 'sf']
__version_major__ = '0'
__version_minor__ = '1'

# start:
print '-'*64
print 'Package ' + __name__ + ', version ' + __version_major__ + '.' + __version_minor__
print 'modules: ' + str(__all__[:5])[1:-1]
print 'types:   ' + str(__all__[5:-4])[1:-1]
print '-'*64

