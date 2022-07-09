#!/usr/bin/env python

from src.analysis import analysis

import numpy

this = analysis()
this.create()

print "... testing exposed likelihood calculated in C++"

seed = numpy.array([1.6, 1.0, 0.0, 2.2, 2.67])
pars = seed 
print "lnprob (seed):", this.get_lnprob(pars)
