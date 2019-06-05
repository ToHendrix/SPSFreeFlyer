# -*- coding: utf-8 -*-
"""
Created on Wed May 15 16:56:24 2019

@author: Willem van Lynden
"""

"Sun-Synchronous Orbit Calculator"

import numpy as np
pi = np.pi
sqrt = np.sqrt
asin = np.arcsin
acos = np.arccos
cos = np.cos

Re = 6378.137*10**3
ph = 350*10**3.
ah = 850*10**3.
mu = 3.986004415*10**14

"Calculations"
perigee = ph + Re
apogee = ah + Re
a = (perigee + apogee)/2
e = 1. - perigee/a

imax = 108.
imin = 98.
incl = (180/pi)*acos(-((a/12352000)**(7./2.)))
print incl