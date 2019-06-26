# -*- coding: utf-8 -*-
"""
Created on Tue May 14 09:54:11 2019

@author: Willem van Lynden
"""
import numpy as np
import matplotlib.pyplot as plt
plt.close('all')
pi = np.pi
sqrt = np.sqrt
asin = np.arcsin
acos = np.arccos
cos = np.cos
R = 6378137
h1 = 850*10**3
#h2 = 800*10**3
mu = 3.986004415*10**14
AU = 1.496*10**11
L = 3.85*10**26

period = []
eperiod = []
dperiod = []
orbit = []
angle = []
intensity = []
velocity = []
incl = []

#for i in h2-h1:
a = float(R)+float(h1)
a = 6978137.0
T = 2*pi*sqrt((a**3)/(mu))
v = sqrt(mu/a)
theta = asin(R/a)
Teclipse = T*theta/(2*pi)
Tday = T-Teclipse
d = a*cos(theta)+AU
S = L/(4*pi*(d**2))
#period.append(T)
#eperiod.append(Teclipse)
#dperiod.append(Tday)
akm = a/1000
#orbit.append(akm)
#angle.append(theta)
#intensity.append(S)
#velocity.append(v)
ic = (180/pi)*-acos(a/12352000)**(7/2)
#incl.append(ic)
#fig, ax1 = plt.subplots()
#
#color = 'red'
#ax1.set_xlabel('semi-major axis [km]')
#ax1.set_ylabel('Orbital Period [s]', color=color)
#ax1.plot(orbit, period, color=color)
#ax1.tick_params(axis='y', labelcolor=color)
#
#ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
#
#color = 'blue'
#ax2.set_ylabel('Eclipse Duration [s]', color=color)  # we already handled the x-label with ax1
#ax2.plot(orbit, eperiod, color=color)
#ax2.tick_params(axis='y', labelcolor=color)
#
#fig.tight_layout()  # otherwise the right y-label is slightly clipped
#
#fig,ax3 = plt.subplots()
#ax3.set_xlabel('semi-major axis [km]')
#ax3.set_ylabel('Solar flux [W/m^2]')
#ax3.plot(orbit,intensity)
#
#fig,ax4 = plt.subplots()
#ax4.set_xlabel('semi-major axis [km]')
#ax4.set_ylabel('inclination [deg]')
#ax4.plot(orbit,incl)
#
#plt.show()
#
#print velocity[0]
#print velocity[350000-1]
print T
print v
print Teclipse
print Tday