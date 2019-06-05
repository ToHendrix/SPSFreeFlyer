# -*- coding: utf-8 -*-
"""
Created on Wed May 15 11:03:17 2019

@author: Willem van Lynden
"""
"Import packages"
import numpy as np
import scipy.special as special
import matplotlib as plt

"Set values"
pi = np.pi
sqrt = np.sqrt
asin = np.arcsin
acos = np.arccos
cos = np.cos

Re = 6378.137*10**3
ph = 350*10**3
ah = 850*10**3
orbitnr = [0]
nr = 0.
mu = 3.986004415*10**14
Cd = 2.2
m = 250.
A = 1.

"Calculate other starting values"
perigee = ph + Re
apogee = ah + Re
a = (perigee + apogee)/2
e = 1. - perigee/a
orbit = [a]
ballisticcoefficient = (m/(Cd*A))
print e
print a

T = 2*pi*np.sqrt((a**3)/(mu))
v = np.sqrt(mu*((2/(a*(1-e)))-1/a))
theta = asin(Re/a)
Teclipse = T*theta/(2*pi)
Tday = T-Teclipse
print 'period is', T
print 'perigee velocity is', v
print 'eclipse duration is', Teclipse

print 'decay after this'
"density"
r = a*(1-e)
altitude = r-Re
if altitude > 400000.:
    print 'out of orbit range of 350 km altitude perigee'
elif altitude > 350000.1:
    rho0 = 6.66*10**-12
    starth = 350000.
    h0 = 54800.
elif altitude > 300000.1:
    rho0 = 1.87*10**-11
    starth = 300000.
    h0 = 50300.
elif altitude > 250000.1:
    rho0 = 5.97*10**-11
    starth = 250000.
    h0 = 44800.
elif altitude > 200000.1:
    rho0 = 2.41*10**-10
    starth = 200000.
    h0 = 37500.
elif altitude > 150000.1:
    rho0 = 1.73*10**-9
    starth = 150000.
    h0 = 25500.
elif altitude > 100000.1:
    rho0 = 5.25*10**-7
    starth = 100000.
    h0 = 5900.
else:
    rho0 = 1.2
    starth = 0.
    h0 = 8.4

dh = altitude - starth
print 'altitude is',altitude
print 'height difference is',dh
rhop = rho0*np.exp(-dh/h0)
print 'rhop is',rhop
c = a*e/h0
print 'a is',a
da = -2*pi*(Cd*A/m)*(a**2)*rhop*np.exp(-c)*(special.iv(0,c) + 2*e*special.iv(1,c))
print 'da is',da
de = -2*pi*(Cd*A/m)*a*rhop*np.exp(-c)*(special.iv(1,c)+e*(special.iv(0,c)+special.iv(2,c))/2)
a = a+da
e = e+de
orbit.append(a)
nr += 1
orbitnr.append(nr)
L = -503000/da
print L
print 5801.23178811*L/(60*60*24*356.25)
#"Steps"
#while a*(1-e) > Re:
#    dh = a*(1-e)-Re
#    rhop = rho0*np.exp(-dh/h0)
#    c = a*e/h0
#    da = -2*pi*(Cd*A/m)*(a**2)*rhop*np.exp(-c)*(special.iv(0,c) + 2*e*special.iv(1,c))
#    de = -2*pi*(Cd*A/m)*a*rhop*np.exp(-c)*(special.iv(1,c)+e*(special.iv(0,c)+special.iv(2,c))/2)
#    a = a+da
#    e = e+de
#    orbit.append(a)
#    nr += 1
#    orbitnr.append(nr)
#print orbit
#plt.plot(orbitnr,orbit)
#plt.show