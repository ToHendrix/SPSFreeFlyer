# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 14:36:49 2019

@author: tom-h_000
"""

from Magnetorquer_calc import get_Js
from ADCS_combined import tau_build

def desaturate():
    H = [1,
          1,
          1]                                                                #Nms 
    H_torquer = [0,0,0]
    second = 0
    cycle = 0
    while (H[0] > H_torquer[0]) or (H[1] > H_torquer[1])  or (H[2] > H_torquer[2]):
        if (cycle < 11):
            torquer_torque = get_Js(15, second)
        
            for i in range(3):
                H_torquer[i] += ((torquer_torque[i] * 60) + tau_build[second][i])
        elif (cycle == 11):
            for i in range(3):
                H_torquer[i] +=  tau_build[second][i] 
            cycle = 0

        second +=1
        cycle += 1
        
    return second


print(desaturate())
    
    