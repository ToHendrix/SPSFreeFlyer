# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 14:44:10 2019

@author: Josh
"""

import numpy as np

#Material Properties
#Aluminimum 7075-T73 | 
#Mass Density
Rho = [2.8*10**3] 
#Tensile Ultimate Strength
F_tu = [460*10**6] 
#Tensile Yield Strength
F_ty = [390*10**6] 
#Young's Modulus
E = [71*10**9] 



def Thickness_Ring(Rho,F_tu,F_ty,E):    
    #Given Requirements
    #Axial frequency limit
    fn_axial = 75 
    #Lateral frequency limit
    fn_lateral = 20 
    #Mass of SPS structure
    m_B = 250 
    #Tip Mass (PPL)
    Tp = 1500 
    #Poisson Ratio
    v = 0.33 
    #Diameter of the ring
    D = 0.937 
    #Height of the Ring  
    L = 0.455 
    
    #Axial Rigidity Eq.(11-48). Assumption that it is a complex beam.
    A_axial = ((fn_axial/0.16)**2*(Tp*L+0.333*m_B*L))/(E)
    t_axial = A_axial/(2*np.pi*(D/2))
    
    
    
    #Lateral Rigidity Eq.(11-48). Assumption that it is a complex beam.    
    
    I_lat = ((fn_lateral/0.276)**2*(Tp*L**3+0.236*m_B*L**3))/(E)
    t_lat = I_lat/(np.pi*(D/2)**3)
    
    return I_lat, t_lat, A_axial, t_axial


Aluminimum = Thickness_Ring(Rho[0],F_tu[0],F_ty[0],E[0])