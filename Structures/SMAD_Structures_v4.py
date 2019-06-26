# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 15:38:13 2019

@author: Josh Pho
"""


import numpy as np
import matplotlib.pyplot as plt

#Material Properties
#Aluminimum 7075-T6 | Titanium Ti-6AI-4V | Epoxy/HS carbon fiber, UD prepreg, QI lay-up
#Mass Density
Rho = [2810,4.43*10**3,1.565*10**3] 
#Tensile Ultimate Strength
F_tu = [572*10**6,1031*10**6,670.5*10**6] 
#Tensile Yield Strength
F_ty = [503*10**6,848*10**6,670.5*10**6] 
#Compression Strength
sigma_compression = [461.5*10**6,964*10**6,599.5*10**6]
#Shear Strength
sigma_shear = [331*10**6,550*10**6,70*10**6]
#Flexural Strength (Bending)
sigma_bending = [444.5*10**6,933*10**6,302.5*10**6]
#Young's Modulus
E = [71.7*10**9,114.5*10**9,54.9*10**9] 



def Thickness_Ring(Rho,F_tu,F_ty,sigma_compression,sigma_shear,sigma_bending,E):    
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
    #Gravity
    g = 9.81
    #Distance c.o.g. PPL and bottom separation plane
    d_cog_PPL = 1.9+L
    
    #Sizing for natural frequency requirement
    #Axial Rigidity Eq.(11-48). Assumption that it is a complex beam.
    A_axial = ((fn_axial/0.16)**2*(Tp*d_cog_PPL+0.333*m_B*d_cog_PPL))/(E)
    t_axial = A_axial/(2*np.pi*(D/2))
    
    #Lateral Rigidity Eq.(11-48). Assumption that it is a complex beam.    
    
    I_lat = ((fn_lateral/0.276)**2*(Tp*d_cog_PPL**3+0.236*m_B*d_cog_PPL**3))/(E)
    t_lat = I_lat/(np.pi*(D/2)**3)
    
    if t_axial > t_lat:
        t = t_axial
    else:
        t = t_lat
    
    print('Frequency Calculation',t)
    #Applied & Equivalant axial loads
    
    #Axial|Lateral 
    Load_Factors = [5.5,1.35]
    #Axial|Lateral|Bending Moment
    Limit_Loads = [(Tp+m_B)*g*Load_Factors[0],(Tp+m_B)*g*Load_Factors[1],Tp*g*Load_Factors[1]*(d_cog_PPL-L/2)+m_B*g*Load_Factors[1]*(L/2)]
    #Factor of Safety, |Ultimate FOS|Yield FOS|
    FOS = [1.25,1.10]
    
    #Equivalent Axial Load Eq.(11-42) (Combined Stress SMAD)
    P_eq = Limit_Loads[0]+(2*Limit_Loads[2])/(D/2)
    #Ultimate Load 
    P_ult = P_eq*FOS[0]
    #Yield load
    P_yield = P_eq*FOS[1]
    
    #Sizing for strength (Combined Stress SMAD)
    #Required thickness (F_tu)
    t_strength_ult = (P_ult)/(2*np.pi*(D/2)*F_tu)
    #Required thickness (F_ty)
    t_strength_yield = (P_yield)/(2*np.pi*(D/2)*F_ty)
    
    if t_strength_ult > t_strength_yield:
        t_strength = t_strength_ult
    else:
        t_strength = t_strength_yield
    
    if t_strength > t:
        t = t_strength
        
    print('Bending + Compression (SMAD)',t)
    
        #Sizing for stability Eq.(11-44)-(11-45) (Buckling stress)
    if (D/2)/t < 1500 and L/(D/2) < 5:
        phi = (1/16)*np.sqrt((D/2)/t)
        gamma = 1.0-0.901*(1-np.e**(-phi))
        #Cylinder Buckling Stress
        sigma_cr = 0.6*gamma*((E*t)/(D/2))
        #Crossectional Area
        A_cross = (np.pi*((D/2)+t)**2)-(np.pi*(D/2)**2)
        #print(A_cross,"A_cross")
        #Critical buckling load
        P_cr = A_cross*sigma_cr
    else:
        print("Inelastic buckling methods.")
        
    MS = (P_cr)/(P_ult)-1        
    
    print('Buckling Stress (Cylinder)',t)
    
    if MS<0.25:
        while MS <0.249 or MS > 0.251:
            if MS < 0.249:
                #Sizing for stability Eq.(11-44)-(11-45) (Buckling stress)
                if (D/2)/t < 1500 and L/(D/2) < 5:
                    phi = (1/16)*np.sqrt((D/2)/t)
                    gamma = 1.0-0.901*(1-np.e**(-phi))
                    #Cylinder Buckling Stress
                    sigma_cr = 0.6*gamma*((E*t)/(D/2))
                    #Crossectional Area
                    A_cross = (np.pi*((D/2)+t)**2)-(np.pi*(D/2)**2)
                    #Critical buckling load
                    P_cr = A_cross*sigma_cr
                else:
                    print("Inelastic buckling methods.")
                
                MS = (P_cr)/(P_ult)-1
                
                #print(P_cr,MS,t)
                t += 0.0000001                 
    
    print('Buckling Stress (Cylinder) - After Iteration',t)
      
    #Quasi-Static Loads, m_PPL
    #Compression | Tension | Lateral
    QSL = [80933,44145,19865]

    #Quasi-Static Loads, m_SPS2
    #Compression | Tension | Lateral
    QSL_SPS2 = [13488.75,7357.5,3310.875]    
    
    #Tension Yield Calculation
    sigma_ty = (QSL[1])/(2*np.pi*(D/2)*t)
    MS_ty = F_ty/(sigma_ty*FOS[1]) - 1
    
    #Tension Ultimate Calculation
    sigma_tu = (QSL[1])/(2*np.pi*(D/2)*t)
    MS_tu = F_tu/(sigma_tu*FOS[0]) - 1
    
    #Compression Yield Calculation
    sigma_cy = (QSL[0])/(2*np.pi*(D/2)*t)
    MS_cy = sigma_compression/(sigma_cy*FOS[1]) - 1
    
    #Bending Yield Calculation
    sigma_b = ((QSL[2]*(d_cog_PPL-L/2)+QSL_SPS2[2]*(L/2))*(D/2))/(np.pi*(D/2)**3*t)
    MS_bending = sigma_bending/(sigma_b*FOS[1]) - 1
    
    #Shear Yield Calculation
    sigma_s = (QSL[2])/(2*np.pi*(D/2)*t)
    MS_shear = sigma_shear/(sigma_s*FOS[1])
    
    #Bending + Tension Yield Calculation
    sigma_comb_bt = (sigma_b+sigma_ty)*FOS[1]
    MS_comb_bt = F_ty/sigma_comb_bt - 1 
    
    #Bending + Compression Yield Calculation
    sigma_comb_bc = (sigma_b+sigma_cy)*FOS[1]
    MS_comb_bc = F_ty/sigma_comb_bc - 1

    if MS_ty < 0.25 or MS_tu < 0.25 or MS_cy < 0.25 or MS_bending < 0.25 or MS_shear < 0.25 or MS_comb_bt < 0.25 or MS_comb_bc < 0.25:
        while MS_ty < 0.249 or MS_tu < 0.249 or MS_cy < 0.249 or MS_bending < 0.249 or MS_shear < 0.249 or MS_comb_bt < 0.249 or MS_comb_bc < 0.249:
            #Tension Yield Calculation
            sigma_ty = (QSL[1])/(2*np.pi*(D/2)*t)
            MS_ty = F_ty/(sigma_ty*FOS[1]) - 1
            
            #Tension Ultimate Calculation
            sigma_tu = (QSL[1])/(2*np.pi*(D/2)*t)
            MS_tu = F_tu/(sigma_tu*FOS[0]) - 1
            
            #Compression Yield Calculation
            sigma_cy = (QSL[0])/(2*np.pi*(D/2)*t)
            MS_cy = sigma_compression/(sigma_cy*FOS[1]) - 1
            
            #Bending Yield Calculation
            sigma_b = ((QSL[2]*(d_cog_PPL-L/2)+QSL_SPS2[2]*(L/2))*(D/2))/(np.pi*(D/2)**3*t)
            MS_bending = sigma_bending/(sigma_b*FOS[1]) - 1
            
            #Shear Yield Calculation
            sigma_s = (QSL[2])/(2*np.pi*(D/2)*t)
            MS_shear = sigma_shear/(sigma_s*FOS[1])
            
            #Bending + Tension Yield Calculation
            sigma_comb_bt = (sigma_b+sigma_ty)*FOS[1]
            MS_comb_bt = F_ty/sigma_comb_bt - 1 
            
            #Bending + Compression Yield Calculation
            sigma_comb_bc = (sigma_b+sigma_cy)*FOS[1]
            MS_comb_bc = F_ty/sigma_comb_bc - 1
            
            #print(MS_buckling,t)
            t += 0.000001
    else:
        #Tension Yield Calculation
        sigma_ty = (QSL[1])/(2*np.pi*(D/2)*t)
        MS_ty = F_ty/(sigma_ty*FOS[1]) - 1
        
        #Tension Ultimate Calculation
        sigma_tu = (QSL[1])/(2*np.pi*(D/2)*t)
        MS_tu = F_tu/(sigma_tu*FOS[0]) - 1
        
        #Compression Yield Calculation
        sigma_cy = (QSL[0])/(2*np.pi*(D/2)*t)
        MS_cy = sigma_compression/(sigma_cy*FOS[1]) - 1
        
        #Bending Yield Calculation
        sigma_b = (QSL[2]*(d_cog_PPL-L/2)+QSL_SPS2[2]*(L/2))/(np.pi*(D/2)**3*t)
        MS_bending = sigma_bending/(sigma_b*FOS[1]) - 1
        
        #Shear Yield Calculation
        sigma_s = (QSL[2])/(2*np.pi*(D/2)*t)
        MS_shear = sigma_shear/(sigma_s*FOS[1])
        
        #Bending + Tension Yield Calculation
        sigma_comb_bt = (sigma_b+sigma_ty)*FOS[1]
        MS_comb_bt = F_ty/sigma_comb_bt - 1 
        
        #Bending + Compression Yield Calculation
        sigma_comb_bc = (sigma_b+sigma_cy)*FOS[1]
        MS_comb_bc = F_ty/sigma_comb_bc - 1
        
    m = 2*np.pi*(D/2)*t*L*Rho
    print('Dutch Space Report (Airbus)',t)
    
    
    #--------Skin-Stringer------
    #Initial Variables Setup
    #ID=Number of stringers
    ID = 4
    ID_max = 12
    
    # |ID|m_skin|m_stiffeners|m_comb_skin_stiffeners|t_skin|MS_buckling_panel|A_stringers|A_skin|A_comb_stiffener_skin|Area
    Final_Stiffeners_Panel = []
    while ID<=ID_max:
        I = I_lat
        t_skin = t_axial
        I_skin = np.pi*(D/2)**3*t_skin
        I_stringers = I - I_skin
        Stiffener_Width = (2*np.pi*(D/2))/ID
        Area = 2*np.pi*(D/2)**3*t
        
        #Distribute Stringers evenly around circumference.
        def circle_points(r, n):
            circles = []
            for r, n in zip(r, n):
                t = np.linspace(0, 2*np.pi, n)
                x = r * np.cos(t)
                y = r * np.sin(t)
                circles.append(np.c_[x, y])
            return circles
    
        r = [D/2]
        n = [ID]
        circles = circle_points(r, n)
        y_loc = np.array(circles[0][:,1])
        #print(y_loc)
#        fig, ax = plt.subplots()
#        for circle in circles:
#            ax.scatter(circle[:, 0], circle[:, 1])
#        ax.set_aspect('equal')
#        plt.show()
        #print(circles)
        #print(y_loc)
        
        #AMOI calculation
        Parallel_Axis_Distance_Squared = y_loc**2
        #print(Parallel_Axis_Distance_Squared)
        Sum_Parallel_Axis_Distance_Squared = np.sum(Parallel_Axis_Distance_Squared)
        #print(Sum_Parallel_Axis_Distance_Squared*10**4)
        A_stringers = (I_stringers)/Sum_Parallel_Axis_Distance_Squared
        #print(A_stringers*10**4)
    
        #-------Panel Stability------
        
        #-------Calculate k in Eq(11-65)----
        x = np.array([1, 4,10,40, 100,400, 1000])
        y = np.array([4, 5,6,30, 70, 300,850])
        c = np.polyfit(x, y, 2)
        
        Z = ((Stiffener_Width**2)/((D/2)*t_skin))*(1-v**2)**0.5
        k = c[0]*Z**2+c[1]*Z**1+c[2]
        
        #Compressive Buckling Stress for Curved skin panel---- 
        sigma_cr_skin = (k*np.pi**2*E)/(12*(1-v**2))*(t_skin/Stiffener_Width)**2
        #print(sigma_cr_skin/10**6,Stiffener_Width)
        MS_buckling_panel = (sigma_cr_skin*Area)/(P_ult) - 1
        print(MS_buckling_panel)
        
        print(t_skin)
        #------Re-iterate-------
        if MS_buckling_panel > 0.25:
            print('test')
            while MS_buckling_panel > 0.251:
                I_skin = np.pi*(D/2)**3*t_skin
                I_stringers = I - I_skin
    
                A_stringers = (I_stringers)/Sum_Parallel_Axis_Distance_Squared
                #print(A_stringers*10**4)
            
                #-------Panel Stability------
                
                #-------Calculate k in Eq(11-65)----
                x = np.array([1, 4,10,40, 100,400, 1000])
                y = np.array([4, 5,6,30, 70, 300,850])
                c = np.polyfit(x, y, 2)
                
                Z = ((Stiffener_Width**2)/((D/2)*t_skin))*(1-v**2)**0.5
                k = c[0]*Z**2+c[1]*Z**1+c[2]
                
                #Compressive Buckling Stress for Curved skin panel----
                sigma_cr_skin = (k*np.pi**2*E)/(12*(1-v**2))*(t_skin/Stiffener_Width)**2
                #print(sigma_cr_skin/10**6)
                MS_buckling_panel = (sigma_cr_skin*Area)/(P_ult) - 1
                #print(MS_buckling_panel)
    
                A_skin = 2*np.pi*(D/2)*t_skin
                #print(A_skin)
                A_comb_stiffener_skin = A_stringers + A_skin
                #print(A_comb_stiffener_skin)
                
                #print(A_stringers*L*Rho*ID)
                if A_comb_stiffener_skin > Area:
                    break
                t_skin = t_skin - 0.00001
        
        elif MS_buckling_panel < 0.249:
            print('hoi')
            while MS_buckling_panel < 0.249:
                    I_skin = np.pi*(D/2)**3*t_skin
                    I_stringers = I - I_skin
        
                    A_stringers = (I_stringers)/Sum_Parallel_Axis_Distance_Squared
                    #print(A_stringers*10**4)
                
                    #-------Panel Stability------
                    
                    #-------Calculate k in Eq(11-65)----
                    x = np.array([1, 4,10,40, 100,400, 1000])
                    y = np.array([4, 5,6,30, 70, 300,850])
                    c = np.polyfit(x, y, 2)
                    
                    Z = ((Stiffener_Width**2)/((D/2)*t_skin))*(1-v**2)**0.5
                    k = c[0]*Z**2+c[1]*Z**1+c[2]
                    
                    #Compressive Buckling Stress for Curved skin panel----
                    sigma_cr_skin = (k*np.pi**2*E)/(12*(1-v**2))*(t_skin/Stiffener_Width)**2
                    #print(sigma_cr_skin/10**6)
                    MS_buckling_panel = (sigma_cr_skin*Area)/(P_ult) - 1
                    #print(MS_buckling_panel)
        
                    A_skin = 2*np.pi*(D/2)*t_skin
                    #print(A_skin)
                    A_comb_stiffener_skin = A_stringers + A_skin
                    #print(A_comb_stiffener_skin)
                    
                    #print(A_stringers*L*Rho*ID)
                    #if A_comb_stiffener_skin > Area:
                    #    break
                    t_skin = t_skin + 0.00001
        
        print(t_skin)
        #print(MS_buckling_panel)
        m_stiffeners = A_stringers*L*Rho*ID
        #print(m_stiffeners)
        m_skin = 2*np.pi*(D/2)*t_skin*L*Rho
        #print(m_skin)
        m_comb_skin_stiffeners = m_stiffeners + m_skin
        #print(m_comb_skin_stiffeners)
    
        Vector_Final_Stiffeners_Panel = [[ID,m_skin,m_stiffeners,m_comb_skin_stiffeners,t_skin,MS_buckling_panel,A_stringers,A_skin,A_comb_stiffener_skin,Area]]
        Final_Stiffeners_Panel.append(Vector_Final_Stiffeners_Panel)
        
        I = I_lat
        t_skin = t_axial
        
        
        ID += 1
    return m,t, MS, MS_ty,MS_tu, MS_cy,MS_bending,MS_shear,MS_comb_bt,MS_comb_bc,Final_Stiffeners_Panel

Aluminimum = Thickness_Ring(Rho[0],F_tu[0],F_ty[0],sigma_compression[0],sigma_shear[0],sigma_bending[0],E[0])
Titanium = Thickness_Ring(Rho[1],F_tu[1],F_ty[1],sigma_compression[1],sigma_shear[1],sigma_bending[1],E[1])
Composite = Thickness_Ring(Rho[2],F_tu[2],F_ty[2],sigma_compression[2],sigma_shear[2],sigma_bending[2],E[2])

    
    
    
    