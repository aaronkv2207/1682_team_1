import numpy as np
import matplotlib.pyplot as plt

#Estimated forces/loads we would need: 

# L_main = ... 
# L_ho = ... 
# L_vert = ... 

# span = ... 
# max_thick = ... 
# width = ...

# I = max_thick * (width**3)/ 12 #rectangle cross section beam 
# E = ... # material

class Beam: 
    """
    input Lift force, I of beam, E of material, span 

    output max stresses (axial and shear), max deflection 

    """
    def __init__(self, L, I, E, span): 
        self.L = L 
        self.I = I 
        self.E = E 
        self.span = span 

    def cant_u_max (self):  
        return (-(self.L*0.5)*(self.span*0.5)**4)/(8*self.E*self.I) #max deflection 
    
    def cant_moment_max (self): 
        return ((self.L*0.5)*(self.span*0.5)**2)/2 #at fixed end 
    
    def cant_stress_max (self): 
        return (self.L*0.5) * (self.span * 0.5) #at fixed end (x = span)
    
    def pointload_root_stress_max (self): #landing case ? 
        return self.L #look into calculating the L required in this landing case 
    
    def thrust_max_stress (self, propnum, propspacing): 


#canteliver beam eqations, evenly distributed load, span and L divided by 2 inside functions

# def cant_u (L, x, I, E, span): #deflection
#     return (-(L*0.5)*x**2)*(x**2 - 4*span*0.5*x+6*(span*0.5)**2)/(24*E*I)

# def cant_u_max (L, E, I, span): #max deflection 
#     return (-(L*0.5)*(span*0.5)**4)/(8*E*I)

# def cant_stress_max (L, span): 
#     return (L*0.5) * (span * 0.5)  #at fixed end (x = span)

# def cant_moment_max (L, span): 
#     return ((L*0.5)*(span*0.5)**2)/2 #at fixed end 

#landing case, weight, normal force, lift 



#takeoff case on wing 







