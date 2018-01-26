#
# Mualem - van Genuchten model, equations from van Genuchten, MT 1980
# 
import math

#
# returns the hydraulic conductivity according to the van genuchten model
#
def hydraulic_conductivity(h,sp):
    h = min(h,0)
    # compute the volumetric moisture content (Eqn 21)
    theta = sp.theta_R + (sp.theta_S-sp.theta_R)/pow(1. + pow(sp.alpha*abs(h),sp.n),sp.m)
    # compute the effective saturation (dimensionless water content, Eqn 2)
    se = (theta-sp.theta_R)/(sp.theta_S-sp.theta_R)
    # Compute the hydraulic conductivity (Eqn 8)
    K = sp.Ksat*math.sqrt(se)*( (1. - pow(1. - pow(se,1./sp.m),sp.m))**2 )
    return K 

#
# class containing the va genuchten parameters
#
class Parameters:
    def __init__(self, R, S, alpha, n, m, Ksat):
        self.theta_R = R
        self.theta_S = S        
        self.alpha = alpha # [1/cm]         
        self.n = n
        self.m = m
        self.Ksat = Ksat
    def __init__(self, R, S, alpha, n, Ksat):
        self.theta_R = R
        self.theta_S = S        
        self.alpha = alpha # [1/cm]         
        self.n = n
        self.m = 1.-1./n
        self.Ksat = Ksat        


