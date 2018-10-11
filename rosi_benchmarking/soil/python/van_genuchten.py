#
# Mualem - van Genuchten model, equations from van Genuchten, MT 1980
# 
import math
import numpy as np

#
# pascal to pressure head, converts a numpy array 
#
def pa2head(pa_, pnref = 1.e5, g = 9.81):
    h = np.zeros(len(pa_))
    for i,p in enumerate(pa_):
        h[i] = (p-pnref) * 100./1000./g
    return h
   
#
# pressure to pascal, converts a numpy array 
#      
def head2pa(h_, pnref = 1.e5, g = 9.81):
    pa = np.zeros(len(h_))
    for i,h in enumerate(h_):
        pa[i] =  pnref + h/100.*1000.*g
    return pa

#
# class containing the van genuchten parameters
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

#
# returns the volumetric water content at a given pressure head  according to the van genuchten model (Eqn 21)
#
def water_content(h, sp):
    return sp.theta_R + (sp.theta_S-sp.theta_R)/pow(1. + pow(sp.alpha*abs(h),sp.n),sp.m)

#
# returns pressure head at a given volumetric water content according to the van genuchten model
#
def pressure_head(theta, sp): 
    theta = min(theta,sp.theta_S) # saturated water conent is the maximum 
    return - pow( pow( (sp.theta_S - sp.theta_R)/(theta - sp.theta_R), (1./sp.m))-1., 1./sp.n) / sp.alpha

#
# returns the effective saturation according to the van genuchten model (dimensionless water content, Eqn 2)
#
def effective_saturation(h,sp):
    h = min(h,0) # pressure head is negative, zero the maximum
    theta = water_content(h,sp)
    se = (theta-sp.theta_R)/(sp.theta_S-sp.theta_R)
    return se

#
# returns the hydraulic conductivity according to the van genuchten model (Eqn 8)
#
def hydraulic_conductivity(h,sp):
    se = effective_saturation(h,sp) 
    K = sp.Ksat*math.sqrt(se)*( (1. - pow(1. - pow(se, 1. / sp.m),sp.m)) ** 2 )
    return K 

#
# returns the specific moisture storage according to the van genuchten model
#
def specific_moisture_storage(h,sp):
    C = -sp.alpha*sp.n*np.sign(h)*(1. / sp.n - 1.) * pow(sp.alpha*abs(h), sp.n-1.) * (sp.theta_R-sp.theta_S) * pow(pow(sp.alpha*abs(h),sp.n) + 1., 1./sp.n-2.)
    return C
#
# returns the water diffusivity (Eqn 11)
#
def water_diffusivity(TH,theta_i,theta_sur, sp):
    theta=TH*(theta_i-theta_sur) + theta_sur
    Se = (theta-sp.theta_R)/(sp.theta_S-sp.theta_R)
    m = sp.m
    D=(1-m)*sp.Ksat/(sp.alpha*m*(sp.theta_S-sp.theta_R)) * pow(Se,0.5-1./m) * (pow(1-pow(Se,1./m),-m) + pow(1-pow(Se,1/m),m)-2)
    return D


