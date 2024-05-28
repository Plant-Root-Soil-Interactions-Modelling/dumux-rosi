
def theta2H(vg,theta):#(-) to cm
    thetar =vg[0]
    thetas = vg[1]
    alpha = vg[2]
    n = vg[3]
    nrev = 1/(1-1/n)
    H =-(((( (thetas - thetar)/(theta - thetar))**nrev) - 1)**(1/n))/alpha
    return(H)#cm


def getCoefhours(t):
    return (np.sin(np.pi*t*2)+1)/2 #( (t%1) < 0.5)#


def qair2rh(qair, press,TK):
    T0 = 273.16
    RH = 26.3 * press * qair  /(exp((17.67*(TK - T0))/(TK- 29.65)))
    return RH


def setShape1D(s,r_in, r_out,length,nCells = 10, doLogarithmic=True):
    
    logbase = 0.5
    s.r_in = r_in
    s.r_out = r_out
    if doLogarithmic:
        s.points = np.logspace(np.log(r_in) / np.log(logbase), np.log(r_out) / np.log(logbase), nCells, base = logbase)
        
    else:
        s.points = [r_in + (r_out - r_in)/(nCells-1) * i for i in range(nCells) ]
        
    s.createGrid1d(s.points)
    nCells -= 1
    s.setParameter( "Soil.Grid.Cells", str(nCells))
    s.setParameter( "Problem.segLength", str(length))
    return s

def setBC1D(s):    
    s.setInnerBC("fluxCyl",  s.win)  # [cm/day] #Y 0 pressure?
    s.setOuterBC("fluxCyl",0.)
    
    s.setParameter( "Soil.BC.Bot.C1Type", str(3))
    s.setParameter( "Soil.BC.Top.C1Type", str(3))
    s.setParameter( "Soil.BC.Bot.C1Value", str(s.exuds_in)) 
    s.setParameter( "Soil.BC.Top.C1Value", str(0.)) 


    s.setParameter( "Soil.BC.Bot.C2Type", str(3))
    s.setParameter( "Soil.BC.Top.C2Type", str(3))
    s.setParameter( "Soil.BC.Bot.C2Value", str(s.exudl_in)) 
    s.setParameter( "Soil.BC.Top.C2Value", str(0.)) 

    for i in range(s.numFluidComp, s.numComp):
        s.setParameter( "Soil.BC.Bot.C"+str(i)+"Type", str(3))
        s.setParameter( "Soil.BC.Top.C"+str(i)+"Type", str(3))
        s.setParameter( "Soil.BC.Bot.C"+str(i)+"Value", str(0)) 
        s.setParameter( "Soil.BC.Top.C"+str(i)+"Value", str(0 )) 
        

def set_all_sd(rs, s):
    """ # sets all standard deviation to a percantage, i.e. value*s """
    for p in rs.getRootRandomParameter():
        p.a_s = p.a * s
        p.lbs = p.lb * s
        p.las = p.la * s
        p.lns = p.ln * s
        p.lmaxs = p.lmax * s
        p.rs = p.r * s
        p.thetas = p.theta * s
        p.rlts = p.rlt * s  # no used
        p.ldelays = p.ldelay * s
    seed = rs.getRootSystemParameter()  # SeedRandomParameter
    seed.firstBs = seed.firstB * s
    seed.delayBs = seed.delayB * s
    seed.maxBs = seed.maxB * s
    seed.firstSBs = seed.firstSB * s
    seed.delaySBs = seed.delaySB * s
    seed.delayRCs = seed.delayRC * s
    seed.nCs = seed.nCs * s
    seed.nzs = seed.nzs * s
    # todo seed position s

    
    
def init_conductivities_const(r, kr_const = 1.8e-4, kx_const = 0.1):
    """ Hydraulic conductivities  kr [1/day], kx [cm3/day] """
    r.setKr([0, kr_const, kr_const, kr_const, kr_const, kr_const])
    r.setKx([1.e3, kx_const, kx_const, kx_const, kx_const, kx_const])
    return r