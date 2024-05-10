
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