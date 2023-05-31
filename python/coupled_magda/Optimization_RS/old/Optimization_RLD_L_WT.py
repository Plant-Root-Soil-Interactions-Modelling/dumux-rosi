""" updated soil core example """

import sys
sys.path.append("../../..")
import plantbox as pb
import vtk_plot as vp
import vtk 
import math
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt

path = "rootsystem_params/"
name = "Optimization_loam1"

def get_result(rs, time :float):
    """ 
    Retrieves a the state of the root systems at a certain time 
    in a single SegmentAnalyser object
    
    @param allRS     list of root systems 
    @param time      of the simulation result (days)
    @return SegmentAnalyser object conaining the segments of 
    all root sytstems at the specific time 
    """
    a = pb.SegmentAnalyser(rs)
    a.filter("creationTime", 0., time)
    a.pack()
    return a


def soil_cores(x :list, y :list, r :float,  h :list, up: list):
    """
    A lsit of soil core geometries with a fixed location in the field  
 
    @param x     x coordinates of the soil cores (cm)
    @param y     y coordinates of the soil cores (cm)
    @param r     radius of the soil core (cm)
    @param h     height of the soil core (cm)
    """
    assert len(x) == len(y), "coordinate length must be equal"
    core = []
    for i in range(0,len(h)): 
        core.append(pb.SDF_PlantContainer(r, r, h[i], False))
    cores = []
    for i in range(0, len(x)):
        cores.append(pb.SDF_RotateTranslate(core[i], 0., pb.SDF_Axis.xaxis, pb.Vector3d(x[i], y[i], up[i])))  # just translate
    return cores;

def set_all_sd(rs):
    for p in rs.getRootRandomParameter():
        p.las = 0
        p.lbs = 0
        p.rs = 0
        p.lmaxs = 0
        p.thetas = 0

def err(fitparams):
    r = fitparams[0]; ln1=fitparams[1]; ln1s=fitparams[2]; ln2 = fitparams[3];
    ln2s = fitparams[4]; maxB=fitparams[5]; maxBs=fitparams[6]; 
    
    rs = pb.RootSystem()
    rs.readParameters(path + name + ".xml")

    #replace the relevant parameters with the fit parameters
    p1 = rs.getRootRandomParameter(1)  # tap root type
    p1.r = r;
    p1.ln = ln1;
    p1.lns = ln1s;
    p4 = rs.getRootRandomParameter(4)  # basal root type
    p4.r = r;
    p4.ln = ln1;
    p4.lns = ln1s;

    p2 = rs.getRootRandomParameter(2)
    p2.ln = ln2
    p2.lns = ln2s
    
    srp = rs.getRootSystemParameter()
    srp.maxB = maxB
    srp.maxBs = maxBs

    rs.setOrganRandomParameter(p1)
    rs.setOrganRandomParameter(p2)
    rs.setOrganRandomParameter(p4)
    rs.setRootSystemParameter(srp)

    #define the cores 
    A = np.zeros((3,4))
    A = [[10,    0,     0,    -20],
        [10,     0,   -20,    -40],
        [10,     0,   -40,    -60]]
    A_ = np.array(A)

    x = A_[:,0]
    y = A_[:,1]
    h = np.subtract(A_[:,2],A_[:,3])
    up = A_[:,2]
    dow = A_[:,3]

    cores = soil_cores(x, y, r1, h, up)
    
    #make simulations
    set_all_sd(rs) #set the std zero
    rs.setGeometry(cube) #set the geometry of the cube

    for i in range(0,sims):
        rs.setSeed(i+1) 
        rs.initialize()
        rs.simulate(times)
        if i == 0: 
            rs.writeParameters("rootsystem_params/p1_Optimization_"+mat+"_test.xml")

        for j in range(0, len(cores)):
            core_analyser = get_result(rs, times)
            core_analyser.crop(cores[j]);
            core_analyser.pack()
            tl1 = core_analyser.distribution("length", up[j], dow[j], 1, True)  # vertical length distribution
            tl1 = np.array(tl1) / ( r1 * r1 * math.pi * h[j])  #RLD within the ind cores in cm/cm³
            rm_[j,i] = tl1 #matrix that contains the root length density

    rm = np.mean(rm_,axis = 1)
    std = np.std(rm_,axis = 1)
    rmstd = np.concatenate((rm, std))
    np.savetxt("RLDs/virt_rm_std_"+soil+"_"+genotype+".txt",rmstd,fmt='%.2f')

    #compare measured and simulated core RLDs 
    with open("RLDs/measured_RLD_std_"+soil+"_"+genotype+".txt") as f: #read in measured RLDs
        real_RLD_ = [[float(xx) for xx in line.split()] for line in f]

    #flatten the vectors
    RLD = np.reshape(rmstd, (corenum*2,1))
    real_RLD = np.reshape(np.array(real_RLD_), (corenum*2,1))

    NRMSE = math.sqrt(sum(((np.subtract(RLD,real_RLD)))**2)/len(RLD))
    err = NRMSE
    print('error= ', err)
    return err

######################################################################################
# Set geometry of soil cube, x =20 cm, y=45 cm, height=75 cm
cube = pb.SDF_PlantContainer(75, 20, 45)

#Parameters to be fitted and initial values
r = 1.6 #cm/d
ln1 = 0.2
ln2 = 0.3 #cm
maxB = 50
ln1std = 0.1
ln2std = 0.1 #cm
maxBstd = 20

p0=[r,ln1,ln1std, ln2, ln2std,maxB,maxBstd] #initial guess
soil = 'loam' #or sand
genotype = 'WT' #or RTH3
sims = 20

corenum = 3 #number of soil cores taken 
times = [42, 63, 98, 154]; #days after planting
r1 = 1.5  # core radius
exportVTP = False  # export single root system cores (for debugging)
rm_ = np.zeros((corenum,sims)) #preallocation of result matrices
rm = np.zeros((corenum)) #preallocation of result matrices
std = np.zeros((corenum)) #preallocation of result matrices
rmstd = np.zeros((corenum*2)) #preallocation of result matrices

res = scipy.optimize.minimize(err, p0, method='BFGS')

print("initial guess= ", p0)
print("optimized params = ", res.x)

#Run final simulation and plot the final root system 
rs = pb.RootSystem()
rs.readParameters(path + name + ".xml")
p1 = rs.getRootRandomParameter(1)
p1.r = res.x[0]
p1.ln = res.x[1]
p1.lns = res.x[2]
rs.setOrganRandomParameter(p1)

p4 = rs.getRootRandomParameter(4)
p4.r = res.x[0]
p4.ln = res.x[1]
p4.lns = res.x[2]
rs.setOrganRandomParameter(p4)

p2 = rs.getRootRandomParameter(2)
p2.ln = res.x[3]
p2.lns = res.x[4]
rs.setOrganRandomParameter(p2)

srp = rs.getRootSystemParameter()
srp.maxB = res.x[5]
srp.maxBs = res.x[6]
rs.setRootSystemParameter(srp)

set_all_sd(rs) #set the std zero
rs.writeParameters("rootsystem_params/Optimization_"+soil+"_"+genotype+".xml")

sys.exit()
rs.setGeometry(cube)
rs.initialize()
rs.simulate(21, False)
rs.write("results/p1_example_Maxime_"+mat+".vtp")
ana = pb.SegmentAnalyser(rs)
vp.plot_roots(ana, "creationTime")






