import sys; sys.path.append("../modules_fpit/"); sys.path.append("../../../CPlantBox/");  sys.path.append("../../build-cmake/cpp/python_binding/")
sys.path.append("../../../CPlantBox/src/python_modules")
sys.path.append("../../../CPlantBox/");
sys.path.append("../../../CPlantBox/src")

import matplotlib; matplotlib.use('agg')

from rosi_richards10c_cyl import Richards10CCylFoam  # C++ part (Dumux binding)

from richards import RichardsWrapper  # Python part
from richards_no_mpi import RichardsNoMPIWrapper  # Python part of cylindrcial 

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank()
# import smallTest_ads_functions as stf
import scenario_setup as stf
from scenario_setup import write_file_array, write_file_float

""" 
Cylindrical 1D model, diffusion only (DuMux), Michaelis Menten

everything scripted, no input file needed, also works parallel with mpiexec
"""

class saveData: 
    pass

obj = saveData() 

def continueLoop(rs,n_iter):
    cL = ((n_iter < 2) or (rs.err1d > rs.max_err) or 
            (rs.err3d > rs.max_err) or  
            (rs.diff1d3dCurrant_rel>rs.max_err))   
    
    print('continue loop?',cL,'errs',rs.err1d, rs.err3d, 'diff',rs.diff1d3dCurrant_rel)

    return cL
    
results_dir="./results/"
if not os.path.exists(results_dir):
    os.makedirs(results_dir)
else:
    test = os.listdir(results_dir)
    for item in test:
        try:
            os.remove(results_dir+item)
        except:
            pass
usemoles = True
#s = RichardsWrapper(Richards10CCylFoam(), usemoles)
year = 2013
min_b =np.array( [-0, -0, -10.] )/2
max_b = np.array([5, 5, 0.] )/2
cell_number = [5,5,5]
soil_type = "loam"
genotype = "WT"
comp = "phenolics"
usemoles = True
soil_ = stf.vg_SPP(0)
mode = "dumux_10c"
p_mean = -100
css1Function_ = 8
paramIdx = 1640
dt = 1./24.
times = np.array([dt*i for i in range(2)])  # days
noAds_ = False

s, soil = stf.create_soil_model(soil_type, year, soil_,#comp, 
                min_b, max_b, cell_number, demoType = mode, 
                times = None, net_inf = None,
                usemoles = usemoles, dirResults = results_dir, 
                p_mean_ = p_mean, css1Function = css1Function_,
                paramIndx=paramIdx,noAds = noAds_)
                         
cell_volumes = comm.bcast(s.getCellVolumes() , root = 0) #cm3               
##### 1D
cell2rhizoId =0
demoType = mode
dirResults = results_dir
p_mean_ = p_mean
css1Function = css1Function_
paramIndx=paramIdx; noAds = noAds_

s1d=RichardsNoMPIWrapper(Richards10CCylFoam(), usemoles)      
s1d.dirResults = dirResults

volCel = cell_volumes[cell2rhizoId]
r_in = 0.02
length_ = 2.
r_out = np.sqrt(volCel/(np.pi*length_)+r_in*r_in)#0.2
# (r_out*r_out-r_in*r_in)*np.pi*length_
s1d.initialize()

#
#print(s1d.dimWorld,s.dimWorld)
stf.setShape1D(s1d,r_in, r_out,length = length_,nCells = 10,doLogarithmic=False)

stf.setDefault(s1d)
#s1d.setParameter("Newton.MaxSteps", "200")
#s1d.setParameter("Newton.MaxTimeStepDivisions", "100")
#s1d.setParameter("Newton.EnableResidualCriterion", "true") 
#s1d.setParameter("Newton.EnableAbsoluteResidualCriterion", "true")
#s1d.setParameter("Newton.SatisfyResidualAndShiftCriterion", "true")
s1d.setParameter("Problem.reactionExclusive", "0")  
#s1d.setParameter("Newton.MaxRelativeShift","1e-10")
stf.setSoilParam(s1d,paramIndx)
stf.getBiochemParam(s1d,paramIndx,noAds_)
stf.setBiochemParam(s1d)
stf.setIC(s1d,paramIndx)
s1d.win = -1.; 
s1d.exudl_in = 0.
s1d.exuds_in =1e-3
stf.setBC1D(s1d)
p_mean_ = comm.bcast(s.getSolutionHead(), root = 0)[cell2rhizoId]
s1d, __ = stf.setupOther(s1d, css1Function, p_mean_)
#s1d.Qexud = s1d.exuds_in * (2 * np.pi * r_in * s1d.segLength)
s1d.proposed_inner_fluxes = s1d.win* (2 * np.pi * r_in * s1d.segLength)
s1d.proposed_inner_fluxes_sol = s1d.exuds_in* (2 * np.pi * r_in * s1d.segLength)
s1d.proposed_inner_fluxes_mucil = s1d.exudl_in* (2 * np.pi * r_in * s1d.segLength)

# s1d.proposed_outer_wat_fluxes = 0.
# s1d.proposed_outer_fluxes_sol = 0.
# s1d.proposed_outer_fluxes_mucil = 0.
s.outer_R_bc_wat = np.array([0. for i in range(len(cell_volumes))] )
s.outer_R_bc_sol = np.array([np.array([0. for i in range(len(cell_volumes))]) for i in range(s.numComp+1)])

cellIds = np.array([i for i in range(len(cell_volumes))])
mol_total_ = [comm.bcast(s.getContent(idComp, idComp <= 2), root = 0)[cell2rhizoId]  for  idComp in range(1,s.numComp+2) ]# solute content [mol].
wat_total_ = comm.bcast(s.getWaterVolumes(), root = 0)
print(wat_total_[cell2rhizoId], sum(s1d.getWaterVolumesCyl()))

if False:
    print(mol_total_)
    print([sum(s1d.getContentCyl( idComp, idComp <= 2)) for idComp in range(1,s.numComp+2) ])
    print(s.getSolutionHead())
    print(s1d.getSolutionHead())
    print(s.getCellVolumes())
    print(s1d.getCellVolumes(), sum(s1d.getCellVolumes()), (r_out*r_out-r_in*r_in)*np.pi*length_ )
    print(s1d.getPoints(), r_out)
    
#######
s.numFluidComp = 2;s1d.numFluidComp = 2


buTotCBeforeEach = comm.bcast(s.getTotCContent_each(), root = 0) 
s1d.buTotCBeforeEach1d = s1d.getTotCContent_each()


#print('buTotCBeforeEach3d',buTotCBeforeEach[[0]][cell2rhizoId])
print('buTotCBeforeEach1d',s1d.buTotCBeforeEach1d[[0,-1]].sum())
         

s.ddt = 1.e-5

#s.base.reportParams()
#print(s.Qexud/ 24. / 3600.)

def prepare1d(s1d,s, dt):
    # proposed_outer_sol_fluxes
    s1d.proposed_outer_wat_fluxes = s.outer_R_bc_wat[cell2rhizoId] / dt
    s1d.proposed_outer_fluxes_sol = s.outer_R_bc_sol[0][cell2rhizoId] / dt
    s1d.proposed_outer_fluxes_mucil = s.outer_R_bc_sol[1][cell2rhizoId] / dt
    s1d.Q_outer_proposed = (s1d.proposed_outer_fluxes_sol + s1d.proposed_outer_fluxes_mucil) 
    
def solve1d(s1d):
    s1d.ddt =min( 1.e-5,s1d.ddt)
    
    QflowIn = s1d.proposed_inner_fluxes 
    qIn = s1d.win#QflowIn/ (2 * np.pi *s1d.r_in * s1d.segLength) # [cm3/day] -> [cm /day]
    QflowOut = s1d.proposed_outer_wat_fluxes
    qOut = s1d.distributeSource(QflowOut, 0,s1d.numFluidComp)
    s1d.setInnerFluxCyl(qIn)  
    
    botVal = s1d.proposed_inner_fluxes_sol
    topVal = s1d.proposed_outer_fluxes_sol
    botVal_mucil = s1d.proposed_inner_fluxes_mucil
    topVal_mucil = s1d.proposed_outer_fluxes_mucil
    topVals = np.array([ topVal, topVal_mucil])
    
    typeBC = np.full(s1d.numComp,3)
    valueTopBC = np.array([np.zeros(s1d.numberOfCellsTot) for _ in range(s1d.numFluidComp)])
    valueBotBC = np.full(s1d.numComp,0.)


    valueBotBC[0] = botVal / (2 * np.pi * s1d.r_in * s1d.segLength) # [mol/day] -> [mol/cm2 /day]
    valueBotBC[1] = botVal_mucil / (2 * np.pi * s1d.r_in * s1d.segLength) # [mol/day] -> [mol/cm2 /day]
    for nc in range(s1d.numFluidComp):
        if (s1d.getContent(nc+1, isDissolved = (nc < s1d.numFluidComp))).all() == 0.:
            cellIndx = 0
        else:
            cellIndx = None
        valueTopBC[nc] = s1d.distributeSource(topVals[nc], nc+1,s1d.numFluidComp, cellIndx)
        
    # valueTopBC = s1d.distributeSources(np.array([ topVal, topVal_mucil]), 
    #                                          np.array([nc+1 for nc in range(s1d.numFluidComp+1)]), 
    #                                                       s1d.numFluidComp)
    
    print(s1d.getContent(1, isDissolved = True),
            valueTopBC[0],  topVal)
                                                       
    s1d.setSoluteBotBC(typeBC, valueBotBC)
    
    s1d.soil_solute_bu = np.array([sum(s1d.getContent(i+1, isDissolved = (i < s1d.numFluidComp))) for i in range(s1d.numComp+1)])
    
    print('solve 1d', s1d.simTime)
    
    s1d.buTotWBeforeEach1d = comm.bcast(s1d.getWaterVolumes(), root = 0) 
    s1d.buTotCBeforeEach1d =  comm.bcast(s1d.getTotCContent_each(), root = 0) 
    
    s1d.solve(dt, maxDt = 250./(24.*3600.))
    Q_in_m, Q_out_m = s1d.getSavedBC(s1d.r_in,s1d.r_out)     
    #QflowIn_limited = Q_in_m[0]
    s1d.seg_fluxes_limited = Q_in_m[0] #QflowIn_limited
    s1d.seg_fluxes_limited_sol_In = Q_in_m[1]
    s1d.seg_fluxes_limited_mucil_In = Q_in_m[2]

def setSource(s, s1d):

    # water_content = comm.bcast( np.array(s.getWaterContent()),root= 0)  # theta per cell [1]
    soil_water = comm.bcast( np.array(s.getWaterVolumes()),root= 0)  # np.multiply(water_content, cell_volumes)  # water per cell [cm3]
    soil_solute_content = comm.bcast(np.array([np.array(
        s.getContent(i+1, isDissolved = (i < s.numFluidComp))) for i in range(s.numComp+1)]),
                                        root = 0) # mol
    soil_contents = np.concatenate((np.array([soil_water]),soil_solute_content ))  # water_content
    assert soil_contents.shape == (s.numComp+s.numFluidComp, s.numberOfCellsTot)

    soil_fluxes_limited = np.zeros(s.numberOfCellsTot)
    soil_fluxes_limited[cell2rhizoId ] = s1d.seg_fluxes_limited
    new_soil_solute = np.array([sum(s1d.getContent(i+1, isDissolved = (i < s.numFluidComp))) for i in range(s.numComp+1)])
    soil_source_sol = np.full(soil_solute_content.shape,0. )
    for nc in range(s.numComp+1):
        soil_source_sol[nc][cell2rhizoId] = np.array( new_soil_solute[nc]  \
                                                - s1d.soil_solute_bu[nc]  \
                                                - s.outer_R_bc_sol[nc][cell2rhizoId] )/dt
    soil_sources_limited = np.concatenate((np.array([soil_fluxes_limited]),soil_source_sol ))
    assert soil_sources_limited.shape == (s.numComp+s.numFluidComp, s.numberOfCellsTot)
    
    
    #print('soil_sources_limited',soil_sources_limited)
    for idComp in range(s.numComp+1):#cm3/day, mol/day
                
        SSL = soil_sources_limited[idComp].copy()
        soil_contents_temp = soil_contents[idComp].copy()
        if idComp == 1:# add source of css1
            SSL += soil_sources_limited[s.numComp+1].copy()
            soil_contents_temp += soil_contents[s.numComp+1].copy()
            
        if (max(abs(SSL)) != 0.):
            SSL = np.maximum(SSL, -soil_contents_temp/dt)
            
            toAdd= np.maximum(0., -(soil_contents_temp/dt + SSL))
            SSL[np.where(toAdd>0.)] += toAdd[np.where(toAdd>0.)] #+ 1e-20
            
            k_limit_source3d = 0
            epsilon_source3d = 1e-25
            while (not (SSL*dt >= -soil_contents_temp).all()) and (k_limit_source3d <= 10):
                SSL[np.where(SSL*dt < -soil_contents_temp)] += epsilon_source3d
                epsilon_source3d *= 10
                k_limit_source3d += 1
            
            try:
                assert min(soil_contents_temp + SSL*dt) >=0.
            except:
                print(soil_sources_limited[idComp], SSL,dt,  min(soil_contents_temp + SSL*dt) )
                write_file_array("setsourceLim2_"+str(idComp), SSL, directory_ =results_dir, fileType =".csv") 
                raise Exception
            
            test_values = list(SSL)
            test_keys = np.array([i for i in range(len(test_values))])
            res = {}
            for key in test_keys:
                for value in test_values:
                    res[key] = value
                    test_values.remove(value)
                    break                        
            #write_file_float("setsource_"+str(idComp), res, directory_ =results_dir) #mol/day
            write_file_array("setsourceLim1_"+str(idComp),  soil_sources_limited[idComp], 
                             directory_ =results_dir, fileType =".csv") 
            write_file_array("setsourceLim2_"+str(idComp), SSL, directory_ =results_dir, fileType =".csv") 
            # print('res',soil_sources_limited[idComp],res)
            s.setSource(res.copy(), eq_idx = idComp)  # [mol/day], in modules/richards.py


buTotWAfterEach = comm.bcast(s.getWaterVolumes(), root = 0) 
buTotWAfterEach_cell = buTotWAfterEach[cell2rhizoId]
buTotWAfterEach1d = comm.bcast(s1d.getWaterVolumes(), root = 0) 
print('diff init',buTotWAfterEach_cell,sum(buTotWAfterEach1d))
print('diff init Vol',comm.bcast(s.getCellVolumes(), root = 0)[cell2rhizoId],sum(s1d.getCellVolumes()))

for i, dt in enumerate(np.diff(times)):
    s1d.n_iter = 0
    obj.k_iter = 4
    obj.max_err = 1.e-10
    obj.err1d = 100.
    obj.err3d = 100.
    obj.diff1d3dCurrant_rel = 100.
    while  continueLoop(obj, s1d.n_iter):

        if (s1d.n_iter > 0) :
            print('resets')            
            s1d.reset() 
            s.reset() 
            
        assert s1d.n_iter < obj.k_iter
        
        s.ddt =min( 1.e-5,s.ddt)

        if rank == 0:
            print("\n\n\n\t*****", "ext time step", dt, 
                  " d, sim time", s.simTime, 
                  "d, n_iter",s1d.n_iter)
                  
        prepare1d(s1d,s, dt) #get the outer BC

        solve1d(s1d) # solve and get teh actual inner (and outer?) BC
        
        setSource(s, s1d) # set source for 3d
        
        print('solve 3d')
            
        
        buTotWBeforeEach = comm.bcast(s.getWaterVolumes(), root = 0)     
        buTotCBeforeEach =  comm.bcast(s.getTotCContent_each(), root = 0) 
        s.solve(  dt, maxDt = 250./(24.*3600.))

        outer_R_bc = s.getFlux_10c()# mol
        bulkSoil_sources = s.getSource_10c() 
        
        if True: #rank == 0:
            s.outer_R_bc_wat = -outer_R_bc[0]
            # print('s.outer_R_bc_wat',s.outer_R_bc_wat)
            s.outer_R_bc_sol = -outer_R_bc[1:]# mol
            s.sources_wat =  bulkSoil_sources[0] #mol
            s.sources_sol =  bulkSoil_sources[1:] #mol
            # print('outer_R_bc_sol',#outer_R_bc_sol,
            #   s.outer_R_bc_sol[0])#,'css1Flow',css1Flow[-1])
            s.outer_R_bc_sol = np.vstack([s.outer_R_bc_sol, s.outer_R_bc_sol[0]*0.]) #add dummy val for css1
            s.sources_sol = np.vstack([s.sources_sol, s.sources_sol[0]*0.]) #add dummy val for css1
            assert s.outer_R_bc_sol.shape == (s.numComp+1, s.numberOfCellsTot)
            assert s.sources_sol.shape == (s.numComp+1, s.numberOfCellsTot)
              
                
        
        ## err3d
        buTotCAfterEach = comm.bcast(s.getTotCContent_each(), root = 0) 
        if True: #rank == 0:
            #print('buTotCAfterEach',buTotCAfterEach[[0,-1]])
            print('buTotCBefore',sum(buTotCBeforeEach.flatten()),
                'buTotCAfter',sum(buTotCAfterEach.flatten()),
                's.sources_sol',sum(s.sources_sol.flatten()),
                'diff',sum(buTotCAfterEach.flatten()) - sum(buTotCBeforeEach.flatten()) - sum(s.sources_sol.flatten()) \
                - sum(s.outer_R_bc_sol.reshape(-1)))
        
            s.bulkMassCError3dsAll_real = buTotCAfterEach - (buTotCBeforeEach +  s.sources_sol + s.outer_R_bc_sol )
            
            print('s.bulkMassCError3dsAll_real_A',#s.bulkMassCError3dsAll_real, 
                    'cs',s.bulkMassCError3dsAll_real[0][cell2rhizoId],
                   'css1', s.bulkMassCError3dsAll_real[-1][cell2rhizoId], 
                   'cstot',s.bulkMassCError3dsAll_real[0][cell2rhizoId] + s.bulkMassCError3dsAll_real[-1][cell2rhizoId])#.flatten())
            s.bulkMassCError3dsAll_real[0] += s.bulkMassCError3dsAll_real[-1]#put together CS and CSS1
            s.bulkMassCError3dsAll_real[-1][:] = 0.
            obj.err3d =abs( sum(s.bulkMassCError3dsAll_real.flatten())/sum(buTotCAfterEach.flatten()))
            print('s.bulkMassCError3dsAll_real_B',s.bulkMassCError3dsAll_real[0][cell2rhizoId],#s.bulkMassCError3dsAll_real,
                    sum(s.bulkMassCError3dsAll_real.flatten()),'rel',obj.err3d)
                    
        buTotWAfterEach = comm.bcast(s.getWaterVolumes(), root = 0) 
        if True: #rank == 0:
            # cm3
            s.bulkMassWError3dsAll_real = buTotWAfterEach - (buTotWBeforeEach +  s.sources_wat + s.outer_R_bc_wat )
            obj.err3dW =abs( sum(s.bulkMassWError3dsAll_real)/sum(buTotWAfterEach.flatten()))
            print('s.bulkMassWError3dsAll',s.bulkMassWError3dsAll_real[cell2rhizoId],'rel',obj.err3dW)
                
        
        ## err1d
                
        buTotCAfterEach1d = s1d.getTotCContent_each()
        obj.err1d = abs((sum(buTotCAfterEach1d.flatten()) - (sum(s1d.buTotCBeforeEach1d.flatten()) \
            + (s1d.Q_outer_proposed + s1d.proposed_inner_fluxes_sol)*dt) )/sum(buTotCAfterEach1d.flatten()))
        
        print('\nbuTotCAfterEach1d_obs',sum(buTotCAfterEach1d.flatten()),'buTotCAfterEach1d_th',(sum(s1d.buTotCBeforeEach1d.flatten()) \
            + (s1d.Q_outer_proposed + s1d.proposed_inner_fluxes_sol)*dt),s1d.Q_outer_proposed , s1d.proposed_inner_fluxes_sol,dt)
        print('buTotCBeforeEach1d',sum(s1d.buTotCBeforeEach1d.flatten()),
            's1d.proposed_inner_fluxes_sol *dt',s1d.proposed_inner_fluxes_sol *dt,
            's1d.proposed_outer_fluxes_sol_mucil *dt',s1d.Q_outer_proposed *dt,
            'diff',abs((sum(buTotCAfterEach1d.flatten()) - (sum(s1d.buTotCBeforeEach1d.flatten()) \
            + (s1d.Q_outer_proposed + s1d.proposed_inner_fluxes_sol)*dt) )),
            'rel',obj.err1d)
        print("\n")
        
        buTotCAfterEach_cell = np.array([cc3d[cell2rhizoId] for cc3d in buTotCAfterEach])
        diff1d3dCurrant = sum(buTotCAfterEach1d.flatten()) - sum(buTotCAfterEach_cell.flatten())
        obj.diff1d3dCurrant_rel = abs(diff1d3dCurrant/sum(buTotCAfterEach1d.flatten()))
        print('diff1d3dCurrant',diff1d3dCurrant,'1d',sum(buTotCAfterEach1d.flatten()),'3d', sum(buTotCAfterEach_cell.flatten()),
                'rel',obj.diff1d3dCurrant_rel)
                
        buTotWAfterEach_cell = buTotWAfterEach[cell2rhizoId]
        buTotWAfterEach1d = comm.bcast(s1d.getWaterVolumes(), root = 0) 
        
        if True: #rank == 0:
            s1d.bulkMassWError1dsAll_real = (sum(buTotWAfterEach1d) - (sum(s1d.buTotWBeforeEach1d) \
            + (s1d.proposed_outer_wat_fluxes + s1d.proposed_inner_fluxes)*dt) )
            
            obj.err1dW = abs(s1d.bulkMassWError1dsAll_real )/sum(buTotWAfterEach1d)
            
            diff1d3dCurrantW = sum(buTotWAfterEach1d) - buTotWAfterEach_cell
            obj.diff1d3dCurrantW_rel = abs(diff1d3dCurrantW/sum(buTotWAfterEach1d))
            print('\n\ns1d.bulkMassWError1dsAll_real',s1d.bulkMassWError1dsAll_real,'rel',obj.err1dW)
            print('sum(buTotWAfterEach1d)',sum(buTotWAfterEach1d), 'sum(s1d.buTotWBeforeEach1d)',sum(s1d.buTotWBeforeEach1d),
            'proposed inner',s1d.proposed_inner_fluxes*dt,'proposed outer',s1d.proposed_outer_wat_fluxes *dt)
        print('\n\ndiff1d3dCurrantW',diff1d3dCurrantW,'1d',sum(buTotWAfterEach1d),'3d', buTotWAfterEach_cell,
                'rel',obj.diff1d3dCurrantW_rel)
        print('s1d.proposed_inner_fluxes',s1d.proposed_inner_fluxes,s1d.seg_fluxes_limited,
                (s.sources_wat[cell2rhizoId])/dt)#, s.getCellVolumes()[cell2rhizoId] )
        #print('buTotWBeforeEach',buTotWBeforeEach, buTotWAfterEach, (buTotWBeforeEach-buTotWAfterEach)/dt)
        s1d.n_iter += 1
        

