"""
functions for the aggregated approach
"""
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt  # for debugging, e.g. to check suf

import plantbox as pb  # CPlantBox
from functional.xylem_flux import *  # root system Python hybrid solver
import visualisation.vtk_plot as vp  # for debugging

import sra


def double_(rsx, rsx2):
    """ inserts dummy values for the artificial segments """
    rsx2[:, 1] = rsx  # rsx2.shape = (ns, 2)
    return np.array(rsx2.flat)  # 0, rsx[0], 0, rsx[1], ...


def get_aggregated_params(r, rs_age, min_b, max_b, cell_number):
    """ returns aggregated root system parameters krs, and per layer: suf_, kr_surf_, surf_, l_, a_, 
        ordered top to bot 
        
        r                     XylemFluxPython
        rs_age                root system age [day]
        min_b                 minimal soil bounds [cm]
        max_b                 maximal soil bounds [cm]
        cell_number           spatial resolution of fv scheme [1] 
    
        returns:
        krs [cm2/day]         root system conductivity 
        suf_ [1]              standard uptake fraction
        kr_surf_ [cm2/day]    kr times surface (summed up per layer)
        surf_ [cm2]           surface (summed up per layer)
        l_ [cm]               length (summed up per layer)
        a_ [cm]               mean layer radius 
    """
    ana = pb.SegmentAnalyser(r.rs)
    krs, _ = r.get_krs(rs_age)
    suf = r.get_suf(rs_age)  # SUF per segment
    ana.addData("suf", suf)
    suf_ = ana.distribution("suf", max_b[2], min_b[2], cell_number[2], False)  # SUF per layer
    segs = r.rs.segments
    nodeCTs = r.rs.nodeCTs
    subTypes = r.rs.subTypes
    lengths = r.rs.segLength()
    radii = r.rs.radii
    kr_surf = np.zeros((len(segs),))  # kr times surface summed up per layer
    for i, s in enumerate(segs):
        age = rs_age - nodeCTs[s.y]
        st = subTypes[i]
        kr_surf[i] = 2 * radii[i] * np.pi * lengths[i] * r.kr_f(age, st)  # [cm2 / day]
    ana.addData("kr_surf", kr_surf)
    kr_surf_ = ana.distribution("kr_surf", max_b[2], min_b[2], cell_number[2], False)
    surf_ = np.array(ana.distribution("surface", max_b[2], min_b[2], cell_number[2], False))
    l_ = np.array(ana.distribution("length", max_b[2], min_b[2], cell_number[2], False))
    a_ = np.divide(surf_, 2 * np.pi * l_)
    return krs, np.array(suf_), np.array(kr_surf_), surf_, l_, a_


def create_aggregated_rs(r, rs_age, min_b, max_b, cell_number):
    """ ceates a root system with one segment per layer, 
        connected to the root collar by an artificial segments,
        with adjsuted conductivities
    
        r                XylemFluxPython
        rs_age 
        min_b
        max_b
        cell_number    
    """
    krs, suf_, kr_surf_, surf_, l_, a_ = get_aggregated_params(r, rs_age, min_b, max_b, cell_number)  # TODO check if it holds for RS

    print("krs", krs)
    print("\nSUF", suf_.shape, np.sum(suf_))
    # print(suf_)
    # print(list(suf_[0:100]))
    # print("\nkr*surf", kr_surf_.shape)
    # print(kr_surf_)
    # print("\nsurface", surf_.shape)
    # print(surf_)
    # print("\nlength", l_.shape)
    # print(l_)
    # print("\nradius", a_.shape)
    # print(a_)
    # print("\n\n")
    # z_ = np.linspace(max_b[2], min_b[2], cell_number[2])  # top to bot
    # plt.plot(suf_, z_)
    # plt.show()

    nz = int(cell_number[2])
    nodes = [pb.Vector3d(0, 0, 0)]
    segs, radii = [], []  # maximal ns segments
    dx = (max_b[2] - min_b[2]) / cell_number[2]
    z_ = np.linspace(max_b[2] - dx / 2, min_b[2] + dx / 2, nz)  # cell centers
    # print(z_)

    c = 0
    for i in range(0, nz):
        if l_[i] > 0:
            nodes.append(pb.Vector3d(0, 0, z_[i]))  # node 2*i+1
            nodes.append(pb.Vector3d(l_[i], 0, z_[i]))  # node 2*i+2
            segs.append(pb.Vector2i(0, 2 * c + 1))  # artificial shoot segment
            radii.append(0.)
            segs.append(pb.Vector2i(2 * c + 1, 2 * c + 2))  # normal segment
            radii.append(a_[i])
            c += 1
    # print("number of segments", len(segs))

    rs = pb.MappedSegments(nodes, segs, radii)
    rs.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]), pb.Vector3d(max_b[0], max_b[1], max_b[2]),
                            pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), cut = False)
    rs.setSoilGrid(r.rs.soil_index)  # copy the picker...
    r2 = XylemFluxPython(rs)  # wrap the root system with the xylem model
    r.test()  # sanity checks

    ll = np.abs(z_)
    # print(ll)

    suf_krs = suf_ * krs
    # print(krs)
    # print(suf_krs)
    # print(kr_surf_)

    kr_up, kx_up = [], []
    for i in range(0, nz):
        if l_[i] > 0:

            kr_up.append(0.)  # artificial segment
            if surf_[i] > 0:
                kr_up.append(kr_surf_[i] / surf_[i])  # regular segment mean, layer kr [1/day]
            else:  # no segments in layer
                kr_up.append(0.)  # regular segment

            if kr_surf_[i] - suf_krs[i] > 0:
                # kx_up.append(ll[i] * kx_up_[i])
                kx_up.append((ll[i] * suf_krs[i] * kr_surf_[i]) / (kr_surf_[i] - suf_krs[i]))  # artificial segment
                # Kxupscale=Krs*SFF*Krupscale/(Krupscale-Krs*SUF));  mit Kxupscale*(Hx-Hcollar)=Q
            else:  # no segments in layer
                raise ValueError('create_aggregated_rs() no segment in layer')
            kx_up.append(1.e3)  # regular segment

    r2.setKrValues(kr_up)
    r2.setKxValues(kx_up)

    kr_ = np.array(kr_up)
    kx_ = np.array(kx_up)

    # print(np.array(kr_up))
    # print(np.array(kx_up))

    print("krs", krs)
    print("suf", np.min(suf_), np.max(suf_), np.sum(suf_))
    print("kr_up", np.min(kr_[1::2]), np.max(kr_[1::2]), np.mean(kr_[1::2]))
    print("kx_up", np.min(kx_[0::2]), np.max(kx_[0::2]), np.mean(kx_[0::2]))
    print("kx_up", kx_.shape)
    print("kx_up")
    print(list(kx_[0::2]))

    print(kr_surf_[0])
    print(suf_krs[0])

    # vp.plot_roots(pb.SegmentAnalyser(rs), "radius")
    return r2


def simulate_const(s, r, sra_table_lookup, trans, sim_time, dt):
    """     
    simulates the coupled scenario       
        root architecture is not gowing  
        conductivities are not changing over time
        
    s
    r
    sra_table_lookup             potentials a root soil interface  
    trans
    sim_time    
    dt
    
    TODO recyle factorisation of left hand side ... 
    """
    wilting_point = -15000  # cm
    skip = 6  # for output and results, skip iteration
    rs_age = 0.  # day
    max_iter = 10  # maximum for fix point iteration

    if isinstance(sra_table_lookup, RegularGridInterpolator):
        root_interface = sra.soil_root_interface_table
    else:
        root_interface = sra.soil_root_interface

    start_time = timeit.default_timer()

    nodes = r.rs.nodes
    segs = r.rs.segments
    ns = len(r.rs.segments)
    mapping = np.array([r.rs.seg2cell[j] for j in range(0, ns)])
    outer_r = r.rs.segOuterRadii()
    inner_r = r.rs.radii
    types = r.rs.subTypes
    rho_ = np.divide(outer_r, np.array(inner_r))
    rho_ = np.minimum(rho_, np.ones(rho_.shape) * 200)  ############################################ (too keep within table)

    psi_x_, psi_s_, sink_ , x_, y_, psi_s2_ = [], [], [], [], [], []  # for post processing

    sx = s.getSolutionHead()  # inital condition, solverbase.py
    cell_centers = s.getCellCenters()
    cell_centers_z = np.array([cell_centers[mapping[2 * j + 1]][2] for j in range(0, int(ns / 2))])
    seg_centers_z = np.array([0.5 * (nodes[segs[2 * j + 1].x].z + nodes[segs[2 * j + 1].y].z)  for j in range(0, int(ns / 2))])
    hsb = np.array([sx[mapping[2 * j + 1]][0] for j in range(0, int(ns / 2))])  # soil bulk matric potential per segment
    # print(list([mapping[2 * j + 1] for j in range(0, int(ns / 2))]))

    kr_ = np.zeros((ns,))
    rsx = hsb.copy()  # initial values for fix point iteration
    rsx2 = np.zeros((rsx.shape[0], 2))

    # r.init_solve_static(rs_age, double_(rsx, rsx2), False, wilting_point, soil_k = [])  # speed up & and forever static...

    rx = r.solve(rs_age, -trans * sinusoidal2(0., dt), 0., double_(rsx, rsx2), False, wilting_point, soil_k = [])
    rx_old = rx.copy()

    kr_ = np.array([r.kr_f(rs_age, types[j], 2, j) for j in range(0, len(outer_r))])
    inner_kr_ = np.multiply(inner_r, kr_)  # multiply for table look up

    N = int(np.ceil(sim_time / dt))  # number of iterations

    """ simualtion loop """
    for i in range(0, N):

        t = i * dt  # current simulation time

        wall_iteration = timeit.default_timer()
        wall_fixpoint = timeit.default_timer()

        err = 1.e6
        c = 0
        while err > 1 and c < max_iter:

            """ interpolation """
            wall_interpolation = timeit.default_timer()
            rx_ = rx[1::2] - seg_centers_z  # from total matric potential to matric potential
            hsb_ = hsb - cell_centers_z  # from total matric potential to matric potential
            rsx = root_interface(rx_, hsb_, inner_kr_[1::2], rho_[1::2], sra_table_lookup)  # [1::2] every second entry, starting from 1
            rsx = rsx + seg_centers_z  # from matric potential to total matric potential
            wall_interpolation = timeit.default_timer() - wall_interpolation

            """ xylem matric potential """
            wall_xylem = timeit.default_timer()
            rx = r.solve(rs_age, -trans * sinusoidal2(t, dt), 0., double_(rsx, rsx2), False, wilting_point, soil_k = [])  # xylem_flux.py, cells = False
            err = np.linalg.norm(rx - rx_old)
            wall_xylem = timeit.default_timer() - wall_xylem

            rx_old = rx.copy()
            c += 1

        wall_fixpoint = timeit.default_timer() - wall_fixpoint

        wall_soil = timeit.default_timer()
        fluxes = r.segFluxes(rs_age, rx, double_(rsx, rsx2), approx = False, cells = False)
        err2 = np.linalg.norm(-trans * sinusoidal2(t, dt) - np.sum(fluxes))
        if r.last == "neumann":
            if err2 > 1.e-6:
                print("error: potential transpiration differs summed radial fluxes in Neumann case" , err2, -trans * sinusoidal2(t, dt), np.sum(fluxes))

        soil_fluxes = r.sumSegFluxes(fluxes)
        s.setSource(soil_fluxes.copy())  # richards.py
        s.solve(dt)
        sx = s.getSolutionHead()[:, 0]  # richards.py
        hsb = np.array([sx[mapping[2 * j + 1]] for j in range(0, int(ns / 2))])
        wall_soil = timeit.default_timer() - wall_soil

        wall_iteration = timeit.default_timer() - wall_iteration

        """ remember results ... """
        if i % skip == 0:
            psi_x_.append(rx.copy()[::2])  # cm (per root node)
            psi_s_.append(rsx.copy())  # cm (per root segment)
            sink = np.zeros(sx.shape)
            for k, v in soil_fluxes.items():
                sink[k] += v
            sink_.append(sink)  # cm3/day (per soil cell)
            x_.append(t)  # day
            y_.append(np.sum(sink))  # cm3/day
            psi_s2_.append(sx.copy())  # cm (per soil cell)
            print("{:g}/{:g} {:g} iterations".format(i, N, c), wall_interpolation / (wall_interpolation + wall_xylem), wall_xylem / (wall_interpolation + wall_xylem))

    print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")

    return psi_x_, psi_s_, sink_, x_, y_, psi_s2_

