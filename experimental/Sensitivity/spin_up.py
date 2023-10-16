""" 
Functions for the spin up simulation to obtain realistic initial conditions
"""

import sys; sys.path.append("../../build-cmake/cpp/python_binding/"); sys.path.append("../modules/");
sys.path.append("../../../CPlantBox/src/python_modules"); sys.path.append("../../../CPlantBox/");

import numpy as np
import timeit
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from rosi_richardsnc import RichardsNCSP  # C++ part (Dumux binding), macroscopic soil model
from rosi_richards import RichardsSP  # C++ part (Dumux binding), macroscopic soil model
from richards import RichardsWrapper  # Python part, macroscopic soil model

import plantbox as pb  # CPlantBox
import van_genuchten as vg
import evapotranspiration as evap
from xylem_flux import *
from datetime import *

SMALL_SIZE = 16
MEDIUM_SIZE = 16
BIGGER_SIZE = 164
plt.rc('font', size = SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize = SMALL_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize = MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize = SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize = BIGGER_SIZE)  # fontsize of the figure title
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']


def vg_enviro_type(i:int):
    """ Van Genuchten parameter for enviro-types, called by maize() and soybean() """
    soil = {}
    soil[0] = [0.0639, 0.3698, 0.0096, 1.4646, 4.47]
    soil[1] = [0.0619, 0.3417, 0.0132, 1.3258, 2.03]
    soil[36] = [0.0760, 0.3863, 0.0091, 1.4430, 2.99]
    soil[5] = [ 0.0451, 0.3721, 0.0325, 1.4393, 34.03]
    soil[59] = [0.0534, 0.3744, 0.0171, 1.4138, 13.09]
    table_name = "envirotype{:s}".format(str(i))
    return soil[i], table_name


def maize(i:int):
    """ Parameters for the maize scenario """
    soil, table_name = vg_enviro_type(i)
    min_b = [-38., -8., -200.]  # data from INARI
    max_b = [38., 8., 0.]
    cell_number = [1, 1, 200]
    area = 76. * 16  # cm2
    Kc_maize = 1.2  # book "crop evapotranspiration" Allen, et al 1998
    # Init. (Lini), Dev. (Ldev), Mid (Lmid), Late (Llate)
    # 30, 40, 50, 50 = 170 (Idaho, USA)
    return soil, table_name, min_b, max_b, cell_number, area, Kc_maize


def plot_netinf(start_date, end_date):
    """ For debugging; to check if the stepwise function is passed correctly """
    start_date2 = datetime.strptime(start_date, '%Y-%m-%d %H:%M:%S')
    end_date = datetime.strptime(end_date, '%Y-%m-%d %H:%M:%S')
    timedelta_ = end_date - start_date2
    sim_time = timedelta_.days  # = 191
    times, net_inf = evap.net_infiltration_table_beers_csv(start_date, sim_time, evap.lai_noroots, Kc = 1.)
    plt.plot(times, net_inf, "*-")
    plt.plot([times[0], times[-1]], [0, 0], 'k:')
    plt.show()


def create_initial_soil(soil_, min_b , max_b , cell_number, area, p_top, p_bot, start_date, end_date):
    """
        Creates initial soil conditions by perfoming a spin up simulation from @param start_date to @param end_date
        
        soil_        van genuchten parameter as a list 
        min_b        minimum of soil domain bounding box [cm]
        max_b        maximum of soil domain bounding box [cm]
        cell_number  resolution in x, y, and z direction [1]
        area         area per plant (for debugging and plausibility) [cm2]
        p_top        soil top matric potential [cm]
        p_bot        soil bottom matric potential [cm]
        start_date   start date (format '2020-10-31 00:00:00')
        end_date     end date (format ''2020-10-31 00:00:00')
        
        returns:
        final soil matric potential [cm] 
        final nitrate concentration [kg/m3]
    """
    dt = 3600. / (24.*3600)
    z_ = [0., -30., -30., -200.]  # initial nitrate: top soil layer of 30 cm
    v_ = 0 * np.array([2.6e-4, 2.6e-4, 0.75 * 2.6e-4, 0.75 * 2.6e-4])  #  initial nitrate concentrations: kg / m3 (~4.e-4)
    f_time = 17  # initial fertilisation time: days before planting
    f1 = 4.08e-4  # initial fertilisation amount g/cm2
    nit_flux = 1.e-7 * (75 * 16 * 1)  # nitrification rate [g/day]

    start_date2 = datetime.strptime(start_date, '%Y-%m-%d %H:%M:%S')
    end_date = datetime.strptime(end_date, '%Y-%m-%d %H:%M:%S')
    timedelta_ = end_date - start_date2
    sim_time = timedelta_.days  # = 191

    soil = vg.Parameters(soil_)
    vg.create_mfp_lookup(soil, -1.e5, 1000)

    s = RichardsWrapper(RichardsNCSP())  # water and one solute
    # s = RichardsWrapper(RichardsSP())  # water only
    s.initialize()
    s.createGrid(min_b, max_b, cell_number, False)  # [cm]
    # IC
    s.setICZ_solute(v_[::-1], z_[::-1])  # ascending order...
    s.setLinearIC(p_top, p_bot)  # cm pressure head, equilibrium
    # BC
    times, net_inf = evap.net_infiltration_table_beers_csv(start_date, sim_time, evap.lai_noroots, Kc = 1.)
    s.setVGParameters([soil_])
    s.setTopBC("atmospheric", 0.5, [times, net_inf])  # 0.5 is dummy value
    # s.setTopBC("noflux")
    # s.setBotBC("freeDrainage")
    s.setBotBC("noflux")  # to approximate deep drainage
    # hard coded initial fertilizer application
    sol_times = np.array([0., sim_time - f_time, sim_time - f_time, sim_time - f_time + 1, sim_time - f_time + 1, 1.e3])
    sol_influx = -np.array([0., 0., f1, f1, 0., 0.])  # g/(cm2 day)
    s.setTopBC_solute("managed", 0.5, [sol_times, sol_influx])
    s.setBotBC_solute("outflow", 0.)
    s.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
    s.setParameter("Component.MolarMass", "6.2e-2")
    s.setParameter("Component.LiquidDiffusionCoefficient", "1.7e-9")
    s.initializeProblem()

    print()
    print("1 m2 / area", 1.e4 / area)
    print("domain water volume", s.getWaterVolume(), "cm3  = ", s.getWaterVolume() / 1000, "l")  # OK
    theta = s.getWaterContent()
    print("water content to water volume", np.sum(theta) * area, "cm3")  # OK
    print("domain water volume", s.getWaterVolume() / area, "cm3/cm2  = ", s.getWaterVolume() / area * 10, "l/m2")  # OK
    print("water content to water volume", np.sum(theta) * 1, "cm3/cm  = ", np.sum(theta) * 1 * 10, "l/m2")  # OK
    print("sum net inf", np.sum(net_inf) / 2 * 10, "mm  = ", np.sum(net_inf) / 2 * 1 * 10, "l/m2")  # (cm m2)/m2 = 10 l / m2, (each value is twice)
    print()

    wilting_point = -10000
    s.setCriticalPressure(wilting_point)  # for boundary conditions constantFlow, constantFlowCyl, and atmospheric
    s.ddt = 1e-4  # [day] initial Dumux time step
    c, h = [], []  # resulting solute concentration

    N = int(np.ceil(sim_time / dt))
    for i in range(0, N):
        t = i * dt  # current simulation time
        print(t, "days")
        soil_sol_fluxes = {}  # empy dict
        evap.add_nitrificatin_source(s, soil_sol_fluxes, nit_flux = nit_flux)  # nitrification debendent on tillage practice
        s.setSource(soil_sol_fluxes.copy(), eq_idx = 1)  # richards.py [g/day]
        s.solve(dt)
        # print("last internal ddt", s.ddt * 24 * 3600, "s")
        c.append(s.getSolution_(1))
        h.append(s.getSolutionHead_())

    print("\nAFTER SIMULATION")
    print("domain water volume", s.getWaterVolume(), "cm3  = ", s.getWaterVolume() / 1000, "l")
    theta = s.getWaterContent()
    print("water content to water volume", np.sum(theta) * area, "cm3")
    print("domain water volume", s.getWaterVolume() / area, "cm  = ", s.getWaterVolume() / area * 10, "l/m2")
    print("water content to water volume", np.sum(theta) * 1, "cm  = ", np.sum(theta) * 1 * 10, "l/m2")

    # write all data (for figures created by plot_initial_soil)
    np.save("data/initial_potential_all.npy", h)
    np.save("data/initial_concentration_all.npy", c)

    return s.getSolutionHead_(), s.getSolution_(1)  # matric potential [cm], concentration [kg/m3]


def plot_initial_soil(start_date, end_date):
    """ creates figures presenting the results of the spin up simulation """

    start_date2 = datetime.strptime(start_date, '%Y-%m-%d %H:%M:%S')
    end_date = datetime.strptime(end_date, '%Y-%m-%d %H:%M:%S')
    timedelta_ = end_date - start_date2
    sim_time = timedelta_.days  # = 191
    times, net_inf = evap.net_infiltration_table_beers_csv(start_date, sim_time, evap.lai_noroots, Kc = 1.)

    h = np.load("data/initial_concentration.npy")
    print("mean initial", np.mean(h))

    h = np.load("data/initial_potential_all.npy")
    c = np.load("data/initial_concentration_all.npy")
    c = np.transpose(c)
    c = c[::-1,:]
    h = np.transpose(h)
    h = h[::-1,:]

    fig, ax = plt.subplots(3, 1, figsize = (18, 10), gridspec_kw = {'height_ratios': [1.5, 3, 3]})
    bar = ax[0].bar(times[::2], np.array(net_inf[::2]) * 10, 0.7)  # cm -> mm
    ax[0].set_ylabel("net inf [mm/day]")
    ax[0].set_xlim(times[0], times[-1])
    divider = make_axes_locatable(ax[0])
    cax0 = divider.append_axes('right', size = '5%', pad = 0.05)
    cax0.axis('off')

    divider = make_axes_locatable(ax[1])
    cax = divider.append_axes('right', size = '5%', pad = 0.05)
    cmap_reversed = matplotlib.cm.get_cmap('jet_r')
    im = ax[1].imshow(h, cmap = cmap_reversed, aspect = 'auto', vmin = -1.e3, extent = [0 , sim_time, -200., 0.])  #  interpolation = 'bicubic', interpolation = 'nearest',
    # vmin = -1000, vmax = 0,
    cb = fig.colorbar(im, cax = cax, orientation = 'vertical')
    cb.ax.get_yaxis().labelpad = 30
    cb.set_label('soil matric potential [cm]', rotation = 270)
    ax[1].set_ylabel("depth [cm]")
    ax[1].set_xlabel("time [days]")

    divider = make_axes_locatable(ax[2])
    cax = divider.append_axes('right', size = '5%', pad = 0.05)
    cmap_ = matplotlib.cm.get_cmap('jet')
    im = ax[2].imshow(c, cmap = cmap_, aspect = 'auto', vmin = 0., vmax = 1.e-3, extent = [0 , sim_time, -200., 0.])  #  interpolation = 'bicubic', interpolation = 'nearest',
    cb = fig.colorbar(im, cax = cax, orientation = 'vertical')
    cb.ax.get_yaxis().labelpad = 30
    cb.set_label('nitrate concentration [g/cm$^3$]', rotation = 270)
    ax[2].set_ylabel("depth [cm]")
    ax[2].set_xlabel("time [days]")

    print()
    print("range", np.min(h), np.max(h), "cm")
    print("range", np.min(c), np.max(c), "g/cm3")
    plt.tight_layout()
    plt.show()

    fig, ax1 = plt.subplots()
    h = h[:, -1]
    color = 'tab:red'
    ax1.plot(h, np.linspace(0., -200, h.shape[0]), color = color)
    ax1.set_xlabel("soil matric potential [cm]", color = color)
    ax1.set_ylabel("depth [cm]")
    ax1.tick_params(axis = 'x', labelcolor = color)
    ax2 = ax1.twiny()
    c = c[:, -1]
    color = 'tab:blue'
    ax2.plot(c, np.linspace(0., -200, c.shape[0]), ':', color = color)
    ax2.set_xlabel("nitrate concentration [g/cm$^3$]", color = color)
    ax2.set_ylabel("depth [cm]")
    ax2.tick_params(axis = 'x', labelcolor = color)
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':

    # spin up period
    start_date = '2020-10-31 00:00:00'
    end_date = '2021-05-10 00:00:00'

    # # growing season
    # start_date = '2021-05-01 00:00:00'
    # end_date = '2021-08-31 00:00:00'

    soil_, table_name, min_b, max_b, cell_number, area, Kc = maize(0)
    s = vg.Parameters(soil_)
    # s.plot_retention_curve()

    print("theta at -100", vg.water_content(-100, s))
    print("theta at -330", vg.water_content(-330, s))
    print("theta at -350", vg.water_content(-350, s))
    print("theta at -400", vg.water_content(-400, s))

    # cell_number = [1, 1, 800] # to check stability regarding spatial resolution
    # h, c = create_initial_soil(soil_, min_b , max_b , cell_number, area, -301., -100., start_date, end_date)
    plot_initial_soil(start_date, end_date)

    # write initial soil conditions
    # np.save("data/initial_potential.npy", h)
    # np.save("data/initial_concentration.npy", c)

