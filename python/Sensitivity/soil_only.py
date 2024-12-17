"""
    Simulates water movement in soil without any roots (see also spin_up.py for water and nitrate)
"""
import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

import timeit
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

import scenario_setup as scenario
import soil_model
import evapotranspiration as evap
import sra_new


def plot_soil(sim_time, times, net_inf, h, soil_times, top, top2):
    """ creates figures presenting the results of the macroscopic water movement 
    
        sim_time       final simulation time [day]
        times          sampling times for net infiltration [day] 
        net_inf        potential net infiltration [cm/day]   
        h              soil matric potential [cm]
        soil_times     sampling times for actual net infiltration [day]
        top            actual net infiltration [cm/day]
    """
    h = np.transpose(h)
    h = h[::-1,:]

    print(net_inf.shape, times.shape)
    print("Mean net infitration", np.mean(net_inf), "ranging from", np.min(net_inf), "to", np.max(net_inf), "[cm/day]")
    fig, ax = plt.subplots(2, 1, figsize = (18, 10), gridspec_kw = {'height_ratios': [1.5, 3]})
    bar = ax[0].bar(times[::2], np.array(net_inf[::2]) * 10, 100 * 0.8 / len(times[::2]))  # cm -> mm
    ax[0].plot(soil_times, -10.*top, "r:")
    ax[0].plot(soil_times, -10.*top2, "g:")
    ax[0].set_ylabel("net inf [mm/day]")
    ax[0].set_xlim(times[0], times[-1])
    divider = make_axes_locatable(ax[0])
    cax0 = divider.append_axes('right', size = '5%', pad = 0.05)
    cax0.axis('off')

    print(h.shape[1], "times")
    print("Mean initial matric potential", np.mean(h), "ranging from", np.min(h), "to", np.max(h), "[cm]")
    print()
    divider = make_axes_locatable(ax[1])
    cax = divider.append_axes('right', size = '5%', pad = 0.05)
    cmap_reversed = matplotlib.cm.get_cmap('jet_r')
    extent = [0 , sim_time, -200., 0.]
    im = ax[1].imshow(h, cmap = cmap_reversed, aspect = 'auto', vmin = -1000, extent = extent)  #  interpolation = 'bicubic', interpolation = 'nearest',

    x = np.linspace(0, sim_time, h.shape[1])
    y = np.linspace(0, -200, h.shape[0])
    X, Y = np.meshgrid(x, y)
    contours = ax[1].contour(X, Y, h, [0.], colors = 'black')
    # vmin = -1000, vmax = 0,
    cb = fig.colorbar(im, cax = cax, orientation = 'vertical')
    cb.ax.get_yaxis().labelpad = 30
    cb.set_label('soil matric potential [cm]', rotation = 270)
    ax[1].set_ylabel("depth [cm]")
    ax[1].set_xlabel("time [days]")

    plt.tight_layout()
    plt.show()


def plot_soil_1d(sim_time, h):
    """ plots 1d soil profiles for some  points in time 
    """
    times_plt = [0, 20, 40, 60, 80]
    h = np.transpose(h)
    h = h[::-1,:]

    fig, ax1 = plt.subplots()
    for t in times_plt:
        ind = int(np.round(h.shape[1] * (t / sim_time)))
        print("ind", ind, "/", h.shape[1])
        ax1.plot(h[:, ind], np.linspace(0., -200, h.shape[0]), label = "day {:g}".format(t))
        print(t, "Mean initial matric potential", np.mean(h[:, ind]), "ranging from", np.min(h[:, ind]), "to", np.max(h[:, ind]), "[cm]")
    ax1.set_xlabel("soil matric potential [cm]")
    ax1.set_ylabel("depth [cm]")
    plt.legend()
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':

    """ parameters   """
    envirotype_number = 0
    sim_time = 21  # 87.5  # [day]
    dt = 360 / (24 * 3600)  # time step [day]

    # file = np.load("data/soil_only{:g}.npz".format(envirotype_number))
    # plot_soil(sim_time, file["times"], file["net_inf"], file["h"], file["soil_times"], file["top"])
    # # plot_soil_1d(sim_time, file["h"])
    # dd

    soil_, table_name, min_b, max_b, cell_number, area, Kc = scenario.soybean(envirotype_number)

    start_date = '2021-05-10 00:00:00'  # INARI csv data
    times, net_inf = evap.net_infiltration_table_beers_csvS(start_date, sim_time, evap.lai_soybean2, Kc, initial_age = 0.)
    # trans_soybean = evap.get_transpiration_beers_csvS(start_date, sim_time, area, evap.lai_soybean2, Kc)

    """ initialize """
    s = soil_model.create_richards(soil_, min_b, max_b, cell_number, times = times, net_inf = net_inf, bot_bc = "potential", bot_value = 80)  # .

    """ numerical solution """
    top_ind = s.pick([0., 0., -0.5])
    bot_ind = s.pick([0., 0., -199.5])
    print("Top layer index", top_ind)
    print("Bottom layer index", bot_ind)

    water0 = s.getWaterVolume()  # total initial water volume in domain

    soil_times, h, top_, net_change = [], [], [], []
    top_new, bot_new = [], []

    N = int(np.ceil(sim_time / dt))
    start_time = timeit.default_timer()

    for i in range(0, N):

        t = i * dt  # current simulation time
        if i % 240 == 0:
            print(t, "days")

        water_ = s.getWaterVolume()
        s.solve(dt)
        net_change.append(s.getWaterVolume() - water_)

        h.append(s.getSolutionHead_())
        soil_times.append(t)

        top_.append(s.getNeumann(top_ind))  # the neumann BC that is described

        velocities = s.getVelocities_()
        top_new.append(velocities[top_ind])
        bot_new.append(velocities[bot_ind])

    water = s.getWaterVolume()

    top_new = np.array(top_new)
    bot_new = np.array(bot_new)

    print("Computation time", timeit.default_timer() - start_time)

    np.savez("data/soil_only{:g}".format(envirotype_number), h = h, times = times, net_inf = net_inf, soil_times = soil_times, top = top_, net_change = net_change)

    topflux_ = np.sum(top_) * dt
    topflux = np.sum(top_new[:, 2] * dt)  # positive value means evaporation
    botflux = np.sum(bot_new[:, 2] * dt)  # positive value means influx

    print()
    print("top", topflux, topflux_, "cm", "g/cm2", topflux * area, "cm3", "with area", area, "cm2")
    print("bot", botflux, "cm", "g/cm2", botflux * area, "cm3", "with area", area, "cm2")
    print("Change in water balance", water - water0, "cm3")
    print("top", topflux * area)
    print("bot", botflux * area)
    print("balances (should be 0)", (water - water0) + topflux * area - botflux * area, "cm3")
    print()

    plot_soil(sim_time, times, net_inf, h, soil_times, np.array(top_), top_new[:, 2])

    # print()
    #
    # plt.plot(soil_times, top_new[:, 2])
    # plt.plot(soil_times, bot_new[:, 2])

    plt.show()

