"""
    tabularized root hydraulic conductivities (from literature) 
    
    spring barley 
    maize 
    lupine 
"""

import numpy as np
import matplotlib.pyplot as plt

SMALL_SIZE = 16
MEDIUM_SIZE = 16
BIGGER_SIZE = 16
plt.rc('font', size = SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize = SMALL_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize = MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize = SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize = BIGGER_SIZE)  # fontsize of the figure title


def springbarley_conductivities(r, skr = 1., skx = 1.):

    kr00 = np.array([[0., 0.]])  # artificial shoot
    kx00 = np.array([[0., 1.e3]])  # artificial shoot

    # Knippe and Frike 2011 (seminal roots, barley)
    dist_kx_sem = np.array([1.332830601107, 3.93064842903628, 2.96315044962497, 2.7857742914885, 4.93006311115263, 5.81859476576794,
                   7.43130539900036, 7.96783617723113, 8.85893584240891, 9.93346483347755, 12.9730721928149, 14.9377837024319,
                   14.9463162678248, 22.9782247735318, 29.9461545761478, 35.4844362127705, 44.9536083032342, 45.4888550761837,
                   74.968515757407, 85.1521450722157])  # mm
    dist_kx_sem = dist_kx_sem / 10.  # mm -> cm

    kx_sem = np.array([2.77122506635745E-11, 3.7398202493309E-11, 1.38415061517059E-10, 1.79567359509706E-10, 1.96306339548237E-10,
              8.89390331220058E-11, 1.96385540541133E-10, 2.28209789046996E-10, 1.38507900546885E-10, 2.55936584710635E-10, 1.20544773836736E-09,
              6.50362310436064E-10, 2.91544497103574E-10, 1.7099886994521E-09, 2.91948868267954E-09, 1.74056180004427E-09, 3.02847418108279E-09,
              7.83723083349083E-10, 3.27282554680027E-09, 3.36494131181473E-09])  # m3/(Mpa s);
    kx_sem = kx_sem * (24 * 3600 / 0.01012)  # cm2 / day; 1 pa = 0.0102 cm pressure head

    # ii = np.argsort(dist_kx_sem)
    # dist_kx_sem = dist_kx_sem[ii]
    # kx_sem = kx_sem[ii]
    """ data too noisy """
    kx_sem_f = np.poly1d(np.polyfit(dist_kx_sem, kx_sem, 1))
    # plt.plot(dist_kx_sem, kx_sem, "*")
    # x = np.linspace(0, np.max(dist_kx_sem), 200)
    # plt.plot(x, kx_sem_f(x))
    # plt.show()
    dist_kx_sem = np.linspace(0, np.max(dist_kx_sem), dist_kx_sem.shape[0])  # resample
    kx_sem = kx_sem_f(dist_kx_sem)

    dist_kr_sem = np.array([3.83499003558418, 5.78832892699593, 8.82224997723071, 14.8044305825092, 22.9202611065502,
                            29.8988298100679, 35.2452443718087, 44.893443484342, 45.6080841380127, 74.9213744206531, 84.926343572044])  # mm
    dist_kr_sem = dist_kr_sem / 10.  # mm -> cm

    kr_sem = np.array([2.86530346707473E-11, 1.98908062845181E-11, 1.83385640585257E-11, 1.71497794027052E-11, 1.32740876949445E-11,
                       1.62939310267509E-11, 1.2716769213763E-11, 1.283388149226E-11, 1.28339781725778E-11, 1.63668761621923E-11, 1.6369077753957E-11])  # m3/(Mpa s)

    kr_sem = kr_sem * (24 * 3600 / 0.01012)  # cm2 / day; 1 pa = 0.0102 cm pressure head
    # print("kr_sem", kr_sem)

    # ii = np.argsort(dist_kr_sem)
    # dist_kr_sem = dist_kr_sem[ii]
    # kr_sem = kr_sem[ii]
    """ data too noisy """
    kr_sem_f = np.poly1d(np.polyfit(dist_kr_sem, kr_sem, 1))
    # plt.plot(dist_kr_sem, kr_sem, "*")
    # x = np.linspace(0, np.max(dist_kr_sem), 200)
    # plt.plot(x, kr_sem_f(x))
    # plt.show()
    dist_kr_sem = np.linspace(0, np.max(dist_kr_sem), dist_kr_sem.shape[0])  # resample
    kr_sem = kr_sem_f(dist_kr_sem)

    growth_rates = [2.722, 0.3, 0.1]  # growth rates cm/day 2.722
    radii = [0.0325, 0.03, 0.02]

    # kr_sem = np.minimum(skr * kr_sem, 1.)  # values
    # kx_sem = np.minimum(skr * kx_sem, 1.)

    age_kx_sem = []
    age_kr_sem = []
    kr_sem_ = []
    for i, lgr in enumerate(growth_rates):
        age_kx_sem.append(dist_kx_sem / lgr)  # added ages [0.]
        age_kr_sem.append(dist_kr_sem / lgr)  # added ages [1.e-4, -0.1]
        kr_sem_.append(kr_sem / (2 * radii[i] * np.pi))

    r.setKrTables([[kr00[:, 1], kr_sem, kr_sem, kr_sem, kr_sem, kr_sem]],
                  [[kr00[:, 0], age_kr_sem[0], age_kr_sem[1], age_kr_sem[2], age_kr_sem[0], age_kr_sem[0]]])

    r.setKxTables([[kx00[:, 1], kx_sem, kx_sem, kx_sem, kx_sem, kx_sem]],
                  [[kx00[:, 0], age_kx_sem[0], age_kx_sem[1], age_kx_sem[2], age_kx_sem[0], age_kx_sem[0]]])

    # r.plot_conductivities(monocot = True, plot_now = True, axes_ind = [1, 4, 5], lateral_ind = [2, 3])


def maize_conductivities(r, skr = 1., skx = 1.):
    """ Hydraulic conductivities for maize modified from Couvreur et al. (2012), originally from Doussan et al. (1998) """

    kr00 = np.array([[0., 0.]])  # artificial shoot
    kx00 = np.array([[0., 1.e3]])  # artificial shoot

    kr0 = np.array([[-1e4, 0.], [-0.1, 0.], [0., 0.000181], [8., 0.000181], [10, 0.0000648], [18, 0.0000648], [25, 0.0000173], [300, 0.0000173]])
    kr1 = np.array([[-1e4, 0.], [-0.1, 0.], [0., 0.000181], [10., 0.000181], [16, 0.0000173], [300, 0.0000173]])

    kx0 = np.array([[0., 0.000864], [5., 0.00173], [12., 0.0295], [15., 0.0295], [20., 0.432], [300., 0.432]])
    kx1 = np.array([[0., 0.0000864], [5., 0.0000864], [10., 0.0000864], [12., 0.0006048], [20., 0.0006048], [23., 0.00173], [300., 0.00173]])

    kr01 = np.minimum(skr * kr0[:, 1], 1.)  # kr0[:, 1] are values
    kr11 = np.minimum(skr * kr1[:, 1], 1.)  # kr1[:, 1] are values
    r.setKrTables([[kr00[:, 1], kr01, kr11, kr11, kr01, kr01]],
                  [[kr00[:, 0], kr0[:, 0], kr1[:, 0], kr1[:, 0], kr0[:, 0], kr0[:, 0]]])
    kx01 = np.minimum(skx * kx0[:, 1], 1.)
    kx11 = np.minimum(skx * kx1[:, 1], 1.)
    r.setKxTables([[kx00[:, 1], kx01, kx11, kx11, kx01, kx01]],
                  [[kx00[:, 0], kx0[:, 0], kx1[:, 0], kx1[:, 0], kx0[:, 0], kx0[:, 0]]])

    # r.plot_conductivities(monocot = True, plot_now = True, axes_ind = [1, 4, 5], lateral_ind = [2, 3])
    # dd


def lupine_conductivities(r, skr = 1., skx = 1.):
    """ Hydraulic conductivities for lupine following Zarebanadkouki et al. (2016) """
    kr00 = np.array([[0., 0.]])  # artificial shoot
    kx00 = np.array([[0., 1.e3]])  # artificial shoot
    kr0 = np.array([[-1.e4, 0.], [-0.1, 0.], [0., 1.14e-03], [2, 1.09e-03], [4, 1.03e-03], [6, 9.83e-04], [8, 9.35e-04], [10, 8.90e-04],
                    [12, 8.47e-04], [14, 8.06e-04], [16, 7.67e-04], [18, 7.30e-04], [20, 6.95e-04], [22, 6.62e-04], [24, 6.30e-04], [26, 5.99e-04],
                    [28, 5.70e-04], [30, 5.43e-04], [32, 5.17e-04]])
    kr1 = np.array([[-1e4, 0.], [-0.1, 0.], [0., 4.11e-03], [1, 3.89e-03], [2, 3.67e-03], [3, 3.47e-03], [4, 3.28e-03],
                    [5, 3.10e-03], [6, 2.93e-03], [7, 2.77e-03], [8, 2.62e-03], [9, 2.48e-03], [10, 2.34e-03], [11, 2.21e-03],
                    [12, 2.09e-03], [13, 1.98e-03], [14, 1.87e-03], [15, 1.77e-03], [16, 1.67e-03], [17, 1.58e-03]])
    kx0 = np.array([[0., 6.74e-02], [2, 7.48e-02], [4, 8.30e-02], [6, 9.21e-02], [8, 1.02e-01], [10, 1.13e-01],
                    [12, 1.26e-01], [14, 1.40e-01], [16, 1.55e-01], [18, 1.72e-01], [20, 1.91e-01], [22, 2.12e-01], [24, 2.35e-01],
                    [26, 2.61e-01], [28, 2.90e-01], [30, 3.21e-01], [32, 3.57e-01]])
    kx1 = np.array([[0., 4.07e-04], [1, 5.00e-04], [2, 6.15e-04], [3, 7.56e-04], [4, 9.30e-04], [5, 1.14e-03],
                    [6, 1.41e-03], [7, 1.73e-03], [8, 2.12e-03], [9, 2.61e-03], [10, 3.21e-03], [11, 3.95e-03], [12, 4.86e-03],
                    [13, 5.97e-03], [14, 7.34e-03], [15, 9.03e-03], [16, 1.11e-02], [17, 1.36e-02]])

    kr01 = np.minimum(skr * kr0[:, 1], 1.)
    kr11 = np.minimum(skr * kr1[:, 1], 1.)
    r.setKrTables([kr00[:, 1], kr01, kr11, kr11, kr01, kr01],
                  [kr00[:, 0], kr0[:, 0], kr1[:, 0], kr1[:, 0], kr0[:, 0], kr0[:, 0]])
    kx01 = np.minimum(skx * kx0[:, 1], 1.)
    kx11 = np.minimum(skx * kx1[:, 1], 1.)
    r.setKxTables([[kx00[:, 1], kx01, kx11, kx11, kx01, kx01]],
                  [[kx00[:, 0], kx0[:, 0], kx1[:, 0], kx1[:, 0], kx0[:, 0], kx0[:, 0]]])

    # r.plot_conductivities(monocot = False, plot_now = True, axes_ind = [1, 4], lateral_ind = [2, 3])


def const_conductivities(r):
    """ const conductivities """
    kr00 = 0.  # artificial shoot
    kx00 = 1.e3  # artificial shoot
    kr_const_ = 1.73e-4  # [1/day]
    kx_const_ = 4.32e-2  # [cm3/day]
    r.setKr([kr00, kr_const_, kr_const_, kr_const_, kr_const_, kr_const_, kr_const_])
    r.setKx([kx00, kx_const_, kx_const_, kx_const_, kx_const_, kx_const_, kx_const_])
