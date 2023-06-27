import sys;
sys.path.append("../../../../CPlantBox/");
sys.path.append("../../../../CPlantBox/src/")
sys.path.append("../")
sys.path.append("../../../build-cmake/cpp/python_binding/");
sys.path.append("../../modules/");
sys.path.append("/../data");
sys.path.append("../../../../CPlantBox/src/functional/");
sys.path.append("../../../../CPlantBox/src/rsml/");
sys.path.append("../../../../CPlantBox/src/visualisation/")
sys.path.append("../../../../CPlantBox/src/structural/")
sys.path.append("../../../../CPlantBox/src/external/")
import numpy as np
import matplotlib.pyplot as plt
from functional.xylem_flux import sinusoidal2
import scenario_setup as scenario
import evapotranspiration as evap
import os

"""scenario"""
year = 2019
soil_type = ["loam", "sand"]
genotype = ["WT", "RTH3"]
name = "maize_exudate_2019"
dt_ = 360 / (24 * 3600)
t = np.linspace(1,154,int(154/dt_))
col = ["darkblue", "red", "cornflowerblue", "darkorange"]

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
dummy = 0
for i in range(0,len(soil_type)):
    for j in range(0,len(genotype)):

        os.chdir("../")
        print(os.getcwd())
        soil_, min_b, max_b, cell_number, area, Kc = scenario.maize_SPP(soil_type[i])
        sim_time = 154   #  [day]
        x_, y_, lai = evap.net_infiltration(year, soil_type[i], genotype[j], sim_time, Kc)
        trans, tpot, evap, t_ = evap.get_transpiration(year, sim_time, area, lai, Kc)
        os.chdir("plotting/")

        ax1.plot(t_, tpot,color = col[dummy], label = soil_type[i] + ', '+genotype[j])  # potential transpiration
        ax2.plot(t_, evap,color = col[dummy], linestyle = '--')  # potential evaporation

        #pt = 10 * np.array([ -tpot(t[i], dt_) / area for i in range(0, t.shape[0]) ])
        #ax1.plot(t, pt,color = col[dummy], label = soil_type[i] + ', '+genotype[j])  # potential transpiration
        #ax2.plot(x_, np.cumsum(y_),color = col[dummy], linestyle = '--')  # potential evaporation
        dummy =dummy+1
        
ax1.set_xlabel("Time (d)")
ax1.set_ylabel(r'Cumulative potential transpiration ($mm$)')
ax2.set_ylabel("cumulative net infiltration [mm]")
#ax1.set_ylim([0, 10])

# ask matplotlib for the plotted objects and their labels
lines, labels = ax1.get_legend_handles_labels()
ax1.legend(lines, labels, loc='best')
plt.show()
