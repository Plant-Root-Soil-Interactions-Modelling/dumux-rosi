import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter


plt.rc('font', size=18)

path = "../../scripts/results_singleroot/Exudate"
dt = 20 / 60 / 24  # days

scenarios = [
    ('loam','sand','low','low','100','100','0','0',
     ['Loam','Sand']),

    ('loam','loam','low','medium','100','100','0','0',
     ['Low sorption','Medium sorption']),

    ('loam','loam','low','low','100','100','0','0.6',
     ['Transpiration = 0 cm/d','Transpiration = 0.6 cm/d']),

    ('loam','loam','low','low','100','1000','0.6','0.6',
     ['Initial SWP = -100 cm','Initial SWP = -1000 cm'])
]

columns = ['Exud_tot', 'Exud_ads', 'Exud_decay', 'Exud_liq']
linestyles = ['-', '--']

fig, axes = plt.subplots(2, 2)
axes = axes.flatten()

color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

for ax, scenario in zip(axes, scenarios):

    soil1, soil2, sorp1, sorp2, swp1, swp2, trans1, trans2, labels = scenario

    name1 = f"{soil1}_res1_sorption{sorp1}_SWP_ini{swp1}_trans{trans1}"
    name2 = f"{soil2}_res1_sorption{sorp2}_SWP_ini{swp2}_trans{trans2}"

    df1 = pd.read_csv(f"{path}/{name1}/exud.csv")
    df2 = pd.read_csv(f"{path}/{name2}/exud.csv")

    x = np.arange(len(df1)) * dt

    variable_lines = []

    # Plot exudate variables
    for idx, col in enumerate(columns):
        color = color_cycle[idx % len(color_cycle)]

        line, = ax.plot(x, df1[col], color=color, linestyle='-')
        ax.plot(x, df2[col], color=color, linestyle='--')

        variable_lines.append(line)

    # --- First legend (colored variables) ---
    if ax == axes[0]:
        legend1 = ax.legend(
            variable_lines,
            ['Total amount of exudates',
             'Sorbed exudates',
             'Decomposed exudates', 
             'Dissolved exudates'],
            loc='upper right'
        )
        ax.add_artist(legend1)

    # --- Second legend (scenario comparison) ---
    scenario_lines = []
    for ls in linestyles:
        line, = ax.plot([], [], color='k', linestyle=ls)
        scenario_lines.append(line)

    ax.legend(scenario_lines, labels, loc='upper left')
    
    for ax in axes:
        formatter = ScalarFormatter(useMathText=True)
        formatter.set_scientific(True)
        formatter.set_powerlimits((0, 0))  # always scientific
        ax.yaxis.set_major_formatter(formatter)

        # Move exponent (×10^x) to top of axis
        ax.yaxis.get_offset_text().set_size(14)
        ax.yaxis.get_offset_text().set_x(-0.1)  # adjust horizontal position if needed

fig.supxlabel('Time (d)')
fig.supylabel('Exudate content in soil domain (mol)', x=0.05)
plt.show()