#!/usr/bin/env python

import os
import subprocess
import multiprocessing as mp
import matplotlib.pyplot as plt
from numpy import genfromtxt
import numpy as np

#data_tpot = np.genfromtxt("a_pot.out")
#data_a_act = np.genfromtxt("a_actual.out")
data_a_cumul = np.genfromtxt("a_cumul.out")
#data_b_act = np.genfromtxt("b_actual.out")
#data_b_cumul = np.genfromtxt("b_cumul.out")

#plot_color = ['r', 'g', 'b', 'r', 'g', 'b']

t = data_a_cumul[:,0]
dt = np.insert((t[1:] - t[:-1]), 0, 0)
#print(dt)
Alpha = 0.04
N = 1.6
Qr = 0.08
Qs = 0.43
M = 1-1/N

initialP = -652.3
NR = (1 + (Alpha * abs(initialP))**N)**M
WC_i = Qr + (Qs - Qr)/NR

print("Initial water content: ", WC_i)

DomainWater = WC_i * 960
print("Amount of water in domain: ", DomainWater)

WC = (124.17 - data_a_cumul[:,1])/960
#print(WC)

head = 1/Alpha*(((Qs-Qr)/(WC - Qr))**(1/M) - 1)**(1/N)

fig, ax = plt.subplots()
ax.plot(t, -head)
plt.show()

"""
# prepare plot
fig, ax = plt.subplots(1, 1, dpi=300, figsize=(8,4))

ax.set_xlabel('Time $[d]$')
ax.set_ylabel('Transpiration rate [$cm^3$ $d^{-1}$]')
ax. set_xlim(0, 3)
ax.set_ylim(0, 13)
ax12 = ax.twinx()
ax12.set_ylabel('Cumulative Transpiration [$cm^3$]')
ax12.set_ylim(0, 4.5)

ax.plot(data_tpot[:,0], data_tpot[:,1], 'k', label=r"$Q_\mathrm{pot}$")
ax.plot(data_tpot[:,0], data_a_act[:,1], "--", color=plot_color[0], linewidth=1, label="C1.2a")
ax.plot(data_tpot[:,0], data_b_act[:,1], color=plot_color[1], linewidth=1, label="C1.2b")
ax12.plot(data_tpot[:,0], data_a_cumul[:,1], "--", color=plot_color[0], linewidth=1)
ax12.plot(data_tpot[:,0], data_b_cumul[:,1], color=plot_color[1], linewidth=1)

ax.legend()
fig.tight_layout(rect=[0.01, 0.01, 0.98, 0.98], pad=0.4, w_pad=0.5, h_pad=1.0)
plt.show()
fig.savefig("c12_transpiration.pdf")
"""
