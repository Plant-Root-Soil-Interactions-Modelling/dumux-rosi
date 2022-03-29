import sys; sys.path.append("../../../modules/"); sys.path.append("../../../../../CPlantBox/");  sys.path.append("../../../../../CPlantBox/src/python_modules")
sys.path.append("../../../../build-cmake/cpp/python_binding/"); sys.path.append("../../../modules/fv/");
sys.path.append("../");

import numpy as np
import matplotlib.pyplot as plt

import plantbox as pb  # CPlantBox

"""
    old, use aggregated_rs.py insead
    
    detaching the root segments and predermining effective conductivities (= krs*suf)
"""


def detached_conductivities(rs:pb.XylemFlux, suf, krs):

    print("\ndetached conductivities")
    old_kx = rs.kx_f(0., 0)
    print("old kx", old_kx)
    print("old kr", rs.kr_f(0., 0))
    print("krs", krs)
#     z_ = np.linspace(-0.25, -49.75, suf.shape[0])
#     plt.plot(suf, z_)
#     print(np.min(suf), np.max(suf), np.mean(suf), np.sum(suf))
#     plt.show()

    a = 0.05
    L = 50
    l = L / suf.shape[0]  #
    print("l", l)
    LL = np.linspace(50, 0.5, suf.shape[0])
    # print(LL)

    # suf = (1. / suf.shape[0]) * np.ones((suf.shape[0],))
    suf_krs = suf * krs
    print("krs/surf = kr", krs / (2 * a * L * np.pi), rs.kr_f(0., 0))
    print("suf*krs", np.min(suf_krs), np.max(suf_krs), np.sum(suf_krs))
    kr_up = rs.kr_f(0, 0) * 2 * a * l * np.pi * np.ones((len(rs.rs.segments),))
    print("kr_up", np.min(kr_up), np.max(kr_up), np.mean(kr_up))
    print("kr_up - suf_krs", np.min(kr_up - suf_krs), np.max(kr_up - suf_krs), np.mean(kr_up - suf_krs))

    kx_up = np.divide(np.multiply(suf_krs, kr_up), kr_up - suf_krs)

    rs.setKrValues(kr_up / (2 * a * l * np.pi))
    rs.setKxValues(np.multiply(kx_up, LL))  # kx_up * l

    kr_ = kr_up / (2 * a * l * np.pi)
    kx_ = np.multiply(kx_up, LL)  # divide by surface, multiply by volume kx_up * l

    print("segment surface", 2 * a * l * np.pi, "1/s", 1 / (2 * a * l * np.pi))

    print("suf", np.min(suf), np.max(suf), np.sum(suf))
    print("kx_up", np.min(kx_), np.max(kx_), np.mean(kx_), kx_[0], kx_[-1], old_kx / np.mean(kx_))
    print(suf)
    print(kx_)
    # print(kx_up)
    print()

