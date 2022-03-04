import sys; sys.path.append("../../../modules/"); sys.path.append("../../../../../CPlantBox/");  sys.path.append("../../../../../CPlantBox/src/python_modules")
sys.path.append("../../../../build-cmake/cpp/python_binding/"); sys.path.append("../../../modules/fv/");
sys.path.append("../");

import numpy as np
import plantbox as pb  # CPlantBox

"""
    Jans upsacling:
    detaching the root segments and predermining effective conductivities (= krs*suf)
"""


def detached_conductivities2(rs:pb.XylemFlux, suf, krs):
    print("old kr", rs.kr_f(0, 0))
    print("old kx", rs.kx_f(0, 0))

    n = len(rs.rs.segments)
    a = 0.05
    l = 0.5
    L = 50
    LL = np.linspace(0.25, L - 0.25, int(n / 2))     
    # print(LL)

    suf_krs = suf * krs
    
    kr_up = []
    kx_up = []
    for i in range(0, int(n / 2)):
        kr_up.append(0.)  # artificial segment
        kr_up.append(rs.kr_f(0, 0) * 2 * a * l * np.pi)  # regular segment
        kx_up.append(LL[i] * (suf_krs[i] * kr_up[2 * i + 1]) / (kr_up[2 * i + 1] - suf_krs[i]))  # artificial segment
        kx_up.append(1.)  # regular segment

    rs.setKrValues(np.array(kr_up) / (2 * a * l * np.pi))
    rs.setKxValues(kx_up) 

    print(np.array(kr_up) / (2 * a * l * np.pi))
    print(kx_up)

    kr_ = np.array(kr_up) / (2 * a * l * np.pi)
    kx_ = np.array(kx_up) 

    print("krs", krs)
    print("suf", np.min(suf), np.max(suf), np.sum(suf))
    print("kr_up", np.min(kr_[1::2]), np.max(kr_[1::2]), np.mean(kr_[1::2]))
    print("kx_up", np.min(kx_[0::2]), np.max(kx_[0::2]), np.mean(kx_[0::2]))
