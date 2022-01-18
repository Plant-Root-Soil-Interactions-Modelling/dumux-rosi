import sys; sys.path.append("../../../modules/"); sys.path.append("../../../../../CPlantBox/");  sys.path.append("../../../../../CPlantBox/src/python_modules")
sys.path.append("../../../../build-cmake/cpp/python_binding/"); sys.path.append("../../../modules/fv/");
sys.path.append("../");

import numpy as np
import plantbox as pb  # CPlantBox

"""
    Jans upsacling:
    detaching the root segments and predermining effective radial conductivity (= krs*suf)
"""


def detached_conductivities(rs:pb.XylemFlux):
    suf = rs.get_suf(sim_time=0.)  # TODO this possibly is time dependent
    krs, _ = rs.get_krs(sim_time=0.)  # TODO this possibly is time dependent

    a = 0.05
    l = 0.5

    suf_krs = suf * krs
    kr_up = rs.kr_f(0, 0) * 2 * a * l * np.pi * np.ones((len(rs.rs.segments),))
    kx_up = np.divide(np.multiply(suf_krs, kr_up), kr_up - suf_krs)

    rs.setKrValues(kr_up / (2 * a * l * np.pi))
    rs.setKxValues(kx_up * l)

    kr_ = kr_up / (2 * a * l * np.pi)
    kx_ = kx_up * l

    print("kr", np.min(kr_), np.max(kr_), np.mean(kr_))
    print("kx", np.min(kx_), np.max(kx_), np.mean(kx_))
