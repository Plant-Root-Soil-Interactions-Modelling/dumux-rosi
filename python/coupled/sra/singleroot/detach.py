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
    suf = rs.get_suf(sim_time = 0.)  # TODO this possibly is time dependent
    krs, _ = rs.get_krs(sim_time = 0.)  # TODO this possibly is time dependent
    rs.setKrValues(suf * krs)
    rs.setKx(np.array([1.]))  # const big

