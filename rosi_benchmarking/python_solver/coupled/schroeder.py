import van_genuchten as vg
from root_conductivities import *


def mfp(h, soil):
#     return vg.matric_flux_potential(h, soil)
    return vg.fast_mfp[soil](h)


def imfp(mfp, soil):
#     return vg.matric_potential_mfp(h, soil)
    return vg.fast_imfp[soil](mfp)


def schroder_nostress(p, rp, q_root, q_out, r_in, r_out, soil):
    rho = r_out / r_in
    mfp = vg.fast_mfp[soil](p) + (q_root * r_in - q_out * r_out) * (r ** 2 / r_in ** 2 / (2 * (1 - rho ** 2)) + rho ** 2 / (1 - rho ** 2) * (np.log(r_out / r) - 0.5)) + q_out * r_out * np.log(r / r_out)
    return vg.fast_imfp[soil](mfp)


def schroder_stress(p, rp, q_root, q_out, r_in, r_out, soil):
    rho = r_out / r_in
    mfp = (vg.fast_mfp[soil](p) + q_out * r_out * np.log(1 / rho)) * ((r ** 2 / r_root ** 2 - 1 + 2 * rho ** 2 * np.log(r_root / r)) / (rho ** 2 - 1 + 2 * rho ** 2 * np.log(1 / rho))) + q_out * r_out * np.log(r / r_root)
    return vg.fast_imfp[soil](mfp)


def getInnerHead(r, p, q_root, q_out, r_in, r_out, soil):
    """ returns the pressure head at the root surface according to Schroeder et al. """
    # print(rp,p)
    if rp < p:  # flux into root
        r = r_in
        rho = r_out / r_in
        mfp = vg.fast_mfp[soil](p) + (q_root * r_in - q_out * r_out) * (r ** 2 / r_in ** 2 / (2 * (1 - rho ** 2)) + rho ** 2 / (1 - rho ** 2) * (np.log(r_out / r) - 0.5)) + q_out * r_out * np.log(r / r_out)
        if mfp > 0:
            h = vg.fast_imfp[soil](mfp)
        else:
            h = -15000.

        if rp < h:  # flux into root
            # print("hello", rp, h, p)
            return h
        else:  # flux into soil
            return rp  # don't use schröder
    else:  # flux into soil
        return p  # don't use schröder
