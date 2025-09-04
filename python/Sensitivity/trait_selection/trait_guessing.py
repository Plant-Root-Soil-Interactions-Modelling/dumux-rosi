

def guess_krs(r, ln, src, a, kr, kx):
    return 0.2 * r - 1.*ln + 1 * src + 0.*a + 0.2 * kr + 0.2 * kx


def guess_SUFz(r, ln, src, a, kr, kx):
    return 1.*r + 1.*ln - 1.*src - 0.2 * a


def guess_transpiration(r, ln, src, a, kr, kx):
    return 0.*r + 0.*ln - 0.2 * src + 0. *a + 1.*kr + 1.* kx


def guess_carbon(r, ln, src, a, kr, kx):
    return -1.*r + 1.*ln + 0.*src - 1.*a


def guess_(r, ln, src, a, kr, kx):
    p = (r, ln, src, a, kr, kx)
    return [guess_krs(*p), guess_SUFz(*p), guess_transpiration(*p), guess_carbon(*p) ]


def guess(primary_root_growth_rate, lateral_root_density, total_inflated_cortical_cell_area_primary_root,
          total_xylem_vessel_Area_primary_root, total_inflated_cortical_cell_area_basal_root, total_xylem_vessel_Area_basal_root, basal_root_number):
    r = primary_root_growth_rate
    ln = 1. / (1. + lateral_root_density)
    src = basal_root_number
    a = 0.25 * (total_inflated_cortical_cell_area_primary_root + total_xylem_vessel_Area_primary_root + total_inflated_cortical_cell_area_basal_root + total_xylem_vessel_Area_basal_root)
    kr = 1. / (1. + a)
    kx = total_xylem_vessel_Area_primary_root + total_xylem_vessel_Area_basal_root
    # print("guess_({:g},{:g},{:g},{:g},{:g},{:g})".format(r, ln, src, a, kr, kx))
    return guess_(r, ln, src, a, kr, kx)
