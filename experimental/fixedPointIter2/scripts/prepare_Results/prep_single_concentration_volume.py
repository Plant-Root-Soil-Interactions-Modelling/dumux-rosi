import numpy as np
import re
import os
import csv
import gzip
import pickle
from pathlib import Path


# =========================================================
# SETTINGS
# =========================================================

path2file = "../results_field_mixedmodel_reduced/"


# =========================================================
# SCENARIO DEFINITIONS
# =========================================================

soiltype = ['loam', 'sand']
diffusion = ['low', 'medium', 'mediumhigh', 'high']
sorption = ['low', 'medium', 'high']

thresh = [1e-10, 1e-9, 1e-8, 1e-7, 5e-7]  # mol/cm^3


# =========================================================
# REMOVE ZONE IDENTIFIERS
# =========================================================

root = Path(path2file)

for file in root.rglob("*Zone.Identifier*"):
    try:
        file.unlink()
        print(f"Deleted: {file}")
    except Exception as e:
        print(f"Error deleting {file}: {e}")


# =========================================================
# FUNCTIONS
# =========================================================

def pad_rows_top(arr, target_rows):
    """
    Add NaN rows at the top if arr has fewer rows than target_rows.
    """
    if arr.shape[0] < target_rows:
        padding = np.full(
            (target_rows - arr.shape[0], arr.shape[1]),
            np.nan,
            dtype=arr.dtype
        )
        return np.vstack([padding, arr])

    return arr


# =========================================================
# MAIN DATA STRUCTURE
# =========================================================

single_conc_volume = {}


# =========================================================
# LOOP THROUGH SCENARIOS
# =========================================================

for i in range(len(soiltype)):

    for m in range(len(sorption)):

        for n in range(len(diffusion)):

            scenario = (
                soiltype[i]
                + '_diffusion'
                + diffusion[n]
                + '_sorption'
                + sorption[m]
                + '/'
            )

            print("\n========================================")
            print("Processing:", scenario)
            print("========================================")


            # -------------------------------------------------
            # TIME
            # -------------------------------------------------

            time = np.loadtxt(
                path2file + scenario + "time.txt",
                delimiter=","
            )[:-1, 0]

            time = np.asarray(time, dtype=np.float32)

            print("shape time:", time.shape)


            # -------------------------------------------------
            # DEPTH
            # -------------------------------------------------

            with open(
                path2file + "/" + scenario + "nodes_Z.txt",
                "r"
            ) as f:

                ndepth_ = f.readlines()[-1]

            ndepth = np.fromstring(
                ndepth_,
                sep=','
            ).astype(np.float32)

            print("shape ndepth:", ndepth.shape)


            # -------------------------------------------------
            # THETA
            # -------------------------------------------------

            theta = np.loadtxt(
                path2file + scenario + "theta.csv",
                delimiter=','
            ).mean(axis=1)

            theta = np.asarray(
                theta,
                dtype=np.float32
            )

            print("shape theta:", theta.shape)

            # -------------------------------------------------
            # PSI XYLEM
            # -------------------------------------------------

            with open(path2file + scenario + "psiXyl.txt", "r") as f:
                psixyl = np.array(
                    [
                        np.mean([float(x) for x in line.strip().split(',') if x.strip()])
                        for line in f
                        if line.strip()
                    ],
                    dtype=np.float32
                )

            print("shape psixyl:", psixyl.shape)

            # -------------------------------------------------
            # MICRO / CYLINDER FILES
            # -------------------------------------------------

            path_cyl = path2file + "/" + scenario + "cyl_val/"

            folder = Path(path_cyl)

            files = list(
                folder.glob("*time*.txt")
            )

            counts = np.sort([
                int(
                    re.search(
                        r'time_(\d+)',
                        str(f)
                    ).group(1)
                )
                for f in files
            ])

            print("Number of time files:", len(counts))


            # -------------------------------------------------
            # ALLOCATE VOLUME ARRAY
            #
            # Shape:
            #   threshold × depth × time
            #
            # float32 instead of float64
            # -------------------------------------------------

            vols = np.zeros(
                (
                    len(thresh),
                    np.max(counts) + 1,
                    len(time)
                ),
                dtype=np.float32
            )


            # -------------------------------------------------
            # DETERMINE MAXIMUM NUMBER OF DEPTH ROWS
            # -------------------------------------------------

            max_rows = None


            # -------------------------------------------------
            # LOOP THROUGH TIME STEPS
            # -------------------------------------------------

            for o in counts:


                # ---------------------------------------------
                # CELL VOLUME
                # ---------------------------------------------

                with open(
                    path_cyl + 'Cyl_cellVol_' + str(o) + ".txt"
                ) as f:

                    cellvol_ = [
                        list(map(float, row))
                        for row in csv.reader(f)
                    ]

                max_len = max(
                    len(row)
                    for row in cellvol_
                )

                cellvol = np.array(
                    [
                        row + [np.nan] * (
                            max_len - len(row)
                        )
                        for row in cellvol_
                    ],
                    dtype=np.float32
                )


                # ---------------------------------------------
                # WATER CONTENT
                # ---------------------------------------------

                with open(
                    path_cyl
                    + 'Cyl_watercontent_'
                    + str(o)
                    + ".txt"
                ) as f:

                    wc_ = [
                        list(map(float, row))
                        for row in csv.reader(f)
                    ]

                wc = np.array(
                    [
                        row + [np.nan] * (
                            max_len - len(row)
                        )
                        for row in wc_
                    ],
                    dtype=np.float32
                )


                # ---------------------------------------------
                # WATER VOLUME
                # ---------------------------------------------

                watvol = (
                    wc[:len(cellvol)]
                    * cellvol[:len(wc)]
                )


                # ---------------------------------------------
                # TOTAL CONCENTRATION
                # ---------------------------------------------

                with open(
                    path_cyl
                    + 'Cyl_content1_'
                    + str(o)
                    + ".txt"
                ) as f:

                    totC_ = [
                        list(map(float, row))
                        for row in csv.reader(f)
                    ]

                totC = np.array(
                    [
                        row + [np.nan] * (
                            max_len - len(row)
                        )
                        for row in totC_
                    ],
                    dtype=np.float32
                )


                # ---------------------------------------------
                # CONCENTRATION
                # ---------------------------------------------

                conc = np.divide(
                    totC[:len(watvol)],
                    watvol[:len(totC)]
                )


                # ---------------------------------------------
                # PAD DEPTH DIMENSION
                # ---------------------------------------------

                if o == 0:

                    max_rows = conc.shape[0]

                else:

                    conc = pad_rows_top(
                        conc,
                        max_rows
                    )

                    watvol = pad_rows_top(
                        watvol,
                        max_rows
                    )


                # ---------------------------------------------
                # APPLY EACH THRESHOLD
                # ---------------------------------------------

                for j in range(len(thresh)):

                    # IMPORTANT:
                    # Work on a copy so one threshold does
                    # not permanently modify the next one.

                    watvol_threshold = watvol.copy()

                    mask = conc < thresh[j]

                    watvol_threshold[mask] = 0

                    dummy_vols = np.sum(
                        watvol_threshold,
                        axis=1
                    )

                    vols[
                        j,
                        o,
                        :len(dummy_vols)
                    ] = dummy_vols.astype(
                        np.float32,
                        copy=False
                    )


            # -------------------------------------------------
            # STORE ONE ENTRY PER SCENARIO
            #
            # All thresholds are contained in vols.
            # -------------------------------------------------

            single_conc_volume[(i, m, n)] = {

                "vols": vols,

                "ndepth": ndepth,

                "time": time,

                "theta": theta, 
                
                "psixyl": psixyl
            }


            print(
                "shape vols:",
                vols.shape
            )


# =========================================================
# SAVE COMPRESSED PICKLE
# =========================================================

output_file = "single_concentration_volume.pkl.gz"

print("\nSaving compressed pickle...")

with gzip.open(
    output_file,
    "wb"
) as f:

    pickle.dump(
        single_conc_volume,
        f,
        protocol=pickle.HIGHEST_PROTOCOL
    )


print("\nFinished!")
print("Output:", output_file)