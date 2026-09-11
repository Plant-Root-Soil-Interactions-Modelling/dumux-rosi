from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
import csv
import os


# ============================================================
# SETTINGS
# ============================================================

# Main results directory
results_dir = Path("../results_field_mixedmodel_reduced/")

# csvfile containing the list of files/prefixes
csv_file = Path("zip_info_maindir.csv")

# Output ZIP file
output_zip = results_dir.parent / "results_field_mixedmodel.zip"


# ============================================================
# READ CSV FILE
# ============================================================

with open(csv_file, newline="", encoding="utf-8") as f:
    items = [row[0] for row in csv.reader(f)]


# ============================================================
# estimate size (zipped and unzipped)
# ============================================================

SCENARIO_DIR = "../results_field_mixedmodel_all/loam_diffusionhigh_sorptionhigh/"

NUMBER_OF_SCENARIOS = 24

# files = [...]  # your existing list of files

scenario_size = 0

for item in items:
    for extension in [".csv", ".txt", ".xlsx"]:
        file_path = os.path.join(SCENARIO_DIR, item + extension)

        if os.path.isfile(file_path):
            scenario_size += os.path.getsize(file_path)
            break

# Rough ZIP estimate for CSV + TXT + XLSX
ESTIMATED_COMPRESSION_RATIO = 0.30

estimated_zip_size = scenario_size * ESTIMATED_COMPRESSION_RATIO

total_uncompressed = scenario_size * NUMBER_OF_SCENARIOS
total_zipped = estimated_zip_size * NUMBER_OF_SCENARIOS

print("One scenario:")
print(f"  Uncompressed: {scenario_size / 1024**3:.2f} GB")
print(f"  Estimated ZIP: {estimated_zip_size / 1024**3:.2f} GB")

print(f"\n{NUMBER_OF_SCENARIOS} scenarios:")
print(f"  Uncompressed: {total_uncompressed / 1024**3:.2f} GB")
print(f"  Estimated ZIP: {total_zipped / 1024**3:.2f} GB")

input("\nPress Enter to continue...")

# ============================================================
# DISPLAY SETTINGS
# ============================================================

print("=" * 60)
print("ITEMS TO INCLUDE")
print("=" * 60)

print("\nMain scenario directory:")


# ============================================================
# FIND SCENARIO DIRECTORIES
# ============================================================

scenario_dirs = [
    p for p in results_dir.iterdir()
    if p.is_dir()
]

print()
print(f"Found {len(scenario_dirs)} scenario directories.")


# ============================================================
# CREATE ZIP
# ============================================================

with ZipFile(output_zip, "w", ZIP_DEFLATED) as zip_file:

    for scenario_dir in scenario_dirs:

        print()
        print("=" * 60)
        print(f"Processing: {scenario_dir.name}")
        print("=" * 60)

        # ====================================================
        # MAIN SCENARIO DIRECTORY
        # ====================================================

        for prefix in items:

            # Find all files whose stem starts with the
            # requested name.
            #
            # Example:
            # Excel: "file1"
            #
            # Matches:
            # file1.txt
            # file1.xlsx
            # file1.csv
            #
            # Does NOT match:
            # myfile1.txt

            matches = [
                p for p in scenario_dir.iterdir()
                if p.is_file()
                and p.stem == prefix
            ]

            if not matches:

                print(
                    f"  WARNING: No file found for '{prefix}' "
                    f"in {scenario_dir.name}"
                )

            for file_path in matches:

                arcname = file_path.relative_to(results_dir)

                zip_file.write(
                    file_path,
                    arcname
                )

                print(f"  Added: {arcname}")


# ============================================================
# DONE
# ============================================================

print()
print("=" * 60)
print("ZIP CREATION COMPLETE")
print("=" * 60)
print(f"Output ZIP: {output_zip}")