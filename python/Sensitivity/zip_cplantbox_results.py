"""
    Macrcoscopic:

    zips the results from result_cplantbox/ folder?
    
    Daniel Leitner, 2025      
"""
import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

import global_optimization_tools as got

""" 1 load everything & merge npz results into input parameter json"""
folder_path = "results_cplantbox/"
exp_name = "soybean_all14_"

all = got.load_json_files(exp_name, folder_path)  # open parameter files
got.merge_results(folder_path, all)  # add results
all_json = {}
for a in all:
    all_json[a["exp_name"]] = a

# Write all JSON data to a zip file
with zipfile.ZipFile(exp_name + ".zip", "w", zipfile.ZIP_DEFLATED) as zipf:
    with zipf.open(exp_name + ".json", "w") as json_file:
        json_file.write(json.dumps(all_json, indent = 4).encode("utf-8"))

