"""
    Dash App - parallel coordinates plot (copy exp names to rerun simulations, or do full dynamic simulations)
"""
import os
import json
import zipfile
import numpy as np
import pandas as pd

from dash import Dash, dcc, Input, Output, Patch, html
import dash_ag_grid as dag
import dash_bootstrap_components as dbc
import plotly.express as px

import global_optimization_tools as got


def create_df(all, dims):
    """ creates a pandas data frame with columns given by dims """
    d = {}
    for key in dims:
        d[key] = np.array([a[key] for a in all])
    return pd.DataFrame(d, columns = dims)


folder_path = "results_cplantbox/"
# exp_name = "soybean_length14"
exp_name = "soybean_all14"
target_names = ["surface", "depth", "-volume", "krs", "SUFz", "RLDz"]  # ["length", "surface",

""" 1 load everything & merge npz results into input parameter json"""
# all = got.load_json_files(exp_name, folder_path)  # open parameter files
# got.merge_results(folder_path, all)  # add results

# Read JSON data from the ZIP file
with zipfile.ZipFile(exp_name + ".zip", "r") as zipf:
    with zipf.open(exp_name + ".json", "r") as json_file:
        all = json.load(json_file)  # Deserialize JSON data
all = list(all.values())

""" 2 filter and add cluster index"""
all = got.filter_list(all, "length", 200., 30000)  # 76 * 3  *100 * 0.6 = 13680 cm;
m_neurons = 3
n_neurons = 4
node2sample, sample2node, som = got.label_clusters(all, n_neurons, m_neurons, target_names, "som")

# # Plot Pareto solutions
# ind = got.pareto_list(all, target_names)
# all_new = []
# for i, a in enumerate(all):
#     if ind[i] == True:
#         all_new.append(a)
# all = all_new
# print("number of pareto solutions: ", np.sum(ind))

""" 3 parameter space regarding nodes """
pbounds = {
    'src_a': (3, 11),
    'src_delay_a': (3, 14),
    'lmax145_a': (50, 150),
    'lmax2_a': (5., 50.),
    'lmax3_a': (5., 50.),
    'ln145_a': (0.5, 10.),
    'ln2_a': (0.5, 10.),
    'r145_a': (0.2, 7.),
    'r2_a': (0.2, 7.),
    'r3_a': (0.2, 7.),
    }
# pbounds = {
#     "conductivity_age1": (1, 21),
#     "conductivity_age2": (1, 21),
#     "conductivity_age3": (1, 21),
#     # 'src_a': (3, 11),
#     # 'src_first_a': (3, 14),
#     # 'src_delay_a': (3, 14),
#     # 'lmax145_a': (50, 150),
#     # 'ln145_a': (0.5, 10.),
#     # 'r145_a': (0.2, 7.),
#     # 'theta145_a': (np.pi / 8., np.pi / 2.),
#     # 'tropismN145_a': (0., 3.5),
#     'hairsLength145_a': (0., 0.1),
#     'hairsZone145_a': (0., 5.),
#     # 'lmax2_a': (5., 50.),
#     # 'ln2_a': (0.5, 10.),
#     # 'r2_a': (0.2, 7.),
#     # 'theta2_a': (np.pi / 8., np.pi / 2.),
#     # 'tropismN2_a': (0., 3.5),
#     'hairsLength2_a': (0., 0.1),
#     'hairsZone2_a': (0., 10.),
#     # 'lmax3_a': (0.5, 50.),
#     # 'r3_a': (0.2, 7.),
#     # 'theta3_a': (np.pi / 8., np.pi / 2.),
#     # 'tropismN3_a': (0., 3.5),
#     'hairsLength3_a': (0., 0.1),
#     'hairsZone3_a': (0., 10.)
# }

dims = list(pbounds.keys())  # for the parallel coordinates
# dims = ['r145_a', 'r2_a', 'r3_a', 'lmax145_a', 'lmax2_a', 'lmax3_a', 'node' ]
target_names = ["node", "surface", "volume", "depth", "RLDz", "krs", "SUFz"]
dims.extend(target_names)

more_dims = dims.copy()
more_dims.extend(["exp_name"])

color_ = 'node'
df = create_df(all, more_dims)

app = Dash(__name__, external_stylesheets = [dbc.themes.BOOTSTRAP])

fig = px.parallel_coordinates(
    df,
    color = color_ ,
    dimensions = dims,
    color_continuous_scale = px.colors.diverging.Spectral,  # Tealrose
    color_continuous_midpoint = 2,
)

app.layout = dbc.Container(
    [
        html.H4("Parallel Coordinates"),
        dcc.Graph(id = "my-graph", figure = fig),
        dag.AgGrid(
            id = "table",
            columnDefs = [{"field": i} for i in df.columns],
            columnSize = "sizeToFit",
            defaultColDef = {"minWidth":120},
            dashGridOptions = {"enableCellTextSelection": True, "ensureDomOrder": True, "animateRows": False}
        ),
        dcc.Store(id = "activefilters", data = {}),
    ]
)


@app.callback(
    Output("table", "rowData"),
    Input("activefilters", "data"),
)
def udpate_table(data):
    if data:
        dff = df.copy()
        for col in data:
            if data[col]:
                rng = data[col][0]
                if isinstance(rng[0], list):
                    # if multiple choices combine df
                    dff3 = pd.DataFrame(columns = df.columns)
                    for i in rng:
                        dff2 = dff[dff[col].between(i[0], i[1])]
                        dff3 = pd.concat([dff3, dff2])
                    dff = dff3
                else:
                    # if one choice
                    dff = dff[dff[col].between(rng[0], rng[1])]
        return dff.to_dict("records")
    return df.to_dict("records")


@app.callback(
    Output("activefilters", "data"),
    Input("my-graph", "restyleData"),
)
def updateFilters(data):
    if data:
        key = list(data[0].keys())[0]
        col = dims[int(key.split("[")[1].split("]")[0])]
        newData = Patch()
        newData[col] = data[0][key]
        return newData
    return {}


if __name__ == "__main__":
    app.run_server(debug = True, dev_tools_hot_reload = False)
