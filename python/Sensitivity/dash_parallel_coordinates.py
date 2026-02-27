"""
Macroscopic:

Dash App - parallel coordinates plot (copy exp names to rerun simulations, or do full dynamic simulations)

pick your simulation based on traits or target properties (filename right)

Daniel Leitner, 2025
"""

import json
import os
import zipfile

import dash_ag_grid as dag
import dash_bootstrap_components as dbc
import global_optimization_tools as got
import numpy as np
import pandas as pd
import plotly.express as px
from dash import Dash, Input, Output, Patch, State, dcc, html
from scipy.cluster.hierarchy import leaves_list, linkage


def create_df(all, dims):
    """creates a pandas data frame with columns given by dims"""
    d = {}
    for key in dims:
        d[key] = np.array([a[key] for a in all])
    return pd.DataFrame(d, columns=dims)


folder_path = "results_cplantbox/"
exp_name = "soybean_all14"

cluster_targets = ["carbon", "surface", "krs", "SUFz", "RLDz"]  # node (=cluster index) is added accordingly
pareto_targets = ["-carbon", "surface", "depth", "krs", "SUFz", "RLDz"]
target_names_parcord = ["node", "carbon", "krs", "SUFz", "RLDz", "depth"]

sort_ = False  # use correlation-based sorting
reduce10_ = False  # take every 10th solution
pareto_ = True  # reduce solutions to pareto solutions (

# ### 1 load everything & merge npz results into input parameter json

# all = got.load_json_files(exp_name, folder_path)  # open parameter files
# got.merge_results(folder_path, all)  # add results

# Read JSON data from the ZIP file
with zipfile.ZipFile(exp_name + ".zip", "r") as zipf:  # zip fuile contains a single json file
    with zipf.open(exp_name + ".json", "r") as json_file:
        all = json.load(json_file)  # Deserialize JSON data
all = list(all.values())  #

# ###  2 filter and cluster index, and (optionally) restrict to pareto solutions
print("data", len(all))
all = got.filter_list(all, "length", 200.0, 30000)  # 76 * 3  *100 * 0.6 = 13680 cm;
all = np.array(all)
print("data after filtering", all.shape)
if reduce10_:
    all = all[::10]
    print("reduce factor 10", all.shape)

m_neurons = 3
n_neurons = 4
node2sample, sample2node, som = got.label_clusters(all, n_neurons, m_neurons, cluster_targets, "som")  # add key "node" to all
print("clustering done")

if pareto_:  # restrict to Pareto solutions
    ind = got.pareto_list(all, pareto_targets)
    all_new = []
    for i, a in enumerate(all):
        if ind[i] == True:
            all_new.append(a)
    all = all_new
    print("number of pareto solutions: ", np.sum(ind))

####  3 parameters

pbounds = {
    "src_a": (3, 11),
    "src_delay_a": (3, 14),
    "lmax145_a": (50, 150),
    "lmax2_a": (5.0, 50.0),
    "lmax3_a": (5.0, 50.0),
    "ln145_a": (0.5, 10.0),
    "ln2_a": (0.5, 10.0),
    "r145_a": (0.2, 7.0),
    "r2_a": (0.2, 7.0),
    "r3_a": (0.2, 7.0),
}

params_conductivity = {
    "anatomy_index145": 0,
    "anatomy_index2": 0,
    "anatomy_index3": 0,
    "kr_young145": 0,
    "kr_young2": 0,
    "kr_young3": 0,
    "kr_old145": 0,
    "kr_old145": 0,
    "kr_old145": 0,
    "kx_young145": 0,
    "kx_young2": 0,
    "kx_young3": 0,
    "kx_old145": 0,
    "kx_old2": 0,
    "kx_old3": 0,
    "conductivity_age1": (1, 21),
    "conductivity_age2": (1, 21),
    "conductivity_age3": (1, 21),
}
params_more = {
    "src_a": (3, 11),
    # 'src_first_a': (3, 14),
    "src_delay_a": (3, 14),
    "lmax145_a": (50, 150),
    "ln145_a": (0.5, 10.0),
    "r145_a": (0.2, 7.0),
    # 'theta145_a': (np.pi / 8., np.pi / 2.),
    # 'tropismN145_a': (0., 3.5),
    # 'hairsLength145_a': (0., 0.1),
    # 'hairsZone145_a': (0., 5.),
    "lmax2_a": (5.0, 50.0),
    # 'ln2_a': (0.5, 10.),
    # 'r2_a': (0.2, 7.),
    # 'theta2_a': (np.pi / 8., np.pi / 2.),
    # 'tropismN2_a': (0., 3.5),
    # 'hairsLength2_a': (0., 0.1),
    # 'hairsZone2_a': (0., 10.),
    "lmax3_a": (0.5, 50.0),
    # 'r3_a': (0.2, 7.),
    # 'theta3_a': (np.pi / 8., np.pi / 2.),
    # 'tropismN3_a': (0., 3.5),
    # 'hairsLength3_a': (0., 0.1),
    # 'hairsZone3_a': (0., 10.)
    "anatomy_index145": (0.0, 0.0),
}

dims = list(params_more.keys())  # for the parallel coordinates
dims.extend(target_names_parcord)

more_dims = dims.copy()
if not sort_:
    more_dims.extend(["exp_name"])  # breaks sorting

color_ = "node"
df2 = create_df(all, more_dims)

if sort_:  # ## Correlation ordering
    corr = df2.corr().values
    dist = 1 - np.abs(corr)  # distance = 1 - |correlation|
    linkage_matrix = linkage(dist, method="average")
    order = leaves_list(linkage_matrix)  # reordering of variables
    ordered_columns = df2.columns[order]
    df = df2[ordered_columns]
else:
    df = df2

# ## START DASH

app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

app.layout = dbc.Container(
    [
        html.H4("Parallel Coordinates"),
        dcc.Graph(
            id="my-graph",
            figure=px.parallel_coordinates(
                df, color=color_, dimensions=df.columns, color_continuous_scale=px.colors.diverging.Spectral, color_continuous_midpoint=2
            ),
        ),
        dbc.Button(
            "Download pick",
            id="download-btn",
            color="primary",
            className="mb-3",
        ),
        dcc.Download(id="download-expnames"),
        dag.AgGrid(
            id="table",
            columnDefs=[{"field": i} for i in df.columns],
            columnSize="sizeToFit",
            defaultColDef={"minWidth": 120},
            dashGridOptions={
                "enableCellTextSelection": True,
                "ensureDomOrder": True,
                "animateRows": False,
            },
        ),
        dcc.Store(id="activefilters", data={}),
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
                    dff3 = pd.DataFrame(columns=df.columns)
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
    State("activefilters", "data"),
)
def updateFilters(data, current_filters):

    if not data:
        return current_filters or {}

    if current_filters is None:
        current_filters = {}

    key = list(data[0].keys())[0]

    # example key: "dimensions[3].constraintrange"
    dim_index = int(key.split("[")[1].split("]")[0])
    col = dims[dim_index]

    current_filters[col] = data[0][key]

    return current_filters


def apply_filters(df, filters):
    dff = df.copy()

    for col, value in filters.items():
        if value:
            if isinstance(value[0], list):
                dff = pd.concat(
                    [dff[dff[col].between(r[0], r[1])] for r in value],
                    ignore_index=True,
                )
            else:
                dff = dff[dff[col].between(value[0], value[1])]

    return dff


@app.callback(
    Output("download-expnames", "data"),
    Input("download-btn", "n_clicks"),
    State("activefilters", "data"),
    prevent_initial_call=True,
)
def download_expnames(n_clicks, filters):

    # Apply current filters
    if filters:
        dff = apply_filters(df, filters)
    else:
        dff = df.copy()

    if "exp_name" not in dff.columns:
        return None

    # One exp_name per line
    content = "\n".join(dff["exp_name"].astype(str))

    return dict(
        content=content,
        filename="exp_name_list.txt",
        type="text/plain",
    )


if __name__ == "__main__":

    app.run(debug=True)
