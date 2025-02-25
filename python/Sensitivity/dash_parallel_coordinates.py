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


""" 1 load everything & merge npz results into input parameter json"""

folder_path = "results_cplantbox/"
exp_name = "soybean_length14"
target_names = ["length", "volume", "depth", "RLDz", "krs", "SUFz"]  # ["length", "surface",

all = got.load_json_files(exp_name, folder_path)  # open parameter files
got.merge_results(folder_path, all)  # add results

print("\nraw data")
d = got.fetch_features(target_names, all)
got.print_info(d, target_names)
print(d.shape)

""" 2 filter """
all = got.filter_list(all, "length", 200., 20000)  # 76 * 3  *100 * 0.6 = 13680 cm;

""" 3 parameter space regarding nodes """
pbounds = {
    'src_a': (3, 11),
    'src_delay_a': (3, 14),
    'lmax145_a': (50, 150),
    'ln145_a': (0.5, 10.),
    'r145_a': (0.2, 7.),
    'lmax2_a': (5., 50.),
    'ln2_a': (0.5, 10.),
    'r2_a': (0.2, 7.),
    'lmax3_a': (5., 50.),
    'r3_a': (0.2, 7.),
    }

dims = list(pbounds.keys())  # for the parallel coordinates

more_dims = dims.copy()
more_dims.extend(["exp_name"])

color_ = 'lmax145_a'
df = create_df(all, more_dims)

app = Dash(__name__, external_stylesheets = [dbc.themes.BOOTSTRAP])

fig = px.parallel_coordinates(
    df,
    color = color_ ,
    dimensions = dims,
    color_continuous_scale = px.colors.diverging.Tealrose,
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
    app.run_server(debug = True)
