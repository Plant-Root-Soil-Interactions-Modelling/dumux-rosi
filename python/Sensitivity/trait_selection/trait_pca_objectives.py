import sys; sys.path.append("../../modules"); sys.path.append("../../../build-cmake/cpp/python_binding/");
sys.path.append("../../../../CPlantBox");  sys.path.append("../../../../CPlantBox/src");
sys.path.append("../")

import pandas as pd
from scipy.cluster.vq import vq, whiten, kmeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler
import plotly.graph_objects as go
import plotly.io as pio
import numpy as np

from trait_guessing import *
import global_optimization_tools as got

# Load the Excel file
data = pd.read_excel("Root_Traits_Of_Interest.xlsx")

# Remove rows with missing values
data = data.dropna()

# Separate the features and the genotype column
genotypes = data['Genotype']
features = data.drop(columns = ['Genotype', 'basal_root_number'])
features2 = data.drop(columns = ['Genotype'])

# Standardize the features
scaler = MinMaxScaler()
scaled_features = scaler.fit_transform(features)
scaled_features2 = scaler.fit_transform(features2)

# guess objectives
sf2 = scaled_features2  # rename
obj = []
for i in range(0, scaled_features2.shape[0]):
     obj.append(guess(sf2[i, 0], sf2[i, 1], sf2[i, 2] , sf2[i, 3], sf2[i, 4], sf2[i, 5], sf2[i, 6]))

obj = np.array(obj)  # krs, SUFz, transpiration, carbon

# add cluser indices as color
cluster_n = 4
centers, dist = kmeans(scaled_features, cluster_n)
cluster_indices, _ = vq(scaled_features, centers)

colors = cluster_indices; title_ = f'PCA Plot {cluster_n} clusters'
# colors = features2['basal_root_number']; title_ = f'PCA Plot - basal roots'
# colors = obj[:, 0]; title_ = f'PCA Plot - krs'
# colors = obj[:, 1]; title_ = f'PCA Plot - SUFz'
# colors = obj[:, 2]; title_ = f'PCA Plot - transpiration'
# colors = obj[:, 3]; title_ = f'PCA Plot - carbon'

print("obj", obj.shape)
pareto = got.pareto_data(-obj, [0, 1, 2, 3])
print("pareto", pareto.shape, np.sum(pareto), np.min(pareto), np.max(pareto))
print(pareto)
colors = 1.*pareto; title_ = f'PCA Plot pareto solutions'

# Perform PCA
pca = PCA()
pca_result = pca.fit_transform(scaled_features)

# Calculate explained variance
explained_var = (pca.explained_variance_ratio_ * 100).round(1)

# Create PCA plot without text labels
fig = go.Figure(data = go.Scatter(
    x = pca_result[:, 0],
    y = pca_result[:, 1],
    mode = 'markers',
    marker = dict(
        size = 10,
        color = colors,  # Color by the scalar feature
        colorscale = 'Viridis',  # Choose from Plotly color scales
        colorbar = dict(title = 'cluster indices')  # Customize colorbar title
    ),
    # Add custom hover information
    hovertemplate =
        'Genotype: %{customdata[0]}<br>' +
        'PC1: %{x:.2f}<br>' +
        'PC2: %{y:.2f}<extra></extra>',
    customdata = pd.DataFrame({
 #       'Row': row_numbers,
        'Genotype': genotypes
    }).values
))

fig.update_layout(
    title = title_,
    xaxis_title = f'PC1 ({explained_var[0]}%)',
    yaxis_title = f'PC2 ({explained_var[1]}%)'
)

# Open in a new browser tab (best cross-platform method)
pio.renderers.default = 'browser'
fig.show()
