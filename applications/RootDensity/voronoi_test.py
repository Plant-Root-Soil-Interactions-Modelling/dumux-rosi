import numpy as np

points = np.array([[0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [1, 2],
                   [2, 0], [2, 1], [2, 2]])

from scipy.spatial import Voronoi, voronoi_plot_2d

vor = Voronoi(points)

print(vor)

print("points", points.shape)
print("pioint_region", vor.point_region.shape, np.min(vor.point_region), np.max(vor.point_region))
