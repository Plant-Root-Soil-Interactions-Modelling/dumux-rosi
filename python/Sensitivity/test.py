import numpy as np


def clip_polygon_by_box(polygon, box):
    # Define the vertices of the box
    vertices = np.array([
        [box[0][0], box[0][1], box[0][2]],
        [box[0][0], box[1][1], box[0][2]],
        [box[1][0], box[1][1], box[0][2]],
        [box[1][0], box[0][1], box[0][2]],
        [box[0][0], box[0][1], box[1][2]],
        [box[0][0], box[1][1], box[1][2]],
        [box[1][0], box[1][1], box[1][2]],
        [box[1][0], box[0][1], box[1][2]]
    ])

    # Define the edges of the box
    edges = [
        (vertices[0], vertices[1]),
        (vertices[1], vertices[2]),
        (vertices[2], vertices[3]),
        (vertices[3], vertices[0]),
        (vertices[4], vertices[5]),
        (vertices[5], vertices[6]),
        (vertices[6], vertices[7]),
        (vertices[7], vertices[4]),
        (vertices[0], vertices[4]),
        (vertices[1], vertices[5]),
        (vertices[2], vertices[6]),
        (vertices[3], vertices[7])
    ]

    # Define the vertices and edges of the polygon
    polygon_vertices = np.array(polygon)
    polygon_edges = [(polygon_vertices[i - 1], polygon_vertices[i]) for i in range(len(polygon_vertices))]

    # Define a list to store the intersection points
    intersection_points = []

    # Iterate over the edges of the box
    for edge in edges:
        # Iterate over the edges of the polygon
        for polygon_edge in polygon_edges:
            # Compute the intersection point of the two edges
            intersection_point = intersect(edge, polygon_edge)
            if intersection_point is not None:
                # Add the intersection point to the list if it is not already in it
                if not any(np.allclose(intersection_point, p) for p in intersection_points):
                    intersection_points.append(intersection_point)

    # Define a list to store the new polygon vertices
    new_vertices = []

    # Iterate over the vertices of the polygon
    for vertex in polygon_vertices:
        # Add the vertex to the new list if it is inside the box
        if is_inside_box(vertex, box):
            new_vertices.append(vertex)

    # Iterate over the intersection points
    for intersection_point in intersection_points:
        # Add the intersection point to the new list if it is inside the box
        if is_inside_box(intersection_point, box):
            new_vertices.append(intersection_point)

    # Return the new polygon
    return new_vertices


def intersect(edge1, edge2):
    # Compute the direction vectors of the edges
    dir1 = edge1[1] - edge1[0]
    dir2 = edge2[1] - edge2[0]

    print(dir1, dir2)

    # Compute the cross product of the direction vectors
    cross = np.cross(dir1, dir2)
    print(cross)

    # If the cross product is close to zero, the edges are parallel and do not intersect.
    if np.allclose(cross, [0, 0, 0]):
        return None

    # Compute the parameters of the line segments where they intersect
    t1 = np.cross(edge2[0] - edge1[0], dir2) / np.dot(dir1, cross)
    t2 = np.cross(edge1[0] - edge2[0], dir1) / np.dot(dir2, cross)

    # If the parameters are between 0 and 1, the line segments intersect
    if 0 <= t1 <= 1 and 0 <= t2 <= 1:
        # Compute the intersection point
        intersection_point = edge1[0] + t1 * dir1
        return intersection_point
    else:
        return None


if __name__ == "__main__":
    polygon = [
    [0, 0, 0],
    [0, 1, 0],
    [1, 1, 0],
    [1, 0, 0],
    [0.5, 0.5, 1]
    ]
    box = [
    [-0.5, -0.5, -0.5],
    [1.5, 1.5, 1.5]
    ]
    new_polygon = clip_polygon_by_box(polygon, box)
    print(new_polygon)

