import numpy as np
import math

""" Daniel Leitner 2018 """


def rotation_matrix(axis, theta):  # from stackoverflow
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])


class Sdf:
    """Base class for all sdf 
    """
    far = 1e9

    def f(self, p):  # signed distance function
        return -self.far

    def bbox(self):  # bounding box of the geometry
        return np.array([-self.far, -self.far, -self.far, self.far, self.far, self.far])

    # def bbox_poitns


class Box(Sdf):
    """Box defined by box = [minx, miny, minz, maxx, maxy, maxz]  
    """

    def __init__(self, box):
        self.b = np.zeros(3,)  # half widths of box
        self.mid = np.zeros(3,)  # mid of box
        for i in range(0, 3):
            self.b[i] = 0.5 * (box[i + 3] - box[i])
            self.mid[i] = 0.5 * (box[i + 3] + box[i])

    def f(self, p):
        m = np.minimum(self.b[0] + (p[:, 0] - self.mid[0]), self.b[0] - (p[:, 0] - self.mid[0]))
        for i in range(1, 3):
            m = np.minimum(m, self.b[i] + (p[:, i] - self.mid[i]))
            m = np.minimum(m, self.b[i] - (p[:, i] - self.mid[i]))
        return -m

    def bbox(self):
        return np.array([self.mid[0] - self.b[0], self.mid[1] - self.b[1], self.mid[2] - self.b[2],
                         self.mid[0] + self.b[0], self.mid[1] + self.b[1], self.mid[2] + self.b[2]])


class Ball(Sdf):
    """ Ball defined by radius and mid point
    """

    def __init__(self, r = 1, pos = np.zeros((1, 3))):
        self.r2 = r * r
        self.pos = pos

    def f(self, p):
        p = np.subtract(p, self.pos)
        return np.sqrt((p ** 2).sum(1) - self.r2 * np.ones((p.shape[0],)))

    def bbox(self):
        return np.array(self.pos[0] - np.sqrt(self.r2), self.pos[1] - np.sqrt(self.r2), self.pos[2] - np.sqrt(self.r2),
                        self.pos[0] + np.sqrt(self.r2), self.pos[1] + np.sqrt(self.r2), self.pos[2] + np.sqrt(self.r2))


class Plant_Container(Sdf):
    """ Square or circular container, given by top radius (rt), bottom radius (rb), height (h), and if it is square or not (square)
    """

    def __init__(self, rt = 0.05, rb = 0.05, h = 1, square = False):
        self.rt = rt
        self.rb = rb
        self.h = h
        self.square = square

    def f(self, p):
        z = p[:, 2] / self.h  # scale 0..1
        r = (1 + z) * self.rb - z * self.rt
        d = 0.
        if self.square:
            d = np.maximum(np.abs(p[:, 0]), np.abs(p[:, 1])) - r
        else:  # round pot
            d = np.sqrt(np.square(p[:, 0]) + np.square(p[:, 1])) - r
        return np.maximum(d, -np.minimum(self.h - p[:, 2], 0. + p[:, 2]))

    def bbox(self):
        m = max(self.rt, self.rb)
        return np.array([-m, -m, 0., m, m, self.h])


class Rotate_Translate(Sdf):
    """ Rotates and or translates a given sdf (sdf). Translates for a vector (p), and rotates for angle (alpha) around an axis (axis)
    """

    def __init__(self, sdf, p = np.zeros((1, 3)), alpha = 0, axis = 0):
        self.sdf = sdf
        self.p = p
        a = np.zeros((3,))
        a[axis] = 1
        self.A = rotation_matrix(a, alpha / 180 * math.pi)

    def f(self, p):
        p1 = np.subtract(p, self.p)
        B = np.transpose(p1)
        p2 = np.matmul(self.A, B)
        return self.sdf.f(np.transpose(p2));

    def bbox(self):
        m = self.sdf.bbox()
        m[0:3] += self.p
        m[3:6] += self.p
        return m


class Union(Sdf):
    """ Union of a list of sdfs
    """

    def __init__(self, sdfs):
        self.sdfs_ = sdfs

    def f(self, p):
        d = self.sdfs_[0].f(p)
        for i in range(1, len(self.sdfs_)):
            d = np.minimum(d, self.sdfs_[i].f(p))
        return d

    def bbox(self):  # TODO test
        m = self.sdfs_[0].bbox()
        for i in range(1, len(self.sdfs_)):
            m2 = self.sdfs_[i].bbox()
            m = np.hstack((np.minimum(m[0:3], m2[0:3]), np.maximum(m[3:6], m2[3:6])))
        return m


class Intersection(Sdf):
    """ Intersection of a list of sdfs
    """

    def __init__(self, sdfs):
        self.sdfs_ = sdfs

    def f(self, p):
        d = self.sdfs_[0].f(p)
        for i in range(1, len(self.sdfs_)):
            d = max(d, self.sdfs_[i].f(p))
        return d

    def bbox(self):  # TODO test
        m = self.sdfs_[0].bbox()
        for i in range(1, len(self.sdfs_)):
            m2 = self.sdfs_[i].bbox()
            m = np.hstack(np.maximum(m[0:3], m2[0:3]), np.minimum(m[3, 6], m2[3, 6]))
        return m


class Difference(Sdf):
    """ Difference between two sdfs (a\b)    
    """

    def __init__(self, a, b):
        self.a = a
        self.b = b

    def f(self, p):
        d = self.a.f(p)
        return  np.maximum(d, -self.b.f(p))

    def bbox(self):  # TODO test
        return self.a.bbox()


class Complement():
    """ Complement of a given sdf  
    """

    def __init__(self, sdf):
        self.sdf = sdf

    def f(self, p):
        return -self.sdf.f(p)

    def bbox(self):
        return self.sdf.bbox()


class Edge_Length(Sdf):
    """ Calculates the edge length function from sdfs
    """

    def __init__(self, sdfs, h_min = 0.5, h_max = 1. , slope = 1.):
        self.sdfs = sdfs
        self.h_min = h_min
        self.h_max = h_max
        self.slope = slope

    def f(self, p):
        d = np.absolute(self.sdfs[0].f(p))
        for i in range(0, len(self.sdfs)):
            d = np.minimum(d, np.absolute(self.sdfs[i].f(p)))
        t = np.minimum(d / self.slope, np.ones(d.shape))
        return t * self.h_max + (1 - t) * self.h_min

    def bbox(self):  # TODO test (same as Unition)
        m = self.sdfs_[0].bbox()
        for i in range(1, len(self.sdfs_)):
            m2 = self.sdfs_[i].bbox()
            m = np.hstack((np.minimum(m[0:3], m2[0:3]), np.maximum(m[3:6], m2[3:6])))
        return m

