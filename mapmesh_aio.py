# -*- coding: utf-8 -*-

import sys
import os
import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
import scipy.io
from numpy.linalg import det
from math import floor
from LoadDatainMatFile import LoadNumericMatrixInMatFile
import warnings

def findintersection(L1, L2):
    """
    find the intersection point of two INFINITE lines
    :param L1: 2d array, 2*2, = [[x11, y11], [x12, y12]], coordinates of two ends of a line
    :param L2: 2d array, 2*2, = [[x21, y21], [x22, y22]], coordinates of two ends of a line
    :return: coordinates of one intersection point (possible to be on the extensions of lines), 1d array
    """

    # coordinates of ends of line 1
    x1, y1 = L1[0, :]
    x2, y2 = L1[1, :]

    # coordinates of ends of line 2
    x3, y3 = L2[0, :]
    x4, y4 = L2[1, :]

    # check if L1 and L2 are parellel
    t1 = np.array([x2 - x1, y2 - y1])
    t2 = np.array([x4 - x3, y4 - y3])

    # extend t1 and t2 to 3d
    t1 = np.concatenate((t1, [0]))
    t2 = np.concatenate((t2, [0]))

    # raise a error if L1 is parellel to L2
    assert np.linalg.norm(np.cross(t1, t2)) > np.finfo(float).eps, 'Parallel lines for findintersection'

    # find intersection point
    x = det([[det([[x1, y1], [x2, y2]]), x1 - x2], [det([[x3, y3], [x4, y4]]), x3 - x4]]) / det(
        [[x1 - x2, y1 - y2], [x3 - x4, y3 - y4]])
    y = det([[det([[x1, y1], [x2, y2]]), y1 - y2], [det([[x3, y3], [x4, y4]]), y3 - y4]]) / det(
        [[x1 - x2, y1 - y2], [x3 - x4, y3 - y4]])

    return np.array([x, y])


def replace(arr, old, new):
    """
    replace old elements in arr by new
    :param arr: 2d array
    :param old: 1d array
    :param new: 1d array
    :return: an other 2d array rather then in place change
    """
    poslist = []
    vallist = []

    for i in range(old.shape[0]):
        oldval = old[i]
        newval = new[i]

        pos = np.where(arr == oldval)

        for (i, j) in zip(pos[0], pos[1]):
            poslist.append([i, j])
            vallist.append(newval)

    # perform replacement
    res = arr.copy()
    for i in range(len(poslist)):
        res[poslist[i][0]][poslist[i][1]] = vallist[i]

    return res


class HOLTRICS:
    """
    class of hollow triangular cross section
    """

    def __init__(self, xc, x0, x1, t):
        """
        constructor
        :param xc:  list, =[x, y], coordintes of centroid
        :param x0:  list, =[x, y], coordinates of vertex 1
        :param x1:  list, =[x, y], coordinates of vertex 2
        :param t:   float, wall thickness
        """
        self.xc = np.array(xc)  # 1d array
        self.t = t

        # coordinates of external triangle
        self.x0 = np.array(x0)
        self.x1 = np.array(x1)
        self.x2 = np.zeros(2)
        self.setx2()            # compute coordinates of vertex 3

        # coordinates of internal triangle, x[N]i is nearest to x[N]
        self.x0i = np.zeros(2)
        self.x1i = np.zeros(2)
        self.x2i = np.zeros(2)
        self.setxi()  # compute vertices of internal edges

    def setx2(self):
        """
        compute coordinates of vertex 3,
        xc = (x0+x1+x2)/3, 1,2,3 are anti-clockwise
        :return:
        """
        self.x2 = 3 * self.xc - self.x0 - self.x1
        return self

    def setxi(self):
        """
        This function moves each edge of a 2d triangle inward by a thickness t and returns a triangular ring
        :return: none
        """
        # = [[x0, y0], [x1, y1], [x2, y2]], vertices of outer triangle
        Po = np.concatenate((self.x0.reshape((1, 2)), self.x1.reshape((1, 2)), self.x2.reshape((1, 2))), axis=0)

        # link head and tail for vertices of outer triangle
        Pout = np.concatenate((Po, Po[0, :].reshape((1, 2))), axis=0)

        lines = {}  # init empty dict
        for i in range(3):
            P1 = Pout[i, :]  # 1d array
            P2 = Pout[i + 1, :]  # 1d array

            dx, dy = P2 - P1

            # normal vector to the left
            n = np.array([-dy, dx])
            n = n / np.linalg.norm(n)  # normalized directional vector

            # translational vector
            disp = n * self.t

            # translate two nodes to obtain translated line
            P1n = P1 + disp  # 1d array
            P2n = P2 + disp  # 1d array

            # store the translated line into dict lines
            lines[i + 1] = np.concatenate((P1n.reshape((1, 2)), P2n.reshape((1, 2))), axis=0)

        # extend lines
        lines[0] = lines[3]

        # Pi stores coordinates of vertices of howllow edge
        Pi = np.zeros((3, 2))
        for i in range(3):
            x, y = findintersection(lines[i], lines[i + 1])
            Pi[i, 0], Pi[i, 1] = (x, y)

        self.x0i = Pi[0, :]     # nearest to x0
        self.x1i = Pi[1, :]     # nearest to x1
        self.x2i = Pi[2, :]     # nearest to x2
        return self

    def draw(self):
        """
        plot triangular cross section
        :return:
        """
        import matplotlib.pyplot as plt

        # plot external boundary
        x = [eval(str, {"self": self}) for str in ['self.x' + str(i) + '[0]' for i in range(0, 3)]]
        x.append(x[0])

        y = [eval(str, {"self": self}) for str in ['self.x' + str(i) + '[1]' for i in range(0, 3)]]
        y.append(y[0])

        plt.plot(x, y, 'go-', label='external', linewidth=2)

        # plot internal boundary
        x = [eval(str, {"self": self}) for str in ['self.x' + str(i) + 'i[0]' for i in range(0, 3)]]
        x.append(x[0])

        y = [eval(str, {"self": self}) for str in ['self.x' + str(i) + 'i[1]' for i in range(0, 3)]]
        y.append(y[0])

        plt.plot(x, y, 'ro-', label='internal', linewidth=2)

        plt.xlabel('Local x axis')
        plt.ylabel('Local y axis')
        plt.grid(True)

    def output(self, filename='shift2centroid.mat'):
        """
        generate mat file storing Pi & Po for cross section mesh
        :return:
        """
        import scipy.io

        Po = np.array([self.x0, self.x1, self.x2]) - np.tile(self.xc, (3, 1))
        Pi = np.array([self.x0i, self.x1i, self.x2i]) - np.tile(self.xc, (3, 1))

        Po = Po.astype('float')
        Pi = Pi.astype('float')

        scipy.io.savemat(filename, {'Po': Po, 'Pi': Pi, 't': self.t})

        return

    def decomp2quads(self):
        """
        decompose this cross section into quadrilateral sub-parts
        :return: dict, {'points': 2d array, 'connectivity': 2d array}
        """
        # projection of x0i
        x0i_in_x0x1 = project2line(self.x0i, np.array([self.x0, self.x1]))
        x0i_in_x2x0 = project2line(self.x0i, np.array([self.x2, self.x0]))

        # projection of x1i
        x1i_in_x0x1 = project2line(self.x1i, np.array([self.x0, self.x1]))
        x1i_in_x1x2 = project2line(self.x1i, np.array([self.x1, self.x2]))

        # projection of x2i
        x2i_in_x1x2 = project2line(self.x2i, np.array([self.x1, self.x2]))
        x2i_in_x2x0 = project2line(self.x2i, np.array([self.x2, self.x0]))

        quads = {}
        quads['points'] = np.array([self.x0,
                                    x0i_in_x0x1,
                                    x1i_in_x0x1,
                                    self.x1,
                                    x1i_in_x1x2,
                                    x2i_in_x1x2,
                                    self.x2,
                                    x2i_in_x2x0,
                                    x0i_in_x2x0,
                                    self.x0i,
                                    self.x1i,
                                    self.x2i
                                    ])
        quads['connectivity'] = np.array([[0, 1, 9, 8],
                                          [1, 2, 10, 9],
                                          [3, 4, 10, 2],
                                          [4, 5, 11, 10],
                                          [6, 7, 11, 5],
                                          [7, 8, 9, 11]
                                          ], dtype=int)

        return quads

    def decomp2tris(self):
        """
        decompose this cross section into triangular sub-parts
        :return: dict, {'points': 2d array, 'connectivity': 2d array}
        """
        # projection of x0i
        x0i_in_x0x1 = project2line(self.x0i, np.array([self.x0, self.x1]))
        x0i_in_x2x0 = project2line(self.x0i, np.array([self.x2, self.x0]))

        # projection of x1i
        x1i_in_x0x1 = project2line(self.x1i, np.array([self.x0, self.x1]))
        x1i_in_x1x2 = project2line(self.x1i, np.array([self.x1, self.x2]))

        # projection of x2i
        x2i_in_x1x2 = project2line(self.x2i, np.array([self.x1, self.x2]))
        x2i_in_x2x0 = project2line(self.x2i, np.array([self.x2, self.x0]))

        tris = {}
        tris['points'] = np.array([self.x0,
                                    x0i_in_x0x1,
                                    x1i_in_x0x1,
                                    self.x1,
                                    x1i_in_x1x2,
                                    x2i_in_x1x2,
                                    self.x2,
                                    x2i_in_x2x0,
                                    x0i_in_x2x0,
                                    self.x0i,
                                    self.x1i,
                                    self.x2i
                                    ])
        tris['connectivity'] = np.array([[0, 1, 9],    # mid point is right angle
                                          [9, 1, 2],
                                          [2, 10, 9],
                                          [10, 2, 3],
                                          [3, 4, 10],
                                          [10, 4, 5],
                                          [5, 11, 10],
                                          [11, 5, 6],
                                          [6, 7, 11],
                                          [11, 7, 8],
                                          [8, 9, 11],
                                          [9, 8, 0]], dtype=int)

        return tris

    def decomp2tris_ar(self):
        """
        decompose this cross section into triangular sub-parts, aspect ratio of sub-parts is optimized
        :return: dict, {'points': 2d array, 'connectivity': 2d array}
        """
        # projection of x0i
        x0i_in_x0x1 = project2line(self.x0i, np.array([self.x0, self.x1]))
        x0i_in_x2x0 = project2line(self.x0i, np.array([self.x2, self.x0]))

        # projection of x1i
        x1i_in_x0x1 = project2line(self.x1i, np.array([self.x0, self.x1]))
        x1i_in_x1x2 = project2line(self.x1i, np.array([self.x1, self.x2]))

        # projection of x2i
        x2i_in_x1x2 = project2line(self.x2i, np.array([self.x1, self.x2]))
        x2i_in_x2x0 = project2line(self.x2i, np.array([self.x2, self.x0]))

        tris = {}
        tris['points'] = np.array([self.x0,
                                   x0i_in_x0x1,
                                   x1i_in_x0x1,
                                   self.x1,
                                   x1i_in_x1x2,
                                   x2i_in_x1x2,
                                   self.x2,
                                   x2i_in_x2x0,
                                   x0i_in_x2x0,
                                   self.x0i,
                                   self.x1i,
                                   self.x2i
                                   ])

        # triangles in the corner of ths cross section
        tris['connectivity'] = np.array([[0, 1, 9],  # mid point is right angle
                                         [10, 2, 3],
                                         [3, 4, 10],
                                         [11, 5, 6],
                                         [6, 7, 11],
                                         [9, 8, 0]], dtype=int)

        # generates seeds on rectangular subparts of the cross section
        rectangles = [[self.x0i, x0i_in_x0x1, x1i_in_x0x1, self.x1i],
                      [self.x1i, x1i_in_x1x2, x2i_in_x1x2, self.x2i],
                      [self.x2i, x2i_in_x2x0, x0i_in_x2x0, self.x0i]]

        for i in range(len(rectangles)):
            P0 = rectangles[i][0]
            P1 = rectangles[i][1]
            P2 = rectangles[i][2]
            P3 = rectangles[i][3]

            # compute aspect ratio
            aspratio = np.linalg.norm(P1 - P2)/self.t

            # number of segments to the long edge of the rectangle
            nsegments = floor(aspratio)

            # generates seeds on long edge
            seed1 = [[x, y] for (x, y) in zip(np.linspace(P1[0], P2[0], nsegments+1),
                                              np.linspace(P1[1], P2[1], nsegments+1))]
            seed2 = [[x, y] for (x, y) in zip(np.linspace(P0[0], P3[0], nsegments+1),
                                              np.linspace(P0[1], P3[1], nsegments+1))]
            # generates points for Delaunay triangulation
            points = np.concatenate((np.array(seed1), np.array(seed2)), axis=0)
            tri = Delaunay(points)

            # append to tris['points'] and tris['connectivity']
            tris['connectivity'] = np.concatenate((tris['connectivity'], tri.simplices + tris['points'].shape[0]), axis=0)
            tris['points'] = np.concatenate((tris['points'], points), axis=0)

        # eliminate repeated nodes in tris['points'] and fix references in tris['connectivity']
        # heuristically determine round precision for node coordinates
        trunc = int(np.format_float_scientific(self.t).split('e')[-1]) - 3

        # extract unique nodes out of allpoints
        unique, unique_indices, unique_inverse = np.unique(tris['points'].round(decimals=abs(trunc)), return_index=True,
                                                           return_inverse=True, axis=0)
        # update references
        tritab = replace(tris['connectivity'], np.arange(unique_inverse.shape[0], dtype=int), unique_inverse)

        tris['points'] = unique
        tris['connectivity'] = tritab

        # # test by plotting
        # plt.triplot(tris['points'][:, 0], tris['points'][:, 1], tris['connectivity'])
        # plt.show()

        return tris

    def mapmesh(self, deltri, vertices, aspratio=False):
        """
        map mesh in deltri onto subparts of this hollow triangular cross section
        :param deltri:  Delaunay instance, deltri.points and deltri.simplices
        :param vertices: 2d array, [[x1, y1], [x2, y2], [x3, y3]] for triangle,
                                   [[x1, y1], [x2, y2], [x3, y3], [x4, y4]] for quadrilateral
        :param aspratio: bool, = True, divide rectangular subparts into more triangles
        :return:  unique: 2d arry, each row is coordinates of a node
                  alltris: 2d array, each row is three node ids (starting from 0) that form a triangle
        """
        # decompose this hollow triangular cross section into sub-parts of triangles/quadrilaterals
        if vertices.shape[0] == 3:
            if aspratio:
                decomp = self.decomp2tris_ar()
            else:
                decomp = self.decomp2tris()
        elif vertices.shape[0] == 4:
            decomp = self.decomp2quads()
        else:
            raise ValueError('incorrect vertices!\n')

        for i in range(decomp['connectivity'].shape[0]):                # for each quadrilateral/triangular subpart
            connec = decomp['connectivity'][i]                          # 1d array of node ids
            polygon = [decomp['points'][j].tolist() for j in connec]    # list of vertex coordinates

            # generate transformation function
            if vertices.shape[0] == 3:
                tranform = create_affine_transform(vertices, polygon)
            elif vertices.shape[0] == 4:
                tranform = create_perspective_transform(vertices, polygon)
            else:
                raise ValueError('incorrect vertices!\n')

            # perform transformation for points in deltri
            points_upd = np.array([tranform(x) for x in deltri.points]) # 2d array, each row is nodal coordinates

            # # plot transformed triangulation
            # import matplotlib.pyplot as plt
            # plt.triplot(points_upd[:, 0], points_upd[:, 1], deltri.simplices)
            # # plt.plot(points_upd[:, 0], points_upd[:, 1], 'o')
            # plt.show()

            if i == 0:
                allpoints = points_upd
                alltris = deltri.simplices
            else:
                allpoints = np.concatenate((allpoints, points_upd), axis=0)
                alltris = np.concatenate((alltris, deltri.simplices + np.max(alltris) + 1), axis=0)

        # find minimum length of edge in alltris
        Lmin = np.inf
        for i in range(alltris.shape[0]):
            for j in range(3):
                nodeid1 = alltris[i][j]
                if j == 2:
                    nodeid2 =  alltris[i][0]
                else:
                    nodeid2 = alltris[i][j+1]

                L = np.linalg.norm(allpoints[nodeid1] - allpoints[nodeid2])

                if L < Lmin:
                    Lmin = L

        # heuristically determine round precision for node coordinates
        trunc = int(np.format_float_scientific(Lmin).split('e')[-1]) - 3

        # extract unique nodes out of allpoints
        unique, unique_indices, unique_inverse = np.unique(allpoints.round(decimals=abs(trunc)), return_index=True,
                                                           return_inverse=True, axis=0)
        # update reference
        alltris = replace(alltris, np.arange(unique_inverse.shape[0], dtype=int), unique_inverse)

        return unique, alltris


def project2line(pt, line):
    """
    project a 2d point onto a 2d INFINITE line
    :param pt: 1d array, =[x, y], coordinates of a point
    :param line: 2d array, =[[x0, y0], [x1, y1]], coordinates of two ends of a line
    :return: 1d array, =[xp, yp], coordinates of pt's projection on the line
    """
    tau = (line[1] - line[0]) / np.linalg.norm(line[1] - line[0])   # line direction
    r = pt - line[0]
    proj = r.dot(tau)
    return line[0] + tau*proj

# -------------------------------------

def plotmesh(allpoints, alltris):
    """
    plot mesh of 2d triangles
    :param allpoints: 2d array of float, each row is coordinates of a point
    :param alltris: 2d array of int, each row is node ids of a triangle
    :return: none
    """
    plt.triplot(allpoints[:, 0], allpoints[:, 1], alltris)
    plt.axis("square")
    plt.show()
    return


def findbound(pts, tris):
    """
    find boundary edges for a Delaunay mesh
    :param pts: 2d array of float, each row is coordinates of a point
    :param tris: 2d array of int, each row is nodes ids of a triangle
    :return: bound_edge: 2d array of int, each row is node ids
             bgp: 1d array of int, ids of nodes on the external boundary
    """
    ntris = tris.shape[0]   # number of triangles
    npts  = pts.shape[0]    # number of points

    edges = []              # record edges

    # check if all triangles are anticlockwise
    for i in range(ntris):
        iP0 = tris[i][0]
        iP1 = tris[i][1]
        iP2 = tris[i][2]

        determinant = np.linalg.det([[pts[iP0][0], pts[iP0][1], 1],
                                     [pts[iP1][0], pts[iP1][1], 1],
                                     [pts[iP2][0], pts[iP2][1], 1]])
        if determinant < 0:
            # adjust node sequence if clockwise
            tmp = tris[i][1]
            tris[i][1] = tris[i][2]
            tris[i][2] = tmp
            print(['triangle i = ', i, ' is reverted.\n'])
        elif determinant == 0:
            raise RuntimeError('colinear nodes for triangle')

        # register edge
        edges.append([tris[i][0], tris[i][1]])
        edges.append([tris[i][1], tris[i][2]])
        edges.append([tris[i][2], tris[i][0]])

    edges = np.array(edges, dtype=int)      # change edges from a list to a 2d array
    edges_sort = np.sort(edges, axis=-1)    # sort each row in ascending order
    # find non-repeated edges in
    unique, unique_indices, unique_inverse, unique_counts = np.unique(edges_sort, return_index=True, return_inverse=True, return_counts=True, axis=0)

    outedges = edges[unique_indices[unique_counts == 1]]

    # # test by plotting outedges
    # for i in range(outedges.shape[0]):
    #     plt.plot([pts[outedges[i][0]][0], pts[outedges[i][1]][0]], [pts[outedges[i][0]][1], pts[outedges[i][1]][1]], '-o')
    # plt.show()

    bound_edge = outedges
    bgp = np.unique(outedges)
    return bound_edge, bgp


def create_perspective_transform_matrix(src, dst):
    """
    Creates a perspective transformation matrix which transforms points
    in quadrilateral "src" to the corresponding points on quadrilateral
    "dst".

    Will raise a "np.linalg.LinAlgError" on invalid input.
    :param src: list of tuple, [(x1,y1), (x2,y2), (x3,y3), (x4,y4)]
    :param dst: list of tuple, [(X1,Y1), (X2,Y2), (X3,Y3), (X4,Y4)]
    :return: 2d array, 3*3
    """
    # See:
    # * http://stackoverflow.com/a/14178717/71522
    in_matrix = []
    for (x, y), (X, Y) in zip(src, dst):
        in_matrix.extend([
            [x, y, 1, 0, 0, 0, -X * x, -X * y],
            [0, 0, 0, x, y, 1, -Y * x, -Y * y],
        ])                                          # extend list of list

    A = np.matrix(in_matrix, dtype=float)        # make list of list into 2d array, 8*8
    B = np.array(dst).reshape(8)                    # flat 2d array to 1d array using row priority, B = [X1,Y1,X2,Y2,X3,Y3,X4,Y4]
    af = np.dot(np.linalg.inv(A.T * A) * A.T, B)    # .T: transpose, af: 1d array, 8
    return np.append(np.array(af).reshape(8), 1).reshape((3, 3))


def create_perspective_transform(src, dst, round=False, splat_args=False):
    """
    Returns a function which will transform points in quadrilateral
    "src" to the corresponding points on quadrilateral "dst":

    e.g. point (5,5) in quadrilateral 1 --> point (74.99999999999639,74.999999999999957) in quadrilateral 2
    >>> transform = create_perspective_transform(
    ...     [(0, 0), (10, 0), (10, 10), (0, 10)],
    ...     [(50, 50), (100, 50), (100, 100), (50, 100)],
    ... )
    >>> transform((5, 5))
    (74.99999999999639, 74.999999999999957)

    If round is True then points will be rounded to the nearest integer and integer values will be returned.

    e.g. point (5,5) in quadrilateral 1 --> point (75,75) in quadrilateral 2
    >>> transform = create_perspective_transform(
    ...     [(0, 0), (10, 0), (10, 10), (0, 10)],
    ...     [(50, 50), (100, 50), (100, 100), (50, 100)],
    ...     round=True,
    ... )
    >>> transform((5, 5))
    (75, 75)

    If splat_args is True the function will accept two arguments instead of a tuple.

    >>> transform = create_perspective_transform(
    ...     [(0, 0), (10, 0), (10, 10), (0, 10)],
    ...     [(50, 50), (100, 50), (100, 100), (50, 100)],
    ...     splat_args=True,
    ... )
    >>> transform(5, 5)
    (74.99999999999639, 74.999999999999957)

    If the input values yield an invalid transformation matrix an identity
    function will be returned and the error attribute will be set to a
    description of the error:

    >>> tranform = create_perspective_transform(
    ...     np.zeros((4, 2)),
    ...     np.zeros((4, 2)),
    ... )
    >>> transform((5, 5))
    (5.0, 5.0)
    >>> transform.error
    'invalid input quads (...): Singular matrix

    :param src: list of tuple, [(x1,y1), (x2,y2), (x3,y3), (x4,y4)]
    :param dst: list of tuple, [(X1,Y1), (X2,Y2), (X3,Y3), (X4,Y4)]
    :param round: bool
    :param splat_args: bool
    :return: function
    """
    try:
        transform_matrix = create_perspective_transform_matrix(src, dst)    # 3*3 2d array
        error = None
    except np.linalg.LinAlgError as e:
        transform_matrix = np.identity(3, dtype=np.float)                   # default transform_matrix to be identity matrix
        error = "invalid input quads (%s and %s): %s" %(src, dst, e)
        error = error.replace("\n", "")
        raise ValueError(error)

    # factory function
    to_eval = "def perspective_transform(%s):\n" %(splat_args and "*pt" or "pt",)   # function declaration, *pt collects all arguments into a list
    to_eval += "  res = np.dot(transform_matrix, ((pt[0], ), (pt[1], ), (1, )))\n"  # ((pt[0], ), (pt[1], ), (1, )) is column vector
    to_eval += "  res = res / res[2]\n"                                             # res is 3*1 2d array, normalize last item in res

    # set return of function
    if round:
        to_eval += "  return (int(round(res[0][0])), int(round(res[1][0])))\n"
    else:
        to_eval += "  return (res[0][0], res[1][0])\n"

    locals = {"transform_matrix": transform_matrix,}
    locals.update(globals())
    exec(to_eval, locals, locals)   # replace exec to_eval in locals, locals
    res = locals["perspective_transform"]
    res.matrix = transform_matrix
    res.error = error
    return res


def create_affine_transform_matrix(src, dst):
    """
    Creates a perspective transformation matrix which transforms points
    in quadrilateral "src" to the corresponding points on quadrilateral
    "dst".

    x = a1 * r + a2 * s + a3
    y = b1 * r + b2 * s + b3

    Will raise a "np.linalg.LinAlgError" on invalid input.
    :param src: list of tuple, [(r1,s1), (r2,s2), (r3,s3)]
    :param dst: list of tuple, [(x1,y1), (x2,y2), (x3,y3)]
    :return: A: 2d array, 2*2
             b: 1d array, 2
    """
    r1, s1 = src[0]
    r2, s2 = src[1]
    r3, s3 = src[2]

    x1, y1 = dst[0]
    x2, y2 = dst[1]
    x3, y3 = dst[2]

    a1, a2, a3 = np.linalg.solve(np.array([[r1, s1, 1], [r2, s2, 1], [r3, s3, 1]], dtype=float), np.array([x1, x2, x3], dtype=float))
    b1, b2, b3 = np.linalg.solve(np.array([[r1, s1, 1], [r2, s2, 1], [r3, s3, 1]], dtype=float), np.array([y1, y2, y3], dtype=float))

    A = np.array([[a1, a2], [b1, b2]], dtype=float)
    b = np.array([a3, b3], dtype=float)
    return A, b         # [x, y]' = A*[r, s]' + b


def create_affine_transform(src, dst, round=False, splat_args=False):
    """
    Returns a function which will transform points in triangle
    "src" to the corresponding points on triangle "dst":

    e.g. point (5,5) in quadrilateral 1 --> point (74.99999999999639,74.999999999999957) in quadrilateral 2
    >>> transform = create_perspective_transform(
    ...     [(0, 0), (10, 0), (10, 10), (0, 10)],
    ...     [(50, 50), (100, 50), (100, 100), (50, 100)],
    ... )
    >>> transform((5, 5))
    (74.99999999999639, 74.999999999999957)

    If round is True then points will be rounded to the nearest integer and integer values will be returned.

    e.g. point (5,5) in quadrilateral 1 --> point (75,75) in quadrilateral 2
    >>> transform = create_perspective_transform(
    ...     [(0, 0), (10, 0), (10, 10), (0, 10)],
    ...     [(50, 50), (100, 50), (100, 100), (50, 100)],
    ...     round=True,
    ... )
    >>> transform((5, 5))
    (75, 75)

    If splat_args is True the function will accept two arguments instead of a tuple.

    >>> transform = create_perspective_transform(
    ...     [(0, 0), (10, 0), (10, 10), (0, 10)],
    ...     [(50, 50), (100, 50), (100, 100), (50, 100)],
    ...     splat_args=True,
    ... )
    >>> transform(5, 5)
    (74.99999999999639, 74.999999999999957)

    If the input values yield an invalid transformation matrix an identity
    function will be returned and the error attribute will be set to a
    description of the error:

    >>> tranform = create_perspective_transform(
    ...     np.zeros((4, 2)),
    ...     np.zeros((4, 2)),
    ... )
    >>> transform((5, 5))
    (5.0, 5.0)
    >>> transform.error
    'invalid input quads (...): Singular matrix

    :param src: list of tuple, [(r1,s1), (r2,s2), (r3,s3)]
    :param dst: list of tuple, [(x1,y1), (x2,y2), (x3,y3)]
    :param round: bool
    :param splat_args: bool
    :return: function
    """
    try:
        transform_matrix1, transform_matrix2 = create_affine_transform_matrix(src, dst)
        error = None
    except np.linalg.LinAlgError as e:
        transform_matrix1 = np.identity(2, dtype=float)
        transform_matrix2 = np.zeros((2, 1))
        error = "invalid input tris (%s and %s): %s" %(src, dst, e)
        error = error.replace("\n", "")
        raise ValueError(error)

    # below generate function
    # function declaration, *pt collects all arguments into a list
    to_eval = "def affine_transform(%s):\n" %(splat_args and "*pt" or "pt",)
    # ((pt[0], ), (pt[1], ), (1, )) is column vector
    to_eval += "  res = np.dot(transform_matrix1, np.array(pt)) + transform_matrix2\n"

    # set return of function
    if round:
        to_eval += "  return (int(round(res[0][0])), int(round(res[1][0])))\n"
    else:
        to_eval += "  return (res[0], res[1])\n"

    locals = {"transform_matrix1": transform_matrix1,
              "transform_matrix2": transform_matrix2}
    locals.update(globals())
    exec(to_eval, locals, locals)   # replace exec to_eval in locals, locals
    res = locals["affine_transform"]
    res.matrix1 = transform_matrix1     # function attribute
    res.matrix2 = transform_matrix2
    res.error = error
    return res


if __name__=='__main__':

    import sys
    try:
        workdir = sys.argv[1]       # designate path where CroSecXX folders exist
    except:
        workdir = '.'


    # generate Delaunay triangulation in a reference triangle
    vertices = np.array([(0, 1), (0, 0), (1, 0)], dtype=float)
    x = np.linspace(0, 1, 5)        # !!! here set the number of seeds on a edge in the reference triangle
    y = x
    xv, yv = np.meshgrid(x, y, sparse=False, indexing='ij')
    xy = [[xcoor, ycoor] for (xcoor, ycoor) in zip(xv.flatten(), yv.flatten()) if (xcoor + ycoor) <= 1]
    points = np.array(xy)       # points in the triangular mesh
    deltri = Delaunay(points)   # triangulation
    # plotmesh(points, deltri)    # plot

    crosecfolders = ['./CroSec' + str(i) for i in range(1, 11)]     # set target folders
    print('Target folders to mesh are: ')
    for folder in crosecfolders:
        print(folder)                       # print info

    for (root, dirs, files) in os.walk(workdir):
        if root in crosecfolders:

            crosecfolders.remove(root)

            print('Meshing ' + root)
            TargetFile = os.path.join(root, 'genTriCros', 'shift2centroid.mat')
            extract = LoadNumericMatrixInMatFile(TargetFile, ['Po', 'thick'])
            Po = extract['Po']
            thick = extract['thick']

            # mesh cross section
            section = HOLTRICS([0, 0], Po[0], Po[1], thick)
            allpoints, alltris = section.mapmesh(deltri, vertices, aspratio=True)
            bound_edge, bgp = findbound(allpoints, alltris)

            try:
                dir = 'tri_mesh'
                os.mkdir(os.path.join(root, dir))    # create subfolder 'tri_mesh'
            except:
                pass                                 # in case folder is existent

            # export to coords_ien_bgp.mat, note that matlab node index starts from 1
            file = 'coords_ien_bgp.mat'
            scipy.io.savemat(os.path.join(root, dir, file),
                             {'coords': allpoints,
                              'ien': alltris.astype('float') + 1,
                              'bound_edge': bound_edge.astype('float') + 1,
                              'bgp': bgp.reshape((-1, 1)).astype('float') + 1})

    if crosecfolders:
        msg = 'Following cross sections are not meshed: '
        for folder in crosecfolders:
            msg = msg + '\n\t' + folder
        raise RuntimeError(msg)
