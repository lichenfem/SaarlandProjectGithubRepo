# -*- coding: utf-8 -*-

import sys
import os
import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
import scipy.io
from numpy.linalg import det, norm
from math import floor
from MiscPyUtilities.line2DIntersect import findintersection
from MiscPyUtilities.triangleGeometry import checkAnticlockwise, TriInscribedCircleRadius
from MiscPyUtilities.affine_perspective_transformation import create_affine_transform
from MiscPyUtilities.project2DPointOnto2DLine import project2DPtOnto2DLine as project2line
from MiscPyUtilities.replaceElementsInArray import replaceIntArray as replace
import warnings


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
        self.setx2()            # compute coordinates of vertex 3
        self.checkValidthicAndAnticlockwise()   # check valid thickness and anticlockwise vertex sequence

        # coordinates of internal triangle, x[N]i is nearest to x[N]
        self.setxi()            # compute vertices of internal edges


    def setx2(self):
        """
        compute coordinates of vertex 3,
        xc = (x0+x1+x2)/3, 1,2,3 are anti-clockwise
        :return:
        """
        self.x2 = 3 * self.xc - self.x0 - self.x1


    def checkValidthicAndAnticlockwise(self):
        Vert = np.array([self.x0, self.x1, self.x2]).flatten()
        r = TriInscribedCircleRadius(Vert)
        if self.t >= r:
            raise RuntimeError('Cross section thickness exceeds radius of inscribed circle!')


        # adjust sequence of node if required
        Anticlockwise = checkAnticlockwise(np.array([self.x0.tolist(),
                                                     self.x1.tolist(),
                                                     self.x2.tolist()
                                                     ]).flatten())
        if not Anticlockwise.all():
            # adjust sequence
            tmp = self.x1
            self.x1 = self.x2
            self.x2 = tmp
            warnings.warn('Node sequence reversed.')


    def setxi(self):
        """
        This function moves each edge of a 2d triangle inward by a thickness t and returns a triangular ring
        :return: none
        """
        # = [[x0, y0], [x1, y1], [x2, y2]], vertices of outer triangle
        Pout = np.array([
            self.x2.tolist(),
            self.x0.tolist(),
            self.x1.tolist(),
            self.x2.tolist(),
                       ])   # 4*2 array

        # translate edges of external bound inward
        edges = []
        for i in range(3):
            P1 = Pout[i, :]
            P2 = Pout[i + 1, :]

            dx, dy = P2 - P1
            n = np.array([-dy, dx])     # normal vector to the left of vector P1->P2
            n = n / np.linalg.norm(n)

            disp = n * self.t           # translational vector

            P1n = P1 + disp             # translate two nodes to obtain translated line
            P2n = P2 + disp

            line = np.array([P1n.tolist(), P2n.tolist()]).flatten().tolist()    # [x1, y1, x2, y2] of a line
            # store the translated line into dict lines
            edges.append(line)

        edges.append(edges[0])

        # Pi stores coordinates of vertices of hollow edge
        Pi = np.zeros((3, 2))
        for i in range(3):
            FlagParallel, intersectPt = findintersection(
                                                        np.array(edges[i]).reshape((2, 2)),
                                                        np.array(edges[i+1]).reshape((2, 2))
                                                        )
            if FlagParallel == 1:
                raise RuntimeError()
            else:
                Pi[i, 0], Pi[i, 1] = intersectPt[0], intersectPt[1]

        # check anticlockwise
        Anticlockwise = checkAnticlockwise(Pi.flatten())
        assert Anticlockwise.all()

        self.x0i = Pi[0, :]
        self.x1i = Pi[1, :]
        self.x2i = Pi[2, :]

        # check if x0i is closest to x0, same for x1, x2
        assert norm(self.x0i - self.x0) < norm(self.x1i - self.x0) and \
               norm(self.x0i - self.x0) < norm(self.x2i - self.x0), 'Among (x0i, x1i, x2i), x0i should be closest to x0'
        assert norm(self.x1i - self.x1) < norm(self.x0i - self.x1) and \
               norm(self.x1i - self.x1) < norm(self.x2i - self.x1), 'Among (x0i, x1i, x2i), x1i should be closest to x1'
        assert norm(self.x2i - self.x2) < norm(self.x0i - self.x2) and \
               norm(self.x2i - self.x2) < norm(self.x1i - self.x2), 'Among (x0i, x1i, x2i), x2i should be closest to x2'


    def decomp2quads(self):
        """
        decompose this cross section into quadrilateral sub-parts
        :return: dict, {'points': 2d array, 'connectivity': 2d array}
        """
        raise DeprecationWarning
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
        # projection of x0i on edges x0x1 & x2x0
        x0i_in_x0x1 = project2line(self.x0i, np.array([self.x0, self.x1]))
        x0i_in_x2x0 = project2line(self.x0i, np.array([self.x2, self.x0]))

        # projection of x1i on edges x0x1 & x1x2
        x1i_in_x0x1 = project2line(self.x1i, np.array([self.x0, self.x1]))
        x1i_in_x1x2 = project2line(self.x1i, np.array([self.x1, self.x2]))

        # projection of x2i on edges x1x2 and x2x0
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

        # triangles in the corner of the cross section, 2 for each corner
        tris['connectivity'] = np.array([
                                        [0, 1, 9],  # mid point is right angle
                                        [10, 2, 3],
                                        [3, 4, 10],
                                        [11, 5, 6],
                                        [6, 7, 11],
                                        [9, 8, 0]
                                         ], dtype=int)

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

            if aspratio < 1:
                nsegments = 1

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
        :param deltri:  Delaunay instance, deltri.points and deltri.simplices, of reference triangle/quadrilateral
        :param vertices: 2d array, [[x1, y1], [x2, y2], [x3, y3]] for triangle,
                                   [[x1, y1], [x2, y2], [x3, y3], [x4, y4]] for quadrilateral, of reference
                                   triangle/quadrilateral
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


def plotmesh(allpoints, alltris):
    """
    plot mesh of 2d triangles
    :param allpoints: 2d array of float, each row is coordinates of a point
    :param alltris: 2d array of int, each row is node ids of a triangle
    :return:
        none
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
            # print(['triangle i = ', i, ' is reverted.\n'])
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
