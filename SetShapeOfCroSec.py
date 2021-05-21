import os
import sys
import numpy as np
from mapmesh_aio import HOLTRICS, findbound
from LoadDatainMatFile import LoadNumericMatrixInMatFile, WriteNumericMatrix2MatFile
import warnings
from scipy.spatial import Delaunay


def cart2pol(x, y):
    """
    cartisian coordinates to polar coordinates
    @param x: float
    @param y: float
    @return: t
        uple of (rho, phi), phi in [0, 2*pi), angle is in rad
    """
    rho = np.sqrt(x**2 + y**2)

    sinval = y/rho
    cosval = x/rho

    if sinval >= 0 and cosval > 0:
        # 1st quarter
        phi = np.arccos(cosval)         # phi in the range of [0, pi]
        assert 0 <= phi < np.pi / 2, "Error"
    elif sinval > 0 and cosval <= 0:
        # 2nd quarter
        phi = np.arccos(cosval)
        assert np.pi / 2 <= phi < np.pi, "Error"
    elif sinval <= 0 and cosval < 0:
        # 3rd quarter
        phi = np.arccos(cosval)  # phi in the range of [0, pi]
        phi = 2*np.pi - phi
        assert np.pi <= phi < 3/2*np.pi, "Error"
    elif sinval < 0 and cosval >= 0:
        # 4th quarter
        phi = np.arccos(cosval)  # phi in the range of [0, pi]
        phi = 2 * np.pi - phi
        assert 3/2*np.pi <= phi < 2*np.pi, "Error"
    else:
        raise ValueError

    assert (phi>=0) and (phi<2*np.pi)
    return (rho, phi)


def pol2cart(rho, phi):
    """
    polar coordinates to cartisian coordinates
    @param rho: float
    @param phi: float
    @return: tuple of (x, y)
    """
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return (x, y)


def checkanticlockwise(XYVertCroSec):
    """
    check of vertices in XYVertCroSec are anticlockwise
    @param XYVertCroSec: 2d array, XYVertCroSec[i] = [x1, y1, x2, y2, x3, y3] for external boundary of triangular cross
                          section centered at (0, 0)
    @return:
    """
    for i in range(XYVertCroSec.shape[0]):
        x1, y1, x2, y2, x3, y3 = XYVertCroSec[i, :]

        if (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1) < 0:
            msg = "Line " + str(i) + " in XYVertCroSec is clockwise!"
            raise ValueError(msg)
        else:   # anti-clockwise
            pass
    return


def CartCroSec2PolarCroSec(CroSecShapeMatFile, t):
    """
    load cartesian coordinates of cross section from output of GenBoundTriOutOfPointCloud.m
    @param CroSecShapeMatFile: str, path to mat file generated by GenBoundTriOutOfPointCloud.m
    @param t: 1d array, thickness of cross section of beams, t[i] for beam i
    @return:
            PolVertCroSec: 2d array, PolVertCroSec[i] = [rho0, phi1, rho1, ohi1, rho3, phi3]
            t: 1d array, t[i] is thickness of hollow triangular cross section of beam i
    """
    extract = LoadNumericMatrixInMatFile(CroSecShapeMatFile, ['trivert'])   # load x y coordinates of vertices of beam cross sections, \
    XYVertCroSec = extract['trivert']                                       # XYVertCroSec[i]=[x1, y1, x2, y2, x3, y3] centered at (0, 0)

    checkanticlockwise(XYVertCroSec)    # check if vertices of triangles are given in anti-clockwise sequence
    nCroSec = XYVertCroSec.shape[0]     # number of cross sections

    PolVertCroSec = np.zeros((nCroSec, 6))  # change cartesian coordinates of vertices to polar coordinates

    for i in range(nCroSec):
        for j in range(3):
            x = XYVertCroSec[i, 2 * j]
            y = XYVertCroSec[i, 2 * j + 1]      # cartesian coordinates
            rho, phi = cart2pol(x, y)           # polar coordinates

            phi = phi + 2*np.pi                 # force phi to be in range of [2*pi, 4*pi]

            PolVertCroSec[i, 2 * j] = rho
            PolVertCroSec[i, 2 * j + 1] = phi

    return PolVertCroSec, t     # note! phi0 is not necessarily smaller than phi1


def PolarCroSec2CartCroSec(rho0, phi0, rho1, phi1, t):
    """
    change polar coordinates of cross section into cartesian coordinates
    @param rho0: 1d array
    @param phi0: 1d array
    @param rho1: 1d array
    @param phi1: 1d array
    @param t: 1d array
    @return:
        XYVertCroSec: 2d array, XYVertCroSec[i] = [x1, y1, x2, y2, x3, y3] for external boundary of triangular cross
                        section centered at (0, 0)
        t: 1d array, t[i] is thickness of hollow triangular cross section of beam i
    """
    nCroSec = len(rho0)
    assert len(phi0) == nCroSec and len(rho1) == nCroSec and len(phi1) == nCroSec and len(t) == nCroSec, \
        "rho0, phi0, rho1, phi1, t should be of identical length!"
    assert (rho0 > 0).all(), "All in rho0 should be positive!"
    assert (rho1 > 0).all(), "All in rho1 should be positive!"
    assert (rho0 >= 2 * np.pi) and (rho0 < 4 * np.pi)
    assert (rho1 >= 2 * np.pi) and (rho1 < 4 * np.pi)


    XYVertCroSec = np.zeros((nCroSec, 6))   # init array to store corner x y coordinates

    for i in range(nCroSec):

        RHO0 = rho0[i]
        PHI0 = phi0[i]  # rho, phi of 1st corners

        RHO1 = rho0[i]
        PHI1 = phi0[i]  # rho, phi of 2nd corners

        x0, y0 = pol2cart(RHO0, PHI0)
        x1, y1 = pol2cart(RHO1, PHI1)
        XYVertCroSec[i][:4] = [x0, y0, x1, y1]

        # compute x2, y2
        x2 = 3*0 - x0 - x1
        y2 = 3*0 - y0 - y1
        XYVertCroSec[i][4:] = [x2, y2]

    return XYVertCroSec, t


def UpdateCroSecFolders(XYVertCroSec, t, FlagOfChange=None):
    """
    update subfolders 'genTriCros' and 'tri_mesh'
    @param XYVertCroSec: 2d array, 10*6
    @param t: 1d array
    @param FlagOfChange: None or 1d array of 1/0, FlagOfChange[i] = 1 means CroSeci is changed
    @return: none
    """
    if FlagOfChange == None:
        index = range(0, XYVertCroSec.shape[0])
    else:
        index = np.nonzero(FlagOfChange.astype(bool))[0]    # 1d array of positions

    CroSecFolders = ['./CroSec' + str(i) for i in [i+1 for i in index]]  # folders of cross sections

    # generate Delaunay triangulation in a reference triangle
    vertices = np.array([(0, 1), (0, 0), (1, 0)], dtype=float)
    x = np.linspace(0, 1, 5)   # here set the number of seeds on a edge in the reference triangle
    y = x
    xv, yv = np.meshgrid(x, y, sparse=False, indexing='ij')
    xy = [[xcoor, ycoor] for (xcoor, ycoor) in zip(xv.flatten(), yv.flatten()) if (xcoor + ycoor) <= 1]
    points = np.array(xy)
    deltri = Delaunay(points)   # triangulation


    # regenerate genTriCros/shift2centroid.mat
    for (intid, CroSecFolder) in zip(index, CroSecFolders):
        Po0 = XYVertCroSec[intid][:2]                    # [x1, y1]
        Po1 = XYVertCroSec[intid][2:4]                   # [x2, y2]
        thick = t[intid]

        section = HOLTRICS([0, 0], Po0, Po1, thick)         # cross section

        Po = np.array([section.x0, section.x1, section.x2])     # cartesian coordinates of external triangle
        Pi = np.array([section.x0i, section.x1i, section.x2i])  # cartesian coordinates of internal triangle

        # create folders if necessary
        if not os.path.exists(CroSecFolder):
            os.mkdir(CroSecFolder)
            warnings.warn('create ' + CroSecFolder, RuntimeWarning)

        genTriCrosFolder = os.path.join(CroSecFolder, 'genTriCros')
        if not os.path.exists(genTriCrosFolder):
            os.mkdir(genTriCrosFolder)
            warnings.warn('create ' + genTriCrosFolder, RuntimeWarning)
        elif os.path.isfile(genTriCrosFolder):
            raise RuntimeError(genTriCrosFolder + ' is a file rather than a folder!')

        # write to genTriCros/shift2centroid.mat
        shift2centroidfile = os.path.join(genTriCrosFolder, 'shift2centroid.mat')
        WriteNumericMatrix2MatFile(shift2centroidfile, {'Po': Po,
                                                       'Pi': Pi,
                                                       'thick': thick
                                                        })

        allpoints, alltris = section.mapmesh(deltri, vertices, aspratio=True)   # map reference triangulation mesh
        bound_edge, bgp = findbound(allpoints, alltris)

        trimeshFolder = os.path.join(CroSecFolder, 'tri_mesh')
        if not os.path.exists(trimeshFolder):
            os.mkdir(trimeshFolder)
            warnings.warn('create ' + trimeshFolder, RuntimeWarning)
        elif os.path.isfile(trimeshFolder):
            raise RuntimeError(trimeshFolder + ' is a file rather than a folder!')

        coords_ien_bgpfile = os.path.join(trimeshFolder, 'coords_ien_bgp.mat')
        WriteNumericMatrix2MatFile(coords_ien_bgpfile, {'coords': allpoints,
                                                       'ien': alltris + 1,
                                                       'bound_edge': bound_edge.astype('float') + 1,
                                                       'bgp': bgp.reshape((-1, 1)).astype('float') + 1})


if __name__=="__main__":

    import numpy as np
    import sys
    import os

    debug = 0   # 1 for debug
    if debug:
        warnings.warn('Running in testing mode! Not for Martin', RuntimeWarning)
        CroSecShapeMatFile = './node_data_trivert_strut1.mat'
        t = 0.01*np.ones(10)
    else:

        if sys.argv[1] == 'SetInitial':
            try:
                CroSecShapeMatFile = sys.argv[2]        # path to mat file generated by GenBoundTriOutOfPointCloud.m
                t = [float(x) for x in sys.argv[3:]]    # thickness of 10 cross sections
                t = np.array(t)
                assert t.shape[0] == 10

                PolVertCroSec, t = CartCroSec2PolarCroSec(CroSecShapeMatFile, t)    # polar coordinates of cross sections

                XYVertCroSec = LoadNumericMatrixInMatFile(CroSecShapeMatFile, ['trivert'])['trivert']   # cartesian coordinates of cross sections

                UpdateCroSecFolders(XYVertCroSec, t, FlagOfChange=None)

                # save rho0, phi0, rho1, phi1 to mat file
                root = os.getcwd()
                rho_phiMatFile = os.path.join(root, 'rho_phi.mat')
                WriteNumericMatrix2MatFile(rho_phiMatFile, {'rho0': PolVertCroSec[:, 0],
                                                            'phi0': PolVertCroSec[:, 1],
                                                            'rho1': PolVertCroSec[:, 2],
                                                            'phi1': PolVertCroSec[:, 3],
                                                            't': t.astype('float')
                                                             })
            except:
                msg = 'Calling SetShapeOfCroSec SetInitial failed.'
                raise RuntimeError(msg)

        elif sys.argv[1] == 'Readin':
            try:
                rho_phiMatFile = sys.argv[2]    # path of rho_phi.mat file

                extract = LoadNumericMatrixInMatFile(rho_phiMatFile, ['rho0', 'phi0', 'rho1', 'phi1', 't'])
                rho0 = extract['rho0']
                phi0 = extract['phi0']
                rho1 = extract['rho1']
                phi1 = extract['phi1']
                t = extract['t']

                XYVertCroSec, t = PolarCroSec2CartCroSec(rho0, phi0, rho1, phi1, t)
                UpdateCroSecFolders(XYVertCroSec, t, FlagOfChange=None)
            except:
                msg = 'Calling SetShapeOfCroSec Readin failed.'
                raise RuntimeError(msg)









