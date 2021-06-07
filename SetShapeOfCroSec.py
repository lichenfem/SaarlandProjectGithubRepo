import os
import sys
import numpy as np
from mapmesh_aio import HOLTRICS, findbound
from MiscPyUtilities.pipeMAT import LoadMatFile2NumpyArray, WriteNumpyArray2MatFile
import warnings
from PlotProfileMesh import PlotSaveHollowTri
from MiscPyUtilities.triangleGeometry import checkAnticlockwise



def cart2pol(x, y, phase=0):
    """
    cartisian coordinates to polar coordinates
    @param x: float
    @param y: float
    #@param phase: float, phase angle to shift phi
    @return: t
        tuple of (rho, phi), phi in [0, 2*pi)+phase, angle is in rad
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

    # adjust phi by phase
    phi = phi + phase

    return (rho, phi)


def pol2cart(rho, phi):
    """
    polar coordinates to cartisian coordinates
    @param rho: float, >0
    @param phi: float, in rad
    @return: tuple of (x, y)
    """
    assert rho > 0
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return (x, y)


def CartCroSec2PolarCroSec(XYVertCroSec, phase=2*np.pi):
    """
    transform cartesian coordinates of cross section vertices to polar coordinates
    @param XYVertCroSec: 2d array, XYVertCroSec[i] = [x1, y1, x2, y2, x3, y3] for external boundary of triangular cross
                          section centered at (0, 0)
    @return:
            PolVertCroSec: 2d array, PolVertCroSec[i] = [rho0, phi1, rho1, ohi1, rho3, phi3], phi is in [2*pi, 4*pi]
    """
    seq = checkAnticlockwise(XYVertCroSec)  # check if vertices of triangles are given in anti-clockwise sequence
    if not seq.all():
        raise RuntimeError()


    nCroSec = XYVertCroSec.shape[0]     # number of cross sections

    PolVertCroSec = np.zeros((nCroSec, 6))  # init array to store polar coordinates for vertices of cross sections

    for i in range(nCroSec):
        for j in range(3):
            x = XYVertCroSec[i, 2 * j]
            y = XYVertCroSec[i, 2 * j + 1]      # cartesian coordinates

            rho, phi = cart2pol(x, y, phase=phase)           # polar coordinates

            PolVertCroSec[i, 2 * j] = rho
            PolVertCroSec[i, 2 * j + 1] = phi   # [0, 2*pi]

    return PolVertCroSec    # note! phi0 is not necessarily smaller than phi1


def PolarCroSec2CartCroSec(PolVertCroSec):
    """
    change polar coordinates of vertices of cross sections into cartesian coordinates
    @param PolVertCroSec: 2d array, PolVertCroSec[i] = [rho0, phi1, rho1, ohi1, rho3/0, phi3/0]
    @return:
        XYVertCroSec: 2d array, XYVertCroSec[i] = [x1, y1, x2, y2, x3, y3] for external boundary of triangular cross
                        section centered at (0, 0)
    """
    nCroSec = PolVertCroSec.shape[0]

    XYVertCroSec = np.zeros((nCroSec, 6))   # init array to store x y coordinates of cross section vertices

    for i in range(nCroSec):

        rho1, phi1, rho2, phi2, rho3, phi3 = PolVertCroSec[i, :].tolist()

        x1, y1 = pol2cart(rho1, phi1)
        x2, y2 = pol2cart(rho2, phi2)

        # compute x3, y3 based on x1, x2, y1, y2
        x3 = 3*0 - x1 - x2
        y3 = 3*0 - y1 - y2

        if rho3 > 0:
            # compute x3, y3 based on rho3, phi3
            x3t, y3t = pol2cart(rho3, phi3)
            if np.allclose(np.array([x3, y3]), np.array([x3t, y3t]), atol=1e-5):
                raise RuntimeError()

        XYVertCroSec[i][:] = [x1, y1, x2, y2, x3, y3]
    return XYVertCroSec


def UpdateCroSecFolders(XYVertCroSec, t, path='.', NUMELEMPEREDGE=5):
    """
    update subfolders 'genTriCros' and 'tri_mesh'
    @param XYVertCroSec: 2d array, 10*6, [x1, y1, x2, y2, x3, y3] of vertices of a triangular cross section
    @param t: 1d array, thickness of cross sections
    @param path: str, full path to the parent folder that include all CroSecXX subfolders
    @param NUMELEMPEREDGE: int, number of element on the Isosceles right triangle
    @return: none
    """
    assert os.path.isdir(path), path + " is not a valid folder!"

    nCroSecs = XYVertCroSec.shape[0]    # number of cross sections
    CroSecFolders = [os.path.join(path, 'CroSec' + str(i)) for i in range(1, nCroSecs+1)]   # folders of cross sections

    # generate Delaunay triangulation in a reference triangle
    from scipy.spatial import Delaunay
    vertices = np.array([(0, 1), (0, 0), (1, 0)], dtype=float)
    x = np.linspace(0, 1, NUMELEMPEREDGE)
    y = x
    xv, yv = np.meshgrid(x, y, sparse=False, indexing='ij')
    xy = [[xcoor, ycoor] for (xcoor, ycoor) in zip(xv.flatten(), yv.flatten()) if (xcoor + ycoor) <= 1]
    points = np.array(xy)
    deltri = Delaunay(points)           # triangulation in reference Isosceles right triangle


    # regenerate genTriCros/shift2centroid.mat
    for (intid, CroSecFolder) in zip(range(0, nCroSecs), CroSecFolders):

        # cartesian coordinates of two vertices of a cross section
        Po0 = XYVertCroSec[intid][:2].tolist()                    # [x1, y1]
        Po1 = XYVertCroSec[intid][2:4].tolist()                   # [x2, y2]
        thick = t[intid]

        section = HOLTRICS([0, 0], Po0, Po1, thick)             # generate cross section object

        # 2d array of cartesian coordinates of external triangle vertices
        Po = np.array([section.x0, section.x1, section.x2])
        # 2d array of cartesian coordinates of internal triangle vertices, x0i is next to x0
        Pi = np.array([section.x0i, section.x1i, section.x2i])

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
        WriteNumpyArray2MatFile(shift2centroidfile, {'Po': Po,       # 2d array
                                                        'Pi': Pi,       # 2d array
                                                        'thick': thick  # float
                                                        })

        # mesh the cross section
        allpoints, alltris = section.mapmesh(deltri, vertices, aspratio=True)
        # generate bounding edge information for cross section mesh
        bound_edge, bgp = findbound(allpoints, alltris)


        trimeshFolder = os.path.join(CroSecFolder, 'tri_mesh')
        if not os.path.exists(trimeshFolder):
            os.mkdir(trimeshFolder)
            warnings.warn('create ' + trimeshFolder, RuntimeWarning)
        elif os.path.isfile(trimeshFolder):
            raise RuntimeError(trimeshFolder + ' is a file rather than a folder!')

        # write to tri_mesh/coords_ien_bgp.mat
        coords_ien_bgpfile = os.path.join(trimeshFolder, 'coords_ien_bgp.mat')
        WriteNumpyArray2MatFile(coords_ien_bgpfile, {'coords': allpoints,
                                                        'ien': alltris + 1,     # index in matlab starts from 1
                                                        'bound_edge': bound_edge.astype('float') + 1,
                                                        'bgp': bgp.reshape((-1, 1)).astype('float') + 1})


def PlotAllCroSecs(workdir):
    """
    plot mesh of CroSec when a cross section is generated
    @param workdir: str, parent folder to all CroSecXX folders
    @return:
    """
    msg = 'Cross sections plotting activated!\nDeactivate for better performance!'
    warnings.warn(msg, RuntimeWarning)

    # find folders of cross sections
    CroSecFolders = []
    for (root, dirs, files) in os.walk(workdir):
        if root == workdir:
            for dir in dirs:
                if dir.startswith('CroSec') and dir[-1].isdigit():
                    CroSecFolders.append(os.path.join(workdir, dir))

    if len(CroSecFolders) == 0:
        raise RuntimeError('No CroSec folders detected in ' + workdir)

    for CroSecFolder in CroSecFolders:
        print('Plotting ' + CroSecFolder + '...')
        PlotSaveHollowTri(CroSecFolder)


if __name__=="__main__":

    import numpy as np
    import sys
    import os

    workdir = os.getcwd()

    debug = 1   # 1 for debug
    if debug:
        warnings.warn('Running in testing mode! Not for Martin', RuntimeWarning)

        # -------------------------------------------------
        # test input of Martin that generates exception of too thick cross section
        file = os.path.join(os.getcwd(), 'design_para_of_excessively_thickness_20210607.txt')
        import pandas as pd
        for line in open(file):
            print(repr(line))
        geomconfig = pd.read_csv(file, header=None, delimiter=r"\s+")   # geometrical configuration of cross sections

        phi0 = np.array(geomconfig[0][:10])
        phi1 = np.array(geomconfig[0][10:20])
        rho0 = np.array(geomconfig[0][20:30])
        rho1 = np.array(geomconfig[0][30:40])
        t = np.array(geomconfig[0][40:50])

        PolVertCroSec = np.zeros((rho0.size, 6))  # init
        PolVertCroSec[:, 0] = rho0
        PolVertCroSec[:, 1] = phi0
        PolVertCroSec[:, 2] = rho1
        PolVertCroSec[:, 3] = phi1

        XYVertCroSec = PolarCroSec2CartCroSec(PolVertCroSec)  # rho3, phi3 = 0, 0
        UpdateCroSecFolders(XYVertCroSec, t, path=workdir)

        PlotAllCroSecs(workdir)  # plotting mesh of cross sections

        # cross section 5 prompts error
        xyvert = np.array([-0.05877377, 0.43822029, 0.0220731, -0.43905166, 0.03670068, 0.00083137])
        thick = t[4]

        from MiscPyUtilities.triangleGeometry import TriInscribedCircleRadius
        tmax = TriInscribedCircleRadius(xyvert)

        # -------------------------------------------------

        # # ---------------------------------------------------------------------------------
        # warnings.warn('Debugging SetShapeOfCroSec Readin ...', RuntimeWarning)
        # rho_phiMatFile = os.path.join(workdir, 'rho_phi.mat')   # path of rho_phi.mat file
        # content = LoadMatFile2NumpyArray(rho_phiMatFile, ['rho0', 'phi0', 'rho1', 'phi1', 't'])
        # rho1 = content['rho0']
        # phi1 = content['phi0']
        # rho2 = content['rho1']
        # phi2 = content['phi1']
        # t = content['t']
        #
        # assert (rho1 > 0).all() and (rho2 > 0).all()
        #
        # PolVertCroSec = np.array([rho1, phi1, rho2, phi2, np.zeros(rho1.size), np.zeros(rho1.size)])# rho3/phi3 set to 0
        # PolVertCroSec = PolVertCroSec.transpose()
        #
        # XYVertCroSec = PolarCroSec2CartCroSec(PolVertCroSec)
        # UpdateCroSecFolders(XYVertCroSec, t)
        # PlotAllCroSecs(workdir)
        # # ---------------------------------------------------------------------------------

    else:
        # non-debug
        if sys.argv[1] == 'SetInitial':
            # read in 'node_data_trivert_strut1.mat'
            CroSecShapeMatFile = sys.argv[2]        # path to mat file generated by GenBoundTriOutOfPointCloud.m, containing cartesian coordinates of cross section vertices
            if not '/' in CroSecShapeMatFile:
                CroSecShapeMatFile = os.path.join(workdir, CroSecShapeMatFile)


            t = [float(x) for x in sys.argv[3:]]    # thickness of 10 cross sections
            t = np.array(t)

            assert t.shape[0] == 10, str(t.shape[0]) + " cross section thickness are given, 10 expected!"

            content = LoadMatFile2NumpyArray(CroSecShapeMatFile, ['trivert', 'node', 'data'])
            XYVertCroSec = content['trivert']       # cartesian coordinates of cross section vertices

            PolVertCroSec = CartCroSec2PolarCroSec(XYVertCroSec, phase=2*np.pi)

            UpdateCroSecFolders(XYVertCroSec, t, path=workdir)

            PlotAllCroSecs(workdir)     # plotting updated mesh of cross sections, name using time flag

            # save rho0, phi0, rho1, phi1 to mat file
            rho_phiMatFile = os.path.join(workdir, 'rho_phi.mat')
            WriteNumpyArray2MatFile(rho_phiMatFile, {'rho0': PolVertCroSec[:, 0],
                                                        'phi0': PolVertCroSec[:, 1],
                                                        'rho1': PolVertCroSec[:, 2],
                                                        'phi1': PolVertCroSec[:, 3],
                                                        't': t
                                                         })
        elif sys.argv[1] == 'Readin':
            # read in 'rho_phi.mat'
            rho_phiMatFile = sys.argv[2]    # path of rho_phi.mat file
            if not '/' in rho_phiMatFile:
                rho_phiMatFile = os.path.join(workdir, rho_phiMatFile)

            content = LoadMatFile2NumpyArray(rho_phiMatFile, ['rho0', 'phi0', 'rho1', 'phi1', 't'])
            rho0 = content['rho0']
            phi0 = content['phi0']
            rho1 = content['rho1']
            phi1 = content['phi1']
            t = content['t']

            PolVertCroSec = np.zeros((rho0.size, 6))    # init
            PolVertCroSec[:, 0] = rho0
            PolVertCroSec[:, 1] = phi0
            PolVertCroSec[:, 2] = rho1
            PolVertCroSec[:, 3] = phi1

            XYVertCroSec = PolarCroSec2CartCroSec(PolVertCroSec)    # rho3, phi3 = 0, 0
            UpdateCroSecFolders(XYVertCroSec, t, path=workdir)

            PlotAllCroSecs(workdir)     # plotting mesh of cross sections

        else:
            raise RuntimeError('Invalid parameter!')










