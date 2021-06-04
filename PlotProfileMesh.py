import os
import matplotlib.pyplot as plt
import numpy as np
from LoadDatainMatFile import LoadNumericMatrixInMatFile
from datetime import datetime

def PlotSaveHollowTri(CroSecFolder):
    """
    plot cross section in a CroSecXX folder
    @param CroSecFolder: str, path to a CroSecXX folder
    @return:
    """
    if CroSecFolder.endswith('/'):
        CroSecFolder = CroSecFolder[:-1]    # remove trailing '/'

    # load profile
    shift2centroidFile = os.path.join(CroSecFolder, 'genTriCros', 'shift2centroid.mat')
    if not os.path.isfile(shift2centroidFile):
        raise FileNotFoundError(shift2centroidFile + 'not exist!')

    extract = LoadNumericMatrixInMatFile(shift2centroidFile, ['Po', 'Pi'])
    Po = extract['Po']          # vertices on external bound
    Pi = extract['Pi']          # vertices on internal bound


    # load mesh
    coords_ien_bgpFile = os.path.join(CroSecFolder, 'tri_mesh', 'coords_ien_bgp.mat')
    extract = LoadNumericMatrixInMatFile(coords_ien_bgpFile, ['ien', 'coords'])
    ien = extract['ien']        # index starts from 1
    coords = extract['coords']  # x y coordinates of mesn nodes

    # start plotting
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.set_aspect(aspect="equal")

    # plot outer bound
    x = Po[:, 0].tolist()
    x.append(x[0])
    y = Po[:, 1].tolist()
    y.append(y[0])
    ax.plot(x, y, 'go--', linewidth=1, label='outer profile')
    for i in range(3):
        ax.text(x[i], y[i], str(i), fontsize=12)    # annotate node index

    # plot inner bound
    x = Pi[:, 0].tolist()
    x.append(x[0])
    y = Pi[:, 1].tolist()
    y.append(y[0])
    ax.plot(x, y, 'ro--', linewidth=1, label='inner profile')

    # plot triangular mesh
    ax.triplot(coords[:, 0], coords[:, 1], ien - 1, linewidth=0.1)

    ax.plot(0, 0, 'rs', markersize=6)          # plot Origin of cross section

    now = datetime.now()
    current_time = now.strftime("%Y-%m-%d-%H-%M-%S")    # year-month-day-hour-min-sec
    pngFile = CroSecFolder.split('/')[-1] + '_' + current_time + '.png' # name of png file
    pngFile = os.path.join(CroSecFolder, pngFile)       # full path
    fig.savefig(pngFile, dpi=300, transparent=True)

    plt.close(fig)
    return


if __name__=='__main__':

    import sys
    import os
    try:
        workdir = sys.argv[1]   # path of parent folder of all CroSecXX
    except:
        workdir = os.getcwd()

    if not os.path.isdir(workdir):
        raise ValueError(workdir + ' is not a valid folders')
    if workdir.endswith('/'):
        workdir = workdir[:-1]

    CroSecFolders = []
    for (root, dirs, files) in os.walk(workdir):
        if root == workdir:
            for dir in dirs:
                if dir.startswith('CroSec'):
                    CroSecFolders.append(dir)

    if len(CroSecFolders) == 0:
        raise RuntimeError('No CroSec folders in ' + workdir)

    for CroSecFolder in CroSecFolders:
        print('Plotting ' + CroSecFolder + '...')
        PlotSaveHollowTri(CroSecFolder)
