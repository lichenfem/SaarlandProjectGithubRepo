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
    Po = extract['Po']
    Pi = extract['Pi']


    # load mesh
    coords_ien_bgpFile = os.path.join(CroSecFolder, 'tri_mesh', 'coords_ien_bgp.mat')
    extract = LoadNumericMatrixInMatFile(coords_ien_bgpFile, ['ien', 'coords'])
    ien = extract['ien']
    coords = extract['coords']


    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.set_aspect(aspect="equal")

    # outer profile
    x = Po[:, 0].tolist(); x.append(x[0])
    y = Po[:, 1].tolist(); y.append(y[0])
    ax.plot(x, y, 'go--', linewidth=1, markersize=12, label='outer profile')
    for i in range(3):
        ax.text(x[i], y[i], str(i), fontsize=12)    # add node index

    # inner profile
    x = Pi[:, 0].tolist(); x.append(x[0])
    y = Pi[:, 1].tolist(); y.append(y[0])
    ax.plot(x, y, 'ro--', linewidth=1, markersize=12, label='inner profile')

    ax.triplot(coords[:, 0], coords[:, 1], ien - 1, linewidth=0.1)

    ax.plot(0, 0, 'rs', markersize=12)          # plot Origin of cross section

    now = datetime.now()
    current_time = now.strftime("%H-%M-%S")
    pngFile = CroSecFolder.split('/')[-1] + '_' + current_time + '.png'
    pngFile = os.path.join(CroSecFolder, pngFile)
    fig.savefig(pngFile,dpi=300, transparent=True)      # save the figure

    plt.close("all")
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
