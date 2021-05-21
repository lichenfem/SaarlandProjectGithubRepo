import h5py
import os
import numpy as np
import scipy.io


class CROSEC:
    """
    class of an irregular cross section
    """
    def __init__(self, xc):
        """

        :param xc: array([x,y,z])
        """
        # global coordinate of centroid
        self.xc = np.array(xc)
        self.outerprofile = []      # coordinates of points on outer profile
        self.innerprofile = []      # coordinates of points on inner profile



# --------------------------------------
if __name__=='__main__':

    file = 'res_tot_strut1_Li.mat'
    matcontents = scipy.io.loadmat(file)

    # import points on the axis
    axisnodes = matcontents['center_loc_center']

    centroidline = []
    for i in range(axisnodes.shape[0]):
        centroidline.append(CROSEC(axisnodes[i, :]))

