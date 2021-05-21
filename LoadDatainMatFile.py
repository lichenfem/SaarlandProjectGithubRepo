import os
import sys
import numpy as np

def LoadNumericMatrixInMatFile(MatFile, variables):
    """
    load numerical arrays in Mat file into python
    @param MatFile: str, path of mat file
    @param variables: list of str, names of variables to be extracted from mat file
    @return: dict
    """
    if not os.path.isfile(MatFile):
        raise FileNotFoundError(MatFile + ' not exist!')

    extract = {}

    try:
        import scipy.io
        mat_contents = scipy.io.loadmat(MatFile)
        # print keys in vars
        for variable in variables:
            if variable not in mat_contents.keys():
                raise ValueError(variable + ' not in ' + MatFile)
            else:
                value = mat_contents[variable]      # array

                # for vectors
                if value.shape[0] == 1 or value.shape[1] == 1:
                    value = value.flatten()

                # for single float
                if value.shape == (1,):
                    value = value[0]
                elif value.shape == (1, 1):
                    value = value[0][0]
                else:
                    pass

            extract[variable] = value
    except:
        import h5py
        with h5py.File(MatFile, 'r') as f:
            for variable in variables:
                if not variable in f.keys():
                    raise ValueError(variable + ' not in ' + MatFile)
                else:
                    value = np.array(f[variable]).transpose()   # array

                    if value.shape[0] == 1 or value.shape[1] == 1:
                        value = value.flatten()

                    # for single float
                    if value.shape == (1,):
                        value = value[0]
                    elif value.shape == (1, 1):
                        value = value[0][0]
                    else:
                        pass
            extract[variable] = value

    return extract

def WriteNumericMatrix2MatFile(MatFile, dict):
    """
    @param MatFile: str, full path to mat file
    @param dict: numerical arrays to be stored
    """
    import scipy.io
    assert MatFile.split('.')[-1] == 'mat', 'File name ends should end with .mat'

    for val in dict.values():
        if isinstance(val, np.ndarray) or isinstance(val, float) or isinstance(val, int):
            continue
        else:
            msg = 'Type ' + str(type(val)) + ' not permitted to write to mat file.'
            raise ValueError(msg)
    scipy.io.savemat(MatFile, dict)
    return

if __name__=="__main":
    # test loading

    MatFile = '/home/lichen/Documents/MATLAB/Martin_Saarland_project/CrossSectionPolarParameterization/PerturbStudyByLiOfXCoordinateOfVertex1OfBeamAtFixedEnd/node_data_trivert_strut1.mat'
    variables = ['data', 'node', 'trivert']

    extract = LoadNumericMatrixInMatFile(MatFile, variables)