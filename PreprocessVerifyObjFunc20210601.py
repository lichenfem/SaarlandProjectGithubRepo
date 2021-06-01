import numpy as np
from LoadDatainMatFile import  WriteNumericMatrix2MatFile
import os

# load cross section shape from file
workdir = os.getcwd()
file = 'DimensionOfOptimalDesign.txt'
path = os.path.join(workdir, file)
assert os.path.isfile(path)

Matrix = []
with open(path, 'r') as f:
    lines = f.readlines()
    for line in lines:
        content = line.split()
        Matrix.append(float(content[0]))

Matrix = np.array(Matrix)

phi0 = Matrix[0:10]
phi1 = Matrix[10:20]
rho0 = Matrix[20:30]
rho1 = Matrix[30:40]
t = Matrix[40:50]

# write to file
outputfile = 'rho_phi.mat'
WriteNumericMatrix2MatFile(os.path.join(workdir, outputfile),
                           {'phi0': phi0,
                            'phi1': phi1,
                            'rho0': rho0,
                            'rho1': rho1,
                            't': t})
print('writing to ' + os.path.join(workdir, outputfile), ' COMPLETED!')
