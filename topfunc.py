import os
import numpy as np
from SetShapeOfCroSec import PolarCroSec2CartCroSec, UpdateCroSecFolders, PlotAllCroSecs
from mapmesh_aio import TooThick
import subprocess
import re
import csv


def ReturnDateTimeStr():
    from datetime import datetime
    now = datetime.now()
    current_time = now.strftime('%Y-%m-%d-%H-%M-%S')
    return current_time


class VertexReverted(Exception):
    """
    exception class when phi0 is larger than phi1
    """
    pass


def counted(f):
    """
    define descriptor
    @param f:
    @return:
    """
    def wrapped(*args, **kwargs):
        wrapped.calls += 1
        return f(*args, **kwargs)
    wrapped.calls = 0
    return wrapped

@counted
def TopFunc(GeomCfg, NumOfPadZeros=4, penalty=1e4):
    """
    top function to bridge with dakota
    @param GeomCfg: dict, {rho0: 1d array, phi0: 1d array, rho1: 1d array, phi1: 1d array}
    @return:
        value of objective function
    """

    workdir = os.getcwd()

    # step 1: mesh cross sections
    try:
        rho0 = GeomCfg['rho0']
        phi0 = GeomCfg['phi0']
        rho1 = GeomCfg['rho1']
        phi1 = GeomCfg['phi1']
        t = GeomCfg['t']

        if not (phi1 > phi0).all():
            raise VertexReverted()

        PolVertCroSec = np.zeros((rho0.size, 6))  # init
        PolVertCroSec[:, 0] = rho0
        PolVertCroSec[:, 1] = phi0
        PolVertCroSec[:, 2] = rho1
        PolVertCroSec[:, 3] = phi1

        XYVertCroSec = PolarCroSec2CartCroSec(PolVertCroSec)  # rho3, phi3 = 0, 0
        UpdateCroSecFolders(XYVertCroSec, t, path=workdir)
        PlotAllCroSecs(workdir)  # plotting mesh of cross sections

    except TooThick as exc:
        # thickness is too large, predefined objective function value is returned
        ObjectiveFuncValue = penalty*exc.gap**2
        WriteObjFunValToFile(TopFunc.calls, ObjectiveFuncValue, NumOfPadZeros)  # generate log of objective function
        return ObjectiveFuncValue
    except VertexReverted:
        # phi0 > phi1 for at least one cross section, predefined objective function value is returned
        ObjectiveFuncValue = penalty*max(phi0 - phi1)**2
        WriteObjFunValToFile(TopFunc.calls, ObjectiveFuncValue, NumOfPadZeros)  # generate log of objective function
        return ObjectiveFuncValue
    else:   # run if no exception raised

        MatlabPath = '/opt/local/MATLAB/R2018a/bin/matlab'  # path of matlab bin

        # --------------------------------------------------
        # evaluate warping function of cross sections
        logfile = 'EvalWarp' + str(TopFunc.calls).zfill(NumOfPadZeros) + '.txt'
        cmd = MatlabPath + ' -nodesktop -nosplash -r "EvalWarp; quit" >  ' + logfile

        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
        (output, err) = p.communicate()
        p_status = p.wait()                     # This makes the wait possible
        # print("Command output: " + output)    # This will give you the output of the command being executed
        if not err == None:
            raise RuntimeError('Error occurred when evaluating warping function of beam cross sections!')
        # --------------------------------------------------

        # --------------------------------------------------
        # run beam simulation
        logfile = 'mtfem' + str(TopFunc.calls).zfill(NumOfPadZeros) + '.txt'
        cmd = MatlabPath + ' -nodesktop -nosplash -r "datf = \'tension\'; mtfem(datf); quit" >  ' + logfile

        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
        (output, err) = p.communicate()
        p_status = p.wait()  # This makes the wait possible
        if not err == None:
            raise RuntimeError('Error occurred when running mtfem simulation!')
        # --------------------------------------------------

        # evaluate objective function
        ObjectiveFuncFolder = os.path.join(workdir)
        cmd = MatlabPath + ' -nodesktop -nosplash -r "type = \'tension\'; strut = \'strut1\'; out = objectivfunction(type, strut); format long; display(out); quit"'
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
        (output, err) = p.communicate()
        p_status = p.wait()  # This makes the wait possible
        if not err == None:
            raise RuntimeError('Error occurred when evaluating objective function!')
        CommandWindowOutput = output.decode('utf-8')
        CommandWindowOutput = CommandWindowOutput.replace('\n', ' ')
        CommandWindowOutput = re.sub(r'\s+', ' ', CommandWindowOutput)
        ObjectiveFuncValueStr = CommandWindowOutput.split(' ')[-2]
        try:
            ObjectiveFuncValue = float(ObjectiveFuncValueStr)

            # ---------------------------------------------
            # generate csv file to store configuration
            logfile = 'geomcfg' + str(TopFunc.calls) + '.csv'
            f = open(logfile, 'w')
            writer = csv.writer(f)
            row = ['phi0', 'rho0', 'phi1', 'rho1', 't']
            writer.writerow(row)                # write header
            for i in range(rho0.size):
                row = [phi0[i], rho0[i], phi1[i], rho1[i], t[i]]
                writer.writerow(row)
            f.close()
            # ---------------------------------------------
            WriteObjFunValToFile(TopFunc.calls, ObjectiveFuncValue, NumOfPadZeros)     # generate log of objective function

            # change name of tensionFU.mat
            cmd = 'mv tensionFU.mat ' + 'tensionFU' + str(TopFunc.calls).zfill(NumOfPadZeros) + '.mat'
            os.popen(cmd)

            return ObjectiveFuncValue
        except:
            raise RuntimeError('ObjectiveFuncValueStr is :' + ObjectiveFuncValueStr)


def WriteObjFunValToFile(counetr, ObjFunVal, NumOfPadZeros):
    logfile = 'objfunval' + str(counetr).zfill(NumOfPadZeros) + '.txt'
    f = open(logfile, 'w')
    f.write(str(ObjFunVal))
    f.close()


def DebugCase1():
    """
    test TooThick exception
    @return:
    """
    GeomCfg = {'phi0': np.array([ 8.43640868, 10.45576289, 10.56926616,  8.30835229,  7.98730527,
                                 10.42666906,  6.29      ,  6.34133918, 10.59525173, 12.1004305 ]),
               'rho0': np.array([0.43018227, 0.29190581, 0.26864673, 0.43729482, 0.44214407,
                                0.51858664, 0.50382642, 0.68494915, 0.57548486, 0.63896565]),
               'phi1': np.array([10.67723593,  6.49803282,  6.50697752, 10.60102591, 11.04580648,
                                 9.69748342,  8.55810619,  8.75830486, 12.24382294,  8.83105787]),
               'rho1': np.array([0.23042822, 0.37032153, 0.34655473, 0.35086148, 0.43960617,
                                0.46288121, 0.69101771, 0.43774986, 0.62543691, 0.6804659 ]),
                't': np.array([0.03004251, 0.02978125, 0.02991518, 0.02978733, 0.02799266,
                                0.02837298, 0.02695488, 0.03035617, 0.02879683, 0.02457201]),
               }

    ObjectiveFuncValue = TopFunc(GeomCfg)
    print('Objective function value of DebugCase1 is %10.5f'%(ObjectiveFuncValue))
    print('TopFunc.calls = ' + str(TopFunc.calls))


def DebugCase2():
    """
    test a normal analysis
    @return:
    """
    GeomCfg = {'phi0': 30/180*np.pi*np.ones(10),
               'rho0': 0.5*np.ones(10),
               'phi1': 150/180*np.pi*np.ones(10),
               'rho1': 0.5*np.ones(10),
                't': 0.05*np.ones(10),
               }                        # a hypothetical cross section shape

    ObjectiveFuncValue = TopFunc(GeomCfg)
    print('Objective function value of DebugCase2 is is %10.5f' % (ObjectiveFuncValue))
    print('TopFunc.calls = ' + str(TopFunc.calls))


def DebugCase3():
    """
    test a normal analysis
    @return:
    """
    GeomCfg = {'phi0': 30/180*np.pi*np.ones(10),
               'rho0': 0.6*np.ones(10),
               'phi1': 150/180*np.pi*np.ones(10),
               'rho1': 0.5*np.ones(10),
                't': 0.05*np.ones(10),
               }                            # a hypothetical cross section shape

    ObjectiveFuncValue = TopFunc(GeomCfg)
    print('Objective function value of DebugCase3 is is %10.5f' % (ObjectiveFuncValue))
    print('TopFunc.calls = ' + str(TopFunc.calls))



if __name__=='__main__':
    # debug
    DebugCase1()
    # DebugCase2()
    # DebugCase3()
