

mapmesh_aio,py
	

LoadDatainMatFile.py
	used to read numeric arrays and floats in .mat file into python

mapmesh_aio.py
	usage: python3 mapmesh_aio.py $WorkDir
	usage: python3 mapmesh_aio.py
	this program meshed $WorkDir/CroSec1-10

GenBoundTriOutOfPointCloud.m
	matlab script to generate minimum area bounding triangles for each beam cross section

EvalWarp.m:
	matlab script tp evaluate warping fuction distribution


SetShapeOfCroSec.py:
	using one of the following format depending on your need:

	format 1:
    		python3 SetShapeOfCroSec.py SetInitial $node_data_trivert_strut1 $t1 $t2 $t3 ... $t10

    		here $node_data_trivert_strut1 is path to 'node_data_trivert_strut1.mat' file in the previous step.
    		$t1 ... $t10 is the thickness of each beam element.

    		after running the above command, a 'rho_phi.mat' file is generated in the folder, which contains rhos and phis of
    		vertices 1 and 2 of each beam cross section, as well as thickness of each beam cross section. Input in this file can
    		be used as input for your DAKOTA program.

    		Also, whenever DAKOTA sets new rhos and phis and thickness, you need to write a glue program so that they are
    		wrapped into 'rho_phi.mat' (used as input for the calling format 2 below)

	format 2:
    		python3 SetShapeOfCroSec.py Readin $rho_phi_mat_file

    		here $rho_phi_mat_file is path of 'rho_phi.mat' file. The rhos, phis and thickness in the file are read so that the cross
    		sections are meshed (for the input of evaluating warping function distribution of the cross section).


PlotProfileMesh.py
	plot triangular mesh of cross sections.