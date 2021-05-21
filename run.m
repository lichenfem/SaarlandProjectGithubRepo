%call function to mesh und simulate

% created by Martin

function f(params,results)

FILE_EXCHANGE = '/home/martin/Documents/sent2martin20210519/Program';
addpath(genpath(FILE_EXCHANGE));

fid = fopen(params,'r');
C = textscan(fid,'%n%s');
fclose(fid);

num_vars = C{1}(1);


for i=1:num_vars
    geo(i)=C{1}(i+1);
end



for j=1:num_vars/5
    phi0(j)=geo(j);
    phi1(j)=geo(num_vars/5+j);
    rho0(j)=geo(2*num_vars/5+j);
    rho1(j)=geo(3*num_vars/5+j);
    t(j)=geo(4*num_vars/5+j);
end

save('rho_phi.mat', 'rho0', 'phi0', 'rho1', 'phi1', 't');


%tic
%MeshCroSec_function(geo);
%toc

%S1_GenBoundTriOutOfPointCloud;
%system('python3 S2_SetShapeOfCroSec.py SetInitial ./node_data_trivert_strut1.mat 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02');

system('python3 SetShapeOfCroSec.py Readin ./rho_phi.mat');
S3_EvalWarp;

tic
%tensile test
datf='tension';
mtfem(datf);
tens=objectivfunction_new('tension', 'strut1');


%bending
%datf='bending0';
%mtfem(datf);
%[bending1]=objectivfunction2('bending0', 'strut1');
%datf='bending1';
%mtfem(datf);
%[bending2]=objectivfunction2('bending1', 'strut1');
%datf='bending2';
%mtfem(datf);
%[bending3]=objectivfunction2('bending2', 'strut1');

%sum for bending
%bend=bending1+bending2+bending3;

%torsion
%datf='torsion'
%mtfem(datf)
%tors=objectivfunction2('torsion', 'strut1');

toc

%total sum

[f]=tens;%+bend+tors;



%------------------------------------------------------------------
% WRITE results.out
%------------------------------------------------------------------
fid = fopen(results,'w');
fprintf(fid,'%20.10e f\r\n', f);


fclose(fid);

