%%objectiv function

function out = objectivfunction(type, strut)
% input:
% type: str, e.g. 'tension'
% strut: str, e.g. 'strut1'


% 10 nodes (id 2-11), 6 dofs per node, each row is for 1 increment
% last 6 columns stores reaction at free end
BEAM=zeros(11,66);  % beam element result
FE=zeros(11,66);    % volume element result

%-------------------------------------------------------------------------%
%load beam data
load([type,'011'], 'reacti');       % load reacti from 'type011.mat'
reacti=full(reacti);    % reac = [2,1;2,2;2,3;2,4;2,5;2,6;....;11,1;11,2;11,3;11,4;11,5;11,6]

%BEAM(:,1..6)= beam1 disp_xyz, rot_xyz ..... BEAM(:,56..60)=beam10
%disp_xyz, rot_xyz

for beam=1:60
    BEAM(:,beam)=reacti(:,beam);
end

% reaction at free end
BEAM(:,61)=reacti(:,115);
BEAM(:,62)=reacti(:,116);
BEAM(:,63)=reacti(:,117);
BEAM(:,64)=reacti(:,118);
BEAM(:,65)=reacti(:,119);
BEAM(:,66)=reacti(:,120);

%-------------------------------------------------------------------------%

%load FE data

load(['../Lars_postprocess/results_postprocess_by_lars/', strut, '_by_lars/res_tot_', strut, '_', type])
% load('res_tot_strut1_tension.mat'); % each row in dis_and_rot is: X, Y, Z, U, V , W, theta, ux, uy, uz

for fe=2:11
    FE(:,1+(fe-2)*6)=dis_and_rot(1:10:101,4+(fe-1)*10); % U
    FE(:,2+(fe-2)*6)=dis_and_rot(1:10:101,5+(fe-1)*10); % V
    FE(:,3+(fe-2)*6)=dis_and_rot(1:10:101,6+(fe-1)*10); % W
    FE(:,4+(fe-2)*6)=dis_and_rot(1:10:101,7+(fe-1)*10).*dis_and_rot(1:10:101,8+(fe-1)*10);  % Rx
    FE(:,5+(fe-2)*6)=dis_and_rot(1:10:101,7+(fe-1)*10).*dis_and_rot(1:10:101,9+(fe-1)*10);  % Ry
    FE(:,6+(fe-2)*6)=dis_and_rot(1:10:101,7+(fe-1)*10).*dis_and_rot(1:10:101,10+(fe-1)*10); % Rz
end

FE(isnan(FE))=0;   % no needed

% allocate reaction at free end
FE(:,61)=-Forcex(1:10:101);
FE(:,62)=-Forcey(1:10:101);
FE(:,63)=-Forcez(1:10:101);
FE(:,64)=-Momentx(1:10:101);
FE(:,65)=-Momenty(1:10:101);
FE(:,66)=-Momentz(1:10:101);

%-------------------------------------------------------------------------%

%calculate output

w1=1;       % w1: weight of displacement error and rotation error
w2=1;       % w2: weight of force error and moment error
f=0;        % f: value of objective function


% only take into account the component that is relevant to the loading case
if strcmp(type, 'tension')
    idof  = 2:6:60;     % column id for Uy of free end in BEAM/FE
    iforc = 62;         % column id for Fy of free end in BEAM/FE
elseif strcmp(type, 'torsion')
    error('not implemented');
elseif strcmp(type, 'bending0')
    error('not implemented');
elseif strcmp(type, 'bending1')
    error('not implemented');
elseif strcmp(type, 'bending2')
    error('not implemented');
else
    error('unidentified type');
end

% determine denominator of dof
% =[maxdofFE_node2, maxdofFE_node3, ..., maxdofFE_node11]
maxdofFE = abs(FE(end, idof));

% determine denominator of forc
maxforcFE = abs(FE(end, iforc));

Dfdisp = 0;     % contribution of displacement/rotation to objective function
Dfforc = 0;     % contribution of force/moment to objective function
for incr=1:11
    df1 = w1*sum(((FE(incr,idof)-BEAM(incr,idof))./maxdofFE).^2);    % rot/disp error
    Dfdisp = Dfdisp + df1;
       
    df2 = w2*((FE(incr,iforc)-BEAM(incr,iforc))/maxforcFE)^2; % force/moment
    Dfforc = Dfforc + df2;
end
f = Dfdisp + Dfforc;

out = f;
end
