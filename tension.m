% inputfile.m

% import coordinates of beam axis

load('./node_data_trivert_strut1.mat', 'node', 'data');

nnode = size(node, 1);
nelem = nnode - 1;

% beam element type
eltp = ['crtw3dp'];

% allocate cross sections
geom = cell(nelem, 2);
for i = 1:nelem
    geom{i,1} = ['./CroSec', num2str(i), '/tri_mesh/coords_ien_bgp.mat'];
    geom{i,2} = ['./CroSec', num2str(i), '/warp_func/wbar.mat'];
end

% set material properties, = [E, v, f_y, Et]
matr = [6.045577e+04, 0.33, 6.2301e+02, 1.39700e+03];

% number of beam FEs per segment
N = 1;

nodetmp = node(1,:);
elem = [];
for i = 1:size(node,1)-1
    
    tau = node(i+1, :) - node(i, :);    % length of current segment
    
    for j = 1:N
    	nodetmp = [nodetmp; node(i, :)+tau/N*j];
        
        % add a beam,=[eltp_id, geom_id, matr_id, data_id, node_id, node_id]
        elem = [elem; ...
                1, i, 1, i, size(nodetmp,1)-1, size(nodetmp,1)];
    end

end

node  = nodetmp;
nnode = size(node,1);
nelem = size(elem,1);


% tension boundary conditions
pdof = [1, 1, 0   % Ux
        1, 2, 0   % Uy
        1, 3, 0   % Uz
        1, 4, 0   % Rx
        1, 5, 0   % Ry
        1, 6, 0   % Rz
        1, 7, 0   % wapring DoF of a node
        nnode, 2, -0.126723367373026      % Uy
        ];                      % prescribed DoFs


% % bending 0 boundary conditions
% pdof = [1, 1, 0     % Ux
%         1, 2, 0     % Uy
%         1, 3, 0     % Uz
%         1, 4, 0     % Rx
%         1, 5, 0     % Ry
%         1, 6, 0     % Rz
%         1, 7, 0     % wapring DoF of a node
%         nnode, 4, pi/2
%         ];                      % prescribed DoFs
    
    
nodf = [];


tolr = [ 1.0e-4 ];    % tolerance of convergence
maxniter = [ 20 ];    % maximum allowable interations in an increment before refinement


% reac = [nnode, 1      %1=Fx
%         nnode, 2      %2=Fy
%         nnode, 3      %3=Fz
%         nnode, 4      %4=Mx
%         nnode, 5      %5=My
%         nnode, 6      %6=Mz
%         nnode, 7];    % output reaction forces

for reaktion=1:6
    for steg=2:11
        reac(reaktion+(steg-2)*6,1)=steg;
        reac(reaktion+(steg-2)*6,2)=reaktion;
    end
end


incr = 0:0.1:1;                                                            % load increments

outi = 1:11;                                                      % load increment steps to be recorded, outi = 0 means to output initial setting

outp = ['displa'
        'rotati'];                                                         % nodal output

restart = 0;                                                               % restart =1 would save latest converged increment to file

%%
% %-----below is the script to generate vtu file for the beams-------------%
% point = [];
% cell = [];
% 
% 
% for ielem = 1:size(elem,1)
%     
%     igeom = elem(ielem, 2);
%     idata = elem(ielem, 4);
%     
%     inode1 = elem(ielem, 5);
%     inode2 = elem(ielem, 6);
%     
%     load(cell2mat( geom(igeom,1) ), 'coords', 'ien');
%     
%     % extend coords
%     coords = [zeros(size(coords,1),1), coords];
%     
%     % rotation matrix
%     r1 = data(idata, 1:3)';
%     r2 = data(idata, 4:6)';
%     r3 = data(idata, 7:9)';
%     R = [r1, r2, r3];
%     %
%     for icoords = 1:size(coords, 1)
%         xyz = coords(icoords,:);
%         xyz = (R*(xyz'))';
%         coords(icoords,:) = xyz;
%     end
%     
%     % two cross sections for two beam ends
%     sec1 = coords + repmat(node(inode1,:), size(coords,1), 1);
%     sec2 = coords + repmat(node(inode2,:), size(coords,1), 1);
%     
%     patch('faces',ien,'vertices',sec1,'FaceColor',[0.8 0.8 1.0],'FaceLighting','gouraud');
%     hold on; axis equal;
%     patch('faces',ien,'vertices',sec2,'FaceColor',[0.8 0.8 1.0],'FaceLighting','gouraud');
%     
%     % add cluster of triangular prisms to point, cell
%     pos1 = size(point,1);
%     point = [point; sec1];  % add points of 1st section
%     pos2 = size(point,1);
%     point = [point; sec2];  % add points of 2nd section
%     pos3 = size(point,1);
%     
%     cell = [cell; replace(ien, 1:size(sec1,1), pos1+1:pos2), ...
%                   replace(ien, 1:size(sec2,1), pos2+1:pos3)];
%     
% end
% 
% save('genvtu.mat', 'point', 'cell');

