% script to postprocess

datf = 'res_tot_strut1_Li.mat';
load(datf);

Nmidcs = size(center_loc_center,1);  % 10 cross sections

% check direction of nodes
center_loc_center(:,2);
% so cross sections are sequenced towards negative y direction

% center node = (end node 1 + end node 2)/2
for i = 1:Nmidcs
    err = center_loc_center(i,:)-(center_loc_ends(i,:)+center_loc_ends(i+1,:))/2;
    if norm(err)>eps
       error() 
    end
end

% global triad
R1 = [1, 0, 0];
R2 = [0, 1, 0];
R3 = [0, 0, 1];


Nppt = 28; % number of profile points of a cross section

node = center_loc_ends;
data = [];
trivert = [];
for i = 1:Nmidcs
    
    % coordinates of center point
    xc = center_loc_center(i, :);
    
    
    %-------------------------------------------------------
    % plot profile of this cross section
    figure;
    axis equal;
    hold on
    
    plot3(xc(1), xc(2), xc(3), '*b');
    
    % lines of Xintersect/Yintersect/Zintersect pertaining to current center point
    idxprofpt = 1+Nppt*(i-1):Nppt*i;
    plot3(Xintersect(idxprofpt), ...
          Yintersect(idxprofpt), ...
          Zintersect(idxprofpt), 'or');
    
    % add text
    for j = 1:numel(idxprofpt)
        idx = idxprofpt(j);     % index in Xintersect
        
        text(Xintersect(idx), ...
             Yintersect(idx), ...
             Zintersect(idx), ...
             [num2str(idx)],'HorizontalAlignment','left','FontSize',10);
    end
    hold off
    %-------------------------------------------------------

    
    
    % compute normal vector using 1st and 13th node in idxprofpt
    v1 = [Xintersect(idxprofpt(1)), ...
            Yintersect(idxprofpt(1)), ...
            Zintersect(idxprofpt(1))] - xc;
    v1 = v1/norm(v1);
        
    v2 = [Xintersect(idxprofpt(13)), ...
            Yintersect(idxprofpt(13)), ...
            Zintersect(idxprofpt(13))] - xc;
    v2 = v2/norm(v2);

    ndir = cross(v1, v2);
    if ndir(2) > 0
        ndir = -ndir;
    end
    
    
    % directional vector linking two end points of current cross section
    % mid point
    vendpt = center_loc_ends(i+1,:)-center_loc_ends(i,:);
    vendpt = vendpt/norm(vendpt);
    
    % ~1 suggesting the consistency of directions
    dot(vendpt, ndir);
    
    
    % local triad for the cross section
    r1 = vendpt;
    r2 = cross([0,0,1], r1);    % for r2 to be in global x-o-y plane
    r2 = r2/norm(r2);
    r3 = cross(r1, r2);
    
    % transformation matrix
    R = [r1*R1', r1*R2', r1*R3'
         r2*R1', r2*R2', r2*R3'
         r3*R1', r3*R2', r3*R3'];
    
    % translate + rotation
    ppt = [Xintersect(idxprofpt), Yintersect(idxprofpt), Zintersect(idxprofpt)];
    ppt = ppt - repmat(xc, Nppt, 1);
    
    ppt_loc = (R*ppt')';    % local coordinates of profile nodes, 
                            % x (w.r.t. r1) ~ 0, no absolute zero because
                            % of misalignment
    
    figure;
    axis equal;
    hold on
    plot(ppt_loc(:,2), ppt_loc(:,3), 'or');
    
    % generate minimum bounding triangle
    [trix,triy] = minboundtri(ppt_loc(:,2),ppt_loc(:,3));
    hold on
    plot(trix, triy, '-r','LineWidth',2);
    hold off
    
    vert = [trix(1:3), triy(1:3)];  % vert = [x1, y1;x2, y2;x3, y3]
    cross([0, vert(2,:)-vert(1,:)], [0, vert(3,:)-vert(1,:)]);  % prove vertices are anti-clockwise
    
    % add initial beam orientation
    data = [data; r1, r2, r3];
    % record vertex
    trivert = [trivert; reshape(vert', 1, [])];
end

save('node_data_trivert_strut1.mat', 'node', 'data', 'trivert');
