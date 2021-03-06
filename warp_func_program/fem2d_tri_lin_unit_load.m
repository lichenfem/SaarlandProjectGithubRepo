function ub = fem2d_tri_lin_unit_load(x, elem_vertex_ids, bound_edge)
% Use variable substitution to transform the f(x, y) and phi_i(x, y) on a
% element to triangle (0,0)-(1,0)-(0,1), then use Gauss quadrature to integrate
% [IN]  x  : 2 * 3 matrix, the geometric coordinates of the element's vertexes
%            1st row is x coordinates, 2nd row is y coordinates, 
%            vertexes shoule be in counter clockwise order\                % x = [x1, x2, x3
%                                                                                 y1, y2, y3]
% [OUT] ub : 3 * 1 vector, integral result
%=========================================================================%
	qx = [1.0/3.0, 0.6, 0.2, 0.2];
	qy = [1.0/3.0, 0.2, 0.2, 0.6];
	w  = [-27.0/96.0, 25.0/96.0, 25.0/96.0, 25.0/96.0];
	n_quadrature = size(w, 2);
	
	ub = zeros(3, 1);
	
	% $\int_{\Omega^e} f * \phi ds$
	for iq = 1 : n_quadrature
		[f_x, f_y] = fem2d_tri_xi_eta_to_x_y(qx(iq), qy(iq), x);
		f_xy = poisson2d_rhs_f(f_x, f_y);
		
		dtm = fem2d_tri_lin_det(qy(iq), qy(iq), x);
		
		N(1) = 1.0 - qx(iq) - qy(iq);
		N(2) = qx(iq);
		N(3) = qy(iq);
		
		ub(1) = ub(1) + dtm * f_xy * N(1) * w(iq);
		ub(2) = ub(2) + dtm * f_xy * N(2) * w(iq);
		ub(3) = ub(3) + dtm * f_xy * N(3) * w(iq);
    end
	
    
	% $\int_{\partial \Omega^e} g * \phi dt $
    if ismember(elem_vertex_ids([1,2]), bound_edge, 'rows')
        %
        tv = [x(1,2)-x(1,1)
              x(2,2)-x(2,1)];                                              % tangential vector
        R  = [cos(-pi/2), -sin(-pi/2)
              sin(-pi/2),  cos(-pi/2)];                                    % rotation matrix
        nv = R*tv; nv = nv/norm(nv);                                       % normal vector
        ny = dot(nv, [1, 0]);
        nz = dot(nv, [0, 1]);
        %
        ub = ub + fem2d_tri_lin_int_g(x, 1, 2, ny, nz);
    end
    
    if ismember(elem_vertex_ids([2,3]), bound_edge, 'rows')
        %
        tv = [x(1,3)-x(1,2)
              x(2,3)-x(2,2)];                                              % tangential vector
        R  = [cos(-pi/2), -sin(-pi/2)
              sin(-pi/2),  cos(-pi/2)];                                    % rotation matrix
        nv = R*tv; nv = nv/norm(nv);                                       % normal vector
        ny = dot(nv, [1, 0]);
        nz = dot(nv, [0, 1]);
        %
        ub = ub + fem2d_tri_lin_int_g(x, 2, 3, ny, nz);
    end
    
    if ismember(elem_vertex_ids([3,1]), bound_edge, 'rows')
        %
        tv = [x(1,1)-x(1,3)
              x(2,1)-x(2,3)];                                              % tangential vector
        R  = [cos(-pi/2), -sin(-pi/2)
              sin(-pi/2),  cos(-pi/2)];                                    % rotation matrix
        nv = R*tv; nv = nv/norm(nv);                                       % normal vector
        ny = dot(nv, [1, 0]);
        nz = dot(nv, [0, 1]);
        %
        ub = ub + fem2d_tri_lin_int_g(x, 3, 1, ny, nz);
    end

end