function k = fem2d_tri_lin_unit_stiff(x, elem_vertex_ids, bound_edge)
% Generate the unit stiffness matrix for a triangular element 
% using linear basis functions
% [IN]  x : 2 * 3 matrix, the geometric coordinates of the element's nodes,
%           first row is x coordinates, second row is y coordinates,
%           points should be in counter clockwise order
% [IN]  elem_vertex_ids: row vector, IDs of three vertices of tri elem in
%           anti-clockwise
% [IN]  bound_edge: each row contains IDs of two nodes, which forms a
%           segment of domain boundary
% [OUT] k : 3 * 3 unit stiffness matrix, k(i, j) = \int_{\Omega^e}
%           phi_{i}^{'}(x, y) phi_{j}^{'}(d, y) dx dy
%           phi_{1, 2, 3} are in counter clockwise order
%=========================================================================%
	k = zeros(3, 3);
	
	% Gauss quadrature points and weight function
	qx = [1.0/3.0, 0.6, 0.2, 0.2];
	qy = [1.0/3.0, 0.2, 0.2, 0.6];
	w  = [-27.0/96.0, 25.0/96.0, 25.0/96.0, 25.0/96.0];
	n_quadrature = size(w, 2);
	
	% In-element 2D integral for $\phi_i^{'} * \phi_j^{'}$
	% Numerical integral using Gauss quadrature
	for iq = 1 : n_quadrature
		[sh, dtm] = fem2d_tri_lin_shape(qx(iq), qy(iq), x);
		d_N = sh(1 : 2, :);
		k = k + d_N' * d_N * dtm * w(iq);
	end
	
	% Boundary line integral for alpha(x, y)
	k2 = zeros(3, 3);
    
    if ismember(elem_vertex_ids([1,2]), bound_edge, 'rows')
        k2 = k2 + fem2d_tri_lin_int_alpha(x, 1, 2);
    end
    %
    if ismember(elem_vertex_ids([2,3]), bound_edge, 'rows')
        k2 = k2 + fem2d_tri_lin_int_alpha(x, 2, 3);
    end
    %
    if ismember(elem_vertex_ids([3,1]), bound_edge, 'rows')
        k2 = k2 + fem2d_tri_lin_int_alpha(x, 3, 1);
    end
	
	k = k + k2;
end