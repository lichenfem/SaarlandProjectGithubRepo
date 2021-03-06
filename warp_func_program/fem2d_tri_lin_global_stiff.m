function K = fem2d_tri_lin_global_stiff(coords, ien, bound_edge)
% Generate the global stiffness matrix for $-\nabla^2 u$ in 2D case
% using linear basis functions and triangular elements
% [IN]  coords : n * 2 matrix, n is the number of grid points,
%                each row has the x and y coordinate of a grid point
% [IN]  ien    : m * 3 matrix, m is the number of triangular elements,
%                each row has 3 vertex point ids of the element, in 
%                counter clockwise order
% [IN]  bgp    : row vector, contains IDs of boundary nodes
% [OUT] K      : n * n Global stiffness matrix
%=========================================================================%
	n = size(coords, 1);                                                   % n: total number of nodes
	m = size(ien, 1);                                                      % m: total number of tri elem
	
	K = sparse(n, n);
	
	for i_elem = 1 : m
		elem_vertex_ids = ien(i_elem, :);                                  % elem_vertex_ids = [i, j, k]
		vertex_coords = coords(elem_vertex_ids, :);                        % vertex_coords = [xi, yi
                                                                           %                  xj, yj
                                                                           %                  xk, yk]
		
		% Transpose vertex_coords to use in local_stiff_quad
		k = fem2d_tri_lin_unit_stiff(vertex_coords', elem_vertex_ids, bound_edge);
		
		% Assemble to global stiffness matrix
		K(elem_vertex_ids, elem_vertex_ids) = K(elem_vertex_ids, elem_vertex_ids) + k;
	end
end