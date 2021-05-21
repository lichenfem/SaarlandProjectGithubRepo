function [wb, uTb] = test_fem2d()
% Test my 2D FEM solver
	
	% import tri mesh
    load('../tri_mesh/coords_ien_bgp.mat', 'coords', 'ien', 'bgp', 'bound_edge');
	
	% Generate the global stiffness matrix
	K = fem2d_tri_lin_global_stiff(coords, ien, bound_edge);
	
	% Generate the right hand side
	b = fem2d_tri_lin_global_load(coords, ien, bound_edge);
	
	% Apply the Dirichlet boundary condition
	% If we want to use Dirichlet boundary condition, set use_dirichlet_bc = 1,
	% and set poisson2d_robin_bc_alpha() = 0, poisson2d_robin_bc_g() = 0
	use_dirichlet_bc = 0;
	if (use_dirichlet_bc == 1)
		n_bgp = max(size(bgp));
		for i = 1 : n_bgp
			bgp_id = bgp(i);
			K(bgp_id, :) = 0;
			K(bgp_id, bgp_id) = 1;
			bv_x = coords(bgp_id, 1);
			bv_y = coords(bgp_id, 2);
			bv_val = poisson2d_dirichlet_bc_g(bv_x, bv_y);
			b(bgp_id) = bv_val;
		end
    end
    
    % suppress u of an arbitrary point
    bgp_id = bgp(1);
    K(bgp_id, :) = 0;
    K(bgp_id, bgp_id) = 1;
    b(bgp_id) = 100;
	
	% Solve the linear system
	u = K \ b;
	
    %% calibrate wb
    wb = u;
    
    A = sum(simpvol(coords,ien));                                              % A: total area of tri elem
    n = size(coords, 1);                                                       % n: total number of nodes
    m = size(ien, 1);                                                          % m: total number of tri elem
    % compute $\bar{\omega} \:= \bar{\omega}-1/A*\int_{A} \bar{\omega} dA$
    int_wb = 0;
    for i_elem = 1 : m
        elem_vertex_ids = ien(i_elem, :);                                      % elem_vertex_ids: node IDs of tri elem, =[i, j, k]
        vertex_coords = coords(elem_vertex_ids, :);                            % vertex_coords: nodal coordinates, [xi, yi
                                                                               %                                    xj, yj
                                                                               %                                    xk, yk]
        % for each tri elem
        % Gauss quadrature points and weight function
        qx = [1.0/3.0, 0.6, 0.2, 0.2];
        qy = [1.0/3.0, 0.2, 0.2, 0.6];
        w  = [-27.0/96.0, 25.0/96.0, 25.0/96.0, 25.0/96.0];
        n_quadrature = size(w, 2);

        % In-element 2D integral for $\bar{\omega}$
        % Numerical integral using Gauss quadrature
        for iq = 1 : n_quadrature
            [sh, dtm] = fem2d_tri_lin_shape(qx(iq), qy(iq), vertex_coords');   % sh = [ dN1/dx, dN2/dx, dN3/dx
                                                                               %        dN1/dy, dN2/dy, dN3/dy
                                                                               %        N1,     N2,     N3     ]
                                                                               % dtm: determinant of jacobian
            int_wb = int_wb + dot(sh(end,:), wb(elem_vertex_ids)) * dtm * w(iq);
        end
    end
    %
    wb = wb - int_wb/A;
    uTb = dot(wb, K*wb);                                                   % for the computation of J(Saint Venant torsion modulus)
    
    save('wbar.mat', 'wb', 'uTb');
    
%     % Plot the result
%     trisurf(ien, coords(:,1),coords(:,2), wb(:)); axis equal;
%     colorbar
%     view([0,0,1]);
%     savefig('warp_func.fig');
    
    close all;
end