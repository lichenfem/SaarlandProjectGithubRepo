function res = poisson2d_robin_bc_g(x, y, nx, ny)
% Function $g(x, y)$ in Robin boundary condition 
% $\frac{\partial u}{\partial \vec{n}} + \alpha * u = g$
% If you want to use Dirichlet boundary condition, set res = 0
	res = x + y;
	%res = 0;
    res = nx*y - ny*x;  % lichen
end