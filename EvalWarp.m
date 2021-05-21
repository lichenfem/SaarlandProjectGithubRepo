
% add necessary folders to path
warp_func_code = './warp_func_program';
addpath(genpath(warp_func_code));
homedir = pwd;		% where CroSecXX folders exist

for ii = 1:10

    dir = ['CroSec', num2str(ii)];

    fprintf(1, 'Evaluating warp function for %s\n', dir);
    
    eval(['cd ', dir]);

    % create and fill subfolder warp_func
    subfolder = 'warp_func';
    eval(['mkdir ', subfolder]);
    eval(['cd ', subfolder]);      % dive into the subfolder
    [wb, uTb] = test_fem2d();
    
    eval(['cd ', homedir]);
end
% approximately 1.4s


