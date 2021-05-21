function [Po,Pi] = genTriRing(Po,t)
% This function moves each edge of a 2d triangle inward by a thickness t
% and returns a triangular ring

% INPUT:
%   Po      : 3*2 matrix, each row is x, y coordinates of a vertex,
%             sequence of vertices are anti-clockwise
%   t       : thickness of triangular ring

% OUTPUT:
%   Po      : 3*2 matrix, each row is x, y coordinates of a vertex of outer 
%             triangle, of anti-clockwise sequence
%   Pi      : 3*2 matrix, each row is x, y coordinates of a vertex of inner 
%             triangle, of anti-clockwise sequence

%=========================================================================%

% % example parameters
% Po = [0, 0
%       1, 0
%       0, 1];
% t = 0.05;


% extend Po
Pout = [Po; Po(1, :)];

lines = cell(4,1);
for i = 1:3
    P1 = Pout(i, :);
    P2 = Pout(i+1, :);

    dx = P2(1) - P1(1);
    dy = P2(2) - P1(2);
    
    % normal vector to the right 
    n = [-dy, dx];
    n = n/norm(n);
    
    disp = n*t;
    
    P1n = P1 + disp;
    P2n = P2 + disp;
    
    lines{i+1} = [P1n;        % x1, y1
                  P2n];       % x2, y2
end

% extend lines
lines{1} = lines{end};

Pi = [];
for i = 1:3

    [x, y] = findintersection(lines{i}, lines{i+1});

    Pi = [Pi; x, y];
end

% % plot
% figure;
% hold on;
% axis equal;
% 
% plot([Po(:,1); Po(1,1)], [Po(:,2); Po(1,2)], '-or');
% plot([Pi(:,1); Pi(1,1)], [Pi(:,2); Pi(1,2)], '-sb');
% 
% hold off

end



