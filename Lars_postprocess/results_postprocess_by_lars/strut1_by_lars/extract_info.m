% extract information

% in dis_and_rot, each line is X,Y,Z,Ux,Uy,Uz,angle, unit rotation axis

% dis_and_rot(:,104:106) are Ux,Uy,Uz of free end
% dis_and_rot(:,107:110) are theta,x,y,z of free end


% for bending 0, rotation at free end: angle = pi/2, axis = (1, 0, 0), other DoFs not fixed
% for bending 1, rotation at free end: angle = pi/2, axis = (sqrt(2)/2, 0, -sqrt(2)/2), other DoFs not fixed
% for bending 2, rotation at free end: angle = pi/2, axis = (0, 0, -1), other DoFs not fixed
% for torsion, rotation at free end: angle = 0.000004961336873, axis = (-1, 0, 0), other DoFs not fixed
% for tension, displacement at free end: Uy = -0.126723367373026, rest of nodes not fixed

%%
% load tension result
load('res_tot_strut1_tension.mat');

node = [dis_and_rot(1, 1:10:end)',dis_and_rot(1, 2:10:end)',dis_and_rot(1, 3:10:end)'];

Ux = dis_and_rot(:, 104);
Uy = dis_and_rot(:, 105);
Uz = dis_and_rot(:, 106);
Rx = [];
Ry = [];
Rz = [];
for i = 1:size(dis_and_rot,1)
    quat = dis_and_rot(i,107:end);
    if any(isnan(quat))
        Rx = [Rx;NaN];
        Ry = [Ry;NaN];
        Rz = [Rz;NaN];
    else
        Rx = [Rx;quat(1)*quat(2)];
        Ry = [Ry;quat(1)*quat(3)];
        Rz = [Rz;quat(1)*quat(4)];
    end
    
end

Fx = Forcex;
Fy = Forcey;
Fz = Forcez;
Mx = Momentx;
My = Momenty;
Mz = Momentz;


figure;
hold on;
U = -Uy;
F = Fy;
plot(U, F, '--sb');

%%

% load bending 0 result
load('res_tot_strut1_bending0.mat');

node = [dis_and_rot(1, 1:10:end)',dis_and_rot(1, 2:10:end)',dis_and_rot(1, 3:10:end)'];

Ux = dis_and_rot(:, 104);
Uy = dis_and_rot(:, 105);
Uz = dis_and_rot(:, 106);
Rx = [];
Ry = [];
Rz = [];
for i = 1:size(dis_and_rot,1)
    quat = dis_and_rot(i,107:end);
    if any(isnan(quat))
        Rx = [Rx;NaN];
        Ry = [Ry;NaN];
        Rz = [Rz;NaN];
    else
        Rx = [Rx;quat(1)*quat(2)];
        Ry = [Ry;quat(1)*quat(3)];
        Rz = [Rz;quat(1)*quat(4)];
    end
    
end

Fx = Forcex;
Fy = Forcey;
Fz = Forcez;
Mx = Momentx;
My = Momenty;
Mz = Momentz;


figure;
hold on;
U = Rx;
F = Mx;
plot(U, F, '--sb');
