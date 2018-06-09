clear all;
close all;

figure;
[x, y, z] = sphere(1000);
% set(gca, 'Color',[0.3 0.1 0.3]);
hold on;

patch(surf2patch(x,y,z,z), 'FaceColor','r', 'EdgeColor','r');
view(3);
% l = 0.5;% completely indadmissible
l = 0.7;% non-connecting
% l = 0.75;% connecting
% l = 1.1;% completely admissible
topsurf = l*[-1 1 1 -1; -1 -1 1 1; 1 1 1 1];
bottomsurf = l*[-1 1 1 -1; -1 -1 1 1; -1 -1 -1 -1];
backsurf = l*[1 1 1 1; -1 -1 1 1; -1 1 1 -1];
frontsurf = l*[-1 -1 -1 -1; -1 -1 1 1; -1 1 1 -1];
leftsurf = l*[-1 -1 1 1; 1 1 1 1; -1 1 1 -1];
rightsurf = l*[-1 -1 1 1; -1 -1 -1 -1; -1 1 1 -1];
patch(topsurf(1,:), topsurf(2,:), topsurf(3,:), 'b', 'EdgeColor','k');
patch(bottomsurf(1,:), bottomsurf(2,:), bottomsurf(3,:), 'b', 'EdgeColor','k');
patch(backsurf(1,:), backsurf(2,:), backsurf(3,:), 'b', 'EdgeColor','k');
patch(frontsurf(1,:), frontsurf(2,:), frontsurf(3,:), 'b', 'EdgeColor','k');
patch(leftsurf(1,:), leftsurf(2,:), leftsurf(3,:), 'b', 'EdgeColor','k');
patch(rightsurf(1,:), rightsurf(2,:), rightsurf(3,:), 'b', 'EdgeColor','k');

set(gca, 'XTick', [], 'YTick', [], 'ZTick', []);
axis equal;
axis off;

% figure;
% 3D Graphics: Cube
% Dr. P.Venkataraman
% 
% set(gcf,'Menubar','none','Name','Cube', ...
%     'NumberTitle','off','Position',[10,350,300,200], ...
%     'Color',[0.3 0.1 0.3]);
% % the cube
% h(1) = axes('Position',[0.2 0.2 0.6 0.6]);
% vert = [1 1 1; 1 2 1; 2 2 1; 2 1 1 ; ...
%         1 1 2;1 2 2; 2 2 2;2 1 2];
% fac = [1 2 3 4; ...
%     2 6 7 3; ...
%     4 3 7 8; ...
%     1 5 8 4; ...
%     1 2 6 5; ...
%     5 6 7 8];
% 
% patch('Faces',fac,'Vertices',vert,'FaceColor','r');  % patch function
% light('Position',[1 3 2]);
% light('Position',[-3 -1 3]);
% material shiny;
% alpha('color');
% alphamap('rampdown');
% camlight(45,45);
% lighting phong
% view(30,30);