function UAV_Patch
clear all;
close all;

figure('color','w');
view(3);
hold on;
%% Rods:
nRod = 100;
hRod = 10;
rRod = .18;
EAnglesRod1 = [0;3*pi/4;0];
EAnglesRod2 = [pi;pi/4;pi/2];
EAnglesRod3 = [pi/4;pi/2;pi/2];
[facesRod verticesRod1] = cylrod(rRod,-hRod,hRod,nRod,EAnglesRod1);
[facesRod verticesRod2] = cylrod(rRod,-hRod,hRod,nRod,EAnglesRod2);
[facesRod verticesRod3] = cylrod(rRod,-hRod,hRod,nRod,EAnglesRod3);
patch('faces',facesRod,'vertices',verticesRod1,'FaceColor',[.5 .5 .5],'EdgeColor',[.5 .5 .5]);
patch('faces',facesRod,'vertices',verticesRod2,'FaceColor',[.5 .5 .5],'EdgeColor',[.5 .5 .5]);
patch('faces',facesRod,'vertices',verticesRod3,'FaceColor',[.5 .5 .5],'EdgeColor',[.5 .5 .5]);
xlabel('x')
ylabel('y')
zlabel('z')
axis(5*[-10 10 -10 10 -10 10]);

%% Propellers
nProp = 100;
hProp = .5;
rProp = 4;
EAnglesProp1 = [0;0;0];
EAnglesProp2 = [0;0;0];
EAnglesProp3 = [0;pi/2;0];
EAnglesProp4 = [0;pi/2;0];
EAnglesProp5 = [0;pi/2;pi/2];
EAnglesProp6 = [0;pi/2;pi/2];
[facesProp verticesProp1] = cyl(rProp,0,hProp,nProp,EAnglesProp1);
verticesProp1(:,1) = verticesProp1(:,1)+hRod/sqrt(2)*ones(2*nProp+2,1);
verticesProp1(:,2) = verticesProp1(:,2)+hRod/sqrt(2)*ones(2*nProp+2,1);
hprop(1)=patch('faces',facesProp,'vertices',verticesProp1,'FaceColor',[0 0 1]);
[facesProp verticesProp2] = cyl(rProp,-hProp,0,nProp,EAnglesProp2);
verticesProp2(:,1) = verticesProp2(:,1)-hRod/sqrt(2)*ones(2*nProp+2,1);
verticesProp2(:,2) = verticesProp2(:,2)-hRod/sqrt(2)*ones(2*nProp+2,1);
hprop(2)=patch('faces',facesProp,'vertices',verticesProp2,'FaceColor',[0 0 1]);
[facesProp verticesProp3] = cyl(rProp,-hProp,0,nProp,EAnglesProp3);
verticesProp3(:,3) = verticesProp3(:,3)+hRod/sqrt(2)*ones(2*nProp+2,1);
verticesProp3(:,2) = verticesProp3(:,2)+hRod/sqrt(2)*ones(2*nProp+2,1);
hprop(3)=patch('faces',facesProp,'vertices',verticesProp3,'FaceColor',[0 0 1]);
[facesProp verticesProp4] = cyl(rProp,0,hProp,nProp,EAnglesProp4);
verticesProp4(:,3) = verticesProp4(:,3)-hRod/sqrt(2)*ones(2*nProp+2,1);
verticesProp4(:,2) = verticesProp4(:,2)-hRod/sqrt(2)*ones(2*nProp+2,1);
hprop(4)=patch('faces',facesProp,'vertices',verticesProp4,'FaceColor',[0 0 1]);
[facesProp verticesProp5] = cyl(rProp,0,hProp,nProp,EAnglesProp5);
verticesProp5(:,1) = verticesProp5(:,1)+hRod/sqrt(2)*ones(2*nProp+2,1);
verticesProp5(:,3) = verticesProp5(:,3)-hRod/sqrt(2)*ones(2*nProp+2,1);
hprop(5)=patch('faces',facesProp,'vertices',verticesProp5,'FaceColor',[0 0 1]);
[facesProp verticesProp6] = cyl(rProp,-hProp,0,nProp,EAnglesProp6);
verticesProp6(:,1) = verticesProp6(:,1)-hRod/sqrt(2)*ones(2*nProp+2,1);
verticesProp6(:,3) = verticesProp6(:,3)+hRod/sqrt(2)*ones(2*nProp+2,1);
hprop(6)=patch('faces',facesProp,'vertices',verticesProp6,'FaceColor',[0 0 1]);

axis equal;
set(hprop,'FaceLighting','phong','AmbientStrength',0.8,'Linestyle','none');
alpha(hprop,0.6);
set(gca,'DataAspectRatioMode','manual','visible','off');
light('Position',[-100 -100 50],'Style','infinite');
set(gca,'box','on','XTick',[],'YTick',[],'ZTick',[]);

view(-158,78);

print('hexrotorDiagram','-dpng');
end

function [faces vertices] = cyl(r,h1,h2,n,EAngles)
theta = linspace(0,2*pi,n);
x = r*cos(theta);
y = r*sin(theta);
top = transpose(EAngles321(EAngles)*[x',y',h2*ones(n,1)]');
bottom = transpose(EAngles321(EAngles)*[x',y',h1*ones(n,1)]');
vertices = [bottom;top;0 0 h1; 0 0 h2];
faces = [(1:n)',(n+1:2*n)',((n+1:2*n)+1)'];
faces(end,end) = 1;
faces = [faces;(1:n)',(2:n+1)',((n+1:2*n)+1)'];
faces(end,2) = 1;
faces(end,3) = n+1;
faces=[faces; (1:n-1)' (2:n)' (2*n+1)*ones(n-1,1)];
faces=[faces; (n+1:2*n-1)' (n+2:2*n)' (2*n+2)*ones(n-1,1)];

end

function [faces vertices] = cylrod(r,h1,h2,n,EAngles)
theta = linspace(0,2*pi,n);
x = r*cos(theta);
y = r*sin(theta);
top = transpose(EAngles321(EAngles)*[x',y',h2*ones(n,1)]');
bottom = transpose(EAngles321(EAngles)*[x',y',h1*ones(n,1)]');
vertices = [bottom;top];
faces = [(1:n)',(n+1:2*n)',((n+1:2*n)+1)'];
faces(end,end) = 1;
faces = [faces;(1:n)',(2:n+1)',((n+1:2*n)+1)'];
faces(end,2) = 1;
faces(end,3) = n+1;

end

