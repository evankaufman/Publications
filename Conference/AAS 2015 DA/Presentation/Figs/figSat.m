clear all;
close all;

a=figure('units','normalized','position',[0.1 0.1 0.7 0.7],'Color','w','renderer','opengl');


load('topo.mat','topo','topomap1');
[x,y,z] = sphere(100);
props.AmbientStrength = 0.1;
props.DiffuseStrength = 1;
props.SpecularColorReflectance = .5;
props.SpecularExponent = 20;
props.SpecularStrength = 1;
props.FaceColor= 'texture';
props.EdgeColor = 'none';
props.FaceLighting = 'phong';
props.Cdata = topo;
radE=8;
surface(radE*x,radE*y,radE*z,props);


N=20;
a=linspace(0,2*pi,N)';
vertice.upper=[sin(a) 1*ones(N,1) cos(a)];
vertice.lower=[sin(a) -1*ones(N,1) cos(a)];
W=1.6;
l=5;
d=0.1;
vertice.plate=[-l -d W/2;
    l -d W/2;
    l d -W/2;
    -l d -W/2];
W=0.2;
l=3;
d=0.5;
vertice.pin=[0 W/2 -l;
    0 W/2 0;
    0 -W/2 0;
    0 -W/2 -l];

vertice=[vertice.upper;vertice.lower;0 1 0;0 -1.0 0;vertice.plate;vertice.pin];
faces=[1:N-1 1:N-1 N+1:2*N-1 2*N+3 2*N+7;
    2:N 2:N N+2:2*N 2*N+4 2*N+8;
    N+2:2*N (2*N+1)*ones(1,N-1) (2*N+2)*ones(1,N-1) 2*N+5 2*N+9;
    N+1:2*N-1 (2*N+1)*ones(1,N-1) (2*N+2)*ones(1,N-1) 2*N+6 2*N+10]';
set(gca,'Position',[0.025 0.025 0.95 0.95],'XColor',[1 1 1],'YColor',[1 1 1],'ZColor',[1 1 1]);
light('Position',[-100 -200 20],'Style','Local');

r0=14;
orbit.x=r0*cos(linspace(0,2*pi,100));
orbit.z=r0*sin(linspace(0,2*pi,100));
line(orbit.x,0*orbit.x,orbit.z,'Color',[0.1 0.1 0.1]);

theta=20;
rc=r0*[cosd(theta); 0; sind(theta)];
rd=rc+[6; 0; 12];

acolor=[0 0 0];
lwidth=1;
alength=1;
awidth=alength*tand(12);


bb=6;
e1=[1 0 0]';
e2=[0 1 0]';
e3=[0 0 1]';
Rc=expmso3(30*pi/180*e1)*expmso3(30*pi/180*e2);
Rd=expmso3(-10*pi/180*e1)*expmso3(60*pi/180*e2)*expmso3(170*pi/180*e3);
for aa=1:size(vertice,1)
    verticec(aa,:)=(Rc*vertice(aa,:)')'+rc';
    verticed(aa,:)=(Rd*vertice(aa,:)')'+rd';
end

hc=patch('Vertices',verticec,'Faces',faces,'FaceColor',[ 0 0 1]);
set(hc,'FaceLighting','phong','AmbientStrength',0.5);
myarrow([0 0 0],[0 0 -r0],acolor,lwidth,alength,awidth)

axis equal;
view(0,0);
axis(r0*[-1.2 1.8 -1.2 1.8 -1.2 1.8]);

print('PF0','-depsc2');

myarrow(rc',rc'+8*[cosd(theta); 0; sind(theta)]',acolor,lwidth,alength,awidth)
myarrow(rc',rc'+8*[-sind(theta); 0; cosd(theta)]',acolor,lwidth,alength,awidth)


print('PF1','-depsc2');


hd=patch('Vertices',verticed,'Faces',faces,'FaceColor',[ 1 0 0]);
set(hd,'FaceLighting','phong','AmbientStrength',0.5);
acolor=[1 0 0];
lwidth=2;
alength=1.5;
[ha haa]=myarrow(rc',0.92*(rd'-rc')+rc',acolor,lwidth,alength,awidth)

print('PF2','-depsc2');

delete([ha haa]);
tmp=[rc rd];
hl=line(tmp(1,:),tmp(2,:),tmp(3,:));
set(hl,'Color','k','LineStyle',':');
[ha haa]=myarrow(rc',0.5*(rd'-rc')+rc',acolor,lwidth,alength,awidth);

print('PF3','-depsc2');