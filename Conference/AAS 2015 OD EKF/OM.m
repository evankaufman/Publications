function OM
% Observability criteria and measure for the nonlinear relative orbit
global mu n a
close all;

%%% Initial conditions and Parameters

% Fixed parameters
mu=398600;      % Earth gravitational parameter
RE=6378;        % Earth radius

% Orbital properties of chief
a=7100;
n=sqrt(mu/a^3); % Mean motion

% Simulation time
T=2*pi/n;       % Orbital period of the chief
tf=2*T;        % Terminal simulation time
N=2*6000+1;    % Total number of time steps
t=linspace(0,tf,N);  % Time vector

% Simulation title
filename='OM_IOD_0A';
disp(['Case: ' filename]);

% Initial condition
xx0=[ -7.100116477545843e-02
    -9.913332955272308e-04
    1.307672689210581e-01
    -2.329550204657181e-06
    -9.913316565799100e-04
    1.307672203750978e-01];


%%% Observability Test for the Linear HCW Model
disp('Observablity test for the linear HCW model... ');
ode_options=odeset('RelTol',1e-6,'AbsTol',1e-6);    % Options for ode45
[t xx_HCW]=ode45(@eom_HCW,t,xx0,ode_options);
xx_HCW=xx_HCW';     % State vectors for the HCW model (6 x N)

for k=1:N
    O_H_svd(:,k)=ObservabilityTest_HCW(xx_HCW(:,k));
end

%%% Observability Test for the Nonlinear Model
disp('Observablity test for the nonlinear model... ');
[t xx_TBP]=ode45(@eom_RelTBP,t,xx0,ode_options);
xx_TBP=xx_TBP';     % State vectors for the nonlinear two-body model (6 x N)

for k=1:N
    [O_svd(:,k) O_criteria(:,k)]=ObservabilityTest_TBP(xx_TBP(:,k));
end

%%% Observability Gramian
disp('Observablity gramian for the nonlinear model... ');
[Wo_HCW Wo_svd_HCW Wo_cond_HCW]=ObservabilityGramian_HCW(t,xx_HCW);
eps=1e-8;
[Wo Wo_svd Wo_cond]=ObservabilityGramian(t,xx0,eps,ode_options);

%%% Post Processing

save(filename);             % save results
plot_results(filename);     % plot figures
evalin('base',['load ' filename ';']);  % load results to workspace
end


function Xdot=eom_TBP(t,X)
% Nonlinear equations of motion for the two-body problem
global mu
r=X(1:3);
r_dot=X(4:6);

r_ddot=-mu/norm(r)^3*r;

Xdot=[r_dot; r_ddot];
end

function Xdot=eom_RelTBP(t,X)
% Nonlinear equations of motion for the relative motion
global mu a n
x=X(1);
y=X(2);
z=X(3);
xdot=X(4);
ydot=X(5);
zdot=X(6);

ra=[x+a; y; z];
nra=norm(ra);
xddot=2*n*ydot+n^2*x+n^2*a-mu*(x+a)/nra^3;
yddot=-2*n*xdot+n^2*y-mu*y/nra^3;
zddot=-mu*z/nra^3;

Xdot=[xdot ydot zdot xddot yddot zddot]';
end

function Xdot=eom_HCW(t,X)
% Linear equation of motion for the HCW model
global mu a n
A=[0 0 0 1 0 0;
    0 0 0 0 1 0;
    0 0 0 0 0 1;
    3*n^2 0 0 0 2*n 0;
    0 0 0 -2*n 0 0;
    0 0 -n^2 0 0 0];

Xdot=A*X;
end

function O_H_svd=ObservabilityTest_HCW(xx)
global mu n a

r=xx(1:3);
rdot=xx(4:6);
nr=norm(r);

ra=r+[a 0 0]';
nra=norm(ra);

Sw=[0 -n 0;n 0 0; 0 0 0];

A=[0 0 0 1 0 0;
    0 0 0 0 1 0;
    0 0 0 0 0 1;
    3*n^2 0 0 0 2*n 0;
    0 0 0 -2*n 0 0;
    0 0 -n^2 0 0 0];
A21=diag([3*n^2 0 -n^2]);

xxdot=A*xx;
rddot_H=xxdot(4:6);

h=r/nr;
y=h;

O_00=1/nr*(eye(3)-y*y');
O_10=-1/nr^2*(rdot*y'+y'*rdot*eye(3)+y*rdot'-3*y*y'*rdot*y');
O_11=O_00;

tmp1=-(2*rdot*rdot.'+rdot.'*rdot*eye(3))/nr^3 ...
    +3*(2*rdot*r.'*rdot*r.'+r*rdot.'*rdot*r.'+r.'*rdot*r.'*rdot*eye(3)+2*r*r.'*rdot*rdot.')/nr^5 ...
    -15*(r*r.'*rdot*r.'*rdot*r.')/nr^7;
tmp3_H=-(rddot_H*r.'+r.'*rddot_H*eye(3)+r*rddot_H.')/nr^3 ...
    +3*r*r.'*rddot_H*r.'/nr^5 ...
    +(eye(3)/nr-r*r'/nr^3)*diag([3*n^2 0 -n^2]);

O_20_H=tmp1+tmp3_H;
O_21=2*O_10-2*O_00*Sw;

O_30_H=(9*rdot.'*rdot*rdot*r.')/nr^5 ...
    -15*(3*(r.'*rdot)^2*rdot*r.'+(rdot.'*rdot)*(r.'*rdot)*r*r.'+2*(r.'*rdot)*(rdot.'*rdot)*r*r.')/nr^7 ...
    +105*(r.'*rdot)^3*r*r.'/nr^9 ...
    +(3*(rddot_H*r.'+r.'*rddot_H*eye(3)+r*rddot_H.')/nr^5 ...
    -15*(r*r.'*rddot_H*r.')/nr^7)*rdot*r.' ...
    +3*(6*r.'*rdot*rdot*rdot.'+rdot.'*rdot*r.'*rdot*eye(3)+3*rdot.'*rdot*r*rdot.'+2*r.'*rdot*rdot.'*rdot*eye(3))/nr^5 ...
    -15*((r.'*rdot)^3*eye(3)+3*(r.'*rdot)^2*r*rdot.')/nr^7 ...
    -(rddot_H.'*rdot*eye(3)+rddot_H*rdot.'+rdot*rddot_H.'+(r.'*rdot*eye(3)+rdot*r.'+r*rdot.')*A21)/nr^3 ...
    +3*(r.'*rddot_H*r.'*rdot*eye(3)+r.'*rdot*r*rddot_H.'+r.'*rdot*r*r.'*A21+r.'*rddot_H*r*rdot.')/nr^5 ...
    -((A21*rdot-2*Sw*rddot_H)*r.')/nr^3 ...
    -(r.'*(A21*rdot-2*Sw*rddot_H)*eye(3)+r*(A21*rdot-2*Sw*rddot_H).')/nr^3 ...
    +3*(r*r.'*(A21*rdot-2*Sw*rddot_H)*r.')/nr^5 ...
    +2*(3*(rdot*r.'+r.'*rdot*eye(3)+r*rdot.')/nr^5 -15*r*r.'*rdot*r.'/nr^7)*rddot_H*r.' ...
    -2*(rdot*rddot_H.'+rddot_H*rdot.'+rdot.'*rddot_H*eye(3))/nr^3 ...
    +6*((r.'*rdot)*(r.'*rddot_H)*eye(3)+r.'*rddot_H*r*rdot.'+r.'*rdot*r*rddot_H.')/nr^5 ...
    +(2*O_10-2*O_00*Sw)*A21;

O_31_H=O_20_H-2*O_21*Sw...
    -(2*rdot.'*rdot*eye(3)+4*rdot*rdot.')/nr^3 ...
    +2*(r'*rdot*eye(3)+rdot*r'+r*rdot.')*Sw/nr^3 ...
    +3*(2*(r.'*rdot)^2*eye(3)+4*(r.'*rdot)*rdot*r.'+4*(r.'*rdot)*r*rdot.'+2*(rdot.'*rdot)*r*r.')/nr^5 ...
    -6*(r'*rdot*r*r.')*Sw/nr^5 ...
    -15*2*(r.'*rdot)^2*r*r.'/nr^7 ...
    -2*(r.'*rddot_H*eye(3)+rddot_H*r.'+r*rddot_H.')/nr^3 ...
    +6*r.'*rddot_H*r*r.'/nr^5;

OO_H=[O_00 zeros(3,3);
    O_10 O_11;
    O_20_H O_21;
    O_30_H O_31_H];

O_H_svd=svd(OO_H);
end

function [O_svd O_criteria]=ObservabilityTest_TBP(xx)
global mu n a

r=xx(1:3);
rdot=xx(4:6);

nr=norm(r);

ra=r+[a 0 0]';
nra=norm(ra);

Sw=[0 -n 0;n 0 0; 0 0 0];

rddot=-2*Sw*rdot-Sw^2*ra-mu*ra/nra^3;

h=r/nr;
y=h;

O_00=1/nr*(eye(3)-y*y');
O_10=-1/nr^2*(rdot*y'+y'*rdot*eye(3)+y*rdot'-3*y*y'*rdot*y');
O_11=O_00;

tmp1=-(2*rdot*rdot.'+rdot.'*rdot*eye(3))/nr^3 ...
    +3*(2*rdot*r.'*rdot*r.'+r*rdot.'*rdot*r.'+r.'*rdot*r.'*rdot*eye(3)+2*r*r.'*rdot*rdot.')/nr^5 ...
    -15*(r*r.'*rdot*r.'*rdot*r.')/nr^7;
tmp3=-(rddot*r.'+r.'*rddot*eye(3)+r*rddot.')/nr^3 ...
    +3*r*r.'*rddot*r.'/nr^5 ...
    +(eye(3)/nr-r*r'/nr^3)*(-Sw^2-mu*eye(3)/nra^3+3*mu*ra*ra.'/nra^5);

O_20=tmp1+tmp3;
O_21=2*O_10-2*O_00*Sw;

OO=[O_00 zeros(3,3);
    O_10 O_11;
    O_20 O_21];

O_svd=svd(OO);

o1=rddot-nr*(-Sw^2-mu*eye(3)/nra^3+3*mu*ra*ra.'/nra^5)*y;
o2=rdot+Sw*r;

O_criteria(1)=norm(cross(r,rdot));
O_criteria(2)=norm(cross(r,o1));
O_criteria(3)=norm(cross(r,o2));
O_criteria(4)=norm(cross(O_20*r+O_21*rdot,O_21*r));
end

function [Wo Wo_svd Wo_cond]=ObservabilityGramian(t,xx0,eps,ode_options)

N=length(t);
delt=t(2)-t(1);

I6=eye(6);
yy_p=zeros(3*6,N);
yy_m=zeros(3*6,N);

Wo=zeros(6,6,N);
for ii=1:6
    ei=I6(:,ii);
    [t xx_TBP]=ode45(@eom_RelTBP,t,xx0+eps*ei,ode_options);
    for k=1:N
        yy_p(3*ii-2:3*ii,k)=xx_TBP(k,1:3)/norm(xx_TBP(k,1:3));
    end
    [t xx_TBP]=ode45(@eom_RelTBP,t,xx0-eps*ei,ode_options);
    for k=1:N
        yy_m(3*ii-2:3*ii,k)=xx_TBP(k,1:3)/norm(xx_TBP(k,1:3));
    end
end

diff_yy=yy_p-yy_m;
delWo=zeros(6,6);
Wo(:,:,1)=zeros(6,6);
for k=1:N-1
    for ii=1:6
        for jj=1:6
            delWo(ii,jj)=diff_yy(3*ii-2:3*ii,k)'*diff_yy(3*jj-2:3*jj,k);
        end
    end
    
    Wo(:,:,k+1)=Wo(:,:,k)+delWo;
end
Wo(:,:,1)=Wo(:,:,2);
Wo=Wo/4/eps^2*delt;

Wo_svd=zeros(6,N);
Wo_cond=zeros(1,N);
for k=1:N
    Wo_svd(:,k)=svd(Wo(:,:,k));
    Wo_cond(k)=max(Wo_svd(:,k))/min(Wo_svd(:,k));
end


end

function plot_results(filename)
close all;
load(filename);

%%% Plot Figures
close all;
figure(1);  % 3D trajectory
plot3(xx_HCW(1,:),xx_HCW(2,:),xx_HCW(3,:),'r');hold on;
plot3(xx_TBP(1,:),xx_TBP(2,:),xx_TBP(3,:),'b--');
axis equal;
grid on;

xlabel('$$x$$','interpreter','latex');
ylabel('$$y$$','interpreter','latex');
zlabel('$$z$$','interpreter','latex');

figure(2);  % x,y,z
plot(t/T,xx_HCW,'r',t/T,xx_TBP,'b--');
grid on;
xlabel('$$t/T$$','interpreter','latex');

figure(3);  % Min singular value
plot(t/T,min(O_H_svd),'r', t/T,min(O_svd),'b--');
xlabel('$$t/T$$','interpreter','latex');
set(gca,'yscale','log');
grid on;

figure(4);  % Observability criteria
plot(t/T,O_criteria);
set(gca,'yscale','log');
xlabel('$$t/T$$','interpreter','latex');
grid on;

figure(5);  % Condition number
plot(t/T,max(O_H_svd)./min(O_H_svd),'r', t/T,max(O_svd)./min(O_svd),'b--');
xlabel('$$t/T$$','interpreter','latex');
set(gca,'yscale','log');
grid on;

figure(6);  % Observability gramian: min singular value
plot(t/T,min(Wo_svd_HCW),'r',t/T,min(Wo_svd),'b--');
set(gca,'yscale','log');
xlabel('$$t/T$$','interpreter','latex');
grid on;

figure(7);  % Observability gramian: condition number
plot(t/T,Wo_cond_HCW,'r',t/T,Wo_cond,'b--');
set(gca,'yscale','log');
xlabel('$$t/T$$','interpreter','latex');
grid on;

% % Uncomment these to save figures
% figure(1);print([filename '_traj'],'-depsc2');
% figure(2);print([filename '_r'],'-depsc2');
% figure(3);print([filename '_svd'],'-depsc2');
% figure(4);print([filename '_cri'],'-depsc2');
% figure(5);print([filename '_cond'],'-depsc2');
% figure(6);print([filename '_W_svd'],'-depsc2');
% figure(7);print([filename '_W_cond'],'-depsc2');

end


function [Wo Wo_svd Wo_cond]=ObservabilityGramian_HCW(t,xx_HCW)
global n

N=length(t);
delt=t(2)-t(1);

A=[0 0 0 1 0 0;
    0 0 0 0 1 0;
    0 0 0 0 0 1;
    3*n^2 0 0 0 2*n 0;
    0 0 0 -2*n 0 0;
    0 0 -n^2 0 0 0];
C=[eye(3) zeros(3)];

Wo(:,:,1)=zeros(6,6);
for k=1:N-1
    r=xx_HCW(1:3,k);
    nr=norm(r);
    
    tmp=(eye(3)/nr-r*r'/nr^3)*C*expm(A*t(k));
    delWo=tmp'*tmp;
    Wo(:,:,k+1)=Wo(:,:,k)+delWo;
end
Wo(:,:,1)=Wo(:,:,2);
Wo=Wo*delt;

Wo_svd=zeros(6,N);
Wo_cond=zeros(1,N);
for k=1:N
    Wo_svd(:,k)=svd(Wo(:,:,k));
    Wo_cond(k)=max(Wo_svd(:,k))/min(Wo_svd(:,k));
end

end


