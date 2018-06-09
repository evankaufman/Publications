function EKF_TBP
% Extended Kalman filter for the nonlinear relative orbit
global mu n a
close all;

%%% Initial conditions and Parameters

% Fixed parameters
mu=398600;  % Earth gravitational parameter
RE=6378;    % Earth radius

% Orbital properties of chief
a=7100;         % Orbital radius of the chief
n=sqrt(mu/a^3); % Mean motion

% Simulation time
T=2*pi/n;     % Orbital period of the chief
tf=2*T;      % Terminal simulation time
N=2*6000+1;     % Total number of time steps
t=linspace(0,tf,N);  % Time vector

% Simulation title
filename='EKF_IOD_0A';
disp(['Case: ' filename]);

% Initial condition
xx0=[ -7.100116477545843e-02
    -9.913332955272308e-04
    1.307672689210581e-01
    -2.329550204657181e-06
    -9.913316565799100e-04
    1.307672203750978e-01];

% Initial estimate
xx0_hat=[ -7.100116477545843e-02
    -9.913332955272308e-04
    1.307672689210581e-01
    -2.329550204657181e-06
    -9.913316565799100e-04
    1.307672203750978e-01];

% Covariance matrices
P0=diag([50^2 50^2 50^2 (50*n)^2 (50*n)^2 (50*n)^2]);   % Initial covariance
Qk=10^0*diag([1e-8 1e-8 1e-8 1e-10 1e-10 1e-10]);    % Covariance of process noise
Rk=tand(2)^2*eye(3);    % Covariance of the measuremet

nu=(randn(N,3)*chol(Rk))'; % Generate measurement errors
mea_interval_steps=1;   % Numer of time steps between measurements

%%% EKF

% Initialize variables

y=zeros(3,N);y_hat=zeros(3,N);y_mea=zeros(3,N);
xx=zeros(6,N);xx_hat=zeros(6,N);
P=zeros(6,6,N);

k=1;
xx(:,k)=xx0;            % true state vector
xx_hat(:,k)=xx0_hat;    % estiamted state vector
P(:,:,k)=P0;            % covariance matrix
rk=xx(1:3,k);
y(:,k)=rk/norm(rk);     % true output
r_hat_k=xx_hat(1:3,k);
y_hat(:,k)=r_hat_k/norm(r_hat_k);   % estimated output
y_mea(:,k)=y(:,k)/norm(y(:,k))+nu(:,k+1);
y_mea(:,k)=y_mea(:,k)/norm(y_mea(:,k));

k_mea=[];   % array for time steps when the output is measured

disp('Filtering...');

% true solution
myoptions=odeset('RelTol',1e-6,'AbsTol',1e-6);
[t xx]=ode45(@eom_RelTBP,t,xx0,myoptions);
xx=xx';

xx_ddot=zeros(6,N);
for k=1:N
    xx_ddot(:,k)=eom_RelTBP(0,xx(:,k));
end

for k=1:N-1
    
    % true measurement
    rkp=xx(1:3,k+1);
    y(:,k+1)=rkp/norm(rkp);
    
    % flow update
    [~, xx_hat_tmp]=ode45(@eom_RelTBP,[t(k) t(k+1)],xx_hat(:,k),myoptions);
    xx_hat(:,k+1)=xx_hat_tmp(end,:)';
    r_hat_kp=xx_hat(1:3,k+1);
    y_hat(:,k+1)=r_hat_kp/norm(r_hat_kp);
    Ak=lin_eom_RelTBP(xx_hat(:,k),t(k+1)-t(k));
    xx_hat_nonlin = xx_hat(:,k+1);
    xx_hat_lin = Ak*xx_hat(:,k);
    change_nonlin = xx_hat(:,k+1)-xx_hat(:,k);
    change_lin = xx_hat_lin-xx_hat(:,k);
    P(:,:,k+1)=Ak*P(:,:,k)*Ak'+Qk;%*0;
    
    % measurement update
    if rem(k-1,mea_interval_steps)==0
        k_mea=[k_mea k+1];
        y_mea(:,k+1)=y(:,k+1)/norm(y(:,k+1))+nu(:,k+1);
        y_mea(:,k+1)=y_mea(:,k+1)/norm(y_mea(:,k+1));
        
        Ckp=[(eye(3)/norm(r_hat_kp)-r_hat_kp*r_hat_kp'/norm(r_hat_kp)^3) zeros(3,3)];
        Kkp=P(:,:,k+1)*Ckp'*inv(Ckp*P(:,:,k+1)*Ckp'+Rk);%*0;
        
        xx_hat(:,k+1)=xx_hat(:,k+1)+Kkp*(y_mea(:,k+1)-y_hat(:,k+1));
        P(:,:,k+1)=(eye(6)-Kkp*Ckp)*P(:,:,k+1);
        J(k) = trace(P(:,:,k+1));
    end
    
    if rem(k,(N-1)/100) == 0
        disp(['   ' num2str(k/(N-1)*100) '% done: ' ...
            num2str(norm(xx_hat(:,k+1))/norm(xx(:,k+1))) ' '...
            num2str(acos(min([1 xx_hat(:,k+1)'*xx(:,k+1)/norm(xx_hat(:,k+1))/norm(xx(:,k+1))])))]);
    end
            
end
plot(t(2:end),J);

%%% Post-processing
disp('Post-processing...');

r=xx(1:3,:);
v=xx(4:6,:);

xx_norm=zeros(1,N);xx_dir=zeros(6,N);
xx_hat_norm=zeros(1,N);xx_hat_dir=zeros(6,N);
eigP=zeros(1,N);v_dir=zeros(3,N);v_hat_dir=zeros(3,N);
err_dir=zeros(1,N);err_mag=zeros(1,N);
err_mag_r=zeros(1,N);err_mag_v=zeros(1,N);
err_residual=zeros(1,N);

for k=1:N
    xx_norm(k)=norm(xx(:,k));           % true magnitude of state vector
    xx_dir(:,k)=xx(:,k)/xx_norm(k);     % true direction of state vector
    
    xx_hat_norm(k)=norm(xx_hat(:,k));               % estimated magnitude of state vector
    xx_hat_dir(:,k)=xx_hat(:,k)/xx_hat_norm(k);     % estimated direction of state vector
    
    eigP(k)=max(eig(P(:,:,k)));                     % max eigenvalue of covariance matrix
    
    v_dir(:,k)=xx(4:6,k)/norm(xx(4:6,k));               % true direction of relative velocity
    v_hat_dir(:,k)=xx_hat(4:6,k)/norm(xx_hat(4:6,k));   % estimated direction of relative velocity
    
    err_dir(k)=acos(min([1, xx_dir(:,k)'*xx_hat_dir(:,k)]));    % direction error in rad
    err_mag(k)=abs((xx_hat_norm(k))/xx_norm(k));                % magnitude error
    err_mag_r(k)=abs(norm(xx_hat(1:3,k)))/norm(xx(1:3,k));      % magnitude error of relative position
    err_mag_v(k)=abs(norm(xx_hat(4:6,k)))/norm(xx(4:6,k));      % magnitude error of relative velocity
    
    err_residual(k)=acos(min([1 y_hat(:,k)'*y_mea(:,k)]));      % redisual error
end

rms_err_dir=sqrt(sum(err_dir.*err_dir)/size(err_dir,2));
rms_err_mag=sqrt(sum((1-err_mag).*(1-err_mag))/size(err_mag,2));
rms_err_mag_r=sqrt(sum((1-err_mag_r).*(1-err_mag_r))/size(err_mag_r,2));
rms_err_mag_v=sqrt(sum((1-err_mag_v).*(1-err_mag_v))/size(err_mag_v,2));

fprintf('%.4f, %.4f, %.4f',rms_err_dir,rms_err_mag_r,rms_err_mag_v);

save(filename);             % save results
plot_results(filename);     % plot results
evalin('base',['load ' filename ';']);  % load results to workspace

end

function Ak=lin_eom_RelTBP(X,delt)
% Linearized equations of motion for the relative motion
global mu a n
r=X(1:3);
ra=r+[a 0 0]';
dis=norm(ra);

Sw=[0 -n 0;
    n 0 0;
    0 0 0];
A=[zeros(3) eye(3);
    -Sw^2-mu*eye(3)/dis^3+3*mu*ra*ra'/dis^5 -2*Sw];

Ak=expm(A*delt);
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

function plot_results(filename)
load(filename);

figure(1);
plot(t/T,xx_norm,'r',t/T,xx_hat_norm,'b--');hold on;
xlabel('$$t/T$$','interpreter','latex');
ylabel('$$\|\mathbf{x}\|$$','interpreter','latex');
grid on;

figure(2);
subplot(2,1,1);
plot(t/T,err_mag,'b--');hold on;
ylabel('$$e_{mag}$$','interpreter','latex');
grid on;
subplot(2,1,2);hold on;
plot(t/T,err_dir,'b--');
ylabel('$$e_{dir}$$','interpreter','latex');
xlabel('$$t/T$$','interpreter','latex');
grid on;

figure(3);
plot(t/T,eigP,'b--');hold on;
ylabel('$$\lambda_{\max}[P]$$','interpreter','latex');
xlabel('$$t/T$$','interpreter','latex');
grid on;

figure(4);
plot(t(k_mea)/T,err_residual(k_mea),'b--');hold on;
ylabel('$$e_{residual}$$','interpreter','latex');
xlabel('$$t/T$$','interpreter','latex');
grid on;

% % Uncomment these to save figures
% figure(1);print([filename '_xx_norm'],'-depsc2');
% figure(2);print([filename '_e'],'-depsc2');
% figure(3);print([filename '_eigP'],'-depsc2');
% figure(4);print([filename '_eres'],'-depsc2');

end

