function [r_vec, v_vec]=oe2rv(h,e,theta,Omega,omega,i) 
xhat=[1;0;0]; 
yhat=[0;1;0]; 
zhat=[0;0;1]; 
mu=398600;

Nhat=cos(Omega)*xhat+sin(Omega)*yhat;
hhat=sin(i)*sin(Omega)*xhat-sin(i)*cos(Omega)*yhat+cos(i)*zhat;
Nthat=-sin(Omega)*cos(i)*xhat+cos(Omega)*cos(i)*yhat+sin(i)*zhat;

urhat=cos(theta+omega)*Nhat+sin(theta+omega)*Nthat;
uthat=-sin(theta+omega)*Nhat+cos(theta+omega)*Nthat;

r=h^2/mu*1/(1+e*cos(theta));
E=-1/2*mu^2/h^2*(1-e^2);
v=sqrt(2*(E+mu/r));
gam=atan2(e*sin(theta),1+e*cos(theta));

r_vec=r*urhat;
v_vec=v*cos(gam)*uthat+v*sin(gam)*urhat;

