clc
clear all
close all

x = 0:30;
xnew = 130*x/30+1300;
y1 = .3*x.^2 -.2*x + 80;
y2 = -.01*x.^3 + .4*x.^2 -10;
ytop = 300*ones(size(x));
ybottom = zeros(size(x));
line = ((220-20)/30)*x + 20;


shadedplot(xnew,y1,y2,'g')
hold on
shadedplot(xnew,ytop,y1,'r')
hold on
shadedplot(xnew,y2,ybottom,'r')
hold on
plot(xnew,line,'k','LineWidth',3)
plot(xnew,y1,'b','LineWidth',2)
plot(xnew,y2,'b','LineWidth',2)
axis([1300 1430 0 220])
text(1320,180,'Low Efficiency','FontSize',16)
text(1380,30,'Low Efficiency','FontSize',16)
text(1360,105,'Line of High Efficiency','FontSize',16)
xlabel('Servo Command','FontSize',16)
ylabel('BLDC Motor Command','FontSize',16)