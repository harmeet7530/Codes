%% Optimum Solution of Performance Surface
% Author : Harmeet Singh
% Email : harmeet.esal@gmail.com
% Date : 02/03/2015

clear all
close all
N = 5;
shift = 60;

%Weight Vectors
w0 = -3:0.5:4;
w1 = -6:0.5:1.5;
[x,y] = meshgrid(w0,w1);

%Equation fo Performance Surface is defined
ps = 0.5*(x.^2 + y.^2) + x.*y.*cos(2*pi/N) + 2*y.*sin(2*pi/N) + 2;

%Calculation of Optimal Weight Values
R = [0.5,0.5*cos(2*pi/5);0.5*cos(2*pi/5),0.5];
P = [0;-sin(2*pi/5)];
Wopt = R\P;

%Calculation of eigen vectors for R
[a,b] = eig(R);

figure(1);
surf(x,y,ps);
title('3-D Surface Plot');
xlabel('Weight 0');
ylabel('Weight 1');
zlabel('Performance Surface');
view(shift,10);

figure(2);
contour(x,y,ps);
title('Contour Plot');
xlabel('Weight 0');
ylabel('Weight 1');

figure(3);
quiver(Wopt(1),Wopt(2),a(1,1),a(2,1),3,'g');
hold on;
quiver(Wopt(1),Wopt(2),-a(1,1),-a(2,1),3,'b');
quiver(Wopt(1),Wopt(2),a(2,1),a(2,1),3,'y');
quiver(Wopt(1),Wopt(2),-a(2,1),-a(2,1),3,'r');
title('V prime Axis');
xlabel('Weight 0');
ylabel('Weight 1');

figure(4);
contour(x,y,ps);
hold on;
quiver(Wopt(1),Wopt(2),a(1,1),a(2,1),3,'g');
quiver(Wopt(1),Wopt(2),-a(1,1),-a(2,1),3,'b');
quiver(Wopt(1),Wopt(2),a(2,1),a(2,1),3,'y');
quiver(Wopt(1),Wopt(2),-a(2,1),-a(2,1),3,'r');
title('Superimpose V prime axis on Contour Plot');
xlabel('Weight 0');
ylabel('Weight 1');

