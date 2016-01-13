%% Newtons Method and Steepest Descent
% Author : Harmeet Singh
% Email : harmeet.esal@gmail.com
% Date : 02/10/2015

clear all
close all

N = 10;
shift = 60;

k=1:1:N;

xk = sin((2*pi*k)/N);

%dk = 2*cos((2*pi*k)/N);

%Weight Vectors
w0 = -3:0.5:4;
w1 = -4:0.5:0;
[x,y] = meshgrid(w0,w1);
W = [x,y];
u =1;

%Equation fo Performance Surface is defined
ps = 0.5*(x.^2 + y.^2) + x.*y.*cos(2*pi/N) + 2*y.*sin(2*pi/N) + 2;

[Gx,Gy] = gradient(ps,0.5,0.5);
G = [Gx,Gy];
l1 = size(Gx,1);
L = l1-1;
%Defining Input Correlation Matrix

Z=xcorr(xk,L);
Z(1:L)=[];
b= zeros(numel(Z));
b(1,:) = Z(1,:);
r1 = zeros(numel(Z));
for i = 1:numel(Z)
     for j = 1:numel(Z)-1
       r1(j+1,i+j)=Z(1,j); 
        
     end
     r2 = r1(1:numel(Z),1:numel(Z));
end
R = r2+r2'+b+b';
I = eye(sqrt(numel(R)));

%Calculation of Optimal Weight Values
Wopt = u*(R\G);
Wnew(1:L+1,1:L+1) = 0;

%Newton's Method
for k = 1:size(W,1)
    for l = 1:size(W,2)
        Wnew(k+1,l+1) = W(k,l) - Wopt(k,l);
    end
end
%Steepest Descent Method
%for k = 2:size(W,1)
   % for l = 2:size(W,2)
        Wnew(k,l) = W(k-1,l-1)+u*(-G(k-1,l-1));
  %  end
%end


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
