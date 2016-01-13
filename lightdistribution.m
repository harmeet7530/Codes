%% Propagation of Light in Free Space
% Author : Harmeet Singh
% Email : harmeet.esal@gmail.com
% Date : 03/15/2015

clear all
close all

%%Define wavelength, number of samples and length
lambda = 1;
N = 500;
L = 40*lambda;

%%PART A (a)
%%Aperture with width D1 and D2
d1 = 1.2*lambda;
D1 = round((round(d1/L*(2*N+1)))/2);

d2 = 0.2*lambda;
D2 = round((round(d2/L*(2*N+1)))/2);

%%Predefining input plane wave propagating in z-direction with aperture D1
%%and D2
Gz0 = [];
G0 = zeros(1,2*N+1);
G0(N+1-D1:N+1+D1) = 1;

Gz1 = [];
G1 = zeros(1,2*N+1);
G1(N+1-D2:N+1+D2) = 1;

%%Fourier Transform of input plane wave
fG0 = fftshift(fft(G0,2*N+1));
fG1 = fftshift(fft(G1,2*N+1));

%%Defining x, vx and z
x = linspace((-L/2),(L/2),(2*N+1));
vx = (-N:N)/L;
z = linspace(0,20*lambda,501);

%%Using the transfer function to multiply by the input plane wave in
%%frequency domain and further obtain the inverse fourier transform of the
%%resulting wave
for i=1:length(z)
        H = exp(-2*pi*z(i)*sqrt(vx.^2-(1/lambda).^2));
        fUG0 = fG0.*H;
        fUG1 = fG1.*H;
        Gz0 = [Gz0;ifft(ifftshift(fUG0))];
        Gz1 = [Gz1;ifft(ifftshift(fUG1))];
end
Gz0 = Gz0';
Gz1 = Gz1';

%% PART A (b)

%%Constants a1 and a2 used in Transparency
a1 = 2*lambda;
a2 = 0.2*lambda;

%%Transparency described for a1 and a2
t1 = 1 + cos(((2*pi)/a1)*x);
t2 = 1 + cos(((2*pi)/a2)*x);

%%Fourier transform of transparency
ft1 = fftshift(fft(t1,2*N+1)); 
ft2 = fftshift(fft(t2,2*N+1));

Gz2 = [];
Gz3 = [];

%%Using the transfer function to multiply by the transparency in
%%frequency domain and further obtain the inverse fourier transform of the
%%resulting wave
for i=1:length(z)
        H = exp(-2*pi*z(i)*sqrt(vx.^2-(1/lambda).^2));
        fGz2 = ft1.*H;
        fGz3 = ft2.*H;
        Gz2 = [Gz2;ifftshift(ifft(fGz2))];
        Gz3 = [Gz3;ifftshift(ifft(fGz3))];
end
Gz2 = Gz2';
Gz3 = Gz3';

%%PART B (a)
%%Defining constants to be used
d = 10*lambda;
a3 = 3*lambda;
a4 = 10*lambda;
U1 = [];
Uz1 = [];
U2 = [];
Uz2 = [];

%%Finding the absolute value of x smaller than or equal to a
r1 = (abs(x)<=a3);
r2 = (abs(x)<=a4);

%%Defining the input wave to be used for a3 and a4
U1 = (exp(j*2*pi*sqrt(x.^2+d.^2)/lambda)./sqrt(x.^2+d.^2)).*r1;   
U2 = (exp(j*2*pi*sqrt(x.^2+d.^2)/lambda)./sqrt(x.^2+d.^2)).*r2; 

%%Fourier transform of the input waves obtained
fU1 = fftshift(fft(U1,2*N+1));
fU2 = fftshift(fft(U2,2*N+1));

%%Using the transfer function to multiply by the input wave in
%%frequency domain and further obtain the inverse fourier transform of the
%%resulting wave
for i=1:length(z)
        H = exp(-2*pi*z(i)*sqrt(vx.^2-(1/lambda).^2));
        fUz1 = fU1.*H;
        fUz2 = fU2.*H;
        Uz1 = [Uz1;ifft(ifftshift(fUz1))];
        Uz2 = [Uz2;ifft(ifftshift(fUz2))];
end
Uz1 = Uz1';
Uz2 = Uz2';

%%Finding the focused spot size for a3 and a4 for comparing the results
I1 = abs(Uz1.*Uz1);
I1a = I1(:,round(501/2));
IIa = I1a/max(I1a);

I2 = abs(Uz2.*Uz2);
I2a = I2(:,round(501/2));
IIb = I2a/max(I2a);

%%PART B (b)
%%Predefining the input wave and refractive medium, n
U3 = [];
Uz3 = [];
n = 2;

%%Input wave
U3 = (exp(j*2*pi*sqrt(x.^2+d.^2)/(lambda/n))./sqrt(x.^2+d.^2)).*r2;

%%Fourier tranform of input wave
fU3 = fftshift(fft(U3,2*N+1));

%%Using the transfer function to multiply by the input wave in
%%frequency domain and further obtain the inverse fourier transform of the
%%resulting wave
for i=1:length(z)
        H = exp(-2*pi*z(i)*sqrt(vx.^2-(1/lambda).^2));
        fUz3 = fU3.*H;
        Uz3 = [Uz3;ifft(ifftshift(fUz3))];       
end
Uz3 = Uz3';

%%PLOTTING OF ALL THE RESULTS OBTAINED

nu=1/(x(end)-x(1))*(floor(-((2*N+1)-1)/2):floor(((2*N+1)-1)/2));

figure(1);
subplot(1,2,1);
imagesc(z,x,abs(Gz0));
colormap hot;
axis equal tight xy;
title (['Wave Distribution with width D = ', num2str(d1)]), xlabel 'z(\lambda)', ylabel 'x(\lambda)';
subplot(1,2,2);
imagesc(z,x,abs(Gz1));
colormap hot;
axis equal tight xy;
title (['Wave Distribution with width D = ', num2str(d2)]), xlabel 'z(\lambda)', ylabel 'x(\lambda)';

figure(2);
subplot(2,1,1);
plot(nu,abs(fG0));
title (['Amplitude of Wave distribution with width D = ', num2str(d1)]), xlabel 'x', ylabel 'Magnitude';
subplot(2,1,2);
plot(nu,abs(fG1));
title (['Amplitude of Wave distribution with width D = ', num2str(d2)]), xlabel 'x', ylabel 'Magnitude';

figure(3);
subplot(1,2,1);
imagesc(z,x,abs(Gz2));
colormap hot;
axis equal tight xy;
title (['Wave distribution with Transparency when a = ', num2str(a1)]), xlabel 'z(\lambda)', ylabel 'x(\lambda)';
subplot(1,2,2);
imagesc(z,x,abs(Gz3));
colormap hot;
axis equal tight xy;
title (['Wave distribution with Transparency when a = ', num2str(a2)]), xlabel 'z(\lambda)', ylabel 'x(\lambda)';

figure(4);
subplot(2,1,1);
plot(nu,abs(ft1));
title (['Amplitude of Transparency when a = ', num2str(a1)]), xlabel 'x', ylabel 'Magnitude';
subplot(2,1,2);
plot(nu,abs(ft2));
title (['Amplitude of Transparency when a = ', num2str(a2)]), xlabel 'x', ylabel 'Magnitude';

figure(5);
subplot(1,2,1);
imagesc(z,x,abs(Uz1));
colormap hot;
axis equal tight xy;
title (['Wave distribution when a = ', num2str(a3)]), xlabel 'z(\lambda)', ylabel 'x(\lambda)';
subplot(1,2,2);
imagesc(z,x,abs(Uz2));
colormap hot;
axis equal tight xy;
title (['Wave distribution when a = ', num2str(a4)]), xlabel 'z(\lambda)', ylabel 'x(\lambda)';

figure(6);
subplot(2,1,1);
plot(nu,abs(fU1));
title (['Amplitude of Wave Distribution when a = ', num2str(a3)]), xlabel 'x', ylabel 'Magnitude';
subplot(2,1,2);
plot(nu,abs(fU2));
title (['Amplitude of Wave Distribution when a = ', num2str(a4)]), xlabel 'x', ylabel 'Magnitude';

figure(7);
subplot(2,1,1);
plot(nu,sign(IIa-0.5));
title (['Size of focused spot size when a = ', num2str(a3)]), xlabel 'x', ylabel 'Magnitude';
subplot(2,1,2);
plot(nu,sign(IIb-0.5));
title (['Size of focused spot size when a = ', num2str(a4)]), xlabel 'x', ylabel 'Magnitude';

figure(8);
imagesc(z,x,abs(Uz3));
colormap hot;
axis equal tight xy;
title (['Wave distribution with Refractive Medium when a = ', num2str(a4)]), xlabel 'z(\lambda)', ylabel 'x(\lambda)';

figure(9);
plot(nu,abs(fU3));
title (['Amplitude of Wave Distribution with Refractive Medium when a = ', num2str(a4)]), xlabel 'x', ylabel 'Magnitude';