%% LMS Algorithm for FIR Filters
% Author : Harmeet Singh
% Email : harmeet.esal@gmail.com
% Date : 04/13/2015

clear all
close all

%% Generating Input Signal

N = 1000;
k = 0:1:N;
sk = sin(2*pi*k/(30-k/50));
rk = sqrt(12)*randn(1,N);
zk = sk+0.1*rk;
zk = zk';

%% Passing the signal through the channel

channel = [1 -1 0.89];
xk = filter(channel,1,zk);
xk = xk(1:N);

%% Delaying the signal by 5 samples

d = padarray(zk,5,'pre');
d = d';
d = d(1:N);

%% Obtaining Input Correlation Matrix

L = 31;
r = xcorr(xk,L,'biased');
r(1:L) = [];
r = r';
R = toeplitz(r,[r(1) conj(r(2:end))]);

umax = 1/max(eig(R)); %Finding the step-size by obtaining maximum eigen values of R
u = 0.04*umax;

%% LMS Algorithm

t = 2; 
w = zeros(t,1);
W = [];
J = [];
E = [];

for n = (t+1):1:N;
    x = xk(n-1:-1:n-t);
    y = w'*x; %Output signal
    e = y - d(n); %Error
    w = w - 2*u*e*x; %Weight Update
    W = [W w]; %Weights
    E = [E e]; %Error
    J = [J;e^2]; %Mean Square Error
end
E = E';
f = ones(32,1)/32;
Jk = filter(f,1,J);

figure(1);
plot(J);
title('Mean Square Error');
xlabel('Iterations');
ylabel('{E(n)}^2');
grid on;

figure(2);
plot(Jk);
title('Learning Curve');
xlabel('Iterations');
ylabel('Error');
grid on;



