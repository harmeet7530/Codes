%% Park-McClellan Filter Design Algorithm
% Author : Harmeet Singh
% Email : harmeet.esal@gmail.com
% Date : 12/08/2014

clear
close all
clc

%% Filter Specifications

Wp = 0.35*pi;
Ws = 0.45*pi;
d1 = 0.034;
d2 = 0.045;
dW = Ws - Wp;

% Order of Filter %
M = ((-20*log10((d1*d2)^(1/2))-13)/(2.285*dW)) + 1;
% Length of Filter %
N = M+1;
L1 = M/2;
L = fix(L1);

%% Finding optimum values iteratively 
% 1st Iteration %
% w = [0.257, 0.428, 0.685, 0.842, 1.099, 1.413, 1.702, 1.99, 2.278, 2.566, 2.854, 2.95 , pi];
% 2nd Iteration %
% w = [0.1 0.35 0.6 0.85 1.099 1.413 1.65 1.9 2.15 2.4 2.65 2.9 pi ];
% 3rd Iteration %
% w = [0.1 0.3624 0.7862 0.9888 1.099 1.413 1.566 1.861 2.15 2.43 2.739
% 3.145 pi ];
% 4th Iteration %
% w = [0.1 0.5098 0.7984 1.013 1.099 1.413 1.529 1.843 2.125 2.653 3.028 3.28 pi ]
% 5th Iteration %
% w = [0.1 0.2764 0.5221 0.823 1.099 1.413 1.529 1.806 2.451 2.831 3.145 3.452 pi ]
% 6th Iteration %
% w = [0.1 0.3194 0.6818 0.823 1.099 1.413 1.529 1.775 2.107 2.5 2.85 3.14 pi ]
% 7th Iteration %
% w = [0.1 0.3194 0.6818 0.9581 1.099 1.413 1.529 1.775 2.033 2.371 2.715 3.157 pi ]
% 8th Iteration %
% w = [0.1 0.3194 0.6818 0.9581 1.099 1.413 1.529 1.775 2.033 2.266 2.575 2.91 pi ]
% 9th Iteration %
w = [0.1 0.3194 0.6818 0.90 1.099 1.413 1.529 1.775 2.033 2.266 2.575 2.91 pi ];
H = [0.966; 1.034; 0.966; 1.034; 0.966; 0.045; -0.045; 0.045; -0.045; 0.045; -0.045; 0.045; -0.045];

%% Defining Chebyshev Polynomial
C = zeros(L+1,L+1);
for i = 1:L+1
    for j = 1:L+1
        C(i,j)= 2*cos((j - 1)*w(i));
    end
    C(i,1) = 1;
end

h = (inv(C))*H(1:L+1);
h = h';
h =[fliplr(h(2:end)) h ];

%% Fourier Transform
hn = fft(h,1024);

f = linspace(0,2*pi,1024);

%% Display the results
figure(1);
plot(f,abs(hn));
grid on;
rads = [0 , (0.2*(pi)) , (0.35*(pi)), (0.45*(pi)), (0.6*(pi)) , (0.8*(pi)) , (pi)];
labels = {'0', '0.2pi', '0.35pi', '0.45pi', '0.6pi', '0.8pi', 'pi'};
xlabel ('Normalized Frequency');
ylabel ('H(n)');
set(gca,'XTick',rads);
set(gca,'XTickLabel',labels);
title('Frequency Response');

figure(2);
stem(hn);
