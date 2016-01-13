%% Frequency Response
% Author : Harmeet Singh
% Email : harmeet.esal@gmail.com
% Date : 11/05/2014

clear
close all
clc

%% Filter Specifications
Fs = 44100;
Ts = 1/Fs;
t = 0:Ts:1;
tmax = length(t);
omega = -Fs/2:Fs/(tmax-1):Fs/2;

H1_zeros = [0.98*exp(i*0.8*pi) 0.98*exp(i*-0.8*pi)];

H1_poles = [0.8*exp(i*0.4*pi) 0.8*exp(i*-0.4*pi)];

k = [1 2 3 4];
temp1 = 0.95*exp(i*(0.15*pi + 0.02*pi*k));
temp2 = 0.95*exp(-(i*(0.15*pi + 0.02*pi*k)));
H2_poles = [temp1 temp2];

H2_poles = [H2_poles H2_poles];


temp1 = 1./(0.95*exp(i*(0.15*pi + 0.02*pi*k)));
temp2 = 1./(0.95*exp(-i*(0.15*pi + 0.02*pi*k)));
H2_zeros = [temp1 temp2];

H2_zeros = [H2_zeros H2_zeros];


H_zeros = [H1_zeros H2_zeros];
H_poles = [H1_poles H2_poles];

% H_zeros = [0 0]
% H_poles = [0.2 0.5]

[n d] = zp2tf(H_zeros(:),H_poles(:),1);

%% Pole - Zero Plot
figure
zplane(H_zeros(:),H_poles(:));


norm_omega = 2*pi*Ts*omega;
[H] = freqz(n,d,norm_omega);

%% Display of Results

figure;subplot(2,2,1);plot(norm_omega,abs(H))
title('Frequency Response (Magnitude)');
rads = [-(pi), (-3*pi/4) , -(pi/(2)) , -((pi)/(4)) , 0 , ((pi/4)) , ((pi/2)) , (3*pi/4) , (pi)];
labels = {'-pi','-3pi/4','-pi/2','-pi/4','0','pi/4','pi/2','3pi/4','pi'};
set(gca,'XTick',rads);
set(gca,'XTickLabel',labels)


% figure
grp_delay = -diff(unwrap(angle(H)));

grp_delay = [0 grp_delay];
subplot(2,2,2);
plot(norm_omega,(grp_delay));
set(gca,'XTick',rads);
set(gca,'XTickLabel',labels)
title('Group Delay');


% figure;
subplot(2,2,3);
plot(norm_omega,angle(H));
set(gca,'XTick',rads);
set(gca,'XTickLabel',labels)
title('Phase of H')

% figure;
subplot(2,2,4);
plot(norm_omega,unwrap(angle(H)));
set(gca,'XTick',rads);
set(gca,'XTickLabel',labels)
title('unwrapped Phase of H')
