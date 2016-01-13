%% Fir Filter Design
% Author : Harmeet Singh
% Email : harmeet.esal@gmail.com
% Date : 11/25/2014

clear
close all
clc

%% Filter Specifications

Fs = 44100;
Ts = 1/Fs;
t = 0:Ts:1;
tmax = length(t);
omega = 0:Fs/(tmax-1):2*Fs;

H1_zeros = [0.9*exp(i*0.6*pi) 0.9*exp(i*-0.6*pi) 1.25*exp(i*0.8*pi) 1.25*exp(i*-0.8*pi)];
Hm_zeros = [0.9*exp(i*0.6*pi) 0.9*exp(i*-0.6*pi) 0.8*exp(i*0.8*pi) 0.8*exp(i*-0.8*pi)];
Ha_zeros = [1/(0.8*exp(i*0.8*pi)) 1/(0.8*exp(i*-0.8*pi))];
Ha_poles = [0.8*exp(i*0.8*pi) 0.8*exp(i*-0.8*pi)];
H1_poles = [0 0 0 0];

[n d] = zp2tf(H1_zeros(:),H1_poles(:),1);
[a b] = zp2tf(Hm_zeros(:),H1_poles(:),1);
[c d] = zp2tf(Ha_zeros(:),Ha_poles(:),1);

%% Zero - Pole Plot

figure
zplane(H1_zeros(:),H1_poles(:));
figure
zplane(Hm_zeros(:),H1_poles(:));
figure
zplane(Ha_zeros(:),Ha_poles(:));

norm_omega = 2*pi*Ts*omega;

[H] = freqz(n,d,norm_omega);
[Hm] = freqz(a,b,norm_omega);
[Ha] = freqz(c,d,norm_omega);

%% Display Results
figure;
subplot(2,2,1);
plot(norm_omega,mag2db(abs(H))),'r';hold on;
plot(norm_omega, angle(H),'m');
title('Magnitude(db) & Phase Response');
%rads = [0 , ((pi/4)) , ((pi/2)) , (3*pi/4) , (pi), (5*(pi/4)), (3*(pi/2)), (7*(pi/4)),(2*pi)];
%labels = {'0','pi/4','pi/2','3pi/4','pi', '5pi/4', '3pi/2', '7pi/4','2pi' };
%set(gca,'XTick',rads);
%set(gca,'XTickLabel',labels)

grp_delay = -diff(unwrap(angle(H)));
grp_delay = [0 grp_delay];
subplot(2,2,2);
plot(norm_omega,(grp_delay));
%set(gca,'XTick',rads);
%set(gca,'XTickLabel',labels)
title('Group Delay');

subplot(2,2,3);
plot(norm_omega,abs(H)),'r';hold on;
plot(norm_omega, angle(H),'m');
title('Magnitude(without Db) & Phase Response');
%set(gca,'XTick',rads);
%set(gca,'XTickLabel',labels)

figure;
subplot(2,2,1);
plot(norm_omega,mag2db(abs(Hm))),'r';hold on;
plot(norm_omega, angle(Hm),'m');
title('Magnitude(db) & Phase Response');
%set(gca,'XTick',rads);
%set(gca,'XTickLabel',labels)

g_d1 = -diff(unwrap(angle(Hm)));
g_d1 = [0 g_d1];
subplot(2,2,2);
plot(norm_omega,(g_d1));
%set(gca,'XTick',rads);
%set(gca,'XTickLabel',labels)
title('Group Delay');

subplot(2,2,3);
plot(norm_omega,abs(Hm)),'r';hold on;
plot(norm_omega, angle(Hm),'m');
title('Magnitude(without Db) & Phase Response');
%set(gca,'XTick',rads);
%set(gca,'XTickLabel',labels)

figure;
subplot(2,2,1);
plot(norm_omega,mag2db(abs(Ha))),'r';hold on;
plot(norm_omega, angle(Ha),'m');
title('Magnitude(db) & Phase Response');
%set(gca,'XTick',rads);
%set(gca,'XTickLabel',labels)

g_d2 = -diff(unwrap(angle(Ha)));
g_d2 = [0 g_d2];
subplot(2,2,2);
plot(norm_omega,(g_d2));
%set(gca,'XTick',rads);
%set(gca,'XTickLabel',labels)
title('Group Delay');

subplot(2,2,3);
plot(norm_omega,abs(Ha)),'r';hold on;
plot(norm_omega, angle(Ha),'m');
title('Magnitude(without Db) & Phase Response');
%set(gca,'XTick',rads);
%set(gca,'XTickLabel',labels)



