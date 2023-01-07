clc
close all;
clear;

%constants
R_earth=6.378E3;
mu=3.986E5;
altitude = 150:1:1000;

r = altitude +R_earth; % define radius

theta_dot = sqrt(mu.*(r.^-3));  % calculate mean motion

T = (2*pi).*(theta_dot.^-1); %calculate periods
v = r.*theta_dot; %calculate velocity


%% Plot Figures
figure(1)
plot(altitude,v)
title("Orbit Altitude vs Speed (Circular Orbit)")
xlabel("Altitude(km)")
ylabel("Inertial speed(km/s)")
print('C:\Users\boaza\OneDrive - Carleton University\AERO 3240\output files\2.5\Orbit_Altitude_vs_Velocity','-depsc')

figure(2)
plot(altitude,T)
title("Orbit Altitude vs period (Circular Orbit)")
xlabel("Altitude(km)")
ylabel("Period(s)")
print('C:\Users\boaza\OneDrive - Carleton University\AERO 3240\output files\2.5\Orbit_Altitude_vs_period','-depsc')
