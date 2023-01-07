clc
clear all;
close all;

% save folder
mydir  = pwd;
idcs   = strfind(string(mydir),'\');
newdir = mydir(1:idcs(end)-1);
save_to = strcat(newdir,'\output files\2.13\');
% Add general Functions
addpath(strcat(newdir,'\general functions'));
clear mydir idcs newdir;
% Gravitational parameter
mu = 398600.4418; % km^3/s^2

% Earth radius 
rE = 6378; % km

% Conversion constants
d2r = pi/180; % rad/deg
r2d = 1/d2r; % deg/rad
 
% Orbital elements
a       = rE + 791; % km
e       = 0.001;
i       = 98.28; % deg
w       = 0; % deg
RAAN    = 270; % deg
tp      = 0; % sec

p = a*(1-e^2);

% Eccentric anomaly at t = 0 sec
eano   = 0; % deg

% True anomaly at t = 0 sec
tano = 0; % deg

% Magnitude of position vector at t = 0 sec
r = p/(1+e*cosd(tano)); % km

% Components of position and velocity vectors in perifocal at t = 0 sec $\phantomsection\label{line: 2.13 part B}$
r_P_ini = r*[1 0 0]';

v_P_ini = [0 sqrt(mu/p)*(e+1) 0]';

% Rotation matrix from perifocal to ECI, i.e., C_IP

C_IP = (C_3(w)*C_1(i)*C_3(RAAN))'; % using functions I made  for rot matrices

% Components of position and velocity vectors in ECI at t = 0 sec
r_I_ini = C_IP*r_P_ini;
v_I_ini = C_IP*v_P_ini;

T = 2*pi*sqrt(a^3/mu);

% Simulation from t = 0 to t = T 
open_system('PROBA2mdl.slx')
set_param('PROBA2mdl', 'StopTime', 'T')
disp('Running Simulation...')
sim('PROBA2mdl')

%% Part A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot components of the position vector in Perifocal $\phantomsection\label{line: 2.13 part A}$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% r_P through conversion of simulation variables to perifocal
r_P =C_IP'*r_I';
% find r_p through equation

% Orbit equation definition with 2d vector componenets
syms r_P_func(theta)
r_P_func(theta) = matlabFunction(p/(1+e*cos(theta))*[cos(theta);sin(theta)]);

%calculation of orbit through orbit equation
r_P_eq = r_P_func(0:0.01:2*pi);
r_P_eq = double(cell2sym(r_P_eq)); % convert to symbolic to number

figure
plot(r_P(1,:),r_P(2,:),'k','LineWidth',2);
title('PROBA 2 Orbital Path in Perifocal Reference Frame')
xlabel('rx_P (km)')
ylabel('ry_P (km)')
hold on
plot(r_P_eq(1,:),r_P_eq(2,:),'r','LineWidth',1);
legend('simulation','orbit equation','Location','best')
pbaspect([1 1 1])
daspect([1 1 1])
hold off

print(strcat(save_to,'perifocal_plot.eps'),'-depsc')

%% Part B
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tabulize r_I_ini and v_I_ini
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varNames = ["Component", "r_I_inital", "v_I_inital"];
varTypes = ["string","double","double"];
sz = [3 3];
initial_conditions = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
initial_conditions.r_I_inital = r_I_ini;
initial_conditions.v_I_inital = v_I_ini;
initial_conditions.Component = ["X";"Y";"Z"];
initial_conditions.Properties.VariableNames = ["Component","r_{I,inital}(km)", "v_{I,inital}(km/s)"];
writetable(initial_conditions,strcat(save_to,'initial_conditions.csv'))



%% Part C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot components of the position vector in ECI in 3D $\phantomsection\label{line: 2.13 part C}$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot3(r_I(:,1),r_I(:,2),r_I(:,3),'k','linewidth', 2)
grid on

hold on
% Adding Earth
% Using Will Campbell (2022). Earth-sized 
% Sphere with Topography 
% (https://www.mathworks.com/matlabcentral
% /fileexchange/27123-earth-sized-sphere-with-topography), 
% MATLAB Central File Exchange. Retrieved November 28, 2022.
earth_sphere(gca,'km') 

set(gca,'FontSize',9,'FontName', 'Times')
title('PROBA 2 Orbital Path in ECI Reference Frame')
xlabel('rx_I (km)')
ylabel('ry_I (km)')
zlabel('rz_I (km)')
axis equal




print(strcat(save_to,'3d_plot.png'),'-dpng','-r600')
%% part D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot components of the position vector as function of time $\phantomsection\label{line: 2.13 part D}$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(3,1,1)
plot(t,r_I(:,1),'k');
title('PROBA 2 Position Vector Components in ECI Reference Frame')
xlabel('Time (sec)')
ylabel('rx_I (km)')
xlim([0 max(t)])

subplot(3,1,2)
plot(t,r_I(:,2),'k');
xlabel('Time (sec)')
ylabel('ry_I (km)')
xlim([0 max(t)])

subplot(3,1,3)
plot(t,r_I(:,3),'k');
xlabel('Time (sec)')
ylabel('rz_I (km)')
xlim([0 max(t)])


print(strcat(save_to,'position_components.eps'),'-depsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot components of the velocity vector as function of time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(3,1,1)
plot(t,v_I(:,1),'k');
title('PROBA 2 Velocity Vector Components in ECI Reference Frame')
xlabel('Time (sec)')
ylabel('vx_I (km/s)')
xlim([0 max(t)])

subplot(3,1,2)
plot(t,v_I(:,2),'k');
xlabel('Time (sec)')
ylabel('vy_I (km/s)')
xlim([0 max(t)])

subplot(3,1,3)
plot(t,v_I(:,3),'k');
xlabel('Time (sec)')
ylabel('vz_I (km/s)')
xlim([0 max(t)])

print(strcat(save_to,'velocity_components.eps'),'-depsc')

%% part E
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot magnitude of position and velocity vectors as function of time $\phantomsection\label{line: 2.13 part E}$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_I_mag = vecnorm(r_I,2,2);
v_I_mag = vecnorm(v_I,2,2);

times = 60*[30 60 90];

figure
subplot(2,1,1)
plot(t,r_I_mag,'k');
title('PROBA 2 Radius from Earth''s Center')
xlabel('Time (sec)')
ylabel('Radius (km)')
xlim([0 max(t)])
xline(times(1),'k',"30 mins")
xline(times(2),'k',"60 mins")
xline(times(3),'k',"90 mins")

subplot(2,1,2)
plot(t,v_I_mag,'k');
title('PROBA 2 Speed in ECI')
xlabel('Time (sec)')
ylabel('Speed (km/s)')
xlim([0 max(t)])
xline(times(1),'k',"30 mins")
xline(times(2),'k',"60 mins")
xline(times(3),'k',"90 mins")

% Vis viva equation:
syms r
v(r) = sqrt(mu*(2/r-1/a));

% Create check table
varNames = ["Time","Radius","vis_viva_calculated_speed","simulation_calculated_speed"];
varTypes = ["double","double","double","double"];
sz = [length(times) length(varNames)];
vis_viva_check = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);

% Calculate and tabulate the vis-viva check table
for loop_var = 1:length(times)
    [useless_var,index] = min(abs(t-times(loop_var))); %find index of time
    radius_at_point = r_I_mag(index); % Get radius at time
    vis_viva_check.Radius(loop_var) = radius_at_point;

    % Calculate velocity based on given radius
    vis_viva_check.vis_viva_calculated_speed(loop_var) = double(v(radius_at_point));
    vis_viva_check.simulation_calculated_speed(loop_var) = v_I_mag(index);
    
    vis_viva_check.Time(loop_var) = times(loop_var)/60;
    clear useless_var;
end
vis_viva_check.Properties.VariableNames = ["Time (minutes)", ...
    "Radius from Simulation(km)", ...
    "Calculated Speed from Vis-Viva equation (km/s)",...
    "Speed from Simulation (km/s)"];

writetable(vis_viva_check,strcat(save_to,'vis_viva_check.csv'))
print(strcat(save_to,'radius_and_speed.eps'),'-depsc')

%% part F
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Ground path $\phantomsection\label{line: 2.13 part F}$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wE = 15.04/3600; %converting 15.04 deg/hr to deg/s;

thetaE = wE.*t; %Array of theta_GMT

% Calculate r_F
r_F = cell2mat(arrayfun(@(i) r_I2r_F(thetaE(i), r_I(i, :)), ...
    1:length(thetaE), 'UniformOutput', false))';


% Calculate lattitude
lattitude = asind(r_F(:,3)./r_I_mag);

% Calculate longitude and place in correct quadrant
negate_angle = (le(r_F(:,2), 0) -0.5).*-2;
longitude = negate_angle.*acosd(r_F(:,1)./((r_I_mag.*cosd(lattitude))));

% Remove horizontal lines
dont_connect = find(diff(longitude) > 250);
for loop_var = 1:length(dont_connect)
    index = dont_connect(loop_var) +loop_var -1;
    longitude = [longitude(1:index);NaN;longitude((index+1):(length(longitude)))];
    lattitude = [lattitude(1:index);NaN;lattitude((index+1):(length(lattitude)))];
end

figure
plot(longitude, lattitude,'k','LineWidth',2);
ylim([-90 90])
xlim([-180 180])
title('PROBA 2 Ground Track')
xlabel('Longitude(^\circ)')
ylabel('Lattitude(^\circ)')
hold on
I = imread('map.png'); 
h = image(xlim,-ylim,I); 
uistack(h,'bottom')

print(strcat(save_to,'ground_tracks.eps'),'-depsc')




