clc
clear;
close all;


% save folder
mydir  = pwd;
idcs   = strfind(string(mydir),'\');
newdir = mydir(1:idcs(end)-1);
save_to = strcat(newdir,'\output files\5.1\');
clear mydir idcs newdir;

global n
ALONG_TRACK = 1;
IN_TRACK = 2;
IN_PLANE_ELIPTICAL = 3;
ORBIT_I = 4;
ORBIT_II = 5;

%%%%%%%%%%%%%%%%%%%%%
formation_type = ORBIT_I;
%%%%%%%%%%%%%%%%%%%%%
% Gravitational parameter
mu = 398600.4418; % km^3/s^2

% Earth's radius
rE = 6378; % km

% Orbital elements of Envisat
a = rE+766; % km
e = 1E-4;
i = 98.4; % deg
RAAN = 27.2; % deg
w = 71; % deg

% Earth angular velocity 
wE = 72.922e-6; % rad/s

% Mean motion 
n = sqrt(mu/a^3); % rad/s

% Orbital period
T = 2*pi/n ; % sec

% Initial conditions of relative motion
if formation_type == ALONG_TRACK 
    % Calculate initial conditions of along track formation $\phantomsection\label{line: 5.1 case 1}$
    name = 'Along Track';
    x0 = 0; % km
    z0 = 0; % km 
    vx0 = 0; % km
    vy0 = 0; % km
    vz0 = 0; % km

    y0 = -10; % km 
elseif formation_type == IN_TRACK
    % Calculate initial conditions of in track formation $\phantomsection\label{line: 5.1 case 2}$
    name = 'In Track';
    x0 = 0; % km
    vx0 = 0; % km
    vy0 = 0; % km
    vz0 = 0; % km

    y0 = -10; % km 
    z0 = -(wE/n*sind(i))*y0; % km 
elseif formation_type == IN_PLANE_ELIPTICAL
    % Calculate initial conditions of in plane eliptical formation $\phantomsection\label{line: 5.1 case 3}$
    name = 'In Plane Eliptical';
    alpha = 0; % deg
    A_0 = 0.5/2; %km

    z0 = 0; % km
    vz0 = 0; % km

    % solve in plane eliptical equations
    syms vx0 x0
    eqns = [
        alpha == atand(-vx0/(n*x0));
        A_0 == sqrt((vx0/n)^2+x0^2)
        ];

    [vx0, x0] = vpasolve(eqns,[vx0 x0]);

    vx0 = round(double(vx0), 10);
    x0 =  round(double(x0), 10);
    vy0 = -2*n*x0; % km
    y0 = 2*vx0/n; % km 
elseif formation_type == ORBIT_I
    % Calculate initial conditions of 5.1) E)I formation $\phantomsection\label{line: 5.1 case 4}$
    name = 'Circular Formation, Radius 20m, Phasing Angle 0^\circ';
    r = 20E-3; % km
    phase = 0;  % deg
    
    x0 = (r/2)*cosd(phase); %km
    vx0 = -(r*n/2)*sind(phase); %km
    z0 = sqrt(3)*x0; % km
    vz0 = sqrt(3)*vx0; % km

    vy0 = -2*n*x0; % km
    y0 = 2*vx0/n; % km 
elseif formation_type == ORBIT_II
    % Calculate initial conditions of 5.1) E)II formation $\phantomsection\label{line: 5.1 case 5}$
    name = 'Projected Circular Formation, Radius 20m, Phasing Angle 270^\circ';
    r = 20E-3; % km
    phase = 270;  % deg
    
    x0 = (r/2)*cosd(phase); %km
    vx0 = -(r*n/2)*sind(phase); %km
    z0 = 2*x0; % km
    vz0 = 2*vx0; % km

    vy0 = -2*n*x0; % km
    y0 = 2*vx0/n; % km 
end
% Naming of output graphs
formation_type = num2str(formation_type);

%creating table of inital conditions
initial_conditions = array2table([vx0 vy0 vz0 x0 y0 z0]);
varNames = ["vx0", "vy0", "vz0", "x0", "y0", "z0"];
initial_conditions.Properties.VariableNames = varNames;
writetable(initial_conditions,strcat(save_to,'initial_conditions_',formation_type,'.csv'))

% Integration
options=odeset('RelTol',1e-9,'AbsTol',1e-9);
[t,X]=ode45(@CWHdyn,[0:10:2*T],[vx0 vy0 vz0 x0 y0 z0]',options);

vx = X(:,1);
vy = X(:,2);
vz = X(:,3);
x = X(:,4);
y = X(:,5);
z = X(:,6);

% Plot components of the relative position vector as function of time
figure
subplot(4,1,1)
plot(t/T,x*1000,'k');
title(name)
xlabel('Number of orbits')
ylabel('x (m)')
grid on
subplot(4,1,2)
plot(t/T,y*1000,'k');
xlabel('Number of orbits')
ylabel('y (m)')
grid on
subplot(4,1,3)
plot(t/T,z*1000,'k');
xlabel('Number of orbits')
ylabel('z (m)')
grid on
subplot(4,1,4)
plot(t/T,sqrt((x*1000).^2+(y*1000).^2+(z*1000).^2),'k');
xlabel('Number of orbits')
ylabel('Distance (m)')
axis([0 2 0 40]);
grid on


print(strcat(save_to,'Position_component_plots_',formation_type,'.eps'), '-depsc') % save figure

% Plot the components of the relative position vector in 3D
figure
plot3(z*1000, y*1000, x*1000,'k');
if str2double(formation_type) == ALONG_TRACK
    hold on
    scatter3(z*1000, y*1000, x*1000,20,'k');
end
hold on
scatter3(0,0,0,20,'k');
xlabel('z (m)')
ylabel('y (m)')
zlabel('x (m)')
title(name)
grid on
set(gca,'XDir','reverse')

print(strcat(save_to,'3D-plot_of_path_around_target_',formation_type,'.eps'), '-depsc') % save figure

% Plot and animate the in-plane and out-of-plane motion
planes = figure
set(planes,'Position',[556 33 454 665])
subplot(3,1,1)
plot(y*1000,x*1000,'k');
ip(1) = line(y(1)*1000, x(1)*1000, 'Marker', '.', 'MarkerSize', 22, 'Color', 'k');
ip(2) = line(0, 0,'Marker', 'o', 'MarkerSize', 6, 'LineWidth', 1.5, 'Color', 'k' );
axis([-70 70 -30 30]);
set(gca,'XDir','reverse')
set(gca,'FontSize',10,'FontName', 'Times')
xlabel('y (m)')
ylabel('x (m)')
title(name)
grid on

subplot(3,1,2)
plot(z*1000,x*1000,'k');
axis([-70 70 -30 30]);
oop(1) = line(z(1)*1000, x(1)*1000, 'Marker', '.', 'MarkerSize', 22, 'Color', 'k');
oop(2) = line(0, 0, 'Marker', 'o', 'MarkerSize', 6, 'LineWidth', 1.5, 'Color', 'k');
set(gca,'XDir','reverse')
set(gca,'FontSize',10,'FontName', 'Times')
xlabel('z (m)')
ylabel('x (m)')
grid on

subplot(3,1,3)
plot(z*1000,y*1000,'k');
axis([-70 70 -30 30]);
ct(1) = line(z(1)*1000, y(1)*1000, 'Marker', '.', 'MarkerSize', 22, 'Color', 'k');
ct(2) = line(0, 0, 'Marker', 'o', 'MarkerSize', 6, 'LineWidth', 1.5, 'Color', 'k');
set(gca,'XDir','reverse')
set(gca,'FontSize',10,'FontName', 'Times')
xlabel('z (m)')
ylabel('y (m)')
grid on
print(strcat(save_to,'Path_in_projected_planes_',formation_type,'.eps'), '-depsc')

for i = 1:length(x)
    set(ip(1),  'XData', y(i)*1000, 'YData', x(i)*1000);
    set(oop(1), 'XData', z(i)*1000, 'YData', x(i)*1000);
    set(ct(1), 'XData', z(i)*1000, 'YData', y(i)*1000);
    drawnow
end
