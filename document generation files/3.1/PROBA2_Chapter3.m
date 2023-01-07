clc
clear;
close all;

% save folder
mydir  = pwd;
idcs   = strfind(string(mydir),'\');
newdir = mydir(1:idcs(end)-1);
save_to = strcat(newdir,'\output files\3.1\');
% Add general Functions
addpath(strcat(newdir,'\general functions'));
clear mydir idcs newdir;
% J2 constant
J_2 = 1.08264E-3;

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
i       = 98.28; % deg $\phantomsection\label{line: 3.1 inclination val}$
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

% Components of position and velocity vectors in perifocal at t = 0 sec
r_P_ini = r*[1 0 0]';

v_P_ini = [0 sqrt(mu/p)*(e+1) 0]';

% Rotation matrix from perifocal to ECI, i.e., C_IP

C_IP = (C_3(w)*C_1(i)*C_3(RAAN))'; % using functions I made  for rot matrices

% Components of position and velocity vectors in ECI at t = 0 sec
r_I_ini = C_IP*r_P_ini;
v_I_ini = C_IP*v_P_ini;

orbits = 1;
T = 2*pi*sqrt(a^3/mu)*orbits;

% Simulation from t = 0 to t = T 
open_system('PROBA2mdl_Chapter3.slx')
set_param('PROBA2mdl_Chapter3', 'StopTime', 'T')
disp('Running Simulation...')
sim('PROBA2mdl_Chapter3')

% Fix for simulink adding a dimension
v_I = squeeze(v_I)';
r_I = squeeze(r_I)';

%% Part A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate secular change in RAAN $\phantomsection\label{line: 3.1 part A}$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RAAN_secular_change = -((3*pi*J_2*rE^2)/(p^2))*cosd(i)*orbits;% rad
RAAN_secular_change = RAAN_secular_change*r2d; % Convert to degrees
RAAN_average_change = RAAN_secular_change/T; % Average change of RAAN (deg/s)
%% Part C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate and plot RAAN $\phantomsection\label{line: 3.1 part C}$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = cross(r_I,v_I); % Calculate angular momentum
N = cross(repmat([0 0 1],length(h),1),h); % Obtain vector pointing at RAAN

RAAN_Array = atan2d(N(:,2),N(:,1)); % find RAAN angle in degrees
RAAN_Array = le(RAAN_Array, 0)*360 + RAAN_Array; % add 360 if less than 0


figure
hold off

% Plot actual RAAN
plot(t, RAAN_Array) 

%plot expected change based on secular change
hold on
plot([0,max(t)],[RAAN, RAAN+RAAN_secular_change])

% Graph formatting and saving
legend('Simulated RAAN', 'expected change in RAAN','Location','best')
xlabel('time (s)')
ylabel('RAAN (^o)')
title(strcat("PROBA-2 RAAN with J_2 Pertubations at ",num2str(i),'^{\circ} Inclined Orbit'))
xlim("tight")
ylim("tight")
print(strcat(save_to,'RAAN_',num2str(i),'_degrees.eps'), '-depsc')


%% Part D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate simulated secular and rate of change in RAAN $\phantomsection\label{line: 3.1 part D}$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate change in RAAN from simulation
simulated_RAAN_secular_change = ...
    (RAAN_Array(length( RAAN_Array)) -  RAAN_Array(1))/orbits;

% Calculated  rate of change change in RAAN from simulation
simulated_RAAN_average_rate_of_change = ...
    simulated_RAAN_secular_change/(T/orbits);

%% Part E
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tabulize expected and simulated secular and rate of change in RAAN $\phantomsection\label{line: 3.1 part E}$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate differences
RAAN_secular_difference = abs(RAAN_secular_change-simulated_RAAN_secular_change);
RAAN_roc_difference = abs(RAAN_average_change-simulated_RAAN_average_rate_of_change);
varNames = ["secular_change", "rate_of_change"];
varTypes = repmat("double", 2,1);
sz = [3 2];
comparison_table = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);

comparison_table.secular_change = ...
[RAAN_secular_change;simulated_RAAN_secular_change;RAAN_secular_difference];

comparison_table.rate_of_change = ...
[RAAN_average_change; simulated_RAAN_average_rate_of_change; RAAN_roc_difference];

comparison_table.Properties.VariableNames = ["secular change (degrees/rev)",... 
    "average rate of change(degrees/s)"];
comparison_table.Properties.RowNames = ["expected","simulation", "difference"];

writetable(comparison_table,strcat(save_to,'RAAN_check_for_inc_',num2str(i),'.csv'),'WriteRowNames',true)

