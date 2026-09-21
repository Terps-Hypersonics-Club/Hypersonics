%% 1-Dimensional Implicit Heat Analysis with conduction, convection, and radiation
clear;

here = fileparts(mfilename('fullpath'));
inputfile = fullfile(here,'Materials','graphiteParaPlane.txt');
[fidmat, msg] = fopen(inputfile,'r');
assert(fidmat > 0, 'Could not open "%s": %s', inputfile, msg);

fig = 1;
fidmat = fopen(inputfile,'r');
data = struct;
while ~feof(fidmat)
line = fgetl(fidmat);
if contains(line, ':')  % Only process lines with key-value pairs
  parts = strsplit(line, ':');
  key = strtrim(parts{1});
  value = strtrim(strjoin(parts(2:end), ':'));  % In case value has ":"
   % Convert numeric values if possible
  num = str2double(value);
  if ~isnan(num)
      data.(matlab.lang.makeValidName(key)) = num;
  else
      data.(matlab.lang.makeValidName(key)) = value;
  end
end
end

% Accesses specific heat data based on temperature
mcp = readmatrix('Air_Mixtures/N2.txt');
tcp1 = mcp;
for i = 1:size(mcp,1)
    tcp1(i,2) = mcp(i,3)*1000/28;
    tcp1(i,3) = mcp(i,3)*1000/28;
end
cpfit = polyfit(tcp1(:,1),tcp1(:,2),10); % Fits a function to specific heat capacity, cp based on temperature
cvfit = polyfit(tcp1(:,1),tcp1(:,2),10); % Fits a function to specific heat capacity, cv based on temperature
alpha = data.ThermalConductivity_W_m__C_/(data.Density_kg_m_3_*data.SpecificHeatCapacity_J_kg__C_); % Thermal diffusivity [m^2/s]
[altitude,rhofit,tempfit,localM,pfit] = altgen('thermal_team_data.csv');
eps = 0.5;        % Emissivity of Material
sig = 5.67e-8;    % Stefan-Boltzmann Constant
% Time and Spatial Step Initialization
L = 0.1;          % Characteristic Length
radius = 0.002;   % Nose Radius
Nx = 100;         % # of Spatial Partitions
dx = L/(Nx-1);    % Spatial Steps
time = altitude(end,1); % Time [s]
Fo = 0.3;         % Fourier Number
dt = dx^2*Fo/alpha; % Time Steps
Nt = ceil(time/dt); % # of Time Partitions

% Parameters
T_inf = altitude(1,7);        % Initial Ambient Temperature [K]
T = zeros(Nx,1) + T_inf;      % Initial Temperature Profile (ITP)
T_new = T;                    % Updated Temperature Profile (Initially the same as T)
T_history = zeros(Nx,Nt);     % Temperature profile across all time
T_history(:,1) = T;           % Initially set the first column to the ITP
Pr = 0.71;                    % Prandtl Number
r = 1;                        % Recovery Factor
RgasAir = 296.8;              % Gas specific constant [J/kg-K]
p_inf = altitude(1,4);        % Freestream Pressure [Pa]
rho_inf = altitude(1,6);      % Freestream Density [kg/m^3]
qddot1 = zeros(Nt,1);         % Heat Flux Check
gamma1 = qddot1;              % Gamma Check

% Temperature Profile Loop
Cf = 0.002;
for p = 1:Nt
    fprintf('Progress: %g%%\n' , round(p/Nt*100,2))
    fprintf('Time: %g\n', p*dt)
    lMach = polyval(localM,p*dt); % Local Mach number
    %gamma = polyval(cpfit,T_inf)/polyval(cvfit,T_inf); % Ratio of Specific Heats
    gamma = polyval(cpfit,T_inf)/(polyval(cpfit,T_inf)-RgasAir); % Ratio of Specific Heats
    T_stag = T_inf*(1+(gamma-1)/2*lMach^2); % Stagnation Temperature [K]
    rhoe = polyval(rhofit,p*dt)*((gamma+1)*(lMach*sind(90))^2/((gamma-1)*(lMach*sind(90))^2+2)); % Edge Density
    pe = p_inf*(1+(2*gamma)/(gamma+1)*((lMach*sind(90))^2-1)); % Edge pressure
    ue = lMach*sqrt(RgasAir*gamma*T_inf)*(1-(2*((lMach*sind(90))^2-1)/((gamma+1)*lMach^2))); % Horizontal component of velocity at the edge
    Te = T_inf*(pe/p_inf)/(rhoe/polyval(rhofit,p*dt)); % Edge Temperature
    Taw = Te*(1+(gamma-1)/2*(ue/sqrt(RgasAir*gamma*Te))^2); % Adiabatic Wall Temperature
    mue = 1.716*10^-5.*(Te/273.1).^(3/2)*383.1./(Te+110); % Viscosity at the edge
    dudx = (1/radius)*sqrt(2*(pe-p_inf)/rhoe); % Gradient of the velocity in the x1 direction
    h_aw = polyval(cpfit,Taw)*Taw; % Adiabatic Wall Enthalpy
    h_w = polyval(cpfit,T(1))*T(1); % Wall Enthalpy
    qe = 0.5*rhoe*ue^2;
    lambdaW = T(1)/T_stag;
    if lambdaW < 0.2
        Fra = 1;
    elseif lambdaW >= 0.2 && lambdaW <=0.65
        Fra = 0.8311 + 0.9675*lambdaW - 0.6142*lambdaW^2;
    else
        Fra = 1.2;
    end
    tau_w = Cf*qe;
    %qconv = Fra*tau_w*polyval(cpfit,T(1))*(Taw-T(1))/ue;
    qconv = 0.57*Pr^(-0.6)*(rhoe*mue)^0.5*(h_aw-h_w)*sqrt(dudx);
    qcond = T(2)-T(1);
    qrad = sig*eps*T(1)^4;
    qddot = qconv-qrad;
    % Update the Edge Node
    T_new(1) = dt/(dx*data.Density_kg_m_3_*data.SpecificHeatCapacity_J_kg__C_)*qddot + (alpha*dt/dx^2)*qcond + T(1);
    for i = 2:Nx-1
        % Update the Internal Nodes
        T_new(i) = Fo*(T(i+1)+T(i-1)) + (1-2*Fo)*T(i);
    end
    % Update conductive end BC
    T_new(Nx) = T(Nx) + Fo*(T(Nx-1)-T(Nx));
    % Update for Next Time Step
    T = T_new;
    T_history(:,p) = T;
    T_inf = polyval(tempfit,p*dt); % Update Ambient Temperature
    T(p,1) = T_inf;
    gamma1(p,1) = gamma;
    qddot1(p,1) = qddot;
end

mode = 1; % mode = 0 for steady state, mode = 1 for transient

% Plot Temperature vs Time
if mode == 1
	[x,t] = meshgrid(linspace(0,L,Nx),linspace(0,time,Nt));
	surf(x,t,T_history'); hold on
	colormap turbo;
	tplane = surf(x, t, data.MaximumServiceTemperature_K_*ones(size(x)), 'FaceAlpha', 0.7, 'EdgeColor', 'black', 'FaceColor', 'cyan');
	xlabel('Position [m]');
	ylabel('Time [s]');
	zlabel('Temperature [K]');
    zlim([0 4000]);
	title(sprintf('Material: %s', data.RecordName));
    subtitle(sprintf('Max Service Temperature [K]: %.2f', data.MaximumServiceTemperature_K_));
	colorbar('eastoutside');
	shading interp;
elseif mode == 0
	x = linspace(0,L,Nx);
	plot(x,T_history(:,end))
	xlabel('Position [m]')
	ylabel('Temperature [K]')
	title('1D Transient Heat Transfer');
end
figure;
tiledlayout(2,3);

% Plot Heat Flux vs Time
nexttile;
t1 = linspace(0,time,Nt);
plot(t1,qddot1)
xlabel('Time [s]')
ylabel('Heat Flux [W/m^2]')
title('Heat Flux over Time')
xlim([-10,time])

% Plot Gamma vs Time
nexttile;
y = linspace(1.397,1.402,100);
plot(t1,gamma1,'r-','LineWidth',2); hold on
title('Gamma over Time')

% Create patch coordinates
x_patch = [min(t1), max(t1), max(t1), min(t1)];
y_patch = [min(y), min(y), max(y), max(y)];

% Add translucent horizontal band
patch(x_patch, y_patch, 'yellow', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
xlabel('Time [s]')
ylabel('Gamma')
ylim([1.39 1.41])

% Plot Altitude vs Time
nexttile;
plot(altitude(:,1),altitude(:,2),'b-');
xlabel('Time [s]')
ylabel('Altitude [m]')
title('Altitude over Time')

% Plot Mach Number vs Time
nexttile;
plot(altitude(:,1),altitude(:,3),'r-','LineWidth',2);
xlabel('Time [s]')
ylabel('Mach Number')
title('Mach Number over Time')

function [altitude,rhofit,tempfit,localM,pfit] = altgen(inputfile)
% Access Altitude Data
data = readtable(inputfile);
% Step 2: Remove rows with NaN values
cleanedData = rmmissing(data);
altime = table2array(cleanedData); % Convert table data to a readable matrix
altitude = zeros(size(altime,1),size(altime,2)); % Creates a new matrix to fill with atmospheric data
% Loads all of the .csv data into the altitude matrix
for i = 1:size(altime,2)
   altitude(:,i) = altime(:,i);
end
% Loops through each row of the altitude matrix
% Using the atmospheric model from NASA we can calculate the atmospheric
% pressure at all altitudes
for i = 1:size(altitude,1)
   if altitude(i,2) <= 11000 && altitude(i,2) >= 0
    %altitude(i,7) = (15.04-0.00649*altitude(i,2))+273.15; % Temperature [K]
    altitude(i,4) = (101.29*((altitude(i,7))/288.08)^5.256)*10^3; % Pressure [kPa]
   elseif altitude(i,2) > 11000 && altitude(i,2) <= 25000
    %altitude(i,7) = 216.69; % Temperature [K]
    altitude(i,4) = (22.65*exp(1.73-0.000157*altitude(i,2)))*10^3; % Pressure [kPa]
   else
    %altitude(i,7) = (-131.21 + 0.00299*altitude(i,2))+273.15; % Temperature [K]
    altitude(i,4) = (2.488*(altitude(i,7))^-11.388)*10^3; % Pressure [kPa]
   end
end
rhofit = polyfit(altitude(:,1),altitude(:,6),10); % Fits a function to the air density with respect to the trajectory timescale
tempfit = polyfit(altitude(:,1),altitude(:,7),22); % Fits a function to the air temperature with respect to the trajectory timescale
localM = polyfit(altitude(:,1),altitude(:,3),14); % Fits a function to the Mach number with respect to the trajectory timescale
pfit = polyfit(altitude(:,1),altitude(:,4),14); % Fits a function to the air pressure with respect to the trajectory timescale
end