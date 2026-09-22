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
        altitude(i,7) = (15.04-0.00649*altitude(i,2))+273.15; % Temperature [K]
        altitude(i,4) = (101.29*((altitude(i,7))/288.08)^5.256)*10^3; % Pressure [kPa]
       elseif altitude(i,2) > 11000 && altitude(i,2) <= 25000
        altitude(i,7) = 216.69; % Temperature [K]
        altitude(i,4) = (22.65*exp(1.73-0.000157*altitude(i,2)))*10^3; % Pressure [kPa]
       else
        altitude(i,7) = (-131.21 + 0.00299*altitude(i,2))+273.15; % Temperature [K]
        altitude(i,4) = (2.488*(altitude(i,7))^-11.388)*10^3; % Pressure [kPa]
       end
       altitude(i,6) = altitude(i,4)/((.2869*altitude(i,7))*10^3);
    end
    rhofit = polyfit(altitude(:,1),altitude(:,6),10); % Fits a function to the air density with respect to the trajectory timescale
    tempfit = polyfit(altitude(:,1),altitude(:,7),22); % Fits a function to the air temperature with respect to the trajectory timescale
    localM = polyfit(altitude(:,1),altitude(:,3),14); % Fits a function to the Mach number with respect to the trajectory timescale
    pfit = polyfit(altitude(:,1),altitude(:,4),18); % Fits a function to the air pressure with respect to the trajectory timescale
end
function [cpfit,cvfit] = air(inputfile)
    mcp = readmatrix(inputfile);
    tcp1 = mcp;
    for i = 1:size(mcp,1)
        tcp1(i,2) = mcp(i,3)*1000/28;
        tcp1(i,3) = mcp(i,3)*1000/28;
    end
    cpfit = polyfit(tcp1(:,1),tcp1(:,2),10); % Fits a function to specific heat capacity, cp based on temperature
    cvfit = polyfit(tcp1(:,1),tcp1(:,3),10); % Fits a function to specific heat capacity, cv based on temperature
end
function data = accessMat(inputfile)
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
end
%% 1-Dimensional Implicit Heat Analysis with conduction, convection, and radiation
function [T_history,dx] = run_1d_solver(inputfile,altitudefile,radius,layer,Nxtot,mode,heat_mode)
warning('off','all');
% Access Material Properties
data = accessMat(inputfile{1});
materials = cell(length(inputfile),1);
alphas = zeros(length(inputfile),1);
for j = 1:length(inputfile)
    data_j = accessMat(inputfile{j});
    materials{j} = data_j;
    alphas(j) = data_j.ThermalConductivity_W_m__C_ /(data_j.Density_kg_m_3_ * data_j.SpecificHeatCapacity_J_kg__C_);
end
alpha = alphas(1);
% Accesses specific heat data based on temperature
cpfit = air('N2.txt');
[altitude,rhofit,tempfit,localM,pfit] = altgen(altitudefile);
%eps = 0.5;        % Emissivity of Material
sig = 5.67e-8;    % Stefan-Boltzmann Constant
% Time and Spatial Step Initialization
dxs = zeros(1,length(layer));
Nx = zeros(1,length(layer));
for i = 1:length(layer)
    Nx(i) = round(layer(i)/sum(layer)*Nxtot);
    dxs(i) = layer(i)/(Nx(i)-1);    % Spatial Steps
end
time = altitude(end,1); % Time [s]
%Fo = 0.4;         % Fourier Number
dt = 0.001;
dx = dxs(1);
Fo = alpha*dt/dx^2; % Time Steps
Nt = ceil(time/dt); % # of Time Partitions

% Parameters
T_inf = altitude(1,7);        % Initial Ambient Temperature [K]
T = zeros(Nxtot,1) + T_inf;      % Initial Temperature Profile (ITP)
T_new = T;                    % Updated Temperature Profile (Initially the same as T)
T_history = zeros(Nxtot,Nt);     % Temperature profile across all time
T_history(:,1) = T;           % Initially set the first column to the ITP
Pr = 0.72;                    % Prandtl Number
r = 1;%(Pr)^(1/3);                        % Recovery Factor
RgasAir = 296.8;              % Gas specific constant [J/kg-K]
ems = 0.9;
p_inf = altitude(1,4);        % Freestream Pressure [Pa]
qddot1 = zeros(Nt,1);         % Heat Flux Check
gamma1 = qddot1;              % Gamma Check

% Temperature Profile Loop
Cf = 4e-4;
for p = 1:Nt
    fprintf('Progress: %g%%\n' , round(p/Nt*100,2))
    fprintf('Time: %g\n', p*dt)
    lMach = polyval(localM,p*dt); % Local Mach number
    p_inf = polyval(pfit,p*dt); % Local ambient pressure
    gamma = polyval(cpfit,T_inf)/(polyval(cpfit,T_inf)-RgasAir); % Ratio of Specific Heats
    T_stag = T_inf*(1+(gamma-1)/2*lMach^2); % Stagnation Temperature [K]
    rhoe = polyval(rhofit,p*dt)*((gamma+1)*(lMach)^2/((gamma-1)*(lMach)^2+2)); % Edge Density
    pe = p_inf*(1+(2*gamma)/(gamma+1)*((lMach*sind(90))^2-1)); % Edge pressure [Pa]
    ue = lMach*sqrt(RgasAir*gamma*T_inf)*(1-(2*((lMach)^2-1)/((gamma+1)*lMach^2))); % Horizontal component of velocity at the edge
    Te = T_inf*(pe/p_inf)/(rhoe/polyval(rhofit,p*dt)); % Edge Temperature [K]
    Taw = Te*(1+r*(gamma-1)/2*(ue/sqrt(RgasAir*gamma*Te))^2); % Adiabatic Wall Temperature [K]
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
    if heat_mode == "surface"
    qconv = Fra*tau_w*polyval(cpfit,T(1))*(Taw-T(1))/ue;
    end
    if heat_mode == "nose"
    qconv = 0.763*Pr^(-0.6)*(rhoe*mue)^0.5*(h_aw-h_w)*sqrt(dudx)*(cosd(0))^2;
    end
    qcond = (T(2)-T(1));
    qrad = sig*ems*(T(1))^4;
    qddot = qconv-qrad;
    % Update the Edge Node
    T_new(1) = dt/(dx*data.Density_kg_m_3_*data.SpecificHeatCapacity_J_kg__C_)*qddot + (alpha*dt/dx^2)*qcond + T(1);
    for i = 2:Nxtot-1
        T_new(i) = (alpha*dt/dx^2)*(T(i+1)+T(i-1)) + (1-2*(alpha*dt/dx^2))*T(i);
        % if i >= Nx(2)
        %     data = materials{2};
        %     alpha = alphas(2);
        %     dx = dxs(2);
        % end
        for j = 2:length(Nx)
            if i>=Nx(j)
                data = materials{j};
                alpha = alphas(j);
                dx = dxs(j);
            end
        end
    end
    % Update conductive end BC
    T_new(Nxtot) = T(Nxtot) + (alpha*dt/dx^2)*(T(Nxtot-1)-T(Nxtot));
    %T_new(Nx) = 288.16;
    % Update for Next Time Step
    T = T_new;
    T_history(:,p) = T;
    T_inf = polyval(tempfit,p*dt); % Update Ambient Temperature
    gamma1(p,1) = gamma;
    qddot1(p,1) = qddot;
end

% mode = 0 for steady state, mode = 1 for transient
data = materials{1};
% Plot Temperature vs Time
if mode == "3D"
	[x,t] = meshgrid(linspace(0,sum(layer),Nxtot),linspace(0,time,Nt));
    %pcolor(x,t,T_history')
	surf(x,t,T_history'); hold on
	colormap turbo;
    if heat_mode == "nose"
        surf(x, t, data.MaximumServiceTemperature_K_*ones(size(x)), 'FaceAlpha', 0.7, 'EdgeColor', 'black', 'FaceColor', 'cyan');
        subtitle(sprintf('Max Service Temperature [K]: %.2f', data.MaximumServiceTemperature_K_));
    else
        surf(x,t,ones(size(x))*(74+273.15), 'FaceAlpha', 0.7, 'EdgeColor', 'black', 'FaceColor', 'cyan');
        subtitle(sprintf('Max Internal Temperature: 74 [C]'));
    end
	xlabel('Position [m]');
	ylabel('Time [s]');
	zlabel('Temperature [K]');
    zlim([0 4000]);
	title(sprintf('Material: %s', data.RecordName));
	colorbar('eastoutside');
	shading interp;
elseif mode == "2D"
	x = linspace(0,L,Nxtot);
    t = linspace(0,time,Nt);
	plot(t,T_history(end,:))
	xlabel('Time [s]')
	ylabel('Temperature [K]')
	title('1D Transient Heat Transfer');
    hold on
end
figure;
% tiledlayout(2,3);
% 
% % Plot Heat Flux vs Time
% nexttile;
t1 = linspace(0,time,Nt);
plot(t1,qddot1)
xlabel('Time [s]')
ylabel('Heat Flux [W/m^2]')
title('Heat Flux over Time')
xlim([-10,time])
T_history1 = run_1d_solver({'HfO2.txt','304L.txt','inconel718.txt','silica.txt','inconel718.txt'},'thermal_trajectory_data.csv',0.003,[0.0005 0.003 0.003 0.006 0.003],50,"3D","surface");
