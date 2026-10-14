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
        altitude(i,4) = (2.488*(altitude(i,7)/216.6)^-11.388)*10^3; % Pressure [kPa]
       end
       altitude(i,6) = altitude(i,4)/((.2869*altitude(i,7))*10^3);
    end
    rhofit = @(t) exp(interp1(altitude(:,1), log(altitude(:,6)), min(t,altitude(end,1)), 'pchip')); % Fits a function to the air density with respect to the trajectory timescale
    tempfit = @(t) interp1(altitude(:,1), altitude(:,7), t, 'spline'); % Fits a function to the air temperature with respect to the trajectory timescale
    localM =  @(t) interp1(altitude(:,1), altitude(:,3), t, 'spline'); % Fits a function to the Mach number with respect to the trajectory timescale
    pfit =  @(t) exp(interp1(altitude(:,1), log(altitude(:,4)), min(t,altitude(end,1)), 'pchip')); % Fits a function to the air pressure with respect to the trajectory timescale
end
function [cpfit,hfit] = air(inputfile)
    mcp = readmatrix(inputfile);

    tcp1 = mcp;
    tcp1(:,2) = mcp(:,3) * 1000 / 28;
    T_grid = min(tcp1(:,1)):0.5:max(tcp1(:,1));
    Cp_grid = interp1(tcp1(:,1), tcp1(:,2), T_grid, 'spline');
    h_grid = cumtrapz(T_grid, Cp_grid);

    % Interpolate C_p values on the fine grid
    cpfit = @(T) interp1(tcp1(:,1), tcp1(:,2), T, 'spline');
    hfit = @(T) interp1(T_grid, h_grid, T, 'spline');
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
function [Taw,T_stag,dudx,rhoe,mue,tau_w,ue] = NormShock(lMach,p_inf,rho_inf,gamma,T_inf,R,radius,Pr)
    r = (Pr)^(1/3);      % Recovery Factor
    Cf = 4e-4;           % Skin friction coefficient
    T_stag = T_inf*(1+(gamma-1)/2*lMach^2); % Stagnation Temperature [K]
    rhoe = rho_inf*((gamma+1)*(lMach)^2/((gamma-1)*(lMach)^2+2)); % Edge Density
    pe = p_inf*(1+(2*gamma)/(gamma+1)*(lMach^2-1)); % Edge pressure [Pa]
    ue = lMach*sqrt(R*gamma*T_inf)*(1-(2*((lMach)^2-1)/((gamma+1)*lMach^2))); % Horizontal component of velocity at the edge
    Te = T_inf*(pe/p_inf)/(rhoe/rho_inf); % Edge Temperature [K]
    Taw = Te*(1+r*(gamma-1)/2*(ue/sqrt(R*gamma*Te))^2); % Adiabatic Wall Temperature [K]
    mue = 1.716*10^-5.*(Te/273.1).^(3/2)*383.1./(Te+110); % Viscosity at the edge
    dudx = (1/radius)*sqrt(2*(pe-p_inf)/rhoe); % Gradient of the velocity in the x1 direction
    qe = 0.5*rhoe*ue^2; % Edge dynamic pressure
    tau_w = Cf*qe; % Shear stress estimate
end
%% 1-Dimensional Explicit Heat Analysis with conduction, convection, and radiation
function T_history = run_1d_solver(inputfile,altitudefile,radius,layer,Nxs,mode,heat_mode)
% Access Material Properties
materials = cell(length(inputfile),1);
alphas = zeros(length(inputfile),1);
for j = 1:length(inputfile)
    data_j = accessMat(inputfile{j});
    materials{j} = data_j;
    alphas(j) = data_j.ThermalConductivity_W_m__C_ /(data_j.Density_kg_m_3_ * data_j.SpecificHeatCapacity_J_kg__C_);
end
% Accesses specific heat data based on temperature
[cpfit,hfit] = air('N2.txt');
[altitude,rhofit,tempfit,localM,pfit] = altgen(altitudefile);

% Time and Spatial Step Initialization
dxs = layer./(Nxs-1); % Spatial Steps
time = altitude(end,1); % Time [s]
Fo = 0.4;         % Fourier Number
dt = min(dxs)^2*Fo/max(alphas);
Nt = ceil(time/dt); % # of Time Partitions
dx = dxs(1);
interfaces = cumsum(Nxs(1:end-1)); % Interface nodes

% Parameters
T_inf = altitude(1,7);        % Initial Ambient Temperature [K]
T = zeros(sum(Nxs),1) + T_inf;      % Initial Temperature Profile (ITP)
T_new = T;                    % Updated Temperature Profile (Initially the same as T)
T_history = zeros(sum(Nxs),Nt);     % Temperature profile across all time
T_history(:,1) = T;     % Initially set the first column to the ITP
Pr = 0.72;              % Prandtl Number
R = 287.05;             % Gas specific constant [J/kg-K]
sige = 5.67e-8*0.9;    % Radiation Constants
qddot1 = zeros(Nt,1);   % Heat Flux Check
gamma1 = qddot1;        % Gamma Check
alpha = alphas(1); % First material

% Temperature Profile Loop
for p = 1:Nt
    fprintf('Progress: %g%%\n' , round(p/Nt*100,2))
    %fprintf('Time: %.4g\n', p*dt)
    lMach = localM(p*dt); % Local Mach number
    p_inf = pfit(p*dt); % Local ambient pressure [Pa]
    rho_inf = rhofit(p*dt); % Local ambient density [kg/m^3]
    gamma = cpfit(T_inf)/(cpfit(T_inf)-R); % Ratio of Specific Heats
    [Taw,T_stag,dudx,rhoe,mue,tau_w,ue] = NormShock(lMach,p_inf,rho_inf,gamma,T_inf,R,radius,Pr);
    h_w = hfit(T(1)); % Wall Enthalpy
    h_aw = hfit(Taw); % Adiabatic Wall Enthalpy
    lambdaW = T(1)/T_stag;
    if lambdaW < 0.2
        Fra = 1; % Reynolds' Analogy Factor
    elseif lambdaW >= 0.2 && lambdaW <=0.65
        Fra = 0.8311 + 0.9675*lambdaW - 0.6142*lambdaW^2;
    else
        Fra = 1.2;
    end
    if heat_mode == "surface"
    qconv = Fra*tau_w*cpfit(T(1))*(Taw-T(1))/ue;
    end
    if heat_mode == "nose"
    qconv = 0.753*Pr^(-0.6)*(rhoe*mue)^0.5*(h_aw-h_w)*sqrt(dudx)*(cosd(0))^2;
    end
    % Update the Edge Node
    data = materials{1};
    a = dt/(dx*data.Density_kg_m_3_*data.SpecificHeatCapacity_J_kg__C_);
    qcond = (alphas(1)*dt/dx^2)*(T(2)-T(1));
    qrad = sige*(T(1))^4;
    qddot = qconv-qrad;
    T_new(1) = a*qddot + qcond + T(1);
    for i = 2:sum(Nxs)-1
        jInterface = find(i == interfaces, 1);
        if ~isempty(jInterface)
            j = jInterface;
            data1 = materials{j};
            data2 = materials{j+1};
            k1 = data1.ThermalConductivity_W_m__C_;
            k2 = data2.ThermalConductivity_W_m__C_;
            rc1 = data1.Density_kg_m_3_*data1.SpecificHeatCapacity_J_kg__C_;
            rc2 = data2.Density_kg_m_3_*data2.SpecificHeatCapacity_J_kg__C_;
            rctot = (rc1*dxs(j)+rc2*dxs(j+1))/(dxs(j)+dxs(j+1));
            dudx1 = (T(i) - T(i-1))/dxs(j);
            dudx2 = (T(i+1) - T(i))/dxs(j+1);
            dkdudxdx = (k2*dudx2-k1*dudx1)/((dxs(j+1)+dxs(j))/2);
            T_new(i) = T(i) + dt/rctot*dkdudxdx;
        else
            j = find(i <= cumsum(Nxs), 1);
            data = materials{j};
            alpha = alphas(j);
            dx2 = dxs(j)^2;
            T_new(i) = (alpha*dt/dx2)*(T(i+1)+T(i-1)) + (1-2*(alpha*dt/dx2))*T(i);
        end
    end
    % Update conductive end BC
    T_new(sum(Nxs)) = T(sum(Nxs)) + (alphas(end)*dt/dxs(end)^2)*(T(sum(Nxs)-1)-T(sum(Nxs)));
    %T_new(sum(Nxs)) = 288.16;
    % Update for Next Time Step
    T = T_new;
    T_history(:,p) = T;
    T_inf = tempfit(p*dt); % Update Ambient Temperature
    gamma1(p,1) = gamma;
    qddot1(p,1) = qddot;
    %[Tmin,imin] = min(T); [Tmax,imax] = max(T);
    %fprintf('h=%d Tinf=%.1f M=%.2f pe=%.3f T1=%.1f cp1=%.1f Tmin=%.1f@%d Tmax=%.1f@%d\n', ...
    %h_w, T_inf, lMach, pe, T(1),cp1, Tmin, imin, Tmax, imax);
end

% mode = 0 for steady state, mode = 1 for transient
data = materials{1};
% Plot Temperature vs Time
if mode == "3D"
	[x,t] = meshgrid(linspace(0,sum(layer),sum(Nxs)),linspace(0,time,Nt));
	surf(x,t,T_history'); hold on
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
    colormap(turbo(256));
	colorbar('eastoutside');
	shading interp;
elseif mode == "2D"
    x = linspace(0, sum(layer), sum(Nxs));
    t = linspace(0, time, Nt);

    imagesc(x, t, T_history');
    set(gca, 'YDir', 'normal');
    set(gca, 'Layer', 'top');

    colormap(turbo(256));
    colorbar('eastoutside');
	h = findobj(gca, 'Type', 'image');
    set(h, 'Interpolation', 'bilinear');
    xlabel('Position [m]','FontSize',1);
	ylabel('Time [s]','FontSize',1);
    if heat_mode == "nose"
        title('Forebody Thermal Simulation','FontSize',1)
    elseif heat_mode == "surface"
        title('Surface Thermal Simulation','FontSize',1)
    end
    cb = gca;
    cb.FontSize = 14;
	% x = linspace(0,L,Nxtot);
    % t = linspace(0,time,Nt);
	% plot(t,T_history(end,:))
	% xlabel('Time [s]')
	% ylabel('Temperature [K]')
	% title('1D Transient Heat Transfer');
    % hold on
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

% % Create patch coordinates
% x_patch = [min(t1), max(t1), max(t1), min(t1)];
% y_patch = [min(y), min(y), max(y), max(y)];
% 
% % Add translucent horizontal band
% patch(x_patch, y_patch, 'yellow', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% xlabel('Time [s]')
% ylabel('Gamma')
% ylim([1.39 1.41])

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
end

%% Function Call
clear; clc;
T_history1 = run_1d_solver({'HfO2.txt','inconel718.txt'},'Trajectory_WAV_RID_8_aero.csv',0.002,[0.0005 0.005],[2 5],"2D","surface");
