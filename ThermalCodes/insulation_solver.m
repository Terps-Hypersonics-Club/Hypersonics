function [T_inshistory] = run_insulation_solver(inputfile,altitudefile,L,Nx,T_history)
    warning('off','all');
    altitude = altgen(altitudefile);
    % Access Material Properties
    %inputfile = 'carbonCeramic.txt';
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
    alpha = data.ThermalConductivity_W_m__C_/(data.Density_kg_m_3_*data.SpecificHeatCapacity_J_kg__C_); % Thermal diffusivity [m^2/s]
    dx = L/(Nx-1);    % Spatial Steps
    time = altitude(end,1); % Time [s]
    Fo = 0.4;         % Fourier Number
    dt = dx^2*Fo/alpha; % Time Steps
    Nt = ceil(time/dt); % # of Time Partitions
    T = zeros(Nx,1) + T_history(end,1);      % Initial Temperature Profile (ITP)
    T(1) = T_history(end,end);
    T_new = T;                    % Updated Temperature Profile (Initially the same as T)
    T_inshistory = zeros(Nx,Nt);     % Temperature profile across all time
    T_inshistory(:,1) = T;           % Initially set the first column to the ITP
    qcond = T(2)-T(1);
    for p = 1:Nt-1
        fprintf('Progress: %g%%\n' , round(p/Nt*100,2))
        fprintf('Time: %g\n', p*dt)
        for i = 2:Nx-1
            % Update the Internal Nodes
            T_new(i) = Fo*(T(i+1)+T(i-1)) + (1-2*Fo)*T(i);
        end
        % Update conductive end BC
        T_new(Nx) = T(Nx) + Fo*(T(Nx-1)-T(Nx));
        % Update the Edge Node
        %T_new(1) = alpha*dt/dx^2*(T(2)-T(1)) + T(1);
        T_new(1) = T(1)-(T_new(2)-T(2));
        % Update for Next Time Step
        T = T_new;
        T_inshistory(:,p+1) = T;
    end
    [x,t] = meshgrid(linspace(0,L,Nx),linspace(0,time,Nt));
	surf(x,t,T_inshistory'); hold on
	colormap turbo;
	surf(x, t, (74+273.15)*ones(size(x)), 'FaceAlpha', 0.7, 'EdgeColor', 'black', 'FaceColor', 'cyan');
	xlabel('Position [m]');
	ylabel('Time [s]');
	zlabel('Temperature [K]');
    zlim([0 1500]);
	title(sprintf('Material: %s', data.RecordName));
    subtitle(sprintf('Max Internal Temperature Requirement [C]: %.2f', 74));
	colorbar('eastoutside');
	shading interp;
end
T_history1 = run_insulation_solver('silica.txt','thermal_team_data.csv',0.01,10,T_history);
