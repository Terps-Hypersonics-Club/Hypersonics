%% Clear Workspace
%% Last updated: 11/13/25 3:49 PM by Will
clear; clc; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'Calls'));

%% File Names (resolved relative to this script's folder)
here      = fileparts(mfilename('fullpath'));
vtkFile   = fullfile(here, 'Calls', 'surf_000034193.vtk');
inputFile = fullfile(here, 'Calls', 'input (1).sdf');
stlFile   = fullfile(here, 'Old',   'Odd_seed_Point_101_round.stl');
outDir    = fullfile(here, 'Exports');
if ~exist(outDir, 'dir'), mkdir(outDir); end
csvPressureOut    = fullfile(outDir, 'pressure_mach5dot5_new.csv');
csvHeatFluxOut    = fullfile(outDir, 'heatflux_mach5dot5_new.csv');

%% Load Flow Conditions
[P_inf, T_inf, M_inf, rho_inf, pran, y, Rgas] = load_input(inputFile);

%% Load CFD Surface Data
% Expecting load_vtk_surf to return centroids, surface normals, areas, etc.
[T_e, P_e, U_e, V_e, W_e, centroids, norms, areas] = load_vtk_surf(vtkFile);

%% Calculate Heat Flux
T_w = zeros(size(T_e)) + 294;
r = zeros(size(T_e));
x = zeros(size(T_e));
TsTe = zeros(size(T_e));
v_mag = zeros(size(T_e));
x = centroids(:, 1);
v_mag = sqrt(U_e.^2+V_e.^2+W_e.^2);
M_e = v_mag./sqrt(y*Rgas*T_e);
T_aw = T_e.*(1+r*((y-1)/2).*M_e.^2); 		% adiabatic wall temperature
T_t = T_inf*(1+(y-1)/2*M_inf);		% total stagnation temperature
rho_e = rho_inf*((y+1)*M_inf^2/(y-1)*M_inf^2+2);
squiggle_w = T_e/T_t;
F_RA = zeros(size(squiggle_w));  % Preallocate
% Case 1: squiggle_w < 0.2
F_RA(squiggle_w < 0.2) = 1;
% Case 2: 0.2 <= squiggle_w <= 0.65
idx_mid = squiggle_w >= 0.2 & squiggle_w <= 0.65;
F_RA(idx_mid) = 0.8311 + 0.9675 .* squiggle_w(idx_mid) - 0.6142 .* squiggle_w(idx_mid).^2;
% Case 3: squiggle_w > 0.65
F_RA(squiggle_w > 0.65) = 1.2;
mu_e = 1.716*10^-5.*(T_e/273.1).^(3/2)*383.1./(T_e+110);
Re = rho_e*v_mag.*x./mu_e;
if Re > 4000 % if turbulent flow - gemini
r = (Pr).^(1/3);
TsTe = 0.5+0.039.*M_e.^2 + 0.5.*(T_w./T_e);
musmue  = (TsTe).^(3/2).*(1+110./T_e)./(TsTe+110./T_e);
OsOe = (musmue).^(1/5)./(TsTe).^(4/5);
cfi = 0.0592./(Re).^(1/5);
cfc = cfi.*OsOe;
else % laminar flow
r = sqrt(pran);
TsTe = 0.5+0.039.*M_e.^2 + 0.5.*(T_w./T_e);
musmue  = (TsTe).^(3/2).*(1+110./T_e)./(TsTe+110./T_e);
OsOe = (musmue).^(1/2)./(TsTe).^(1/2);
cfi = 0.664./(Re).^(1/2);
cfc = cfi.*OsOe;
end
q_e = 0.5*rho_e.*v_mag.^2;
Tau_w = cfc.*q_e;					% Wall shear stress
Cp = y*Rgas/(y-1);
q_w = Cp./v_mag.*(T_aw - T_w).*F_RA.*Tau_w;
% Define invalid mask
invalid_idx = (norms(:,1) == 1) | (T_e == 0);
% Set invalid entries to 0 or NaN (your choice)
q_w(invalid_idx) = 0;  % or NaN if you want to clearly flag them

%% Load Structural Mesh
TR = stlread(stlFile);
stlFaces     = TR.ConnectivityList;
stlVertices  = TR.Points;

% Recompute centroids with rotated vertices
stlCentroids = (stlVertices(stlFaces(:,1),:) + ...
                stlVertices(stlFaces(:,2),:) + ...
                stlVertices(stlFaces(:,3),:)) / 3;

%% --- Alignment Step (Scale + Translate) ---
% Compute ranges
rangeCFD = max(centroids) - min(centroids);
rangeSTL = max(stlCentroids) - min(stlCentroids);

% Compute scaling factors (per axis)
scaleFactors = rangeCFD ./ rangeSTL;

% Apply scaling
stlCentroidsScaled = stlCentroids .* scaleFactors;

% Compute means
meanCFD = mean(centroids,1);
meanSTL = mean(stlCentroidsScaled,1);

% Apply translation
stlCentroidsAligned = stlCentroidsScaled - meanSTL + meanCFD;

%% Debug Visualization: Check alignment
figure;
scatter3(centroids(:,1), centroids(:,2), centroids(:,3), 8, 'r', 'filled'); hold on;
scatter3(stlCentroidsAligned(:,1), stlCentroidsAligned(:,2), stlCentroidsAligned(:,3), 8, 'b');
axis equal; grid on;
legend('CFD centroids','Aligned STL centroids');
title('CFD vs STL alignment check');

%% --- Mirror Across Z-axis (symmetry fill) ---

% Positive-Z side assumed to have CFD data
posMask = centroids(:,1) >= 0;
posNodes = centroids(posMask,:);
posP     = P_e(posMask);
posHF    = q_w(posMask);

% For negative-Z side: mirror across Z, then assign via nearest neighbor
negMask = centroids(:,1) < 0;
negNodes = centroids(negMask,:);

if any(negMask)
    mirroredNegNodes = negNodes;
    mirroredNegNodes(:,1) = -mirroredNegNodes(:,1);

    % Build KD-tree from positive-Z nodes
    Mdl = KDTreeSearcher(posNodes);
    idxMirror = knnsearch(Mdl, mirroredNegNodes);

    % Replace pressures on negative side
    P_e(negMask) = posP(idxMirror);
    q_w(negMask) = posHF(idxMirror);
end

% Re-map pressures to STL using the updated symmetric CFD pressures
idx = knnsearch(centroids, stlCentroidsAligned);
mappedPressures = P_e(idx);
mappedHeatFlux = q_w(idx);

%% Diagnostics
fprintf('Unique CFD centroids: %d\n', length(unique(round(centroids,6),'rows')));
fprintf('Unique STL centroids: %d\n', length(unique(round(stlCentroidsAligned,6),'rows')));
fprintf('Mapped pressure range: [%.2f, %.2f]\n', min(mappedPressures), max(mappedPressures));

%% Export to CSV for ANSYS
FaceID = (1:size(stlCentroids,1))';
csvPressureData = [FaceID, stlCentroids, mappedPressures];
csvHeatFluxData = [FaceID, stlCentroids, mappedHeatFlux];

% Export pressure values to csv
headers = {'FaceID','X','Y','Z','Pressure'};
writecell(headers, csvPressureOut);
writematrix(csvPressureData, csvPressureOut, 'WriteMode','append');

fprintf('Export complete. Data saved to %s\n', csvPressureOut);

% Export heat flux values to csv
headers = {'FaceID','X','Y','Z','Heat Flux'};
writecell(headers, csvHeatFluxOut);
writematrix(csvHeatFluxData, csvHeatFluxOut, 'WriteMode','append');

fprintf('Export complete. Data saved to %s\n', csvHeatFluxOut);

%% Visualization: Pressure field on STL

% Pressure Plotting
figure;
trisurf(stlFaces, stlVertices(:,1), stlVertices(:,2), stlVertices(:,3), ...
        'FaceVertexCData', mappedPressures,'FaceColor','flat','EdgeColor','none');
axis equal; colorbar;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Pressure Field on Mesh');

% Convective Heat Flux Plotting
figure;
trisurf(stlFaces, stlVertices(:,1), stlVertices(:,2), stlVertices(:,3), ...
        'FaceVertexCData', mappedHeatFlux,'FaceColor','flat','EdgeColor','none');
axis equal; colorbar;

xlabel('X'); ylabel('Y'); zlabel('Z');
title('Convective Heat Flux Field on Mesh');
% Set scale
clim(prctile(mappedHeatFlux, [5 95]));
colormap(turbo)
