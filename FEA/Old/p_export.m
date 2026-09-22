%% Clear Workspace
clear; clc; close all;

%% File Names
vtkFile   = 'surf_000034193.vtk';
inputFile = 'input (1).sdf';
stlFile   = 'Odd_seed_Point_101_round.stl';
csvPressureOut    = 'pressure_mach5dot5.csv';
csvHeatFluxOut    = 'heatflux_mach5dot5.csv';

%% Load Flow Conditions
[P_inf, T_inf, M_inf, rho_inf, pran, y, Rgas] = load_input(inputFile);

%% Load CFD Surface Data
% Expecting load_vtk_surf to return centroids, surface normals, areas, etc.
[T_e, P_e, U_e, V_e, W_e, centroids, norms, areas] = load_vtk_surf(vtkFile);

%% Calculate Heat Flux
% Speeds and Mach
v_mag      = sqrt(U_e.^2 + V_e.^2 + W_e.^2);
v_mag_safe = max(v_mag, 1e-6);
a_e        = sqrt(y*Rgas.*T_e);
M_e        = v_mag_safe ./ a_e;

% Local state and freestream total T
rho_e = P_e ./ (Rgas .* T_e);
T_t   = T_inf * (1 + (y-1)/2 * M_inf^2);

% Viscosity and Re_x
mu_e = 1.716e-5 .* (T_e/273.15).^(3/2) .* (383.15) ./ (T_e + 110.4);
x_le = min(centroids(:,1));
x    = max(centroids(:,1) - x_le, 1e-6);
Re_x = rho_e .* v_mag_safe .* x ./ mu_e;

% Laminar vs turbulent, recovery factor, T_aw
isTurb = Re_x > 4e3;          % pick a better criterion if you have one
r = zeros(size(Re_x));
r(~isTurb) = sqrt(pran);
r( isTurb) = pran.^(1/3);
T_w  = 294 + zeros(size(T_e));
T_aw = T_e .* (1 + r .* ((y-1)/2) .* M_e.^2);

% Compressible Cf with transformation
Cf_lam = 0.664 ./ sqrt(Re_x);
Cf_tur = 0.0592 ./ (Re_x).^(1/5);
TsTe   = 0.5 + 0.039.*M_e.^2 + 0.5.*(T_w./T_e);
musmue = (TsTe).^(3/2) .* (1 + 110.4./T_e) ./ (TsTe + 110.4./T_e);
OsOe_lam = (musmue).^(1/2) ./ (TsTe).^(1/2);
OsOe_tur = (musmue).^(1/5) ./ (TsTe).^(4/5);
Cf = zeros(size(Re_x));
Cf(~isTurb) = Cf_lam(~isTurb) .* OsOe_lam(~isTurb);
Cf( isTurb) = Cf_tur( isTurb) .* OsOe_tur( isTurb);

% Wall shear and heat flux via Chilton–Colburn
q_dyn     = 0.5 .* rho_e .* v_mag_safe.^2;
Tau_w     = Cf .* q_dyn;
Cp        = y*Rgas/(y-1);
Pr_factor = pran.^(-2/3);
q_w = Tau_w .* (Cp .* Pr_factor) .* (T_aw - T_w) ./ v_mag_safe;

% Mask bad cells
%invalid_idx = ~isfinite(q_w) | (areas < 1e-12);
invalid_idx = (T_e <= 0) | ~isfinite(q_w);
q_w(invalid_idx) = 0;

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
csvPressureData = [FaceID, stlCentroids/1000, mappedPressures];
csvHeatFluxData = [FaceID, stlCentroids/1000, mappedHeatFlux];

% Export pressure values to csv
headers = {'FaceID','X','Y','Z','Pressure'};
writecell(headers, csvPressureOut);
writematrix(csvPressureData, csvPressureOut, 'WriteMode','append');

fprintf('Export complete. Data saved to %s\n', csvPressureOut);

% Export heat flux values to csv
headers = {'Heat Flux','FaceID'};
writecell(headers, csvHeatFluxOut);
writematrix(csvHeatFluxData(:,[5,1]), csvHeatFluxOut, 'WriteMode','append');

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
clim(prctile(mappedHeatFlux, [1 99]));
colormap(turbo)
