%% export_cfd_heating.m
% Export CFD (CHAMPS) surface loads to ANSYS using the solver's OWN wall heat
% flux, instead of the Reynolds-analogy model in exportpressureandheat_working.m.
%
% Requires a surface VTK with fields: P, T, U, V, W, qw, wall_temperature.
% Produces, per structural face (or per CFD cell if no STL is given):
%   <case>_pressure.csv     FaceID, X, Y, Z, Pressure [Pa]
%   <case>_convection.csv   FaceID, X, Y, Z, h [W/m^2K], T_aw [K]      <- use as ANSYS CONV load
%   <case>_heatflux_cfd.csv FaceID, X, Y, Z, qw [W/m^2], T_wall_cfd [K] <- reference only
%
% Why h and T_aw rather than qw: CHAMPS's qw is the flux into the wall AT THE
% WALL TEMPERATURE IT ASSUMED (here radiative equilibrium, ~2000-3500 K). The
% structure starts cold and heats up, so the flux must follow the wall
% temperature. q = h*(T_aw - T_wall) does that automatically if ANSYS is given
% a convection load (film coefficient h, bulk temperature T_aw).
clear; clc; close all;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'Calls'));

%% ---------------- User inputs ----------------
vtkFile   = fullfile(here, 'Calls', 'surf_out_000000001 (1).vtk');
stlFile   = fullfile(here, 'WAV_RID.stl');   % from WAV_RID_STEP.step via Calls/step2stl.py; '' -> export at CFD cell centroids
stlUnits  = 'mm';          % units of the STL file: 'mm' or 'm'  (CFD is assumed to be in m)
caseName  = 'mach8_sealevel';
symAxis   = 'y';           % axis normal to the CFD symmetry plane ('x','y','z'); '' = none
mirrorCFD = 'auto';        % 'auto' | true | false : reflect half-model CFD to cover a full STL
alignMode = 'shift';       % 'none' | 'shift' : translate only (same geometry, different origin)
                           % | 'bbox' : per-axis scale+shift STL onto CFD bounding box
gamma     = 1.4;
Rgas      = 287.15;
pran      = 0.71;          % laminar Prandtl number
recovery  = 'turbulent';   % 'turbulent' (r = Pr^(1/3)) | 'laminar' (r = sqrt(Pr)) | numeric value
dT_min    = 20;            % [K] minimum (T_aw - T_wall) for a trustworthy h; below -> masked
T_cold    = 300;           % [K] reference wall temp for the "cold-wall equivalent" diagnostic
aftMask   = 0.005;         % [m] mask heating on faces within this distance of the base (aft face + rim); 0 = off
makePlots = true;

outDir = fullfile(here, 'Exports');
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% ---------------- Map CFD onto the STL ----------------
opt = struct('stlUnits', stlUnits, 'symAxis', symAxis, 'mirrorCFD', mirrorCFD, 'alignMode', alignMode, ...
    'gamma', gamma, 'Rgas', Rgas, 'pran', pran, 'recovery', recovery, 'dT_min', dT_min, ...
    'T_cold', T_cold, 'aftMask', aftMask, 'verbose', true);
m = map_cfd_to_stl(vtkFile, stlFile, opt);
useSTL = m.useSTL;  stlFaces = m.stlFaces;  stlVertices = m.stlVertices;  stlCentroids = m.stlCentroids;
mappedP = m.P;  mappedH = m.h;  mappedTaw = m.Taw;  mappedQ = m.qw;  mappedTw = m.Tw;

%% ---------------- Export ----------------
FaceID = (1:size(stlCentroids,1))';
XYZ = stlCentroids;                              % metres
csvP = fullfile(outDir, [caseName '_pressure.csv']);
csvC = fullfile(outDir, [caseName '_convection.csv']);
csvQ = fullfile(outDir, [caseName '_heatflux_cfd.csv']);
% ANSYS tables cannot take NaN; h = 0 -> no heating. Zero copies only, so the
% plots below still show masked faces as NaN (pink) instead of 0.
Hcsv = mappedH;      Hcsv(isnan(Hcsv)) = 0;
TawCsv = mappedTaw;  TawCsv(isnan(TawCsv)) = 0;
Qcsv = mappedQ;      Qcsv(isnan(Qcsv)) = 0;
writecell({'FaceID','X','Y','Z','Pressure'}, csvP);
writematrix([FaceID, XYZ, mappedP], csvP, 'WriteMode','append');
writecell({'FaceID','X','Y','Z','h','T_aw'}, csvC);
writematrix([FaceID, XYZ, Hcsv, TawCsv], csvC, 'WriteMode','append');
writecell({'FaceID','X','Y','Z','qw_cfd','T_wall_cfd'}, csvQ);
writematrix([FaceID, XYZ, Qcsv, mappedTw], csvQ, 'WriteMode','append');
fprintf('\nWrote:\n  %s\n  %s\n  %s\n', csvP, csvC, csvQ);
fprintf('ANSYS: apply %s as a CONVECTION load (film coefficient h, bulk temperature T_aw), not HFLUX.\n', [caseName '_convection.csv']);

%% ---------------- Plots ----------------
if makePlots
    if useSTL
        F = stlFaces; Vv = stlVertices;
    else
        % build a triangle list from the CFD file for plotting
        F = []; Vv = [];
    end
    fields = {mappedQ, 'CFD wall heat flux qw [W/m^2]'; mappedH, 'Film coefficient h [W/m^2K]'; ...
              mappedTaw, 'Adiabatic wall temperature [K]'; mappedP, 'Pressure [Pa]'};
    for k = 1:size(fields,1)
        figure('Name', fields{k,2});
        if useSTL
            trisurf(F, Vv(:,1), Vv(:,2), Vv(:,3), 'FaceVertexCData', fields{k,1}, 'FaceColor','flat', 'EdgeColor','none');
            hold on; mk = isnan(fields{k,1});
            if any(mk), trisurf(F(mk,:), Vv(:,1), Vv(:,2), Vv(:,3), 'FaceColor',[1 0.3 0.8], 'EdgeColor','none'); end
            hold off;
        else
            scatter3(stlCentroids(:,1), stlCentroids(:,2), stlCentroids(:,3), 4, fields{k,1}, 'filled');
            hold on; mk = isnan(fields{k,1});
            if any(mk), scatter3(stlCentroids(mk,1), stlCentroids(mk,2), stlCentroids(mk,3), 12, [1 0.3 0.8], 'filled'); end
            hold off;
        end
        axis equal; colorbar; colormap(turbo);
        v = fields{k,1}(~isnan(fields{k,1})); if ~isempty(v), clim(prctile(v,[2 98])); end
        xlabel('X'); ylabel('Y'); zlabel('Z');
        title(sprintf('%s  (pink = masked, %d)', fields{k,2}, nnz(mk)));
    end
end
