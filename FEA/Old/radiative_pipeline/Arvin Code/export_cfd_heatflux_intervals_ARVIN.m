%% export_cfd_heatflux_intervals.m
% Raw CFD heat flux, one CHAMPS case per time interval. No scaling, no film
% coefficient: each case's own wall heat flux qw is applied unchanged for its
% window, then the next case takes over.
%
% (build_ansys_loads.m is the other route: h and T_aw per case, scaled to the
% trajectory at every time. This script is the simpler one.)
%
% Outputs (Exports/<runName>_*):
%   _heatflux_faces.csv      FaceID, X, Y, Z [m], then q_<case> [W/m^2] per case
%   _heatflux_intervals.csv  Case, qColumn, t_start_s, t_end_s (one row per case)
%
% In ANSYS Mechanical run ANSYS_SCRIPTS/apply_heatflux_intervals.py. It imports
% _heatflux_faces.csv as External Data and adds an Imported Heat Flux on SKIN
% whose value switches from one case to the next at each interval boundary.
%
% Note what qw is: CHAMPS's convective heat flux INTO the wall at the wall
% temperature CHAMPS assumed (radiative equilibrium if radiative = true). Applied
% as a fixed flux it does not drop as the structure heats up; apply surface
% radiation in ANSYS as well (the ANSYS script does) so the skin can shed heat.
clear; clc; close all;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'Calls'));

%% ---------------- User inputs (EDIT ALL FOUR) ----------------
runName = 'WAV_RID_9B';
% STL of the vehicle's outer surface, exported from the SAME geometry you mesh in ANSYS
% (half model is fine). Folder and file name relative to this script.
stlFile = fullfile(here, 'geometry', 'WAV_RID_9B.stl');
stlUnits = 'mm';          % units of that STL: 'mm' (Discovery/SolidWorks default) or 'm'
vtkDir  = fullfile(here, 'cfd');   % folder holding your CHAMPS surface VTKs

% One row per CFD run, in flight-time order. Intervals must abut (each starts where the
% previous one ends) and cover the whole flight. EXAMPLE ROWS: replace with your VTKs
% and the time window each one should cover.
caseList = {
%   name           vtk file                          t_start [s]  t_end [s]
    'M8_SL',       'Mach8_aoa1_SeaLevel.vtk',        0,           40
    'M8_24km',     'Mach8_aoa2_24km.vtk',            40,          75
    'M7_30km',     'Mach7_aoa4_30km.vtk',            75,          240
    'M6_25km',     'Mach6_aoa4_25km.vtk',            240,         730
    'M3_SL',       'Mach3_aoa-1_SeaLevel.vtk',       730,         824
};

% mapping options (same meaning as in build_ansys_loads.m; only the geometry
% handling and the aft mask matter here, the h/T_aw settings are unused)
mapOpt = struct('stlUnits', stlUnits, 'symAxis', 'y', 'mirrorCFD', 'auto', 'alignMode', 'shift', ...
    'gamma', 1.4, 'Rgas', 287.15, 'pran', 0.71, 'recovery', 'turbulent', 'dT_min', 20, ...
    'T_cold', 300, 'aftMask', 0.005, 'verbose', false);
makePlots = true;

outDir = fullfile(here, 'Exports');
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% ---------------- Cases ----------------
assert(exist(stlFile, 'file') == 2, 'STL not found: %s  (set stlFile at the top)', stlFile);
assert(exist(vtkDir, 'dir') == 7, 'VTK folder not found: %s  (set vtkDir at the top)', vtkDir);
cs = cell2table(caseList, 'VariableNames', {'name','vtk','t_start','t_end'});
nC = height(cs);
cs.vtkPath = fullfile(vtkDir, cs.vtk);
for c = 1:nC
    assert(exist(cs.vtkPath{c}, 'file') == 2, 'Missing VTK: %s', cs.vtkPath{c});
end
assert(all(cs.t_end > cs.t_start), 'Every interval needs t_end > t_start.');
assert(all(abs(cs.t_start(2:end) - cs.t_end(1:end-1)) < 1e-9), ...
    'Intervals must abut: each t_start must equal the previous t_end.');

%% ---------------- Read qw at the CFD cell centroids ----------------
nF = [];
for c = 1:nC
    fprintf('Reading %s (%s) ...\n', cs.name{c}, cs.vtk{c});
    m = map_cfd_to_stl(cs.vtkPath{c}, '', mapOpt);   % no STL target -> the CFD cells themselves
    if isempty(nF)
        nF = numel(m.qw);  XYZ = m.stlCentroids;  Q = zeros(nF, nC);
    elseif numel(m.qw) ~= nF || max(abs(m.stlCentroids - XYZ), [], 'all') > 1e-9
        error(['%s is on a different CFD surface mesh than %s. All cases must share one ' ...
               'surface mesh for this script.'], cs.vtk{c}, cs.vtk{1});
    end
    Q(:,c) = m.qw;
    fprintf('  qw median %.3g W/m^2, max %.3g W/m^2, masked/zero cells %d\n', ...
        median(m.qw, 'omitnan'), max(m.qw), nnz(isnan(m.qw) | m.qw == 0));
end
Q(isnan(Q)) = 0;   % ANSYS cannot take NaN; 0 = no heating (solver holes, aft mask)

% CFD points are in the solver frame; ANSYS uses the geometry (STL/STEP) frame
shift = cfdToGeometryShift(XYZ, stlFile, stlUnits);
XYZ = XYZ - shift;
fprintf('\nExport points: %d CFD cell centroids, moved into the geometry frame by -[%.5f %.5f %.5f] m\n', nF, shift);

%% ---------------- Export ----------------
csvF = fullfile(outDir, [runName '_heatflux_faces.csv']);
csvI = fullfile(outDir, [runName '_heatflux_intervals.csv']);
qCols = strcat('q_', cs.name)';
writeCompact(csvF, [{'FaceID','X','Y','Z'}, qCols], [(1:nF)', XYZ, Q], ...
    ['%d,%.6f,%.6f,%.6f', repmat(',%.6g', 1, nC), '\n']);
writetable(table(cs.name, strcat('q_', cs.name), cs.t_start, cs.t_end, ...
    'VariableNames', {'Case','qColumn','t_start_s','t_end_s'}), csvI);
fprintf('\nWrote:\n  %s  (%d points x %d cases)\n  %s\n', csvF, nF, nC, csvI);
fprintf('\nIntervals:\n');
for c = 1:nC
    fprintf('  %-12s %7.1f - %7.1f s   qw max %.3g W/m^2\n', cs.name{c}, cs.t_start(c), cs.t_end(c), max(Q(:,c)));
end

%% ---------------- Plot ----------------
if makePlots
    figure('Name', 'Heat flux per interval', 'Position', [100 100 1200 300*ceil(nC/2)]);
    tiledlayout(ceil(nC/2), 2, 'TileSpacing', 'compact');
    lim = [0 prctile(Q(Q > 0), 98)];
    for c = 1:nC
        nexttile;
        scatter3(XYZ(:,1), XYZ(:,2), XYZ(:,3), 2, Q(:,c), 'filled');
        axis equal; view(0, 90); colormap(turbo); clim(lim); colorbar;
        title(sprintf('%s: %.0f-%.0f s', cs.name{c}, cs.t_start(c), cs.t_end(c)), 'Interpreter', 'none');
    end
end

%% ---------------- Local functions ----------------
function shift = cfdToGeometryShift(cfdPts, stlFile, units)
% Rigid translation from the CFD solver frame to the STL/STEP frame (same body,
% different origin): bbox centres, then translation-only point-to-plane ICP.
% (Same as build_ansys_loads.m.)
TR = stlread(stlFile);
V = TR.Points;  F = TR.ConnectivityList;
if strcmpi(units, 'mm'), V = V / 1000; end
C = (V(F(:,1),:) + V(F(:,2),:) + V(F(:,3),:)) / 3;
N = cross(V(F(:,2),:) - V(F(:,1),:), V(F(:,3),:) - V(F(:,1),:), 2);
N = N ./ max(vecnorm(N, 2, 2), eps);
% Half-model STL (cut at y = 0) vs full-body CFD: mirror the STL across y = 0 for the
% alignment only. Without this the loads land ~57 mm off in y on a half model.
if (max(C(:,2)) < 1e-4 || min(C(:,2)) > -1e-4) && range(cfdPts(:,2)) > 1.5*range(C(:,2))
    fprintf('STL is a half model (one side of y = 0): mirroring it across y = 0 for the alignment.\n');
    C = [C; C .* [1 -1 1]];  N = [N; N .* [1 -1 1]];
end
tree = KDTreeSearcher(C);
shift = (max(cfdPts) + min(cfdPts))/2 - (max(C) + min(C))/2;
for it = 1:30
    p = cfdPts - shift;
    j = knnsearch(tree, p);
    dn = sum((p - C(j,:)) .* N(j,:), 2);
    step = median(dn .* N(j,:), 1);
    shift = shift + step;
    if norm(step) < 1e-7, break; end
end
p = cfdPts - shift;  j = knnsearch(tree, p);
gap = abs(sum((p - C(j,:)) .* N(j,:), 2));
fprintf('CFD surface vs %s after alignment: gap median %.3f mm, 95%% %.3f mm, 99%% %.3f mm\n', ...
    stlFile, 1e3*median(gap), 1e3*prctile(gap, 95), 1e3*prctile(gap, 99));
if prctile(gap, 95) > 1e-3
    warning('CFD surface and %s differ by more than 1 mm on 5%% of the points: wrong geometry or units?', stlFile);
end
end

function writeCompact(file, hdr, dat, fmt)
fid = fopen(file, 'w');
assert(fid > 0, 'cannot open %s for writing (open in another program?)', file);
fprintf(fid, '%s\n', strjoin(hdr, ','));
fprintf(fid, fmt, dat.');
fclose(fid);
end
