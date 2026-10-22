%% step2_build_ansys_loads.m   (THERMAL_PIPELINE, see README.md)
% Turn the isothermal-wall CHAMPS runs (one per case from step1) into ANSYS
% convection loads along the trajectory.
%
% Needs:
%   outputs/<traj>_cfd_cases.csv, outputs/<traj>_load_schedule.csv   (from step1)
%   cfd/case01.vtk, cfd/case02.vtk, ...   one CHAMPS surface output per case,
%                                         run with an ISOTHERMAL 315 K wall
%   geometry/<body>.stl                   outer surface of the body in the SAME
%                                         frame/units as the geometry in ANSYS
%
% For every CFD surface cell and case:
%   T_aw = T_e * (1 + r*(gamma-1)/2 * M_e^2)       (local edge state, turbulent r)
%   h    = q_w / (T_aw - T_wall),  T_wall = 315 K  (the isothermal wall CHAMPS used)
% ANSYS then applies q = h * (T_aw - T_structure), so the heating drops as the
% structure heats up. h is the COLD-wall value: exact at 315 K, high by up to
% ~1.5x on the body once the wall nears 1500 K (measured on the Mach 8 sea-level
% 315/1500 K pair). Conservative.
%
% Between cases the loads are scaled to the trajectory's freestream at each time
% (scaleToTrajectory = true): h * (rho V)^0.8 ratio, T_aw * recovery-temperature
% ratio, both relative to the point the case was run at (columns hScaleAcreage
% and TawScale of the step1 schedule). Set it false to hold each case's loads
% constant over its window.
%
% Outputs (outputs/<runName>_*):
%   _faces.csv       FaceID, X, Y, Z [m], then per case: h_caseNN [W/m^2K], Taw_caseNN [K]
%   _time_table.csv  rows for the ANSYS data table: time, case column, scale factors
% Then run step3_apply_loads_mechanical.py inside ANSYS Mechanical.
clear; clc; close all;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'lib'));

%% ---------------- User inputs ----------------
runName  = 'WAV_RID_9_FINS';
trajName = 'Trajectory_WAV_RID_COMP_aero';   % must match the trajectory step1 used
stlFile  = fullfile(here, 'geometry', 'WAV_RID_9_FINS.stl');   % '' = keep the CFD frame (only if ANSYS uses it too)
stlUnits = 'm';                              % units of the STL: 'm' or 'mm'
T_wall_cfd = 315;          % [K] isothermal wall temperature the CHAMPS runs used
scaleToTrajectory = true;  % scale each case to the trajectory freestream (see header)
tableTol = 0.01;           % thin the time table: keep rows needed to reproduce the
                           % scale factors within 1% by linear interpolation
makePlots = true;

% conversion settings (same as the rest of the FEA scripts)
mapOpt = struct('stlUnits', stlUnits, 'symAxis', 'y', 'mirrorCFD', 'auto', 'alignMode', 'shift', ...
    'gamma', 1.4, 'Rgas', 287.15, 'pran', 0.71, 'recovery', 'turbulent', 'dT_min', 20, ...
    'T_cold', T_wall_cfd, 'aftMask', 0.005, 'verbose', false);

outDir = fullfile(here, 'outputs');
cfdDir = fullfile(here, 'cfd');

%% ---------------- Cases and schedule from step1 ----------------
casesCsv = fullfile(outDir, [trajName '_cfd_cases.csv']);
schedCsv = fullfile(outDir, [trajName '_load_schedule.csv']);
assert(exist(casesCsv, 'file') == 2, 'Missing %s: run step1_plan_cfd_cases.m first.', casesCsv);
assert(exist(schedCsv, 'file') == 2, 'Missing %s: run step1_plan_cfd_cases.m first.', schedCsv);
C = readtable(casesCsv, 'TextType', 'char');
S = readtable(schedCsv);
nC = height(C);
names = compose('case%02d', (1:nC)');
missing = C.VTK_file(cellfun(@(f) exist(fullfile(cfdDir, f), 'file') ~= 2, C.VTK_file));
assert(isempty(missing), 'Missing CFD files in %s: %s', cfdDir, strjoin(missing, ', '));

%% ---------------- h and T_aw per case ----------------
nF = [];
for c = 1:nC
    vtk = fullfile(cfdDir, C.VTK_file{c});
    fprintf('Case %d (%s): Mach %.2f, AoA %+.2f, %.1f km, window %.1f-%.1f s\n', c, C.VTK_file{c}, ...
        C.Mach(c), C.Alpha_deg(c), C.Altitude_m(c)/1e3, C.t_start_s(c), C.t_end_s(c));
    m = map_cfd_to_stl(vtk, '', mapOpt);      % no STL target -> the CFD cells themselves
    if isempty(nF)
        nF = numel(m.h);  XYZ = m.stlCentroids;  H = zeros(nF, nC);  TAW = H;
    elseif numel(m.h) ~= nF || max(abs(m.stlCentroids - XYZ), [], 'all') > 1e-9
        error(['%s is on a different CFD surface mesh than %s. All cases must be run on ' ...
               'the same CHAMPS geometry file.'], C.VTK_file{c}, C.VTK_file{1});
    end
    % the runs must be isothermal at T_wall_cfd: check the wall temperature CHAMPS reports
    tw = m.Tw(m.Tw > 0);
    if isempty(tw) || abs(median(tw) - T_wall_cfd) > 1 || (prctile(tw, 95) - prctile(tw, 5)) > 1
        error(['%s does not look like an isothermal %g K run (wall temperature median %.0f K, ' ...
               '5-95%% %.0f-%.0f K). Rerun it with the isothermal wall, or set T_wall_cfd.'], ...
               C.VTK_file{c}, T_wall_cfd, median(tw), prctile(tw, 5), prctile(tw, 95));
    end
    H(:,c) = m.h;  TAW(:,c) = m.Taw;
    fprintf('  h median %.0f W/m^2K (max %.0f), T_aw median %.0f K, masked cells %d of %d\n', ...
        median(m.h, 'omitnan'), max(m.h), median(m.Taw, 'omitnan'), nnz(isnan(m.h)), nF);
end
H(isnan(H)) = 0;  TAW(isnan(TAW)) = 0;   % ANSYS cannot take NaN; h = 0 -> no heating there

%% ---------------- Into the geometry frame ----------------
if isempty(stlFile)
    warning(['No STL given: points stay in the CHAMPS frame. Only correct if the ANSYS ' ...
             'geometry uses exactly that frame.']);
else
    assert(exist(stlFile, 'file') == 2, 'Missing STL %s (see README step 0).', stlFile);
    shift = cfdToGeometryShift(XYZ, stlFile, stlUnits);
    XYZ = XYZ - shift;
    fprintf('\nPoints moved into the geometry frame by -[%.5f %.5f %.5f] m\n', shift);
end

%% ---------------- Time table ----------------
ok = S.Case > 0;                         % points step1 left without a case (below machMin)
t  = S.Time_s(ok);  a = S.Case(ok);
if scaleToTrajectory
    hS = S.hScaleAcreage(ok);  tS = S.TawScale(ok);
else
    hS = ones(size(t));  tS = ones(size(t));
end
nFull = numel(t);
keep = thinTable(t, a, log([hS tS]), log(1 + tableTol));
t = t(keep);  a = a(keep);  hS = hS(keep);  tS = tS(keep);
fprintf('\nTime table: %d of %d trajectory points kept (%d case switches)\n', numel(t), nFull, nnz(diff(a) ~= 0));

%% ---------------- Export ----------------
if ~exist(outDir, 'dir'), mkdir(outDir); end
csvF = fullfile(outDir, [runName '_faces.csv']);
csvT = fullfile(outDir, [runName '_time_table.csv']);
hdrF = [{'FaceID','X','Y','Z'}, reshape([strcat('h_', names), strcat('Taw_', names)]', 1, [])];
datF = zeros(nF, 4 + 2*nC);
datF(:, 1:4) = [(1:nF)', XYZ];
datF(:, 5:2:end) = H;  datF(:, 6:2:end) = TAW;
writeCompact(csvF, hdrF, datF, ['%d,%.6f,%.6f,%.6f', repmat(',%.6g,%.2f', 1, nC), '\n']);
writetable(table(t, a, names(a), strcat('h_', names(a)), hS, strcat('Taw_', names(a)), tS, ...
    'VariableNames', {'Time_s','CaseIdx','Case','hColumn','hScale','TawColumn','TawScale'}), csvT);
fprintf('\nWrote:\n  %s  (%d points x %d cases)\n  %s  (%d rows)\n', csvF, nF, nC, csvT, numel(t));
fprintf('Next: step3_apply_loads_mechanical.py in ANSYS Mechanical (RUN_NAME = "%s").\n', runName);

%% ---------------- Plot ----------------
if makePlots
    figure('Name', 'Load schedule', 'Position', [100 100 900 500]);
    tiledlayout(2, 1, 'TileSpacing', 'compact');
    nexttile; stairs(t, a, 'LineWidth', 1.2); ylabel('case'); grid on;
    title(sprintf('%s: %d cases', runName, nC), 'Interpreter', 'none');
    nexttile; semilogy(t, hS, t, tS, 'LineWidth', 1.2); ylabel('scale vs case'); xlabel('Time [s]');
    legend({'hScale', 'TawScale'}, 'Location', 'best'); grid on;
end

%% ---------------- Local functions ----------------
function shift = cfdToGeometryShift(cfdPts, stlFile, units)
% Rigid translation from the CFD solver frame to the STL frame (same body,
% different origin): bbox centres, then translation-only point-to-plane ICP.
TR = stlread(stlFile);
V = TR.Points;  F = TR.ConnectivityList;
if strcmpi(units, 'mm'), V = V / 1000; end
Cc = (V(F(:,1),:) + V(F(:,2),:) + V(F(:,3),:)) / 3;
N = cross(V(F(:,2),:) - V(F(:,1),:), V(F(:,3),:) - V(F(:,1),:), 2);
N = N ./ max(vecnorm(N, 2, 2), eps);
tree = KDTreeSearcher(Cc);
shift = (max(cfdPts) + min(cfdPts))/2 - (max(Cc) + min(Cc))/2;
for it = 1:30
    p = cfdPts - shift;
    j = knnsearch(tree, p);
    dn = sum((p - Cc(j,:)) .* N(j,:), 2);
    step = median(dn .* N(j,:), 1);
    shift = shift + step;
    if norm(step) < 1e-7, break; end
end
p = cfdPts - shift;  j = knnsearch(tree, p);
gap = abs(sum((p - Cc(j,:)) .* N(j,:), 2));
fprintf('CFD surface vs STL after alignment: gap median %.3f mm, 95%% %.3f mm\n', 1e3*median(gap), 1e3*prctile(gap, 95));
if prctile(gap, 95) > 1e-3
    warning('CFD surface and %s differ by more than 1 mm on 5%% of the points: wrong geometry or units?', stlFile);
end
end

function keep = thinTable(t, a, L, tolL)
% Fewest rows such that linear interpolation between kept rows reproduces every
% dropped row's log scale factors within tolL. Rows never span a case switch, so
% the last row of one case and the first row of the next are always kept.
n = numel(t);
keep = 1;  i = 1;
while i < n
    j = i + 1;
    while j < n && a(j+1) == a(i) && segOk(t, L, i, j + 1, tolL)
        j = j + 1;
    end
    keep(end+1, 1) = j; %#ok<AGROW>
    i = j;
end
end

function ok = segOk(t, L, i, j, tolL)
m = (i+1:j-1)';
if isempty(m), ok = true; return; end
w = (t(m) - t(i)) / (t(j) - t(i));
ok = all(abs(L(i,:) + w .* (L(j,:) - L(i,:)) - L(m,:)) <= tolL, 'all');
end

function writeCompact(file, hdr, dat, fmt)
fid = fopen(file, 'w');
assert(fid > 0, 'cannot open %s for writing (open in another program?)', file);
fprintf(fid, '%s\n', strjoin(hdr, ','));
fprintf(fid, fmt, dat.');
fclose(fid);
end
