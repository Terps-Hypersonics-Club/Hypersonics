%% build_ansys_loads.m
% Turn a set of CHAMPS cases (one surface VTK each) into time-dependent ANSYS
% heating and pressure loads along the trajectory.
%
%   1. Converts every case's VTK to h = qw/(T_aw - T_wall), T_aw, P per CFD cell
%      (map_cfd_to_stl, the same conversion as export_cfd_heating.m).
%      exportAt = 'cfd' writes those values at the CFD cell centroids, moved into
%      the geometry frame: no interpolation, so nothing is lost before ANSYS maps
%      them. Needs every case on the same CFD surface mesh (checked).
%      exportAt = 'stl' maps them onto the STL faces instead (one extra
%      nearest-neighbour step, and a much bigger file for a fine STL).
%   2. Assigns each time to one case: nearest in Mach / AoA / density, or fixed windows.
%   3. Scales that case's loads from the freestream it was run at to the
%      trajectory's freestream at that time.
%   4. Thins the time table: keeps only the rows needed to reproduce every
%      trajectory point's scale factors within tableTol by linear interpolation
%      (ANSYS ramps imported loads between rows), plus both rows at every case switch.
%
% Outputs (Exports/<runName>_*):
%   _faces.csv       FaceID, X, Y, Z [m], then per case: h_<case> [W/m^2K], Taw_<case> [K]   (thermal)
%   _pressure.csv    FaceID, X, Y, Z [m], then per case: P_<case> [Pa]                        (structural)
%   _time_table.csv  per time: which case, and the scale/offset for its columns
%   _coverage.csv    per trajectory point: actual vs case Mach/AoA and the mismatch
%
% In ANSYS Mechanical: run ANSYS_SCRIPTS/apply_trajectory_loads.py from
% Automation > Scripting. It imports _faces.csv as External Data (X/Y/Z in metres)
% and adds one Imported Convection row per _time_table.csv row:
%   film coefficient = h_<case>   * hScale
%   ambient temp     = Taw_<case> * TawScale
% (Imported Pressure, for a structural run: P_<case> * pScale + pOffset from _pressure.csv.)
% Also apply surface radiation on the same faces.
%
% Scaling (same as plan_cfd_cases.m): h ~ (rho V)^0.8 (turbulent acreage) or
% sqrt(rho) V (stagnation), T_aw from the freestream recovery temperature,
% surface pressure rise ~ dynamic pressure at fixed Mach and AoA. Mach and AoA
% differences between a case and the trajectory are NOT corrected; the coverage
% report shows how large they are.
clear; clc; close all;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'Calls'));
addpath(fullfile(here, 'WAV_RID_9'));

%% ---------------- User inputs ----------------
runName  = 'WAV_RID_traj9';
trajFile = fullfile(here, 'WAV_RID_9', 'Trajectory_WAV_RID_9_aero.csv');
stlFile  = fullfile(here, 'WAV_RID_9', 'WAV_RID_9.STL');   % geometry frame (and target faces for exportAt = 'stl')
vtkDir   = fullfile(here, 'WAV_RID_9');      % folder holding the case VTKs
exportAt = 'cfd';       % 'cfd': values at CFD cell centroids (lossless, ~135k rows) | 'stl': mapped onto stlFile faces

% One row per CHAMPS run. VTK files are looked up in vtkDir.
% p_inf / T_inf: NaN -> standard atmosphere at alt_km; otherwise the run's Pref [Pa] / Tref [K].
caseList = {
%   name           vtk file                        Mach   AoA[deg]  alt[km]  p_inf  T_inf
    'M8_a0',       'Mach7.98_aoa0.38_11.8km.vtk',             7.98,    0.38,    11.8,   NaN,   NaN
    'M7_a3p5',     'Mach6.76_aoa3.61_25km.vtk',           6.76,    3.61,    25,   NaN,   NaN
    'M4_a7',       'Mach3.87_aoa6.96_21km.vtk',             3.87,     6.96,    21,   NaN,   NaN
    'M3_a1p5',     'Mach3.25_aoa1.68_26.5km.vtk',           3.25,    1.68,    26.5,   NaN,   NaN
   'M3_aneg1p5',   'Mach3.56_aoa-1.48_15.3km.vtk',        3.56,   -1.48,    15.3,   NaN,   NaN
%   'M8_a0_SL',    'surf_out_000000001 (1).vtk',   8.0,    0.00,     0.00,   NaN,   NaN   % existing sea-level run: covers launch without 2.8x scaling
};
missingVTK = 'skip';    % 'skip': build loads from the cases that exist (warns) | 'error'

assignMode = 'nearest'; % 'nearest': each time -> closest case in (Mach, AoA, density)
                        % 'windows': case c applies from caseStart_s(c) until the next one starts
caseStart_s = [];       % [s] one per case (in caseList order), for 'windows'
dMach   = 0.5;          % distance weights for 'nearest' and the coverage flags:
dAlpha  = 1.0;          %   one unit = this much Mach / AoA [deg] /
dLogRho = log(4);       %   a factor-4 density change (scaling corrects density, so weighted lightly)

hScaleType = 'acreage'; % 'acreage': (rho V)^0.8 | 'stag': sqrt(rho) V  (nose / leading edges)
pRef     = 101325;      % [Pa] subtracted from pressure (gauge vs sealed interior at sea level); 0 = absolute
tStep    = [];          % [s] uniform time-table spacing; [] = adaptive (every trajectory point, thinned to tableTol)
tableTol = 0.01;        % adaptive table: max relative error of hScale / TawScale / pScale between kept rows
pOffTol  = 100;         % [Pa] adaptive table: max error of pOffset between kept rows
gamma    = 1.4;
Rgas     = 287.05;
recovery = 0.71^(1/3);  % for the T_aw scale (same as export_cfd_heating.m)

% mapping options (same meaning as in export_cfd_heating.m)
mapOpt = struct('stlUnits', 'm',   'symAxis', 'y', 'mirrorCFD', 'auto', 'alignMode', 'shift', ...
    'gamma', gamma, 'Rgas', 287.15, 'pran', 0.71, 'recovery', 'turbulent', 'dT_min', 20, ...
    'T_cold', 300, 'aftMask', 0.005, 'verbose', false);
makePlots = true;

outDir = fullfile(here, 'Exports');
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% ---------------- Cases: freestream and VTKs ----------------
cs = cell2table(caseList, 'VariableNames', {'name','vtk','Mach','AoA','alt_km','p_inf','T_inf'});
[Ti, pI, ~] = isa1976(1e3*cs.alt_km);
useIsaP = isnan(cs.p_inf);  cs.p_inf(useIsaP) = pI(useIsaP);
useIsaT = isnan(cs.T_inf);  cs.T_inf(useIsaT) = Ti(useIsaT);
cs.rho  = cs.p_inf ./ (Rgas*cs.T_inf);
cs.V    = cs.Mach .* sqrt(gamma*Rgas*cs.T_inf);
cs.vtkPath = fullfile(vtkDir, cs.vtk);
have = cellfun(@(f) exist(f, 'file') == 2, cs.vtkPath);
fprintf('Cases (%d of %d VTKs found):\n', nnz(have), height(cs));
for c = 1:height(cs)
    fprintf('  %-10s M %.2f  AoA %+5.2f  alt %5.2f km  p %8.1f Pa  T %5.1f K   %s\n', cs.name{c}, cs.Mach(c), cs.AoA(c), ...
        cs.alt_km(c), cs.p_inf(c), cs.T_inf(c), iff(have(c), 'ok', ['MISSING ' cs.vtk{c}]));
end
if ~all(have)
    if strcmp(missingVTK, 'error') || ~any(have)
        error('Missing VTKs: %s', strjoin(cs.vtk(~have), ', '));
    end
    warning('Building loads from %d of %d cases; the trajectory is assigned to the available ones only.', nnz(have), height(cs));
    if strcmp(assignMode, 'windows'), caseStart_s = caseStart_s(have); end
    cs = cs(have, :);
end
nC = height(cs);

%% ---------------- Trajectory ----------------
tr = readtable(trajFile);
traj.t = tr.Time_s;  traj.alt = max(tr.Altitude_m, 0);  traj.M = tr.Mach;  traj.AoA = tr.Alpha_deg;
if isempty(tStep)
    tt = traj.t;
else
    tt = (traj.t(1):tStep:traj.t(end))';
    if tt(end) < traj.t(end), tt(end+1) = traj.t(end); end
end
fs   = freestream(traj, traj.t, gamma, Rgas, recovery);   % at trajectory points (coverage)
fsT  = freestream(traj, tt,     gamma, Rgas, recovery);   % at time-table points

%% ---------------- Assign times to cases ----------------
caseOf = @(f) assignCases(f, cs, assignMode, caseStart_s, dMach, dAlpha, dLogRho);
aTraj  = caseOf(fs);
aTab   = caseOf(fsT);

%% ---------------- Map each VTK onto the export points ----------------
switch exportAt
    case 'cfd', tgtFile = '';        % map_cfd_to_stl with no STL returns the CFD cells themselves
    case 'stl', tgtFile = stlFile;
    otherwise, error('exportAt must be ''cfd'' or ''stl''');
end
nF = [];
for c = 1:nC
    fprintf('\nMapping %s (%s) ...\n', cs.name{c}, cs.vtk{c});
    m = map_cfd_to_stl(cs.vtkPath{c}, tgtFile, mapOpt);
    if isempty(nF)
        nF = numel(m.P);  XYZ = m.stlCentroids;
        H = zeros(nF, nC);  TAW = H;  P = H;
    elseif numel(m.P) ~= nF || max(abs(m.stlCentroids - XYZ), [], 'all') > 1e-9
        error(['%s is on a different CFD surface mesh than %s, so the cases cannot share one point set.\n' ...
               'Set exportAt = ''stl''.'], cs.vtk{c}, cs.vtk{1});
    end
    H(:,c) = m.h;  TAW(:,c) = m.Taw;  P(:,c) = m.P;
    fprintf('  h median %.0f W/m^2K (max %.0f), T_aw median %.0f K, P %.3g-%.3g Pa, masked faces %d, mapping distance 95%% %.2f mm\n', ...
        median(m.h, 'omitnan'), max(m.h), median(m.Taw, 'omitnan'), min(m.P), max(m.P), nnz(isnan(m.h)), 1e3*prctile(m.dNN, 95));
end
H(isnan(H)) = 0;  TAW(isnan(TAW)) = 0;   % ANSYS cannot take NaN; h = 0 -> no heating
if strcmp(exportAt, 'cfd')
    % CFD points are in the solver frame; ANSYS uses the geometry (STL/STEP) frame
    shift = cfdToGeometryShift(XYZ, stlFile, mapOpt.stlUnits);
    XYZ = XYZ - shift;
    fprintf('\nExport points: %d CFD cell centroids, moved into the geometry frame by -[%.5f %.5f %.5f] m\n', nF, shift);
end

%% ---------------- Scale factors ----------------
sc = scales(fsT, cs, aTab, hScaleType, recovery, gamma);
pOffset = fsT.p - sc.p .* cs.p_inf(aTab) - pRef;

%% ---------------- Thin the time table ----------------
if isempty(tStep)
    nFull = numel(tt);
    [keep, err] = thinTable(tt, aTab, log([sc.h sc.Taw sc.p]), log(1 + tableTol), pOffset, pOffTol);
    tt = tt(keep);  aTab = aTab(keep);  pOffset = pOffset(keep);
    sc.h = sc.h(keep);  sc.Taw = sc.Taw(keep);  sc.p = sc.p(keep);
    fsT.M = fsT.M(keep);  fsT.AoA = fsT.AoA(keep);  fsT.alt = fsT.alt(keep);
    fprintf('\nTime table: %d of %d trajectory points kept (tolerance %.1f%%, %d case switches)\n', ...
        numel(keep), nFull, 100*tableTol, nnz(diff(aTab) ~= 0));
    fprintf('  worst linear-interpolation error at a dropped point: hScale %.2f%%, TawScale %.2f%%, pScale %.2f%%, pOffset %.0f Pa\n', ...
        100*(exp(err(1:3)) - 1), err(4));
end

%% ---------------- Coverage report ----------------
dM = cs.Mach(aTraj) - fs.M;  dA = cs.AoA(aTraj) - fs.AoA;
scTraj = scales(fs, cs, aTraj, hScaleType, recovery, gamma);
fprintf('\n--- Coverage (trajectory vs assigned case) ---\n');
fprintf('%-10s %16s %8s | %13s %15s | %9s %9s | %s\n', 'case', 'time windows [s]', 'tot [s]', 'Mach covered', 'AoA covered', 'max|dM|', 'max|dAoA|', 'h scale');
for c = 1:nC
    k = aTraj == c;
    if ~any(k), fprintf('%-10s   (never used)\n', cs.name{c}); continue; end
    [r0, r1] = runs(k);
    win = strjoin(arrayfun(@(a, b) sprintf('%.0f-%.0f', traj.t(a), traj.t(min(b+1, end))), r0, r1, 'UniformOutput', false), ', ');
    dur = sum(traj.t(min(r1+1, end)) - traj.t(r0));
    fprintf('%-10s %16s %8.0f | %5.2f - %5.2f %+6.2f - %+6.2f | %9.2f %9.2f | %.2f-%.2f\n', cs.name{c}, win, dur, ...
        min(fs.M(k)), max(fs.M(k)), min(fs.AoA(k)), max(fs.AoA(k)), max(abs(dM(k))), max(abs(dA(k))), min(scTraj.h(k)), max(scTraj.h(k)));
end
gap = abs(dM) > dMach | abs(dA) > dAlpha;
if any(gap)
    [g0, g1] = runs(gap);
    fprintf('\nCoverage gaps (|dMach| > %.2f or |dAoA| > %.1f deg): %.0f s of %.0f s\n', dMach, dAlpha, ...
        sum(traj.t(min(g1+1, end)) - traj.t(g0)), traj.t(end) - traj.t(1));
    for i = 1:numel(g0)
        k = g0(i):g1(i);
        fprintf('  t %6.0f-%6.0f s  Mach %.2f-%.2f  AoA %+.2f to %+.2f  -> %-10s  max|dM| %.2f  max|dAoA| %.2f\n', traj.t(g0(i)), traj.t(min(g1(i)+1, end)), ...
            min(fs.M(k)), max(fs.M(k)), min(fs.AoA(k)), max(fs.AoA(k)), cs.name{mode(aTraj(k))}, max(abs(dM(k))), max(abs(dA(k))));
    end
else
    fprintf('\nNo coverage gaps beyond dMach %.2f / dAoA %.1f deg.\n', dMach, dAlpha);
end

%% ---------------- Export ----------------
csvF = fullfile(outDir, [runName '_faces.csv']);
csvP = fullfile(outDir, [runName '_pressure.csv']);
csvT = fullfile(outDir, [runName '_time_table.csv']);
csvC = fullfile(outDir, [runName '_coverage.csv']);
% Fixed precision instead of writematrix's 15 digits: 1 um in X/Y/Z, 6 significant
% figures in h and P (1e-6 relative), 0.01 K in T_aw. All far below the mapping
% and model errors, and the files come out ~3x smaller.
hdrF = [{'FaceID','X','Y','Z'}, reshape([strcat('h_', cs.name), strcat('Taw_', cs.name)]', 1, [])];
datF = zeros(nF, 4 + 2*nC);
datF(:, 1:4) = [(1:nF)', XYZ];
datF(:, 5:2:end) = H;  datF(:, 6:2:end) = TAW;
writeCompact(csvF, hdrF, datF, ['%d,%.6f,%.6f,%.6f', repmat(',%.6g,%.2f', 1, nC), '\n']);
writeCompact(csvP, [{'FaceID','X','Y','Z'}, strcat('P_', cs.name)'], [(1:nF)', XYZ, P], ...
    ['%d,%.6f,%.6f,%.6f', repmat(',%.6g', 1, nC), '\n']);
tab = table(tt, aTab, cs.name(aTab), strcat('h_', cs.name(aTab)), sc.h, strcat('Taw_', cs.name(aTab)), sc.Taw, ...
    strcat('P_', cs.name(aTab)), sc.p, pOffset, fsT.M, fsT.AoA, fsT.alt, ...
    'VariableNames', {'Time_s','CaseIdx','Case','hColumn','hScale','TawColumn','TawScale','PColumn','pScale','pOffset_Pa', ...
    'Mach','Alpha_deg','Altitude_m'});
writetable(tab, csvT);
cov = table(traj.t, cs.name(aTraj), fs.M, cs.Mach(aTraj), dM, fs.AoA, cs.AoA(aTraj), dA, fs.alt, 1e3*cs.alt_km(aTraj), scTraj.h, scTraj.Taw, scTraj.p, gap, ...
    'VariableNames', {'Time_s','Case','Mach','Mach_case','dMach','Alpha_deg','Alpha_case','dAlpha','Altitude_m','Altitude_case_m', ...
    'hScale','TawScale','pScale','gap'});
writetable(cov, csvC);
fprintf('\nWrote:\n  %s  (%d points x %d cases, h and T_aw)\n  %s  (pressure)\n  %s  (%d rows)\n  %s\n', ...
    csvF, nF, nC, csvP, csvT, numel(tt), csvC);

%% ---------------- Plots ----------------
if makePlots
    figure('Name', 'ANSYS load schedule', 'Position', [100 100 950 900]);
    tl = tiledlayout(4, 1, 'TileSpacing', 'compact');
    title(tl, sprintf('%s: %d cases along the trajectory', runName, nC), 'Interpreter', 'none');
    col = lines(nC);
    nexttile; hold on;
    plot(traj.t, fs.M, 'k-', 'LineWidth', 1.3);
    for c = 1:nC, k = aTraj == c; plot(traj.t(k), cs.Mach(c)*ones(nnz(k),1), '.', 'Color', col(c,:), 'MarkerSize', 8); end
    ylabel('Mach'); grid on; box on;
    legend([{'trajectory'}; cs.name], 'Location', 'eastoutside', 'Interpreter', 'none');
    nexttile; hold on;
    plot(traj.t, fs.AoA, 'k-', 'LineWidth', 1.3);
    for c = 1:nC, k = aTraj == c; plot(traj.t(k), cs.AoA(c)*ones(nnz(k),1), '.', 'Color', col(c,:), 'MarkerSize', 8); end
    ylabel('AoA [deg]'); grid on; box on;
    legend([{'trajectory'}; cs.name], 'Location', 'eastoutside', 'Interpreter', 'none');
    nexttile; hold on;
    plot(traj.t, scTraj.h, traj.t, scTraj.Taw, traj.t, scTraj.p, 'LineWidth', 1.1);
    set(gca, 'YScale', 'log'); ylabel('scale vs case'); grid on; box on;
    legend({'hScale', 'TawScale', 'pScale'}, 'Location', 'eastoutside');
    nexttile; hold on;
    plot(traj.t, dM/dMach, traj.t, dA/dAlpha, 'LineWidth', 1.1);
    yline([-1 1], 'k--');
    patch([traj.t; flipud(traj.t)], [3*gap - 1.5; -1.5*ones(size(gap))], 'r', 'FaceAlpha', 0.08, 'EdgeColor', 'none');
    ylabel('mismatch / tolerance'); xlabel('Time [s]'); grid on; box on;
    legend({sprintf('\\DeltaMach / %.2f', dMach), sprintf('\\DeltaAoA / %.1f deg', dAlpha)}, 'Location', 'eastoutside');
end

%% ---------------- Local functions ----------------
function f = freestream(traj, t, g, R, r)
f.t   = t;
f.alt = interp1(traj.t, traj.alt, t);
f.M   = interp1(traj.t, traj.M, t);
f.AoA = interp1(traj.t, traj.AoA, t);
[f.T, f.p, ~] = isa1976(f.alt);
f.rho = f.p ./ (R*f.T);
f.V   = f.M .* sqrt(g*R*f.T);
f.Taw = f.T .* (1 + r*(g-1)/2 .* f.M.^2);
f.q   = 0.5*f.rho.*f.V.^2;
end

function a = assignCases(f, cs, mode, tStart, dM, dA, dLR)
switch mode
    case 'nearest'
        D = ((f.M - cs.Mach')/dM).^2 + ((f.AoA - cs.AoA')/dA).^2 + ((log(f.rho) - log(cs.rho'))/dLR).^2;
        [~, a] = min(D, [], 2);
    case 'windows'
        assert(numel(tStart) == height(cs), 'caseStart_s needs one start time per case');
        [ts, ord] = sort(tStart(:));
        k = discretize(f.t, [ts; inf]);
        k(isnan(k)) = 1;                      % before the first start -> first case
        a = ord(k);
    otherwise
        error('assignMode must be ''nearest'' or ''windows''');
end
end

function s = scales(f, cs, a, hType, r, g)
switch hType
    case 'acreage', s.h = (f.rho.*f.V).^0.8 ./ (cs.rho(a).*cs.V(a)).^0.8;
    case 'stag',    s.h = sqrt(f.rho).*f.V ./ (sqrt(cs.rho(a)).*cs.V(a));
    otherwise, error('hScaleType must be ''acreage'' or ''stag''');
end
TawCase = cs.T_inf(a) .* (1 + r*(g-1)/2 .* cs.Mach(a).^2);
s.Taw = f.Taw ./ TawCase;
s.p   = f.q ./ (0.5*cs.rho(a).*cs.V(a).^2);
end

function shift = cfdToGeometryShift(cfdPts, stlFile, units)
% Rigid translation from the CFD solver frame to the STL/STEP frame (same body,
% different origin): bbox centres, then translation-only point-to-plane ICP.
TR = stlread(stlFile);
V = TR.Points;  F = TR.ConnectivityList;
if strcmpi(units, 'mm'), V = V / 1000; end
C = (V(F(:,1),:) + V(F(:,2),:) + V(F(:,3),:)) / 3;
N = cross(V(F(:,2),:) - V(F(:,1),:), V(F(:,3),:) - V(F(:,1),:), 2);
N = N ./ max(vecnorm(N, 2, 2), eps);
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

function [keep, err] = thinTable(t, a, L, tolL, pOff, pOffTol)
% Fewest rows such that linear interpolation between kept rows reproduces every
% dropped row within tolerance: log scale factors within tolL, pOffset within
% pOffTol. Rows never span a case switch, so the last row of one case and the
% first row of the next are always kept (same ramp between them as the full table).
n = numel(t);
keep = 1;  i = 1;
err = zeros(1, 4);
while i < n
    j = i + 1;                             % always reachable
    while j < n && a(j+1) == a(i)
        [ok, ~] = segErr(t, L, pOff, i, j + 1, tolL, pOffTol);
        if ~ok, break; end
        j = j + 1;
    end
    [~, e] = segErr(t, L, pOff, i, j, tolL, pOffTol);
    err = max(err, e);
    keep(end+1, 1) = j; %#ok<AGROW>
    i = j;
end
end

function [ok, e] = segErr(t, L, pOff, i, j, tolL, pOffTol)
m = (i+1:j-1)';
if isempty(m), ok = true; e = zeros(1, 4); return; end
w = (t(m) - t(i)) / (t(j) - t(i));
eL = max(abs(L(i,:) + w .* (L(j,:) - L(i,:)) - L(m,:)), [], 1);
eP = max(abs(pOff(i) + w .* (pOff(j) - pOff(i)) - pOff(m)));
e  = [eL, eP];
ok = all(eL <= tolL) && eP <= pOffTol;
end

function writeCompact(file, hdr, dat, fmt)
fid = fopen(file, 'w');
assert(fid > 0, 'cannot open %s for writing (open in another program?)', file);
fprintf(fid, '%s\n', strjoin(hdr, ','));
fprintf(fid, fmt, dat.');
fclose(fid);
end

function [r0, r1] = runs(k)
% start/end indices of consecutive true runs
d  = diff([false; k(:); false]);
r0 = find(d == 1);  r1 = find(d == -1) - 1;
end

function out = iff(c, a, b)
if c, out = a; else, out = b; end
end
