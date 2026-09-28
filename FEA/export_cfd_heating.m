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
vtkFile   = 'C:\Users\kkdar\Downloads\surf_out_000000001 (1).vtk';
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

%% ---------------- Load CFD surface ----------------
[T_e, P_e, U_e, V_e, W_e, centroids, norms, areas, ex] = load_vtk_surf(vtkFile);
nCFD = numel(T_e);
assert(isfield(ex, 'qw'), 'VTK has no qw field. Use exportpressureandheat_working.m for files without solver heating.');
assert(any(ex.qw ~= 0), 'qw is all zeros in this VTK. Use exportpressureandheat_working.m, or rerun CHAMPS with heating output.');
qw_cfd = ex.qw;
if isfield(ex, 'wall_temperature')
    Tw_cfd = ex.wall_temperature;
else
    error('VTK has no wall_temperature field; cannot convert qw to a film coefficient.');
end
patch_id = [];
if isfield(ex, 'patch_id'), patch_id = ex.patch_id; end

%% ---------------- Adiabatic wall temperature and film coefficient ----------------
v_mag = sqrt(U_e.^2 + V_e.^2 + W_e.^2);
M_e   = v_mag ./ sqrt(gamma*Rgas*max(T_e, 1));
if ischar(recovery) || isstring(recovery)
    switch lower(string(recovery))
        case "turbulent", r = pran^(1/3);
        case "laminar",   r = sqrt(pran);
        otherwise, error('recovery must be ''turbulent'', ''laminar'' or a number');
    end
else
    r = recovery;
end
T_aw = T_e .* (1 + r*(gamma-1)/2 .* M_e.^2);
dT   = T_aw - Tw_cfd;

% Masks: solver holes (qw==0 & Tw==0), and cells where T_aw ~ T_wall (h ill-conditioned)
hole    = (qw_cfd == 0) & (Tw_cfd == 0);
illcond = ~hole & (dT < dT_min);
invalid = hole | illcond;

h = qw_cfd ./ dT;
h(invalid) = NaN;
q_cold = h .* (T_aw - T_cold);        % what a cold structure would see at t=0 (diagnostic)

%% ---------------- Diagnostics ----------------
fprintf('CFD surface: %d cells, extent x[%.3f,%.3f] y[%.3f,%.3f] z[%.3f,%.3f] m\n', nCFD, ...
    min(centroids(:,1)), max(centroids(:,1)), min(centroids(:,2)), max(centroids(:,2)), min(centroids(:,3)), max(centroids(:,3)));
fprintf('Flow direction (mean velocity): U=%+.0f V=%+.0f W=%+.0f m/s\n', mean(U_e), mean(V_e), mean(W_e));
fprintf('Edge state: T_e [%.0f, %.0f] K, M_e [%.2f, %.2f], P [%.3g, %.3g] Pa\n', min(T_e), max(T_e), min(M_e), max(M_e), min(P_e), max(P_e));
fprintf('Recovery factor r = %.3f (%s)\n', r, string(recovery));
fprintf('T_aw   [%.0f, %.0f] K;  CFD wall temperature [%.0f, %.0f] K (mean %.0f)\n', min(T_aw), max(T_aw), min(Tw_cfd(~hole)), max(Tw_cfd), mean(Tw_cfd(~hole)));
sig = 5.670374e-8; eps_impl = qw_cfd(~invalid) ./ (sig*Tw_cfd(~invalid).^4);
if (prctile(eps_impl,95) - prctile(eps_impl,5)) < 0.02
    fprintf('Wall BC looks like RADIATIVE EQUILIBRIUM with emissivity %.3f (qw = eps*sigma*Tw^4 in every cell)\n', median(eps_impl));
else
    fprintf('Wall BC: not radiative equilibrium (implied emissivity spread %.3f..%.3f)\n', prctile(eps_impl,5), prctile(eps_impl,95));
end
fprintf('\n--- Heating ---\n');
fprintf('qw (CFD, at CFD wall temp): [%.3g, %.3g] W/m^2, mean %.3g\n', min(qw_cfd), max(qw_cfd), mean(qw_cfd));
fprintf('h (film coeff):             median %.0f, 5-95%% [%.0f, %.0f] W/m^2K, max %.0f\n', ...
    median(h,'omitnan'), prctile(h,5), prctile(h,95), max(h));
fprintf('Cold-wall equivalent q at %d K: max %.3g W/m^2, mean %.3g  (this is the t=0 load)\n', T_cold, max(q_cold), mean(q_cold,'omitnan'));
fprintf('\n--- Mask accounting (CFD cells) ---\n');
fprintf('solver holes (qw=0 & Tw=0):      %d\n', nnz(hole));
fprintf('T_aw - T_wall < %d K (h masked): %d\n', dT_min, nnz(illcond));
fprintf('total invalid:                   %d of %d (%.2f%%)\n', nnz(invalid), nCFD, 100*nnz(invalid)/nCFD);
if ~isempty(patch_id)
    ids = unique(patch_id);
    for k = 1:numel(ids)
        m = patch_id == ids(k);
        fprintf('patch %d: %6d cells, qw mean %.3g, Tw mean %.0f, x[%.3f,%.3f]\n', ids(k), nnz(m), mean(qw_cfd(m)), mean(Tw_cfd(m)), min(centroids(m,1)), max(centroids(m,1)));
    end
end
[~, kpk] = max(qw_cfd);
fprintf('peak qw at [%.4f %.4f %.4f] m: qw=%.3g, Tw=%.0f, T_aw=%.0f, h=%.0f\n', centroids(kpk,:), qw_cfd(kpk), Tw_cfd(kpk), T_aw(kpk), h(kpk));

%% ---------------- Target faces: STL or CFD centroids ----------------
useSTL = ~isempty(stlFile);
if useSTL
    TR = stlread(stlFile);
    stlFaces    = TR.ConnectivityList;
    stlVertices = TR.Points;
    if strcmpi(stlUnits, 'mm'), stlVertices = stlVertices / 1000; end
    stlCentroids = (stlVertices(stlFaces(:,1),:) + stlVertices(stlFaces(:,2),:) + stlVertices(stlFaces(:,3),:)) / 3;
    fprintf('\nSTL: %d faces, extent (m) x[%.3f,%.3f] y[%.3f,%.3f] z[%.3f,%.3f]\n', size(stlFaces,1), ...
        min(stlCentroids(:,1)), max(stlCentroids(:,1)), min(stlCentroids(:,2)), max(stlCentroids(:,2)), min(stlCentroids(:,3)), max(stlCentroids(:,3)));
else
    stlCentroids = centroids;
    fprintf('\nNo STL given: exporting at %d CFD cell centroids.\n', nCFD);
end

%% ---------------- Mirror half-model CFD if needed ----------------
src   = centroids;  srcP = P_e;  srcH = h;  srcTaw = T_aw;  srcQ = qw_cfd;  srcTw = Tw_cfd;  srcN = norms;
doMirror = false;
if ~isempty(symAxis)
    ax = find('xyz' == lower(symAxis));
    cfdOneSided = min(centroids(:,ax)) >= -1e-6 || max(centroids(:,ax)) <= 1e-6;
    tgtTwoSided = min(stlCentroids(:,ax)) < -1e-3 && max(stlCentroids(:,ax)) > 1e-3;
    % a half model spans ~half the target's width on this axis; a shifted full model spans all of it
    cfdHalfWidth = range(centroids(:,ax)) < 0.75*range(stlCentroids(:,ax));
    if islogical(mirrorCFD), doMirror = mirrorCFD; else, doMirror = cfdOneSided && tgtTwoSided && cfdHalfWidth; end
    if doMirror
        mir = centroids; mir(:,ax) = -mir(:,ax);
        nmir = norms; nmir(:,ax) = -nmir(:,ax); srcN = [norms; nmir];
        src = [centroids; mir]; srcP = [P_e; P_e]; srcH = [h; h]; srcTaw = [T_aw; T_aw]; srcQ = [qw_cfd; qw_cfd]; srcTw = [Tw_cfd; Tw_cfd];
        fprintf('Mirrored CFD half-model across %s = 0 (%d -> %d source cells)\n', symAxis, nCFD, size(src,1));
    end
end

%% ---------------- Align STL to CFD (optional) and map ----------------
tgt = stlCentroids;
if useSTL && strcmpi(alignMode, 'bbox')
    rangeCFD = max(src) - min(src);  rangeSTL = max(tgt) - min(tgt);
    sf = rangeCFD ./ rangeSTL;
    tgt = tgt .* sf;
    tgt = tgt - mean(tgt,1) + mean(src,1);
    fprintf('bbox alignment scale factors [%.3f %.3f %.3f]  (far from 1 => unit or geometry mismatch)\n', sf);
    if any(abs(sf-1) > 0.05), warning('STL and CFD extents differ by >5%% on some axis. Check stlUnits / geometry.'); end
elseif useSTL && strcmpi(alignMode, 'shift')
    % rigid translation only (same geometry, different origin): start from bbox
    % centres, then refine with translation-only ICP against the CFD centroids
    shift = (max(src)+min(src))/2 - (max(tgt)+min(tgt))/2;
    for it = 1:15
        j = knnsearch(src, tgt + shift);
        step = median(src(j,:) - (tgt + shift), 1);
        shift = shift + step;
        if norm(step) < 1e-6, break; end
    end
    tgt = tgt + shift;
    fprintf('shift alignment: STL translated by [%.4f %.4f %.4f] m (%d ICP iterations)\n', shift, it);
end
idx = knnsearch(src, tgt);
if useSTL
    % Wrong-side guard: near the thin leading edge the nearest CFD centroid can lie
    % on the opposite surface. Require the source cell normal to agree with the
    % target face normal; otherwise take the nearest of 32 candidates that does.
    tgtN = cross(stlVertices(stlFaces(:,2),:) - stlVertices(stlFaces(:,1),:), ...
                 stlVertices(stlFaces(:,3),:) - stlVertices(stlFaces(:,1),:), 2);
    tgtN = tgtN ./ max(vecnorm(tgtN, 2, 2), eps);
    if median(sum(srcN(idx,:).*tgtN, 2)) < 0, srcN = -srcN; end   % CFD normals point inward
    bad = sum(srcN(idx,:).*tgtN, 2) < 0.2;
    nBad0 = nnz(bad);
    if nBad0
        cand = knnsearch(src, tgt(bad,:), 'K', 32);
        dots = sum(reshape(srcN(cand',:), 32, [], 3) .* permute(tgtN(bad,:), [3 1 2]), 3)';   % [nBad x 32]
        ok = dots >= 0.2;
        [hasOk, first] = max(ok, [], 2);
        bi = find(bad);
        fix = hasOk == 1;
        idx(bi(fix)) = cand(sub2ind(size(cand), find(fix), first(fix)));
        bad(bi(fix)) = false;
    end
    fprintf('normal check: %d faces had an opposite-facing nearest cell, %d remapped, %d unresolved (kept nearest)\n', ...
        nBad0, nBad0 - nnz(bad), nnz(bad));
end
dNN = sqrt(sum((src(idx,:) - tgt).^2, 2));
mappedP   = srcP(idx);
mappedH   = srcH(idx);
mappedTaw = srcTaw(idx);
mappedQ   = srcQ(idx);
mappedTw  = srcTw(idx);
fprintf('nearest-neighbour distance: median %.2f mm, 95%% %.2f mm, max %.2f mm\n', 1e3*median(dNN), 1e3*prctile(dNN,95), 1e3*max(dNN));
if aftMask > 0
    % base face and its rim: h = qw/(T_aw - T_wall) is ill-conditioned there
    % (T_aw barely above T_wall in the base expansion), and the base is not a design load
    aft = stlCentroids(:,1) > max(stlCentroids(:,1)) - aftMask;
    mappedH(aft) = NaN;  mappedTaw(aft) = NaN;  mappedQ(aft) = NaN;
    fprintf('aft mask: %d faces within %.0f mm of the base excluded from heating\n', nnz(aft), 1e3*aftMask);
end
fprintf('target faces with NaN h: %d of %d\n', nnz(isnan(mappedH)), numel(mappedH));

%% ---------------- Export ----------------
FaceID = (1:size(stlCentroids,1))';
XYZ = stlCentroids;                              % metres
csvP = fullfile(outDir, [caseName '_pressure.csv']);
csvC = fullfile(outDir, [caseName '_convection.csv']);
csvQ = fullfile(outDir, [caseName '_heatflux_cfd.csv']);
Hcsv = mappedH;  Hcsv(isnan(Hcsv)) = 0;           % ANSYS tables cannot take NaN; h = 0 -> no heating
mappedTaw(isnan(mappedTaw)) = 0;  mappedQ(isnan(mappedQ)) = 0;
writecell({'FaceID','X','Y','Z','Pressure'}, csvP);
writematrix([FaceID, XYZ, mappedP], csvP, 'WriteMode','append');
writecell({'FaceID','X','Y','Z','h','T_aw'}, csvC);
writematrix([FaceID, XYZ, Hcsv, mappedTaw], csvC, 'WriteMode','append');
writecell({'FaceID','X','Y','Z','qw_cfd','T_wall_cfd'}, csvQ);
writematrix([FaceID, XYZ, mappedQ, mappedTw], csvQ, 'WriteMode','append');
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
