%% plan_cfd_cases.m
% Pick the CFD cases (one surface VTK each) needed to cover a trajectory for
% transient thermal sims, and the time window each case covers.
%
% The trajectory is split into consecutive time segments inside which Mach,
% angle of attack and the heating level each stay within a tolerance band.
% Each segment gets ONE steady CFD run at its anchor point (Mach, AoA,
% altitude -> freestream p, T, rho). The thermal sim applies that case's loads
% for the segment's time window, scaled in time by the factors in the schedule.
%
% Outputs (in Exports/):
%   <traj>_cfd_cases.csv      one row per CFD run: conditions to run + time window it covers
%   <traj>_load_schedule.csv  per trajectory point: case ID and load scale factors vs the case anchor
%
% Using the schedule in ANSYS (convection load from export_cfd_heating.m):
%   h(x,t)    = h_case(x)    * hScale(t)      (hScaleAcreage on the body, hScaleStag on nose/leading edges)
%   T_aw(x,t) = T_aw_case(x) * TawScale(t)
%   p(x,t)    = p_inf(t) + (p_case(x) - p_inf_case) * pScale(t)
% These are first-order scalings (h ~ (rho V)^0.8 turbulent, ~ sqrt(rho) V stagnation;
% surface pressure rise ~ dynamic pressure at fixed Mach and AoA). The tolerances
% keep them small corrections; tighten the tolerances if they matter.
clear; clc; close all;
here = fileparts(mfilename('fullpath'));

%% ---------------- User inputs ----------------
trajFile  = fullfile(here, 'Calls', 'Trajectory_WAV_RID_9_aero.csv');
nCases    = [5];        % [] -> as many cases as the tolerances need; N -> best split into at most N cases
dMach     = 0.5;       % max Mach spread within one case
dAlpha    = 1.0;       % [deg] max angle-of-attack spread within one case
heatTol   = 1.0;       % max heating spread within one case (1.0 -> max/min <= 2). The schedule's
                       % scale factors correct the heating LEVEL, so this can be looser than the
                       % Mach/AoA bands, which change the flow pattern and need a new CFD run
machMin   = 0;         % skip trajectory points below this Mach (no CFD case assigned)
anchor    = 'center';  % 'center': point nearest the middle of the segment's band
                       % 'peak'  : point of highest heating (conservative if loads are not scaled)
gamma     = 1.4;
Rgas      = 287.05;
recovery  = 0.71^(1/3);% turbulent recovery factor, for T_aw (same as export_cfd_heating.m)
makePlots = true;

outDir = fullfile(here, 'Exports');
if ~exist(outDir, 'dir'), mkdir(outDir); end
[~, trajName] = fileparts(trajFile);

%% ---------------- Load trajectory and freestream ----------------
tr    = readtable(trajFile);
t     = tr.Time_s;
alt   = max(tr.Altitude_m, 0);
M     = tr.Mach;
alpha = tr.Alpha_deg;
V     = tr.Velocity_ms;
n     = numel(t);

[T_inf, p_inf, rho] = isa1976(alt);
mu    = 1.458e-6 * T_inf.^1.5 ./ (T_inf + 110.4);   % Sutherland
Re_m  = rho .* V ./ mu;
qdyn  = 0.5 * rho .* V.^2;
T_aw  = T_inf .* (1 + recovery*(gamma-1)/2 .* M.^2);

% heating indicators (radius and wall temperature cancel in ratios)
qStag = sqrt(rho) .* V.^3;          % stagnation point / leading edge (Sutton-Graves form)
qAcre = (rho .* V).^0.8 .* V.^2;    % turbulent acreage (St ~ Re^-0.2, driving enthalpy ~ V^2)

% sanity: trajectory Mach vs V and the standard atmosphere
M_isa = V ./ sqrt(gamma*Rgas*T_inf);
fprintf('Trajectory: %d points, t = %.1f-%.1f s (dt %.3f s)\n', n, t(1), t(end), median(diff(t)));
fprintf('Mach %.2f-%.2f, AoA %+.2f to %+.2f deg, altitude %.0f-%.0f m\n', min(M), max(M), min(alpha), max(alpha), min(alt), max(alt));
fprintf('Mach check vs V/a(ISA): max difference %.3f\n', max(abs(M_isa - M)));

%% ---------------- Segment the trajectory ----------------
keep = M >= machMin;
lnTol = log(1 + heatTol);
X = [M/dMach, alpha/dAlpha, log(qStag)/lnTol, log(qAcre)/lnTol];   % 1 unit = one full tolerance band

if isempty(nCases)
    s = 1;
    seg = greedySegments(X, keep, s);
else
    % smallest uniform scaling of all tolerances that needs at most nCases segments
    lo = 0;  hi = max(max(X(keep,:)) - min(X(keep,:))) + 1;
    for it = 1:60
        mid = (lo + hi)/2;
        if size(greedySegments(X, keep, mid), 1) <= nCases, hi = mid; else, lo = mid; end
    end
    s = hi;
    seg = greedySegments(X, keep, s);
    fprintf('\nTarget %d cases -> effective tolerances: dMach %.2f, dAlpha %.2f deg, heating x%.2f\n', ...
        nCases, s*dMach, s*dAlpha, (1+heatTol)^s);
end
nSeg = size(seg, 1);

%% ---------------- Anchors, time windows, scale factors ----------------
caseID = zeros(n, 1);
k_a    = zeros(nSeg, 1);
tStart = zeros(nSeg, 1);  tEnd = zeros(nSeg, 1);
for c = 1:nSeg
    i = seg(c,1);  j = seg(c,2);
    caseID(i:j) = c;
    Xs = X(i:j, :);
    switch anchor
        case 'center'
            mid = (max(Xs) + min(Xs))/2;
            [~, k] = min(max(abs(Xs - mid), [], 2));
        case 'peak'
            [~, k] = max(qStag(i:j)/max(qStag) + qAcre(i:j)/max(qAcre));
        otherwise
            error('anchor must be ''center'' or ''peak''');
    end
    k_a(c) = i + k - 1;
    % windows abut: each runs until the next segment starts (or to its own last point)
    tStart(c) = t(i);
    if j < n && keep(j+1), tEnd(c) = t(j+1); else, tEnd(c) = t(j); end
end
dur = tEnd - tStart;

ka = zeros(n, 1);  ka(caseID > 0) = k_a(caseID(caseID > 0));
sc = @(f) iff(caseID > 0, f ./ f(max(ka, 1)), NaN);
hScaleAcre = sc((rho .* V).^0.8);
hScaleStag = sc(sqrt(rho) .* V);
TawScale   = sc(T_aw);
pScale     = sc(qdyn);

%% ---------------- Report ----------------
fprintf('\n%d CFD cases (anchor = %s)\n', nSeg, anchor);
fprintf('%4s %8s %8s %7s | %7s %6s %7s %8s | %9s %7s %9s %9s %7s | %s\n', 'case', 't0[s]', 't1[s]', 'dur[s]', ...
    'Mach', 'AoA', 'alt[km]', 'V[m/s]', 'p_inf[Pa]', 'T_inf', 'rho', 'Re/m', 'T_aw', 'h scale in window (acreage)');
for c = 1:nSeg
    k = k_a(c);  r = seg(c,1):seg(c,2);
    fprintf('%4d %8.1f %8.1f %7.1f | %7.2f %+6.2f %7.2f %8.0f | %9.4g %7.1f %9.4g %9.3g %7.0f | %.2f-%.2f\n', c, tStart(c), tEnd(c), dur(c), ...
        M(k), alpha(k), alt(k)/1e3, V(k), p_inf(k), T_inf(k), rho(k), Re_m(k), T_aw(k), min(hScaleAcre(r)), max(hScaleAcre(r)));
end
if any(~keep), fprintf('%d trajectory points below Mach %.2f have no case.\n', nnz(~keep), machMin); end
fprintf('Covered time: %.1f s of %.1f s\n', sum(dur), t(end) - t(1));

%% ---------------- Export ----------------
csvCases = fullfile(outDir, [trajName '_cfd_cases.csv']);
csvSched = fullfile(outDir, [trajName '_load_schedule.csv']);
cases = table((1:nSeg)', tStart, tEnd, dur, t(k_a), M(k_a), alpha(k_a), alt(k_a), V(k_a), p_inf(k_a), T_inf(k_a), rho(k_a), Re_m(k_a), qdyn(k_a), T_aw(k_a), ...
    arrayfun(@(c) min(M(seg(c,1):seg(c,2))), (1:nSeg)'), arrayfun(@(c) max(M(seg(c,1):seg(c,2))), (1:nSeg)'), ...
    arrayfun(@(c) min(alpha(seg(c,1):seg(c,2))), (1:nSeg)'), arrayfun(@(c) max(alpha(seg(c,1):seg(c,2))), (1:nSeg)'), ...
    'VariableNames', {'Case','t_start_s','t_end_s','duration_s','t_anchor_s','Mach','Alpha_deg','Altitude_m','Velocity_ms', ...
    'p_inf_Pa','T_inf_K','rho_kgm3','Re_per_m','q_dyn_Pa','T_aw_K','Mach_min','Mach_max','Alpha_min_deg','Alpha_max_deg'});
writetable(cases, csvCases);
sched = table(t, caseID, M, alpha, alt, p_inf, T_inf, rho, hScaleAcre, hScaleStag, TawScale, pScale, ...
    'VariableNames', {'Time_s','Case','Mach','Alpha_deg','Altitude_m','p_inf_Pa','T_inf_K','rho_kgm3', ...
    'hScaleAcreage','hScaleStag','TawScale','pScale'});
writetable(sched, csvSched);
fprintf('\nWrote:\n  %s\n  %s\n', csvCases, csvSched);

%% ---------------- Plots ----------------
if makePlots
    figure('Name', 'CFD case plan', 'Position', [100 100 900 900]);
    ys = {M, 'Mach'; alpha, 'AoA [deg]'; alt/1e3, 'Altitude [km]'; [qStag/max(qStag), qAcre/max(qAcre)], 'Heating (norm.)'};
    for p = 1:4
        ax = subplot(4, 1, p); hold on;
        y = ys{p,1};
        yl = [min(y(:)) max(y(:))];  yl = yl + [-1 1]*0.05*max(diff(yl), eps);
        for c = 1:nSeg
            patch([tStart(c) tEnd(c) tEnd(c) tStart(c)], yl([1 1 2 2]), 0.85 + 0.1*mod(c,2)*[1 1 1], 'EdgeColor', 'none');
        end
        hl = plot(t, y, 'LineWidth', 1.2);
        plot(t(k_a), y(k_a,:), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
        ylim(yl); ylabel(ys{p,2}); grid on; box on; ax.Layer = 'top';
        if p == 4
            % launch heating is ~50x cruise; log scale keeps cruise visible
            yl = [min(y(:))/1.5 1.5];
            set(ax, 'YScale', 'log');  ylim(yl);
            set(findobj(ax, 'Type', 'patch'), 'YData', yl([1 1 2 2])');
            legend(hl, {'stagnation', 'acreage'}, 'Location', 'northeast'); xlabel('Time [s]');
        end
        if p == 1, title(sprintf('%d CFD cases (shaded = time window, dots = run conditions)', nSeg), 'Interpreter', 'none'); end
    end
end

%% ---------------- Local functions ----------------
function seg = greedySegments(X, keep, s)
% Longest-first split: extend each segment while every column of X stays within
% a band of width s. Greedy gives the fewest segments for this constraint.
n = size(X, 1);
seg = zeros(0, 2);
i = 1;
while i <= n
    if ~keep(i), i = i + 1; continue; end
    lo = X(i,:);  hi = X(i,:);  j = i;
    while j < n && keep(j+1)
        lo2 = min(lo, X(j+1,:));  hi2 = max(hi, X(j+1,:));
        if any(hi2 - lo2 > s), break; end
        lo = lo2;  hi = hi2;  j = j + 1;
    end
    seg(end+1, :) = [i j]; %#ok<AGROW>
    i = j + 1;
end
end

function [T, p, rho] = isa1976(z)
% US Standard Atmosphere 1976, 0-71 km. z = geometric altitude [m].
r0 = 6356766;  H = r0*z ./ (r0 + z);                     % geopotential altitude
Hb = [0 11000 20000 32000 47000 51000 71000];
Tb = [288.15 216.65 216.65 228.65 270.65 270.65 214.65];
Lb = [-0.0065 0 0.001 0.0028 0 -0.0028 -0.002];
g0 = 9.80665;  R = 287.05287;
pb = zeros(size(Hb));  pb(1) = 101325;
for b = 1:numel(Hb)-1
    dH = Hb(b+1) - Hb(b);
    if Lb(b) == 0, pb(b+1) = pb(b)*exp(-g0*dH/(R*Tb(b)));
    else,          pb(b+1) = pb(b)*(Tb(b+1)/Tb(b))^(-g0/(R*Lb(b))); end
end
b = discretize(H, [Hb inf]);
assert(all(~isnan(b)), 'altitude outside 0-71 km');
T = Tb(b)' + Lb(b)'.*(H - Hb(b)');
p = zeros(size(H));
iso = Lb(b)' == 0;
p(iso)  = pb(b(iso))' .* exp(-g0*(H(iso) - Hb(b(iso))')./(R*Tb(b(iso))'));
p(~iso) = pb(b(~iso))' .* (T(~iso)./Tb(b(~iso))').^(-g0./(R*Lb(b(~iso))'));
rho = p ./ (R*T);
end

function out = iff(cond, a, b)
out = repmat(b, size(cond));
out(cond) = a(cond);
end
