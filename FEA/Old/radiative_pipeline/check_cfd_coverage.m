%% check_cfd_coverage.m
% How much does a CFD case plan (from plan_cfd_cases.m) miss compared with
% flying the real trajectory?
%
% CFD cannot be run at every trajectory point, so a cheap engineering model
% stands in for it: each surface is a wedge (oblique shock, or Prandtl-Meyer
% expansion if alpha turns it leeward) with Eckert reference-temperature
% turbulent flat-plate heating. The model is evaluated
%   truth : at every trajectory point (actual Mach, AoA, altitude)
%   plan  : at each case's anchor, then carried over its time window exactly as
%           ANSYS would use it (h * hScale(t), T_aw * TawScale(t), p rise * pScale(t))
%   raw   : at each case's anchor, held constant over the window (no scaling)
% The plan-vs-truth gap is the discretization error of the case plan. It says
% nothing about the absolute accuracy of the CFD or of this model.
clear; clc; close all;
here = fileparts(mfilename('fullpath'));

%% ---------------- User inputs ----------------
trajName = 'Trajectory_WAV_RID_8_aero';
casesCsv = fullfile(here, 'Exports', [trajName '_cfd_cases.csv']);
schedCsv = fullfile(here, 'Exports', [trajName '_load_schedule.csv']);
% Surface ramp angles at alpha = 0 [deg], from the WAV_RID CFD surface
% (centerline slopes); signAoA = +1 if positive alpha steepens the ramp.
surfaces = struct( ...
    'name',    {'lower',  'upper fwd', 'upper aft'}, ...
    'delta0',  { 2.6,      4.4,         0.0       }, ...
    'signAoA', {+1,       -1,          -1         });
xRun     = 1.0;        % [m] running length for the flat-plate heating
Tw       = 300;        % [K] wall temperature (cold structure = t=0 worst case)
gamma    = 1.4;
Rgas     = 287.05;
cp       = gamma*Rgas/(gamma-1);
Pr       = 0.71;
recovery = Pr^(1/3);
makePlots = true;

%% ---------------- Load plan ----------------
C  = readtable(casesCsv);
S  = readtable(schedCsv);
ok = S.Case > 0;
S  = S(ok, :);
t  = S.Time_s;
a  = S.Case;                                  % case index per time
V  = S.Mach .* sqrt(gamma*Rgas*S.T_inf_K);
Va = C.Velocity_ms;
fprintf('Plan: %d CFD cases covering %.0f s (%d trajectory points)\n', height(C), t(end) - t(1), numel(t));

%% ---------------- Model: truth, plan, raw ----------------
nS = numel(surfaces);  n = numel(t);
q_true = zeros(n, nS);  q_plan = q_true;  q_raw = q_true;
p_true = q_true;        p_plan = q_true;  p_raw = q_true;
for k = 1:nS
    sf = surfaces(k);
    % truth at every point
    for i = 1:n
        [q_true(i,k), ~, ~, p_true(i,k)] = wedgeHeating(S.Mach(i), V(i), S.p_inf_Pa(i), S.T_inf_K(i), ...
            sf.delta0 + sf.signAoA*S.Alpha_deg(i), xRun, Tw, gamma, Rgas, cp, Pr, recovery);
    end
    % at each case anchor
    nC = height(C);  qa = zeros(nC,1);  ha = qa;  Tawa = qa;  pa = qa;
    for c = 1:nC
        [qa(c), ha(c), Tawa(c), pa(c)] = wedgeHeating(C.Mach(c), Va(c), C.p_inf_Pa(c), C.T_inf_K(c), ...
            sf.delta0 + sf.signAoA*C.Alpha_deg(c), xRun, Tw, gamma, Rgas, cp, Pr, recovery);
    end
    % carried over the window the way ANSYS applies it
    q_plan(:,k) = ha(a) .* S.hScaleAcreage .* (Tawa(a) .* S.TawScale - Tw);
    p_plan(:,k) = S.p_inf_Pa + (pa(a) - C.p_inf_Pa(a)) .* S.pScale;
    q_raw(:,k)  = qa(a);
    p_raw(:,k)  = pa(a);
end

%% ---------------- Report ----------------
pct = @(x, ref) 100*(x - ref)./ref;
fprintf('\nMach/AoA mismatch (case anchor vs actual): max |dMach| %.2f, max |dAoA| %.2f deg\n', ...
    max(abs(C.Mach(a) - S.Mach)), max(abs(C.Alpha_deg(a) - S.Alpha_deg)));
fprintf('Model: wedge + Eckert turbulent flat plate, x = %.2f m, T_wall = %.0f K\n\n', xRun, Tw);
fprintf('%-10s | %10s %10s | %12s %16s %16s | %14s %14s\n', 'surface', 'peak q', 'peak q', 'heat load', 'heat load plan', 'heat load raw', 'max |q err|', 'max |p err|');
fprintf('%-10s | %10s %10s | %12s %16s %16s | %14s %14s\n', '', 'true', 'plan', 'true [MJ/m2]', '(err)', '(err, unscaled)', 'plan (at t)', 'plan (at t)');
for k = 1:nS
    Qt = trapz(t, q_true(:,k))/1e6;  Qp = trapz(t, q_plan(:,k))/1e6;  Qr = trapz(t, q_raw(:,k))/1e6;
    [eq, iq] = max(abs(pct(q_plan(:,k), q_true(:,k))));
    [ep, ip] = max(abs(pct(p_plan(:,k), p_true(:,k))));
    fprintf('%-10s | %10.3g %10.3g | %12.2f %9.2f (%+5.1f%%) %9.2f (%+5.1f%%) | %6.1f%% @%5.0fs %6.1f%% @%5.0fs\n', surfaces(k).name, ...
        max(q_true(:,k)), max(q_plan(:,k)), Qt, Qp, pct(Qp,Qt), Qr, pct(Qr,Qt), eq, t(iq), ep, t(ip));
end
fprintf('\nheat load = time-integrated heat flux into a %.0f K wall; err = plan vs truth.\n', Tw);
fprintf('"raw" = anchor loads held constant over each window (what you get without the schedule scaling).\n');

%% ---------------- Plots ----------------
if makePlots
    figure('Name', 'CFD plan coverage', 'Position', [100 100 950 950]);
    tl = tiledlayout(4, 1, 'TileSpacing', 'compact');
    title(tl, sprintf('%d-case plan vs trajectory', height(C)));

    nexttile; hold on;
    yyaxis left;  plot(t, S.Mach, '-');  stairs(t, C.Mach(a), '--');  ylabel('Mach');
    yyaxis right; plot(t, S.Alpha_deg, '-');  stairs(t, C.Alpha_deg(a), '--');  ylabel('AoA [deg]');
    legend({'Mach actual', 'Mach case', 'AoA actual', 'AoA case'}, 'Location', 'northeast');
    grid on; box on;

    nexttile; hold on;
    k = 1;
    plot(t, q_true(:,k), 'k-', 'LineWidth', 1.4);
    plot(t, q_plan(:,k), '-', 'LineWidth', 1.1);
    stairs(t, q_raw(:,k), ':', 'LineWidth', 1.1);
    set(gca, 'YScale', 'log'); ylabel('q [W/m^2]'); grid on; box on;
    legend({'truth', 'plan (scaled)', 'raw (unscaled)'}, 'Location', 'northeast');
    title(sprintf('Heat flux, %s surface (%.0f K wall)', surfaces(k).name, Tw));

    nexttile; hold on;
    for k = 1:nS, plot(t, pct(q_plan(:,k), q_true(:,k)), 'LineWidth', 1.1); end
    yline(0, 'k-'); ylabel('q error [%]'); grid on; box on;
    legend({surfaces.name}, 'Location', 'best'); title('Heat flux error, plan vs truth');

    nexttile; hold on;
    for k = 1:nS, plot(t, pct(p_plan(:,k), p_true(:,k)), 'LineWidth', 1.1); end
    yline(0, 'k-'); ylabel('p error [%]'); xlabel('Time [s]'); grid on; box on;
    legend({surfaces.name}, 'Location', 'best'); title('Surface pressure error, plan vs truth');
end

%% ---------------- Local functions ----------------
function [q, h, Taw, pe] = wedgeHeating(M, V, p, T, delta, x, Tw, g, R, cp, Pr, r)
% Edge state behind a wedge of deflection delta [deg] (expansion if negative),
% then Eckert reference-temperature turbulent flat-plate heating at length x.
d = abs(delta)*pi/180;
if delta > 0
    % weak oblique shock
    f = @(b) tan(d) - 2*cot(b).*(M^2*sin(b).^2 - 1)./(M^2*(g + cos(2*b)) + 2);
    mu = asin(1/M);
    bmax = fminbnd(@(b) -(2*cot(b).*(M^2*sin(b).^2 - 1)./(M^2*(g + cos(2*b)) + 2)), mu, pi/2);
    b = fzero(f, [mu + 1e-9, bmax]);
    Mn = M*sin(b);
    pe = p*(1 + 2*g/(g+1)*(Mn^2 - 1));
    Te = T*pe/p*((g-1)*Mn^2 + 2)/((g+1)*Mn^2);
    Mn2 = sqrt((1 + (g-1)/2*Mn^2)/(g*Mn^2 - (g-1)/2));
    Me = Mn2/sin(b - d);
elseif delta < 0
    % Prandtl-Meyer expansion
    nu = @(m) sqrt((g+1)/(g-1))*atan(sqrt((g-1)/(g+1)*(m.^2 - 1))) - atan(sqrt(m.^2 - 1));
    Me = fzero(@(m) nu(m) - nu(M) - d, [M, 50]);
    Te = T*(1 + (g-1)/2*M^2)/(1 + (g-1)/2*Me^2);
    pe = p*(Te/T)^(g/(g-1));
else
    Me = M;  Te = T;  pe = p;
end
ue  = Me*sqrt(g*R*Te);
Taw = Te*(1 + r*(g-1)/2*Me^2);
Ts  = Te*(1 + 0.032*Me^2 + 0.58*(Tw/Te - 1));            % Eckert reference temperature
rs  = pe/(R*Ts);
mus = 1.458e-6*Ts^1.5/(Ts + 110.4);
Rex = rs*ue*x/mus;
St  = 0.0296*Rex^-0.2*Pr^(-2/3);                          % turbulent flat plate, Reynolds analogy
h   = rs*ue*cp*St;
q   = h*(Taw - Tw);
end
