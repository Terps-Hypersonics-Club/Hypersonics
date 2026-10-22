function [T, p, rho] = isa1976(z)
% US Standard Atmosphere 1976, 0-71 km. z = geometric altitude [m] (column).
% Returns temperature [K], pressure [Pa], density [kg/m^3].
z  = z(:);
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
