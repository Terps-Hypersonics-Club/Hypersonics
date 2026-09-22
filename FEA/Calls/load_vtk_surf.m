function [T, P, U, V, W, centroids, norms, areas, extra] = load_vtk_surf(vtkFile)
% Loads an ASCII VTK POLYDATA surface file into MATLAB arrays.
% Each row is a triangle with geometry and all cell data fields.
%
% Cell-data SCALARS blocks are read BY NAME, in whatever order they appear
% in the file, so the loader does not depend on field ordering.
%
% Outputs:
%   T, P, U, V, W  - cell temperature, pressure, velocity components
%   centroids      - [num_tri x 3] triangle centroids
%   norms          - [num_tri x 3] unit normals
%   areas          - [num_tri x 1] triangle areas
%   extra          - struct with EVERY cell-data field found (e.g. patch_id,
%                    tau, qw), keyed by field name

fid = fopen(vtkFile,'r');
assert(fid > 0, 'Cannot open file: %s', vtkFile);
cleanup = onCleanup(@() fclose(fid));

%% --- Skip header lines until POINTS section
line = fgetl(fid);
while ischar(line) && ~startsWith(strtrim(line), 'POINTS')
    line = fgetl(fid);
end
assert(ischar(line), 'POINTS section not found in %s', vtkFile);
tokens = split(strtrim(line));
npts   = str2double(tokens{2});

%% --- Read points (3 values per point, whitespace separated)
pointarr = fscanf(fid, '%f', [3, npts])';
assert(size(pointarr,1) == npts, 'Expected %d points, read %d', npts, size(pointarr,1));

%% --- Skip until POLYGONS section
line = fgetl(fid);
while ischar(line) && ~startsWith(strtrim(line), 'POLYGONS')
    line = fgetl(fid);
end
assert(ischar(line), 'POLYGONS section not found in %s', vtkFile);
tokens  = split(strtrim(line));
num_tri = str2double(tokens{2});

%% --- Read POLYGONS (each line: 3 i j k, zero-based indices)
poly = fscanf(fid, '%d', [4, num_tri])';
assert(size(poly,1) == num_tri, 'Expected %d polygons, read %d', num_tri, size(poly,1));
assert(all(poly(:,1) == 3), 'Non-triangular polygons found; loader expects triangles');
idx = poly(:,2:4) + 1;

v1 = pointarr(idx(:,1),:);
v2 = pointarr(idx(:,2),:);
v3 = pointarr(idx(:,3),:);
trinorm   = cross(v2 - v1, v3 - v2, 2);
nmag      = sqrt(sum(trinorm.^2, 2));
areas     = 0.5 * nmag;
centroids = (v1 + v2 + v3) / 3;
norms     = trinorm ./ (nmag + eps);

%% --- Read every cell-data SCALARS block by name
extra = struct();
line = fgetl(fid);
while ischar(line)
    s = strtrim(line);
    if startsWith(s, 'SCALARS')
        tok   = split(s);
        fname = tok{2};
        ncomp = 1;
        if numel(tok) >= 4
            ncomp = str2double(tok{4});
        end
        % Next line is LOOKUP_TABLE <name>
        lt = fgetl(fid);
        assert(ischar(lt) && startsWith(strtrim(lt), 'LOOKUP_TABLE'), ...
            'Expected LOOKUP_TABLE after SCALARS %s', fname);
        vals = fscanf(fid, '%f', [ncomp, num_tri])';
        assert(size(vals,1) == num_tri, ...
            'Field %s: expected %d values, read %d', fname, num_tri, size(vals,1));
        extra.(matlab.lang.makeValidName(fname)) = vals;
    end
    line = fgetl(fid);
end

%% --- Map named fields to the legacy outputs
req = {'T','P','U','V','W'};
for k = 1:numel(req)
    assert(isfield(extra, req{k}), 'Required cell-data field "%s" not found in %s', req{k}, vtkFile);
end
T = extra.T;
P = extra.P;
U = extra.U;
V = extra.V;
W = extra.W;
end
