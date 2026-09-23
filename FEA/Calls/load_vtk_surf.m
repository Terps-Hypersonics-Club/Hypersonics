function [T, P, U, V, W, centroids, norms, areas, extra] = load_vtk_surf(vtkFile)
% Loads a legacy VTK POLYDATA surface file (ASCII or BINARY) into MATLAB arrays.
% Each row is a triangle with geometry and all cell data fields.
%
% Cell-data SCALARS blocks are read BY NAME, in whatever order they appear
% in the file, so the loader does not depend on field ordering.
%
% Outputs:
%   T, P, U, V, W  - cell temperature, pressure, velocity components
%   centroids      - [num_tri x 3] triangle centroids
%   norms          - [num_tri x 3] unit normals (from winding order)
%   areas          - [num_tri x 1] triangle areas
%   extra          - struct with EVERY cell-data field found (e.g. patch_id,
%                    tau, qw, wall_temperature), keyed by field name

fid = fopen(vtkFile,'r');
assert(fid > 0, 'Cannot open file: %s', vtkFile);
cleanup = onCleanup(@() fclose(fid));

% --- Detect ASCII vs BINARY from the header (line 3)
fgetl(fid); fgetl(fid);
fmt = strtrim(fgetl(fid));
frewind(fid);

if strcmpi(fmt, 'BINARY')
    [pointarr, tri, extra] = read_binary(fid, vtkFile);
elseif strcmpi(fmt, 'ASCII')
    [pointarr, tri, extra] = read_ascii(fid, vtkFile);
else
    error('Unrecognised VTK format line "%s" in %s', fmt, vtkFile);
end

%% --- Geometry from triangles
v1 = pointarr(tri(:,1),:);
v2 = pointarr(tri(:,2),:);
v3 = pointarr(tri(:,3),:);
trinorm   = cross(v2 - v1, v3 - v2, 2);
nmag      = sqrt(sum(trinorm.^2, 2));
areas     = 0.5 * nmag;
centroids = (v1 + v2 + v3) / 3;
norms     = trinorm ./ (nmag + eps);

%% --- Map named fields to the legacy outputs
req = {'T','P','U','V','W'};
for k = 1:numel(req)
    assert(isfield(extra, req{k}), 'Required cell-data field "%s" not found in %s', req{k}, vtkFile);
end
T = extra.T;  P = extra.P;  U = extra.U;  V = extra.V;  W = extra.W;
end

%% =====================================================================
function [pointarr, tri, extra] = read_ascii(fid, vtkFile)
line = fgetl(fid);
while ischar(line) && ~startsWith(strtrim(line), 'POINTS')
    line = fgetl(fid);
end
assert(ischar(line), 'POINTS section not found in %s', vtkFile);
tokens = split(strtrim(line));
npts   = str2double(tokens{2});
pointarr = fscanf(fid, '%f', [3, npts])';
assert(size(pointarr,1) == npts, 'Expected %d points, read %d', npts, size(pointarr,1));

line = fgetl(fid);
while ischar(line) && ~startsWith(strtrim(line), 'POLYGONS')
    line = fgetl(fid);
end
assert(ischar(line), 'POLYGONS section not found in %s', vtkFile);
tokens  = split(strtrim(line));
num_tri = str2double(tokens{2});
poly = fscanf(fid, '%d', [4, num_tri])';
assert(size(poly,1) == num_tri, 'Expected %d polygons, read %d', num_tri, size(poly,1));
assert(all(poly(:,1) == 3), 'Non-triangular polygons found; loader expects triangles');
tri = poly(:,2:4) + 1;

extra = struct();
line = fgetl(fid);
while ischar(line)
    s = strtrim(line);
    if startsWith(s, 'SCALARS')
        tok   = split(s);
        fname = tok{2};
        ncomp = 1;
        if numel(tok) >= 4, ncomp = str2double(tok{4}); end
        lt = fgetl(fid);
        assert(ischar(lt) && startsWith(strtrim(lt), 'LOOKUP_TABLE'), ...
            'Expected LOOKUP_TABLE after SCALARS %s', fname);
        vals = fscanf(fid, '%f', [ncomp, num_tri])';
        assert(size(vals,1) == num_tri, 'Field %s: expected %d values, read %d', fname, num_tri, size(vals,1));
        extra.(matlab.lang.makeValidName(fname)) = vals;
    end
    line = fgetl(fid);
end
end

%% =====================================================================
function [pointarr, tri, extra] = read_binary(fid, vtkFile)
% Legacy binary VTK stores all numeric data BIG-ENDIAN.
raw = fread(fid, '*uint8')';
txt = char(raw);                % for header searching only
pos = 1;

    function [line, dataStart] = header_after(keyword)
        % find "\n<keyword>" at or after pos, return the header line and index of first data byte
        idx = strfind(txt(pos:end), [newline keyword]);
        assert(~isempty(idx), '%s section not found in %s', keyword, vtkFile);
        s = pos + idx(1);                       % start of keyword (after the newline)
        e = s + find(raw(s:end) == 10, 1) - 1;  % index of terminating newline
        line = strtrim(txt(s:e-1));
        dataStart = e + 1;
    end

    function vals = take(dataStart, n, vtype)
        switch vtype
            case 'double', nb = 8; cls = 'double';
            case 'float',  nb = 4; cls = 'single';
            case 'int',    nb = 4; cls = 'int32';
            case 'long',   nb = 8; cls = 'int64';
            otherwise, error('Unsupported VTK type "%s"', vtype);
        end
        assert(dataStart + n*nb - 1 <= numel(raw), 'Unexpected end of file in %s', vtkFile);
        vals = double(swapbytes(typecast(raw(dataStart:dataStart+n*nb-1), cls)));
        pos = dataStart + n*nb;
    end

% POINTS
[line, ds] = header_after('POINTS');
tok = split(line); npts = str2double(tok{2}); ptype = tok{3};
pointarr = reshape(take(ds, npts*3, ptype), 3, npts)';

% POLYGONS
[line, ds] = header_after('POLYGONS');
tok = split(line); num_tri = str2double(tok{2}); nint = str2double(tok{3});
assert(nint == num_tri*4, 'Non-triangular polygons found; loader expects triangles');
conn = reshape(take(ds, nint, 'int'), 4, num_tri)';
assert(all(conn(:,1) == 3), 'Non-triangular polygons found; loader expects triangles');
tri = conn(:,2:4) + 1;

% CELL_DATA
[line, ~] = header_after('CELL_DATA');
tok = split(line); assert(str2double(tok{2}) == num_tri, 'CELL_DATA count mismatch');

% SCALARS blocks
extra = struct();
while true
    idx = strfind(txt(pos:end), [newline 'SCALARS']);
    if isempty(idx), break; end
    [line, ds] = header_after('SCALARS');
    tok = split(line); fname = tok{2}; vtype = tok{3};
    ncomp = 1; if numel(tok) >= 4, ncomp = str2double(tok{4}); end
    % LOOKUP_TABLE line follows immediately
    e = ds + find(raw(ds:end) == 10, 1) - 1;
    assert(startsWith(strtrim(txt(ds:e-1)), 'LOOKUP_TABLE'), 'Expected LOOKUP_TABLE after SCALARS %s', fname);
    vals = take(e + 1, num_tri*ncomp, vtype);
    if ncomp > 1, vals = reshape(vals, ncomp, num_tri)'; end
    extra.(matlab.lang.makeValidName(fname)) = vals;
end
end
