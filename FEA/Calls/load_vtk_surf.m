function [T, P, U, V, W, centroids, norms, areas] = load_vtk_surf(vtkFile)
% Loads an ASCII VTK POLYDATA surface file into MATLAB arrays.
% Each row is a triangle with geometry and all cell data fields.

fid = fopen(vtkFile,'r');
assert(fid>0, 'Cannot open file.');
%% --- Skip header lines until POINTS section
line = fgetl(fid);
while ischar(line)
    if contains(line, 'POINTS')
        break;
    end
    line = fgetl(fid);
end

%% --- Read points
tokens = split(strtrim(line));
points = str2double(tokens(2));

pointarr = zeros(points,3);
for i = 1:points
    line = fgetl(fid);
    pointarr(i,:) = str2double(split(strtrim(line)));
end

%% --- Skip header lines until POLYGONS section
line = fgetl(fid);
while ischar(line)
    if contains(line, 'POLYGONS')
        break;
    end
    line = fgetl(fid);
end

%% --- Read POLYGONS
tokens   = split(strtrim(line));
num_tri  = str2double(tokens(2));

% Table declarations
norms     = zeros(num_tri, 3);
areas     = zeros(num_tri, 1);
P         = zeros(num_tri, 1);
T         = zeros(num_tri, 1);
U         = zeros(num_tri, 1);
V         = zeros(num_tri, 1);
W         = zeros(num_tri, 1);
centroids = zeros(num_tri, 3);

for i = 1:num_tri
    line   = fgetl(fid);
    tokens = str2double(split(strtrim(line)));
    v1     = pointarr(tokens(2)+1,:);
    v2     = pointarr(tokens(3)+1,:);
    v3     = pointarr(tokens(4)+1,:);
    trinorm      = cross(v2-v1, v3-v2);
    areas(i)     = 0.5*norm(trinorm);
    centroids(i,:) = (v1+v2+v3)/3;
    norms(i,:)   = trinorm / norm(trinorm + eps);
end

%% Helper: function to skip until keyword
    function skipTo(keyword)
        l = fgetl(fid);
        while ischar(l)
            if contains(l, keyword)
                break;
            end
            l = fgetl(fid);
        end
    end

%% --- Pressures
skipTo('SCALARS P');
fgetl(fid); % skip one line after header
for i = 1:num_tri
    P(i) = str2double(fgetl(fid));
end
%% --- Temperatures
fgetl(fid);
for i = 1:num_tri
    T(i) = str2double(fgetl(fid));
end
%% --- U velocities
fgetl(fid);
for i = 1:num_tri
    U(i) = str2double(fgetl(fid));
end

%% --- V velocities
fgetl(fid);
for i = 1:num_tri
    V(i) = str2double(fgetl(fid));
end

%% --- W velocities
fgetl(fid);
for i = 1:num_tri
    W(i) = str2double(fgetl(fid));
end

fclose(fid);
end
