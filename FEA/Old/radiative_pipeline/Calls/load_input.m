function [Pinf, Tinf, mach, rhoinf, pran, gamma, Rgas] = load_input(inputFile)
% Parse top-level numeric keys from your Dictionary-style file.

txt = fileread(inputFile);
lines = regexp(txt, '\r\n|\n|\r', 'split');

% Keys we want: map output var -> key in file
keys = struct( ...
    'Tinf', 'Tref', ...
    'gamma', 'gamma', ...
    'Rgas', 'Rgas', ...
    'mach', 'mach', ...
    'aoa',  'aoa', ...
    'Pinf', 'Pref', ...
    'pran', 'pran');

% Initialize as NaN to catch missing values
Tinf = NaN; gamma = NaN; Rgas = NaN; mach = NaN; aoa = NaN; Pinf = NaN; pran = NaN;

% Helper regex builder: exact key at start, equals, numeric literal
getrx = @(k) ['^\s*' k '\s*=\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*(?:$|//|#|;|%|})'];

% Scan once and fill matches
fields = fieldnames(keys);
found = false(size(fields));
for i = 1:numel(lines)
    L = strtrim(lines{i});
    if isempty(L) || startsWith(L, {'//','#',';','%'})
        continue
    end
    for f = 1:numel(fields)
        if found(f), continue, end
        keyname = keys.(fields{f});
        rx = getrx(keyname);
        m = regexp(L, rx, 'tokens', 'once');
        if ~isempty(m)
            val = str2double(m{1});
            switch fields{f}
                case 'Tinf', Tinf = val;
                case 'gamma', gamma = val;
                case 'Rgas',  Rgas = val;
                case 'mach',  mach = val;
                case 'aoa',   aoa  = val; %#ok<NASGU> % parsed but not returned
                case 'Pinf',  Pinf = val;
                case 'pran',  pran = val;
            end
            found(f) = true;
        end
    end
    if all(found), break, end
end

% Basic checks with clear errors
req = {'Tinf','gamma','Rgas','mach','Pinf','pran'};
for r = 1:numel(req)
    if isnan(eval(req{r}))
        error('Missing or nonnumeric value for %s in %s', req{r}, inputFile);
    end
end

% Derived
rhoinf = Pinf/(Rgas*Tinf);
end
