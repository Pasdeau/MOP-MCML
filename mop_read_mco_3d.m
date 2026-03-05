function S = mop_read_mco_3d(fname, varargin)
% mop_read_mco_3d  Parse 3D MCML (Cartesian) .mco file into a structured MATLAB object.
%
%   S = mop_read_mco_3d('test_3d.mco','Verbose',true);
%
% Key outputs in S:
%   S.grid:     dz, dr, dx, dy, Nz, Nr, Na, Nx, Ny, z(:), r(:), x(:), y(:)
%   S.blocks:   Contains OP_3D of size [Nz, Ny, Nx]
%
% Notes:
%   Custom parser designed to read 3D Cartesian coordinates.

p = inputParser;
p.addParameter('Verbose', true, @(x)islogical(x)||ismember(x,[0 1]));
p.parse(varargin{:});
vb = p.Results.Verbose;

fid = fopen(fname,'rb');
assert(fid>0, 'Cannot open file: %s', fname);
c = onCleanup(@() fclose(fid));

% Helper to read next non-empty line
    function L = nextline()
        L = fgetl(fid);
        if ~ischar(L), error('Unexpected EOF while parsing %s', fname); end
    end

% Helper to strip trailing comments and parse floats
    function v = parsefloats(L)
        Lc = regexprep(L,'#.*$',''); 
        v = sscanf(Lc, '%f');
    end

S = struct();
S.notes = {};

line = '';
for i=1:15
    line = nextline();
    if i==1
        if startsWith(line, 'MCML') || startsWith(line, 'A')
            S.format = 'A';
        elseif startsWith(line, 'B')
            S.format = 'B';
        end
    end
    if vb && i<=3, S.notes{end+1} = sprintf('hdr%02d: %s',i,line); end 
end
u = parsefloats(line);
assert(numel(u)>=4, 'Header parse failed: expected dz dr dx dy on line 15.');
S.grid.dz = u(1); S.grid.dr = u(2); S.grid.dx = u(3); S.grid.dy = u(4);

line = nextline();             % Nz Nr Na Nx Ny
u = parsefloats(line);
assert(numel(u)>=5, 'Expected Nz Nr Na Nx Ny.');
S.grid.Nz = u(1); S.grid.Nr = u(2); S.grid.Na = u(3);
S.grid.Nx = u(4); S.grid.Ny = u(5);

line = nextline();             % skip
line = nextline();             % Nlayers
S.layers.Nlayers = sscanf(line, '%d');
assert(~isempty(S.layers.Nlayers),'Expected Nlayers.');

line = nextline();             % skip
line = nextline();             % nabove
S.layers.nabove = sscanf(line,'%f');
assert(~isempty(S.layers.nabove),'Expected nabove.');

% Layer table
L = S.layers.Nlayers;
S.layers.n   = zeros(L,1);
S.layers.mua = zeros(L,1);
S.layers.mus = zeros(L,1);
S.layers.g   = zeros(L,1);
S.layers.d   = zeros(L,1);
if vb, fprintf('\t#:\tn\tmua\tmus\tg\td\n'); end
for i=1:L
    line = nextline();
    u = parsefloats(line);
    assert(numel(u)>=5, 'Layer line %d malformed.', i);
    S.layers.n(i)   = u(1);  S.layers.mua(i) = u(2);
    S.layers.mus(i) = u(3);  S.layers.g(i)   = u(4);
    S.layers.d(i)   = u(5);
    if vb
        fprintf('\t%d:\t%.3g\t%.3g\t%.3g\t%.3g\t%.4g\n', i, ...
            S.layers.n(i), S.layers.mua(i), S.layers.mus(i), S.layers.g(i), S.layers.d(i));
    end
end

line = nextline();                       % nbelow
S.layers.nbelow = sscanf(line,'%f');
assert(~isempty(S.layers.nbelow),'Expected nbelow.');
line = nextline();                       % likely blank or section header

% PD / Light
pos = ftell(fid);
linePD = nextline(); vPD = parsefloats(linePD);
if numel(vPD)>=6
    S.pd.Rx = vPD(1); S.pd.Ry = vPD(2); S.pd.Rl = vPD(3);
    S.pd.Tx = vPD(4); S.pd.Ty = vPD(5); S.pd.Tl = vPD(6);
    lineLT = nextline(); vLT = parsefloats(lineLT);
    assert(numel(vLT)>=4, 'Light line present but malformed.');
    S.light.type = round(vLT(1));
    S.light.x    = vLT(2); S.light.y = vLT(3); S.light.l = vLT(4);
else
    fseek(fid, pos, 'bof');
    S.pd = struct('Rx',NaN,'Ry',NaN,'Rl',NaN,'Tx',NaN,'Ty',NaN,'Tl',NaN);
    S.light = struct('type',NaN,'x',NaN,'y',NaN,'l',NaN);
end

% RAT block
hdr = nextline(); 
line = nextline(); S.rat.Rsp = sscanf(line,'%f'); 
line = nextline(); S.rat.Rd  = sscanf(line,'%f'); 
line = nextline(); S.rat.A   = sscanf(line,'%f'); 
line = nextline(); S.rat.Td  = sscanf(line,'%f'); 
line = nextline(); % skip

% Photon counters
line = nextline(); % header
line = nextline(); S.nphot.NphR  = sscanf(line,'%d'); 
line = nextline(); S.nphot.NphT  = sscanf(line,'%d'); 
line = nextline(); S.nphot.NphOP = sscanf(line,'%d'); 

% Al
line = nextline(); % header
line = nextline(); % subheader
S.blocks.Al = zeros(L,1);
for i=1:L
    line = nextline();
    S.blocks.Al(i) = sscanf(line,'%f');
end

% A_z
line = nextline(); % header
line = nextline(); % subheader
S.blocks.Az = fscanf(fid, '%f', S.grid.Nz);

% Rd_r
line = nextline();
while all(isspace(line)), line = nextline(); end % seek header
S.blocks.Rr = fscanf(fid,'%f', S.grid.Nr);

% Rd_a
line = nextline();
while all(isspace(line)), line = nextline(); end % seek header
S.blocks.Ra = zeros(S.grid.Na,1);
for i=1:S.grid.Na
    line = nextline();
    S.blocks.Ra(i) = sscanf(line,'%f');
end

% Tt_r
line = nextline();
while all(isspace(line)), line = nextline(); end % seek header
S.blocks.Tr = fscanf(fid, '%f', S.grid.Nr);

% Tt_a
line = nextline();
while all(isspace(line)), line = nextline(); end % seek header
S.blocks.Ta = zeros(S.grid.Na,1);
for i=1:S.grid.Na
    line = nextline();
    S.blocks.Ta(i) = sscanf(line,'%f');
end

% Azr
line = nextline();
while all(isspace(line)), line = nextline(); end % seek header
for i=1:5, line = nextline(); end % skip strings
u = fscanf(fid,'%f', S.grid.Nz * S.grid.Nr);
S.blocks.Azr = reshape(u, [S.grid.Nz S.grid.Nr]);

% Rra
line = nextline();
while all(isspace(line)), line = nextline(); end % seek header
for i=1:5, line = nextline(); end % skip strings
u = fscanf(fid, '%f', S.grid.Nr * S.grid.Na);
S.blocks.Rra = reshape(u, [S.grid.Nr S.grid.Na]);

% Tra
line = nextline();
while all(isspace(line)), line = nextline(); end % seek header
for i=1:5, line = nextline(); end % skip strings
u = fscanf(fid, '%f', S.grid.Nr * S.grid.Na);
S.blocks.Tra = reshape(u, [S.grid.Nr S.grid.Na]);

% OP_3D -> OP_3D (Nz x Ny x Nx)
line = nextline();
while all(isspace(line)), line = nextline(); end % seek header
for i=1:5, line = nextline(); end % skip strings

nx = S.grid.Nx; ny = S.grid.Ny; nz = S.grid.Nz;
total_op = nx * ny * nz;

if isfield(S, 'format') && strcmp(S.format, 'B')
    u = fread(fid, total_op, 'double');
    fgetl(fid); % empty newline
else
    u = fscanf(fid,'%f'); % Read the entire rest of the file
end

need = total_op - numel(u);
if need>0 && need<=10
    u = [u(:); zeros(need,1)]; % padding at the end
end
S.blocks.OP_3D = reshape(u(1:total_op), [nz ny nx]);

if ~isempty(S.nphot.NphOP) && S.nphot.NphOP>0
    S.blocks.OP_3D = S.blocks.OP_3D / double(S.nphot.NphOP);
end

% Derived axes
S.grid.z = ((1:S.grid.Nz)' - 0.5).*S.grid.dz;
S.grid.r = ((1:S.grid.Nr)' - 0.5).*S.grid.dr;

S.grid.x = (((1:S.grid.Nx)') - S.grid.Nx/2 - 0.5) * S.grid.dx;
S.grid.y = (((1:S.grid.Ny)') - S.grid.Ny/2 - 0.5) * S.grid.dy;

end
