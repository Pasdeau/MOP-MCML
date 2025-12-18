function S = mop_read_mco(fname, varargin)
% mop_read_mco  Parse MCML/MOP-MCML .mco file into a structured MATLAB object.\
% 
%   S = mop_read_mco('fat800.mco','Verbose',true);
%
% Key outputs in S:
%   S.grid:     dz, dr, Nz, Nr, Na, z(:), r(:)
%   S.layers:   n(:), mua(:), mus(:), g(:), d(:), nabove, nbelow
%   S.pd:       Rx,Ry,Rl, Tx,Ty,Tl (NaN if not present)
%   S.light:    type, x, y, l      (NaN if not present)
%   S.rat:      Rsp, Rd, A, Td
%   S.nphot:    NphR, NphT, NphOP
%   S.blocks:   Al(:), Az(:), Rr(:), Ra(:), Tr(:), Ta(:), Azr(Nz,Nr), Rra(Nr,Na), Tra(Nr,Na), OP(Nz,Nr)
%   S.aux:      Azrm, rm (mirror-extended), notes
%
% Notes:
%   Robust to presence/absence of the two appended metadata lines (PD, light).
%   Minimizes assumptions but keeps MCML block-skipping behavior where needed.

p = inputParser;
p.addParameter('Verbose', true, @(x)islogical(x)||ismember(x,[0 1]));
p.parse(varargin{:});
vb = p.Results.Verbose;

fid = fopen(fname,'r');
assert(fid>0, 'Cannot open file: %s', fname);
c = onCleanup(@() fclose(fid));

% Helper to read next non-empty line (keep raw line for headers)
    function L = nextline()
        L = fgetl(fid);
        if ~ischar(L), error('Unexpected EOF while parsing %s', fname); end
    end

% Helper to strip trailing comments (after #) and parse floats
    function v = parsefloats(L)
        Lc = regexprep(L,'#.*$',''); % remove trailing comments
        v = sscanf(Lc, '%f');
    end

S = struct();
S.notes = {};

% Read line by line, accumulate total floats; skip comments/empty lines automatically
function u = local_read_block(fid, total, parsefloats)
    u = zeros(total,1); k = 0;
    while k < total
        pos = ftell(fid);
        L = fgetl(fid);
        if ~ischar(L)
            break  % Unexpected EOF
        end
        v = parsefloats(L);   % Value after removing # comment
        if isempty(v)
            % Empty line or header line, continue
            continue
        end
        need = total - k;
        if numel(v) <= need
            u(k+1:k+numel(v)) = v;
            k = k + numel(v);
        else
            % Rare: Line provided more numbers than needed. Take only needed part.
            u(k+1:total) = v(1:need);
            k = total;
            % If worried about discarding extra numbers, can fseek back here.
            % But common MCML output doesn't trigger this branch.
        end
    end
    if k < total
        % Conservative padding strategy matching your original OP section script
        u(k+1:total) = 0;
        warning('local_read_block:short','Insufficient data: padded %d elements with zeros.', total-k);
    end
end

% Legacy header skipping: MCML writes a textual header before grid line.
% Keep behavior but validate the final line has 2 floats: dz dr.
line = '';
for i=1:15
    line = nextline();
    if vb && i<=3, S.notes{end+1} = sprintf('hdr%02d: %s',i,line); end %#ok<AGROW>
end
u = parsefloats(line);
assert(numel(u)>=2, 'Header parse failed: expected dz dr on line 15.');
S.grid.dz = u(1); S.grid.dr = u(2);

line = nextline();             % Nz Nr Na
u = parsefloats(line);
assert(numel(u)>=3, 'Expected Nz Nr Na.');
S.grid.Nz = u(1); S.grid.Nr = u(2); S.grid.Na = u(3);

line = nextline();             % skip
line = nextline();             % Nlayers
S.layers.Nlayers = sscanf(line, '%d');
assert(~isempty(S.layers.Nlayers),'Expected Nlayers.');

line = nextline();             % skip
line = nextline();             % nabove
S.layers.nabove = sscanf(line,'%f');
assert(~isempty(S.layers.nabove),'Expected nabove.');

% Layer table: n, mua, mus, g, d (per layer)
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
    S.layers.n(i)   = u(1);
    S.layers.mua(i) = u(2);
    S.layers.mus(i) = u(3);
    S.layers.g(i)   = u(4);
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

% Optional appended PD line (6 floats) and light line (4 floats)
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
    % Not appended; rewind and mark NaN fields
    fseek(fid, pos, 'bof');
    S.pd = struct('Rx',NaN,'Ry',NaN,'Rl',NaN,'Tx',NaN,'Ty',NaN,'Tl',NaN);
    S.light = struct('type',NaN,'x',NaN,'y',NaN,'l',NaN);
end

% RAT block: header then Rsp, Rd, A, Td as scalars
hdr = nextline(); %#ok<NASGU>
line = nextline(); S.rat.Rsp = sscanf(line,'%f'); if vb, fprintf('Rsp = %g\n',S.rat.Rsp); end
line = nextline(); S.rat.Rd  = sscanf(line,'%f'); if vb, fprintf('Rd  = %g\n',S.rat.Rd ); end
line = nextline(); S.rat.A   = sscanf(line,'%f'); if vb, fprintf('A   = %g\n',S.rat.A  ); end
line = nextline(); S.rat.Td  = sscanf(line,'%f'); if vb, fprintf('Td  = %g\n',S.rat.Td ); end
line = nextline(); % skip

% Photon counters
line = nextline(); % header
line = nextline(); S.nphot.NphR  = sscanf(line,'%d'); if vb, fprintf('NphR  = %d\n',S.nphot.NphR ); end
line = nextline(); S.nphot.NphT  = sscanf(line,'%d'); if vb, fprintf('NphT  = %d\n',S.nphot.NphT ); end
line = nextline(); S.nphot.NphOP = sscanf(line,'%d'); if vb, fprintf('NphOP = %d\n',S.nphot.NphOP); end

% Absorption per layer Al
line = nextline(); % header
line = nextline(); % subheader
S.blocks.Al = zeros(L,1);
for i=1:L
    line = nextline();
    S.blocks.Al(i) = sscanf(line,'%f');
end

% A_z (vector of length Nz)
line = nextline(); % header
line = nextline(); % subheader
S.blocks.Az = fscanf(fid, '%f', S.grid.Nz);

% Rd_r -> Rr (Nr)
line = nextline(); % header or newline
S.blocks.Rr = fscanf(fid,'%f', S.grid.Nr);

% Rd_a -> Ra (Na)
line = nextline(); % header
S.blocks.Ra = zeros(S.grid.Na,1);
for i=1:S.grid.Na
    line = nextline();
    S.blocks.Ra(i) = sscanf(line,'%f');
end

% Tt_r -> Tr (Nr)
line = nextline(); % header
line = nextline(); % subheader
S.blocks.Tr = fscanf(fid, '%f', S.grid.Nr);

% Tt_a -> Ta (Na)
line = nextline(); % header
S.blocks.Ta = zeros(S.grid.Na,1);
for i=1:S.grid.Na
    line = nextline();
    S.blocks.Ta(i) = sscanf(line,'%f');
end

% A_rz -> Azr (Nz x Nr)
line = nextline(); % header
line = nextline(); % subheader
for i=1:5, line = nextline(); end % legacy skip
u = fscanf(fid,'%f', S.grid.Nz * S.grid.Nr);
S.blocks.Azr = reshape(u, [S.grid.Nz S.grid.Nr]);

% Rd_ra -> Rra (Nr x Na)
line = nextline(); % header
for i=1:5, line = nextline(); end % legacy skip
u = fscanf(fid, '%f', S.grid.Nr * S.grid.Na);
S.blocks.Rra = reshape(u, [S.grid.Nr S.grid.Na]);

% Tt_ra -> Tra (Nr x Na)
line = nextline(); % header
for i=1:5, line = nextline(); end % legacy skip
u = fscanf(fid, '%f', S.grid.Nr * S.grid.Na);
S.blocks.Tra = reshape(u, [S.grid.Nr S.grid.Na]);

% OP -> OP (Nz x Nr), normalize by NphOP if present
line = nextline(); % header
line = nextline(); % subheader
for i=1:5, line = nextline(); end % legacy skip
u = fscanf(fid,'%f'); % OP block size may vary across MCML variants
% Some legacy files omit the first few bins; pad if needed to Nz*Nr
nz = S.grid.Nz; nr = S.grid.Nr;
need = nz*nr - numel(u);
if need>0 && need<=10
    u = [zeros(need,1); u]; % conservative padding as in legacy scripts
end
S.blocks.OP = reshape(u(1:nz*nr), [nz nr]);
if ~isempty(S.nphot.NphOP) && S.nphot.NphOP>0
    S.blocks.OP = S.blocks.OP / double(S.nphot.NphOP);
end

% Derived axes
S.grid.z = ((1:S.grid.Nz)' - 0.5).*S.grid.dz;
S.grid.r = ((1:S.grid.Nr)' - 0.5).*S.grid.dr;

% Symmetrized Azr for visualization
S.aux.Azrm = zeros(S.grid.Nz, 2*S.grid.Nr-1);
S.aux.Azrm(:, S.grid.Nr:end) = S.blocks.Azr;
S.aux.Azrm(:, 1:S.grid.Nr)   = S.blocks.Azr(:, S.grid.Nr:-1:1);
S.aux.rm = [((-S.grid.Nr+1:-1) - 0.5)*S.grid.dr, (S.grid.r' - S.grid.dr)]' + S.grid.dr/2;

if vb
    fprintf('Parsed %s successfully. Nz=%d, Nr=%d, Na=%d, layers=%d.\n', ...
        fname, S.grid.Nz, S.grid.Nr, S.grid.Na, S.layers.Nlayers);
end
end
