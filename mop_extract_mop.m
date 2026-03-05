function M = mop_extract_mop(S, varargin)
% mop_extract_mop  Automatically extract the Mean Optical Path (MOP) from
%                  an MCML simulation result structure.
%
% Usage:
%   M = mop_extract_mop(S);
%   M = mop_extract_mop(S, 'PlotOn', true, 'FigureHandle', gca);
%
% Input:
%   S — struct returned by mop_read_mco()
%       Must contain: S.blocks.OP, S.grid.r, S.grid.z,
%                     S.nphot.NphOP, S.pd.Rx, S.pd.Rl,
%                     S.light.x, S.light.type, S.light.l
%
% Output:
%   M — struct with fields:
%         .r_mop     [cm]  r-coords of MOP centre line
%         .z_mop     [cm]  z-coords of MOP centre line
%         .r_upper   [cm]  r-coords of upper envelope
%         .z_upper   [cm]  depth of upper envelope
%         .r_lower   [cm]  r-coords of lower envelope
%         .z_lower   [cm]  depth of lower envelope
%         .max_depth [cm]  maximum penetration depth of MOP
%         .NphOP          number of PD-detected photons
%
% Options (name-value):
%   'PlotOn'        (logical, default false)
%   'FigureHandle'  (axes handle, default gca)
%   'PolyOrder'     (int, default 4)     polynomial order for smoothing
%   'SmoothSigma'   (scalar, default 2)  Gaussian blur sigma (grid cells)
%   'HalfPeakFrac'  (scalar, default 0.3) fraction of peak defining envelope
%   'ColourMOP'     (colour spec, default 'g')
%   'ColourUpper'   (colour spec, default 'c')
%   'ColourLower'   (colour spec, default 'r')
%   'LineWidth'     (scalar, default 2)

%% --- Parse inputs -------------------------------------------------------
p = inputParser;
p.addParameter('PlotOn',       false,  @(x) islogical(x)||ismember(x,[0 1]));
p.addParameter('FigureHandle', [],     @(x) isgraphics(x)||isempty(x));
p.addParameter('PolyOrder',    4,      @(x) isnumeric(x)&&x>=1&&x<=9);
p.addParameter('SmoothSigma',  2,      @isnumeric);
p.addParameter('HalfPeakFrac', 0.3,    @(x) isnumeric(x)&&x>0&&x<1);
p.addParameter('ColourMOP',    'g',    @(x) true);
p.addParameter('ColourUpper',  'c',    @(x) true);
p.addParameter('ColourLower',  'r',    @(x) true);
p.addParameter('LineWidth',    2,      @isnumeric);
p.parse(varargin{:});
opt = p.Results;

%% --- Extract data -------------------------------------------------------
OP    = S.blocks.OP;           % [Nz x Nr]  already normalised by NphOP
r     = S.grid.r(:);
z     = S.grid.z(:);
Nr    = numel(r);
Nz    = numel(z);
dr    = r(2) - r(1);
dz    = z(2) - z(1);

NphOP = S.nphot.NphOP;
if NphOP == 0
    warning('mop_extract_mop: NphOP=0, no PD photons detected.');
    M = struct('r_mop',[],'z_mop',[],'r_upper',[],'z_upper',[],...
               'r_lower',[],'z_lower',[],'max_depth',NaN,'NphOP',0);
    return
end

%% --- Column range [i_src..i_pd] ----------------------------------------
src_r = abs(S.light.x);
pd_r  = abs(S.pd.Rx);
i_src = max(1,  round(src_r / dr));
i_pd  = min(Nr, round(pd_r  / dr));
if i_src >= i_pd
    i_src = max(1, i_pd - 2);
end
margin = max(2, round(0.03 * (i_pd - i_src)));
i1 = max(1, i_src - margin);
i2 = min(Nr, i_pd + margin);
n_cols = i2 - i1 + 1;

%% --- Step 1: 2-D Gaussian smooth to suppress noise ---------------------
sig = opt.SmoothSigma;
if sig > 0
    hw  = ceil(3 * sig);
    [gx, gz] = meshgrid(-hw:hw, -hw:hw);
    kernel   = exp(-(gx.^2 + gz.^2) / (2*sig^2));
    kernel   = kernel / sum(kernel(:));
    OP_sm    = conv2(OP, kernel, 'same');
else
    OP_sm = OP;
end

%% --- Step 2: MOP = depth of maximum OP in each column -----------------
OP_sub = OP_sm(:, i1:i2);   % [Nz x n_cols]
col_max = max(OP_sub, [], 1); % [1 x n_cols]

z_mop_raw   = zeros(1, n_cols);
z_upper_raw = zeros(1, n_cols);
z_lower_raw = zeros(1, n_cols);

for c = 1:n_cols
    col = OP_sub(:, c);
    if col_max(c) == 0; continue; end

    % MOP: index of maximum
    [~, ipeak] = max(col);
    z_mop_raw(c) = z(ipeak);

    % Envelope: find the range where col >= frac * peak
    thresh = opt.HalfPeakFrac * col_max(c);
    above  = find(col >= thresh);
    if isempty(above)
        z_upper_raw(c) = z(ipeak);
        z_lower_raw(c) = z(ipeak);
    else
        z_upper_raw(c) = z(above(1));    % shallowest point above threshold
        z_lower_raw(c) = z(above(end));  % deepest point above threshold
    end
end

%% --- Step 3: Polynomial smooth -----------------------------------------
r_cols = r(i1:i2)';   % [n_cols x 1]
valid  = col_max > 0;
if sum(valid) < opt.PolyOrder + 2
    warning('mop_extract_mop: Not enough valid columns.');
    M = struct('r_mop',[],'z_mop',[],'r_upper',[],'z_upper',[],...
               'r_lower',[],'z_lower',[],'max_depth',NaN,'NphOP',NphOP);
    return
end
rv = r_cols(valid);
poly_ord = min(opt.PolyOrder, sum(valid)-1);

% Trim 5% of endpoints to avoid edge artefacts
ne = max(1, round(0.05 * numel(rv)));
z_mop_v   = z_mop_raw(valid);   z_mop_v(1:ne)=z_mop_v(ne+1);   z_mop_v(end-ne+1:end)=z_mop_v(end-ne);
z_upper_v = z_upper_raw(valid); z_upper_v(1:ne)=z_upper_v(ne+1); z_upper_v(end-ne+1:end)=z_upper_v(end-ne);
z_lower_v = z_lower_raw(valid); z_lower_v(1:ne)=z_lower_v(ne+1); z_lower_v(end-ne+1:end)=z_lower_v(end-ne);

p_mop   = polyfit(rv, z_mop_v,   poly_ord);
p_upper = polyfit(rv, z_upper_v, poly_ord);
p_lower = polyfit(rv, z_lower_v, poly_ord);

r_smooth  = linspace(rv(1), rv(end), 512)';
clip = @(x) max(z(1), min(z(end), x));
z_mop_s   = clip(polyval(p_mop,   r_smooth));
z_upper_s = clip(polyval(p_upper, r_smooth));
z_lower_s = clip(polyval(p_lower, r_smooth));

% Compute MOP as MIDPOINT of upper and lower envelope (true banana centre)
z_mop_s = (z_upper_s + z_lower_s) / 2;

%% --- Pack output --------------------------------------------------------
M.r_mop    = r_smooth;   M.z_mop    = z_mop_s;
M.r_upper  = r_smooth;   M.z_upper  = z_upper_s;
M.r_lower  = r_smooth;   M.z_lower  = z_lower_s;
M.max_depth = max(z_mop_s);
M.NphOP    = NphOP;

fprintf('[mop_extract_mop] NphOP=%d  MOP max depth = %.3f cm\n', NphOP, M.max_depth);

%% --- Optional overlay plot ----------------------------------------------
if opt.PlotOn
    ax = [];
    if ~isempty(opt.FigureHandle)
        if isa(opt.FigureHandle,'matlab.graphics.axis.Axes')
            ax = opt.FigureHandle; axes(ax);
        else
            figure(opt.FigureHandle); ax = gca;
        end
    else
        ax = gca;
    end
    hold(ax, 'on');
    plot(ax, M.r_upper, M.z_upper, '-', 'Color', opt.ColourUpper, ...
         'LineWidth', opt.LineWidth,     'DisplayName', 'Upper envelope');
    plot(ax, M.r_lower, M.z_lower, '-', 'Color', opt.ColourLower, ...
         'LineWidth', opt.LineWidth,     'DisplayName', 'Lower envelope');
    plot(ax, M.r_mop,   M.z_mop,   '-', 'Color', opt.ColourMOP,   ...
         'LineWidth', opt.LineWidth+0.5, 'DisplayName', 'MOP');
    hold(ax, 'off');
end
end
