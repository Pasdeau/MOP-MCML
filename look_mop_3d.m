% look_mop_3d.m  —  3D MCO slice viewer (X/Y/Z draggable planes)
%
% This script mirrors the "edit parameters at the top, then run" workflow of
% `look_mop.m`.
%
% Requires: `mop_read_mco_3d.m`

% ========== USER PARAMETERS ==============================================
mco_file = 'validate_RT/3d_gpu_Ronly.mco';  % 3D .mco file to load

% --- Visual-only options (does not affect .mco) ---
SliceFactor = 1;        % >1: upsample displayed slices by interpolation
Alpha       = 0.85;     % slice transparency (0–1)
CLim        = [-5 -2];  % colour axis for log10(OP)
VerboseRead = false;    % pass to mop_read_mco_3d(...,'Verbose',VerboseRead)
% =========================================================================

opt = struct( ...
    'SliceFactor', SliceFactor, ...
    'Alpha',       Alpha, ...
    'CLim',        CLim, ...
    'VerboseRead', VerboseRead);

local_look_mop_3d(mco_file, opt);

% ========================================================================
function local_look_mop_3d(mco_file, opt)
% Implementation kept as a function so UI callbacks can share state.
% This is still "script-style" usage: you edit parameters above and run.

if nargin < 1 || isempty(mco_file)
    mco_file = 'validate_RT/3d_gpu_Ronly.mco';
end
if nargin < 2 || isempty(opt)
    opt = struct();
end
if ~isfield(opt,'SliceFactor'), opt.SliceFactor = 1; end
if ~isfield(opt,'Alpha'),       opt.Alpha       = 0.85; end
if ~isfield(opt,'CLim'),        opt.CLim        = [-5 -2]; end
if ~isfield(opt,'VerboseRead'), opt.VerboseRead = false; end

fprintf('Loading 3D MCO file: %s\n', mco_file);
S = mop_read_mco_3d(mco_file, 'Verbose', opt.VerboseRead);

% Get the 3D OP block (Nz x Ny x Nx)
OP = S.blocks.OP_3D;

% Avoid log(0)
min_val = min(OP(OP > 0));
if isempty(min_val), min_val = 1e-12; end
OP(OP <= 0) = min_val;
logOP = log10(OP);

% Permute logOP from [Nz, Ny, Nx] to [Ny, Nx, Nz] to match meshgrid ordering
logOP_plot = permute(logOP, [2, 3, 1]);

% Build an interpolant for smooth, high-res slice display
% logOP_plot is [Ny, Nx, Nz] => grid order {y, x, z}
if exist('griddedInterpolant', 'file')
    F = griddedInterpolant({S.grid.y, S.grid.x, S.grid.z}, logOP_plot, 'linear', 'nearest');
    evalF = @(Y, X, Z) F(Y, X, Z);
else
    evalF = @(Y, X, Z) interp3(S.grid.x, S.grid.y, S.grid.z, logOP_plot, X, Y, Z, 'linear');
end

% Precompute query vectors for display resolution
nxq = max(2, round(numel(S.grid.x) * opt.SliceFactor));
nyq = max(2, round(numel(S.grid.y) * opt.SliceFactor));
nzq = max(2, round(numel(S.grid.z) * opt.SliceFactor));
xq = linspace(min(S.grid.x), max(S.grid.x), nxq);
yq = linspace(min(S.grid.y), max(S.grid.y), nyq);
zq = linspace(min(S.grid.z), max(S.grid.z), nzq);

% Create Visualization
fig = figure('Color', 'w', 'Name', '3D MOP Visualization', 'Position', [100 100 900 700]);
ax = axes('Parent', fig);

colormap(ax, jet);
caxis(ax, opt.CLim);
cb = colorbar(ax);
ylabel(cb, 'log_{10}(OP) [1/cm^3]');

xlabel(ax, 'X (cm)', 'FontWeight', 'bold');
ylabel(ax, 'Y (cm)', 'FontWeight', 'bold');
zlabel(ax, 'Depth Z (cm)', 'FontWeight', 'bold');
set(ax, 'ZDir', 'reverse', 'FontSize', 12);
title(ax, '3D Monte Carlo Optical Path (Cartesian Grid)', 'FontSize', 14);

view(ax, -35, 30);
grid(ax, 'on');
camlight(ax);
lighting(ax, 'gouraud');
daspect(ax, [1 1 1]);
xlim(ax, [min(S.grid.x) max(S.grid.x)]);
ylim(ax, [min(S.grid.y) max(S.grid.y)]);
zlim(ax, [min(S.grid.z) max(S.grid.z)]);

% Initial slice positions (can be changed via sliders)
x_slice = clamp_to_range(0, min(S.grid.x), max(S.grid.x));
y_slice = clamp_to_range(0, min(S.grid.y), max(S.grid.y));
z_slice = S.grid.z(round(end/2));

% Draw initial slices
hs = [];
draw_slices();

% ---------------- UI: draggable sliders ----------------
% Use uicontrol sliders for compatibility across MATLAB versions.
uih = 0.04;        % control height (normalized)
margin = 0.012;    % margin
labelw = 0.08;     % label width
valuew = 0.08;     % value text width
sliderw = 1 - 2*margin - labelw - valuew;

% Keep some room at the bottom for UI
ax.Position = [0.08, 0.18, 0.84, 0.78];

% X slider
uicontrol(fig, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [margin, 0.10, labelw, uih], 'String', 'X', 'BackgroundColor', 'w');
sx = uicontrol(fig, 'Style', 'slider', 'Units', 'normalized', ...
    'Position', [margin+labelw, 0.10, sliderw, uih], ...
    'Min', min(S.grid.x), 'Max', max(S.grid.x), 'Value', x_slice, ...
    'Callback', @on_slider);
tx = uicontrol(fig, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [margin+labelw+sliderw, 0.10, valuew, uih], ...
    'String', sprintf('%.3g', x_slice), 'BackgroundColor', 'w');

% Y slider
uicontrol(fig, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [margin, 0.06, labelw, uih], 'String', 'Y', 'BackgroundColor', 'w');
sy = uicontrol(fig, 'Style', 'slider', 'Units', 'normalized', ...
    'Position', [margin+labelw, 0.06, sliderw, uih], ...
    'Min', min(S.grid.y), 'Max', max(S.grid.y), 'Value', y_slice, ...
    'Callback', @on_slider);
ty = uicontrol(fig, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [margin+labelw+sliderw, 0.06, valuew, uih], ...
    'String', sprintf('%.3g', y_slice), 'BackgroundColor', 'w');

% Z slider
uicontrol(fig, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [margin, 0.02, labelw, uih], 'String', 'Z', 'BackgroundColor', 'w');
sz = uicontrol(fig, 'Style', 'slider', 'Units', 'normalized', ...
    'Position', [margin+labelw, 0.02, sliderw, uih], ...
    'Min', min(S.grid.z), 'Max', max(S.grid.z), 'Value', z_slice, ...
    'Callback', @on_slider);
tz = uicontrol(fig, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [margin+labelw+sliderw, 0.02, valuew, uih], ...
    'String', sprintf('%.3g', z_slice), 'BackgroundColor', 'w');

fprintf('Drag the sliders (X/Y/Z) to move the slicing planes.\n');

% ---------------- nested callbacks ----------------
    function on_slider(~, ~)
        x_slice = get(sx, 'Value');
        y_slice = get(sy, 'Value');
        z_slice = get(sz, 'Value');
        set(tx, 'String', sprintf('%.3g', x_slice));
        set(ty, 'String', sprintf('%.3g', y_slice));
        set(tz, 'String', sprintf('%.3g', z_slice));
        draw_slices();
    end

    function draw_slices()
        if ~isempty(hs)
            try
                delete(hs);
            catch
            end
        end
        hold(ax, 'on');

        % X-slice (constant x): sample on a (z,y) grid for smooth display
        [Y1, Z1] = meshgrid(yq, zq);
        X1 = x_slice * ones(size(Y1));
        V1 = evalF(Y1, X1, Z1);
        hsx = surf(ax, X1, Y1, Z1, V1);

        % Y-slice (constant y): sample on a (z,x) grid for smooth display
        [X2, Z2] = meshgrid(xq, zq);
        Y2 = y_slice * ones(size(X2));
        V2 = evalF(Y2, X2, Z2);
        hsy = surf(ax, X2, Y2, Z2, V2);

        % Z-slice (constant z): sample on a (y,x) grid for smooth display
        [X3, Y3] = meshgrid(xq, yq);
        Z3 = z_slice * ones(size(X3));
        V3 = evalF(Y3, X3, Z3);
        hsz = surf(ax, X3, Y3, Z3, V3);

        hs = [hsx; hsy; hsz];
        set(hs, 'EdgeColor', 'none', 'FaceAlpha', opt.Alpha, 'FaceColor', 'interp');
        hold(ax, 'off');
        try
            drawnow limitrate;
        catch
            drawnow;
        end
    end

    function v = clamp_to_range(v, vmin, vmax)
        if v < vmin, v = vmin; end
        if v > vmax, v = vmax; end
    end
end
