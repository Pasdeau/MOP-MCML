%   look_mop.m  — Reads an MCML .mco file, plots the OP banana map,
%                 and optionally auto-extracts / overlays the MOP path.
%
%   Requires: mop_read_mco.m, mop_plot_mco.m, mop_extract_mop.m (if MOPon)

% ========== USER PARAMETERS ==============================================
mco_file = 'validate_RT/test_cpu_Ronly.mco';   % .mco file to load
MOPon    = false;              % true: overlay MOP | false: plain OP plot only
VisMode  = '2D';              % '2D' (cylindrical) | '3D' (Cartesian, future)

% --- Display parameters ---
CLim          = [-6 -1];   % colour axis for log10(OP)
r_margin      = 0.3;       % [cm] display margin around source/PD
ZoomOn        = false;      % true: zoom to banana region | false: show full grid

% --- MOP algorithm parameters ---
SmoothSigma   = 3;         % Gaussian blur radius (grid cells); increase for noisy data
HalfPeakFrac  = 0.25;      % envelope threshold: fraction of column peak (0-1)
PolyOrder_MOP = 4;         % polynomial order for MOP smoothing
% =========================================================================

%% 1. Read .mco file
S = mop_read_mco(mco_file, 'Verbose', true);

%% 2. Plot depending on mode
if strcmp(VisMode, '2D')

    %% --- 2D cylindrical OP banana plot ----------------------------------
    figs = mop_plot_mco(S, 'ShowRPD', true, 'ShowTPD', true, 'CLim', CLim);

    % Auto-zoom to banana region (disable with ZoomOn = false)
    if ZoomOn && isfield(S,'light') && isfield(S,'pd')
        r_src = abs(S.light.x);
        r_pd  = abs(S.pd.Rx);
        figure(figs.OP);
        xlim([max(0, min(r_src,r_pd) - r_margin), max(r_src,r_pd) + r_margin]);
        ylim([0, min(1.0, S.grid.z(end))]);
    end

    %% --- MOP overlay (only if MOPon) ------------------------------------
    if MOPon
        if S.nphot.NphOP > 0
            figure(figs.OP);
            ax = gca;
            M = mop_extract_mop(S, ...
                'PlotOn',       true, ...
                'FigureHandle', ax, ...
                'SmoothSigma',  SmoothSigma, ...
                'HalfPeakFrac', HalfPeakFrac, ...
                'PolyOrder',    PolyOrder_MOP, ...
                'ColourMOP',    [1 0.85 0], ...   % yellow — MOP centre
                'ColourUpper',  'w', ...           % white  — upper (shallow) boundary
                'ColourLower',  'c', ...           % cyan   — lower (deep)  boundary
                'LineWidth',    2);
            legend({'Upper boundary (shallow)', ...
                    'Lower boundary (deep)', ...
                    'MOP (centre)'}, ...
                   'Location','southwest','TextColor','w');
            fprintf('\nMOP max depth: %.3f cm\n', M.max_depth);
        else
            fprintf('NphOP = 0 — cannot extract MOP.\n');
        end
    end

elseif strcmp(VisMode, '3D')
    %% --- 3D Cartesian placeholder (future) ------------------------------
    % The 3D mode will require a .mco file generated with 3D grid output.
    % Planned implementation:
    %   S3 = mop_read_mco_3d(mco_file);       % reads x,y,z grid + OP_xyz
    %   mop_plot_mco_3d(S3, 'CLim', CLim);    % isosurface / slice plot
    %   mop_extract_mop_3d(S3, 'PlotOn', MOPon);
    fprintf('[VisMode=3D] 3D visualisation not yet implemented.\n');
    fprintf('Run CPU/GPU with 3D grid output first (job pending).\n');

else
    error('look_mop: VisMode must be ''2D'' or ''3D''. Got: %s', VisMode);
end

