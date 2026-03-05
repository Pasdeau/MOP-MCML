function figs = mop_plot_mco(S, varargin)
% mop_plot_mco  Visualize OP map and annotate source/PD geometry.
% Usage:
%   figs = mop_plot_mco(S,'ShowRPD',true,'ShowTPD',true,'CLim',[-6 -1]);
%
% Options:
%   ShowRPD, ShowTPD : whether to render PD rectangles at z=0 and z=Zmax
%   CLim             : color limits in log10 space for OP map
%   Colormap         : function handle or name (default: makec2f if exists, else parula)

p = inputParser;
p.addParameter('ShowRPD', true,  @(x)islogical(x)||ismember(x,[0 1]));
p.addParameter('ShowTPD', true,  @(x)islogical(x)||ismember(x,[0 1]));
p.addParameter('CLim',    [],    @(x)isnumeric(x) && numel(x)==2);
p.addParameter('Colormap',[],    @(x)ischar(x)||isa(x,'function_handle')||isempty(x));
p.parse(varargin{:});
opt = p.Results;

r = S.grid.r; z = S.grid.z; OP = S.blocks.OP;
figs = struct();

% xy = 32; % Reference font size (Unused)
figs.OP = figure();
set(figs.OP, ...
    'Color',       'w', ...
    'Name',        '2D MOP Visualization', ...
    'NumberTitle', 'on', ...
    'Position',    [100 100 900 700]);
clf(figs.OP);

% Use reference plotting style (but default fonts)
imagesc(r, z, log10(OP));
axis tight ij
xlabel('r [cm]');
ylabel('z [cm]');
% title('log_{10}(OP)'); % Reverted as requested (Title removed)

% set(gca, 'FontSize', xy, 'FontWeight', 'Bold'); % Reverted as requested

if isempty(opt.Colormap)
    if exist('makec2f','file'), colormap(makec2f);
    else, colormap(parula); % Fallback
    end
else
    colormap(opt.Colormap);
end

cb = colorbar('eastoutside');
% Reference style for colorbar
set(cb, 'TickLabels', {' ', ' '});
set(cb, 'FontSize', 20);
% ylabel(cb,'log_{10}(OP)'); % Reference code commented out/removed title or label

if ~isempty(opt.CLim), caxis(opt.CLim); end
hold on

% Rectangles (PD)
rect_h = S.grid.dz*5;  % Uniform thickness for PD and Source
z_top = z(1); z_bot = z(end); gap = 0;

if opt.ShowRPD && all(isfinite([S.pd.Rx S.pd.Rl]))
    rectangle('Position',[S.pd.Rx - S.pd.Rl/2, z_top - gap - rect_h, S.pd.Rl, rect_h], ...
        'EdgeColor',[0 0 1],'LineWidth',2,'FaceColor',[0 0 1 0.15],'Clipping','off');
    % Marker removed as requested
end
if opt.ShowTPD && all(isfinite([S.pd.Tx S.pd.Tl]))
    rectangle('Position',[S.pd.Tx - S.pd.Tl/2, z_bot + gap, S.pd.Tl, rect_h], ...
        'EdgeColor',[0 0 1],'LineWidth',2,'FaceColor',[0 0 1 0.15],'Clipping','off');
    % Marker removed as requested
end

% Source glyph
if isfinite(S.light.type) && isfinite(S.light.x)
    switch S.light.type
        case 1 % point
            % Increased triangle size
            tri_h = 0.05; tri_w = 0.05;
            P = [ S.light.x,             z_top; ...
                S.light.x - tri_w/2,   z_top - gap - tri_h; ...
                S.light.x + tri_w/2,   z_top - gap - tri_h ];
            patch('XData',P(:,1),'YData',P(:,2), ...
                'FaceColor',[1 0 0], 'EdgeColor','none', 'Clipping','off');
        case {2,3} % gaussian/flat
            % Match thickness to PD (rect_h)
            sideL = rect_h;
            % Width is S.light.l, Height is sideL
            rectangle('Position',[S.light.x - S.light.l/2, z_top - gap - sideL, S.light.l, sideL], ...
                'EdgeColor',[1 0 0],'LineWidth',2,'FaceColor',[1 0 0 0.15], 'Clipping','off');
    end
end
hold off

% Strict limits to exclude external annotations
ylim([z(1) z(end)]);
% xlim([r(1) r(end)]); % automatic xlim usually fine
end
