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

figs.OP = figure('Color','w'); clf;
imagesc(r, z, log10(max(OP, eps)));
axis tight ij
xlabel('r [cm]'); ylabel('z [cm]');
title('log_{10}(OP)');

if isempty(opt.Colormap)
    if exist('makec2f','file'), colormap(makec2f);
    else, colormap(parula);
    end
else
    colormap(opt.Colormap);
end
cb = colorbar('eastoutside'); ylabel(cb,'log_{10}(OP)');

if ~isempty(opt.CLim), caxis(opt.CLim); end
hold on

% PD rectangles
rect_h = S.grid.dz*5;  % plot height outside domain for visibility
z_top = z(1); z_bot = z(end); gap = 0;

if opt.ShowRPD && all(isfinite([S.pd.Rx S.pd.Rl]))
    rectangle('Position',[S.pd.Rx - S.pd.Rl/2, z_top - gap - rect_h, S.pd.Rl, rect_h], ...
        'EdgeColor',[0 0 1],'LineWidth',2,'FaceColor',[0 0 1 0.15],'Clipping','off');
    plot(S.pd.Rx, z_top, 'v','MarkerFaceColor',[0 0 1],'MarkerEdgeColor','none','Clipping','off');
end
if opt.ShowTPD && all(isfinite([S.pd.Tx S.pd.Tl]))
    rectangle('Position',[S.pd.Tx - S.pd.Tl/2, z_bot + gap, S.pd.Tl, rect_h], ...
        'EdgeColor',[0 0 1],'LineWidth',2,'FaceColor',[0 0 1 0.15],'Clipping','off');
    plot(S.pd.Tx, z_bot, '^','MarkerFaceColor',[0 0 1],'MarkerEdgeColor','none','Clipping','off');
end

% Source glyph
if isfinite(S.light.type) && isfinite(S.light.x)
    switch S.light.type
        case 1 % point
            tri_h = 0.02; tri_w = 0.02;
            P = [ S.light.x,             z_top; ...
                  S.light.x - tri_w/2,   z_top - gap - tri_h; ...
                  S.light.x + tri_w/2,   z_top - gap - tri_h ];
            patch('XData',P(:,1),'YData',P(:,2), ...
                  'FaceColor',[1 0 0], 'EdgeColor','none', 'Clipping','off');
        case {2,3} % gaussian/flat
            sideL = max(rect_h, S.light.l);
            rectangle('Position',[S.light.x - sideL/2, z_top - gap - sideL, sideL, sideL], ...
                      'EdgeColor',[1 0 0],'LineWidth',2,'FaceColor',[1 0 0 0.15], 'Clipping','off');
    end
end
hold off
end
