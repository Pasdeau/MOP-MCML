% function get_mop(name)%   get_mop.m reads the MCML output file <name> and returns MOP.% 	If PRINTON == 1 or PLOTON == 1, output is turned on.%   USES makec2f.m which creates a colormap.%   by Steven L. Jacques, Jan. 2008%   Oregon Health & Sciences University, Portland, OR, USAglobal A Al Az Azr OP Fzr Fz Na Nc Nlayers NphOP NphR NphTglobal Nr Nz Rd Ra Rr Rra Rsp T Ta Td Tr Tra d dr dz global g mua mus n nabove nbelow r zglobal rm Azrm Fzrmglobal PLOTON PRINTONPLOTON = 1; PRINTON = 1;SHOW_R_PD = true;    % R mode PDSHOW_T_PD = true;    % T mode PDfid = fopen('960.mco','r');%%% read dr, dz, Nz, Nr, Na, Nlayer %%%%for i=1:15; line = fgetl(fid); disp(line); endu = sscanf(line, '%f %f');dz = u(1);dr = u(2);line = fgetl(fid);u = sscanf(line, '%d %d %d');Nz = u(1);Nr = u(2);Na = u(3);line = fgetl(fid);line = fgetl(fid);Nlayers = sscanf(line, '%d');if PRINTON    fprintf('dr = %0.4f\n',dr)    fprintf('dz = %0.4f\n',dz)    fprintf('Nr = %d\n',Nr)    fprintf('Nz = %d\n',Nz)    fprintf('Na = %d\n',Na)    fprintf('Nlayer = %d\n', Nlayers)endline = fgetl(fid);line = fgetl(fid);nabove = sscanf(line, '%f');if PRINTON    fprintf('nabove = %0.3f\n', nabove)endn = zeros(Nlayers,1); mua = n; mus = n; g = n; d = n;if PRINTON    fprintf('\t#:\tn   \tmua  \tmus  \tg    \td\n')endfor i=1:Nlayers	line = fgetl(fid);	u = sscanf(line, '%f %f %f %f %f');	n(i)   = u(1);	mua(i) = u(2);	mus(i) = u(3);	g(i)   = u(4);	d(i)   = u(5);	if PRINTON        fprintf('\t%d:\t%0.2f\t%0.2f\t%0.1f\t%0.3f\t%0.4f\n', i, n(i), mua(i), mus(i), g(i), d(i) )    endendline = fgetl(fid);nbelow = sscanf(line,'%f');if PRINTON    fprintf('nbelow = %0.3f\n', nbelow)endline = fgetl(fid);% Read PD line (expect 6 numbers)
line = fgetl(fid);
line_clean = regexprep(line, '#.*$', '');   % remove trailing comment for safety
u = sscanf(line_clean, '%f');
if numel(u) < 6
    error('PD line expects 6 numbers, actual got %d: %s', numel(u), line);
end
[Rx, Ry, Rl, Tx, Ty, Tl] = deal(u(1), u(2), u(3), u(4), u(5), u(6));

% Read light line (expect 4 numbers)
line = fgetl(fid);
line_clean = regexprep(line, '#.*$', '');
u = sscanf(line_clean, '%f');   % %f works, light_type can be rounded later
if numel(u) < 4
    error('light line expects 4 numbers, actual got %d: %s', numel(u), line);
end
light_type = round(u(1));
[light_x, light_y, light_l] = deal(u(2), u(3), u(4));%%%% read RAT ---> Rsp, Rd, A, T %%%%line = fgetl(fid); if PRINTON    disp(line)endline = fgetl(fid);Rsp  = sscanf(line,'%f'); fprintf('\tRsp = %f\n', Rsp)line = fgetl(fid);Rd    = sscanf(line,'%f'); fprintf('\tRd = %f\n', Rd)line = fgetl(fid);A    = sscanf(line,'%f'); fprintf('\tA = %f\n', A)line = fgetl(fid);Td   = sscanf(line,'%f'); fprintf('\tT = %f\n', Td)line=fgetl(fid);% %%%% read nph% %%%line=fgetl(fid);line=fgetl(fid);NphR  = sscanf(line,'%d'); fprintf('\tNphR = %f\n', NphR)line = fgetl(fid);NphT    = sscanf(line,'%d'); fprintf('\tNphT = %f\n', NphT)line = fgetl(fid);NphOP    = sscanf(line,'%d'); fprintf('\tNphOP = %f\n', NphOP)%%%% Read A_l --> Al%%%%line = fgetl(fid);line = fgetl(fid); if PRINTON    disp(line)endfor i=1:Nlayers	line = fgetl(fid);	Al(i) = sscanf(line,'%f');	if PRINTON        fprintf('\tAl(%d) = %f\n',i,  Rsp)    endend%%%% Read A_z --> Az%%%%line = fgetl(fid); line = fgetl(fid); if PRINTON    disp(line)endAz = fscanf(fid, '%f');%%%% read Rd_r ---> Rr[Nr]%%%%line = fgetl(fid); if PRINTON    disp(line)endRr = fscanf(fid,'%f');%%%% read Rd_a ---> Ra[Na]%%%%line = fgetl(fid); disp(line)for i=1:Na	line = fgetl(fid);	Ra(i,1) = sscanf(line,'%f');end%%%% read Tt_r ---> Tr[Nr]%%%%line = fgetl(fid); line = fgetl(fid); if PRINTON    disp(line)endTr = fscanf(fid, '%f');%%%% read Tt_a ---> Ta[Na]%%%%line = fgetl(fid); if PRINTON    disp(line)endfor i=1:Na	line = fgetl(fid);	Ta(i,1) = sscanf(line,'%f');end%%%% read A_rz ---> Azr[Nz,Nr]%%%%line = fgetl(fid); line = fgetl(fid); if PRINTON    disp(line)endfor i=1:5; line = fgetl(fid); endu     = fscanf(fid,'%f');Azr   = reshape(u, [Nz Nr]);%%%% read Rd_ra ---> Rra[Nr, Na]%%%%line = fgetl(fid); if PRINTON    disp(line)endfor i=1:5; line = fgetl(fid); endu = fscanf(fid, '%f');Rra = reshape(u, [Nr Na]);%%%% read Tt_ra ---> Tra[Nr, Na]%%%%line = fgetl(fid); if PRINTON    disp(line)endfor i=1:5; line = fgetl(fid); endu = fscanf(fid, '%f');Tra = reshape(u, [Nr Na]);%%%% read OP ---> OP[Nr,Nz]%%%%line = fgetl(fid); line = fgetl(fid); if PRINTON    disp(line)endfor i = 1 : 5; line = fgetl(fid); endu     = fscanf(fid,'%f');OP = u;OP = [0;0;0;0;0;OP];OP   = reshape(OP, [Nz Nr]);OP = OP / NphOP;%%%% done%%%fclose(fid);clear fid%%%%% r, z%%%z     = ((1:Nz)' - 0.5)*dz;r     = ((1:Nr)' - 0.5)*dr;%%%% Azrm%%%Azrm  = zeros(Nz, 2*Nr-1);Azrm(:,Nr:2*Nr-1) = Azr;Azrm(:,1:Nr) = Azr(:,Nr:-1:1);rm = [((-Nr+1:-1) - 0.5)*dr (r'- dr)]' + dr/2;% rm was checked by plot(diff(rm)) and imagesc(rm, z, log10(Azrm))% % find ma(iz)% j = 1;% D = d(j);% for iz=1:Nz+1% 	if z(iz) < D% 		ma(iz) = mua(j);% 	else% 		j = j+1;% 		D = D + d(j);% 		ma(iz) = mua(j);% 	end% end% % Fzrm = Azrm*0; % set size of Fzrm% Fz = zeros(Nz,1);% for iz=1:Nz% 	Fzrm(iz,:) = Azrm(iz,:)/ma(iz);%     Fzr(iz,:) = Azr(iz,:)/ma(iz);% 	Fz(iz) = Az(iz)/ma(iz);% endif PLOTON % 0 = plotting turned OFF%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot Azr, Fzr %%%%%%%%%%%%%%%%% figure(2); clf;% figure(3); clf;% figure(4); clf;% % set(figure(1),'position', [10    15   700   800])% set(figure(2),'position', [735   476   600   360])% set(figure(3),'position', [10    15   600   360])% set(figure(4),'position', [735    15   600   360])figure(1); clf;u = (2:length(rm)-2);v = (1:Nz-1);xy = 32;figure(1)imagesc(r, z, log10(OP));set(gca, 'fontsize', xy)xlabel('r [cm]', 'fontsize', xy)ylabel('z [cm]', 'fontsize', xy)hold onrect_h = dz*5;        % keep your initial thickness
z_top  = z(1);
z_bot  = z(end);
gap    = 0;            % gap with surface, set to 0

% ===== Upper Surface R_PD: Draw outwards (up), size same as initial: width=Rl, height=rect_h =====
if SHOW_R_PD && exist('Rx','var') && exist('Rl','var')
    rectangle('Position', [Rx - Rl/2,  z_top - gap - rect_h,  Rl, rect_h], ...
        'EdgeColor','b','LineWidth',2, ...
        'FaceColor','b', ...
        'Clipping','off');
    plot(Rx, z_top, 'MarkerSize', 18, 'Clipping','off');  % Contact point
end

% ===== Lower Surface T_PD: Draw outwards (down), size same as initial: width=Tl, height=rect_h =====
if SHOW_T_PD && exist('Tx','var') && exist('Tl','var')
    rectangle('Position', [Tx - Tl/2,  z_bot + gap,  Tl, rect_h], ...
        'EdgeColor','b','LineWidth',2, ...
        'FaceColor','b', ...
        'Clipping','off');
    plot(Tx, z_bot, 'MarkerSize', 18, 'Clipping','off');
end

% ===== Light Source: Outside label; type=1 Inverse Triangle; type=2 w=0.05cm; type=3 w=0.10cm =====
if exist('light_type','var') && exist('light_x','var') && ~isempty(z)
    switch light_type
        case 1   % point -> solid inverse triangle; height same as PD thickness (rect_h)
            tri_h = rect_h;       % Same as PD thickness
            tri_w = 0.02;         % Horizontal visual width, adjustable
            P = [ light_x,               z_top; ...
                light_x - tri_w/2,     z_top - gap - tri_h; ...
                light_x + tri_w/2,     z_top - gap - tri_h ];
            patch('XData',P(:,1),'YData',P(:,2), ...
                'FaceColor','r','EdgeColor','r', ...
                'Clipping','off');
        case 2   % gaussian -> fixed width 0.05 cm, vertical thickness same as PD (rect_h)
            w = 0.05;  % cm
            rectangle('Position', [light_x - w/2,  z_top - gap - rect_h,  w, rect_h], ...
                'EdgeColor','r','LineWidth',2, ...
                'FaceColor','r', ...
                'Clipping','off');
        case 3   % flat -> fixed width 0.10 cm, vertical thickness same as PD (rect_h)
            w = 0.10;  % cm
            rectangle('Position', [light_x - w/2,  z_top - gap - rect_h,  w, rect_h], ...
                'EdgeColor','r','LineWidth',2, ...
                'FaceColor','r', ...
                'Clipping','off');
    end
    endhold offc = colorbar('eastoutside');set(c, 'TickLabels', {' ', ' '});colormap(makec2f)set(c, 'fontsize', 20)% title(c, 'Poids des photons', 'fontsize', 18);set(gca, 'FontWeight', 'Bold', 'FontSize', xy);%% tran% xlim([0.4 1.2]);% xticks(0.4:0.2:1.2);% xticklabels({'0.4', '0.6', '0.8', '1.0', '1.2'});% % % xticks(0:0.2:2.0);% % yticks(0:0.2:1.0);% yticklabels({'0', '0.2', '0.4', '0.6', '0.8', '1.0'});%% ref% xlim([0.5 1.5])% xticks(0.5:0.2:1.5);% xticklabels({'0.6', '0.8', '1.0', '1.2', '1.4'});% xticklabels({'0.5', '0.7', '0.9', '1.1', '1.3', '1.5'});% % % % xticks(0:0.2:2.0);% % % ylim([0 0.565]);% yticks(0:0.1:0.565)% yticklabels({'0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.565'});%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot Rr%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% figure(2)% area(r, Tra, 'DisplayName', 'Tra')% xlabel('r [cm]')% ylabel('T [ ]')% title('Transmittance en fonction de la répartition énergétique du signal à la surface du photodétecteur')% legend('Tr')% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % plot Az% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% figure(3)% area(r, Rra, 'DisplayName', 'Rra')% xlabel('r [cm]')% ylabel('R [ ]')% title('Réflectance en fonction de la répartition énergétique du signal à la surface du photodétecteur')% legend('Rr')% % figure(4)% imagesc(rm(u),z(v),log10(Azrm(v,u)));% set(gca,'fontsize',xy)% xlabel('r [cm]')% ylabel('z [cm]')% title(sprintf('log_1_0( Azr [J/cm^3] ),       Rd = %0.5f', Rd))% colorbar% colormap(makec2f)% set(colorbar,'fontsize',xy)% axis('equal')% xlim([0 inf])% ylim([0 inf])endclear D ans a b c2 h h2 i ii j line xy u v w w2 ifile