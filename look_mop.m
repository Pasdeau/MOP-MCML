%   look_mop.m Reads the MCML output file specified (eg., example.mco)
%	and uses the subroutine <get_mop> to read the output file.
%	and return values and MOP.
%	If PRINTON == 1 or PLOTON == 1, then getmcml also provides output. 
%   Uses getmcml.m, which in turn uses makec2f.m
%   by Steven L. Jacques, Jan. 2008
%   Oregon Health & Sciences University, Portland, OR, USA


% PLOTON = 1; PRINTON = 1;    % 0 = off, 1 = on.  Controls plotting option
% 
% get_mop('output.mco')

S = mop_read_mco('450.mco','Verbose',true);
figs = mop_plot_mco(S,'ShowRPD',true,'ShowTPD',true,'CLim',[-6 -1]);
