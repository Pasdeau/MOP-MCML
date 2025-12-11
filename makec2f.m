function c2 = makec2f

% function c2 = makec2f
% 	Creates a colormap.
% USAGE:   colormap(makec2f)
% by Steven L. Jacques, Jan. 2008
% Oregon Health & Sciences University, Portland, OR, USA

% ----> c2(64,3)
%  red 
%  yellow white
%  green yellow
%  blue  red 
%  black black
c2 = zeros(64,3);

for ii=1:20   % black to red
	c2(ii,1) = (ii-1)./19;    % red
	c2(ii,2) = 0;    % green
	c2(ii,3) = 0;     % blue
end
for ii=21:48   % red to yellow
	c2(ii,1) = 1;    % red
	c2(ii,2) = (ii-21)./27;    % green
	c2(ii,3) = 0; % blue
end
for ii=49:64 % yellow to white
	c2(ii,1) = 1;   
	c2(ii,2) = 1;   
	c2(ii,3) = (ii-49)./15; 
end
