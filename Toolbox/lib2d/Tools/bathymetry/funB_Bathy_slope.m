function bathyslope = funB_Bathy_slope(dom,bath)
 parslope = bath.par;

D_left      = parslope(1);              % depth_left 
D_right     = parslope(2);              % depth_right 
sf          = parslope(3);              % steepness fraction 
Xss         = parslope(4);              % start slope 
length_s    = abs(D_right -D_left)/sf;  % length of slope
Xse         = Xss + length_s;
%Xse         = Xss + length_s + 1;		% nunu
D_max = max(D_left,D_right);
D_min = min(D_left,D_right);

if strcmp(bath.type,'Slope (in x-axis)')
% bathyslope = ( D_left*ones(1,dom.Nx) + ...
            % -sign(D_left-D_right)*sf*(dom.X-Xss).*funC_Heaviside(dom.X-Xss) ).*funC_Heaviside((Xse+dom.dx)-dom.X)+ ...
            % + D_right.*funC_Heaviside(dom.X-(Xse+dom.dx));
bathyslope  =   D_left*ones(1,dom.Nx) - sign(D_left-D_right)*sf*(dom.X-Xss); % nunu
bathyslope  = min(D_max,max(D_min,bathyslope)); % cut min-max depth
% bathyslope(bathyslope==2*D_right)=D_right;
bathyslope  = -repmat(bathyslope,dom.Ny,1);        
else
% bathyslope = ( D_left*ones(1,dom.Ny) + ...
            % -sign(D_left-D_right)*sf*(dom.Y-Xss).*funC_Heaviside(dom.Y-Xss) ).*funC_Heaviside((Xse-dom.dy)-dom.Y)+ ...
            % + D_right.*funC_Heaviside(dom.Y-(Xse-dom.dy));   
bathyslope  =   D_left*ones(1,dom.Ny) - sign(D_left-D_right)*sf*(dom.Y-Xss); % nunu
bathyslope  = min(D_max,max(D_min,bathyslope)); % cut min-max depth
% bathyslope(bathyslope==2*D_right)=D_right;
bathyslope =-repmat(bathyslope',1,dom.Nx);      
end   
%[X,Y]=meshgrid(dom.X,dom.Y);
% figure
% surf(X,Y,-bathyslope,'edgecolor','none')   
% pause;

