% mainABflat, version May 2011
% freqspace
% % % % Software developed by 
% % % % LabMath-Indonesia & University of Twente
%%
function [k,N] = funC_freqspace(x,N,kind)
if nargin ==1 
    N=length(x); 
    kind = 'default';
elseif nargin ==2
    kind = 'default';
end

dx = x(2)-x(1);%diff(x);
% if std(dx)>1e-10,
%     error('First input should be a linear spaced vector');
% else
%     dx= mean(dx); 
% end
L = x(end)-x(1);%
dk = 2*pi/L;
fmin = -0.5*2*pi/dx;
if rem(N,2)==0
    fmax = -fmin-dk;
elseif rem(N,2)==1
    fmax = -fmin;
end

switch kind
    case 'default'
        k = ifftshift(linspace(fmin,fmax,N)');%
        k(1)=1e-10;
    case 'single sided'
        k = linspace(0,-fmin,(floor(N/2)+1)*2*pi);%factor 2pi
    case 'unshifted'
        k = linspace(fmin,fmax,N)';
end
       
  