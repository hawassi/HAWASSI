function bathy=bathy_data(bath,x,cfSA)
x=x';
cfSA=cfSA';

xdat=bath.data(:,1);
bathdat=bath.data(:,2);
bathy=interp1(xdat,bathdat,x,'linear','extrap'); 

if any(bathdat>0)
Nx=length(x);
if bathy(Nx,1)>0  
bathy(floor(Nx/2):Nx,1)=bathy(floor(Nx/2):Nx,1); 
elseif bathy(1)>0 
bathy(1:floor(Nx/2),1)=bathy(1:floor(Nx/2):Nx,1);     
end
cfSA(bathy<0)=1;
bathy=bathy.*cfSA;
end
