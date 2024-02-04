function IVP=IVP_data(data,x)
IVP.eta=zeros(length(x),1);
IVP.u  =zeros(length(x),1);
x=x';

xdat=data(:,1);
etadat=data(:,2);
Ndat  = length(data(1,:));
if Ndat>=3
    udat  =data(:,3);
end
max_xmin=max(x(1),xdat(1));
indti=closest(x,max_xmin)+1;

min_xmax=min(x(end),xdat(end));
indtf=closest(x,min_xmax)-1;

IVP.eta(indti:indtf,1)=interp1(xdat,etadat,x(indti:indtf));
if Ndat>=3
IVP.u(indti:indtf,1)=interp1(xdat,udat,x(indti:indtf)); 
end
