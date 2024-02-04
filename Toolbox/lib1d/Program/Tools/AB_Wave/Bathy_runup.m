function [Depthrunup,parslope] = Bathy_runup(x,xinflux,parslope,cfSA)
bath    =zeros(length(x),1);
Depth   =parslope(1);
slope   =parslope(2);              % steepness fraction
xshore  =parslope(3);
cfSAnow=cfSA.';

if xinflux>xshore, slope=-abs(slope);end
xb=xshore-Depth/slope;


indxb=closest(x,xb);

if sign(slope)>0
    if xb>x(1)
        bath(1:indxb)=-Depth;
        bath(indxb+1:end)=(slope.*(x(indxb+1:end)-xb)-Depth);
    else
        bath(1:end)=(slope.*(x(1:end)-xb)-Depth);
    end
    
else
    if xb<x(end)
        bath(indxb:end)=-Depth;
        bath(1:indxb+1)=(slope.*(x(1:indxb+1)-xb)-Depth);
    else
        bath(1:end)=(slope.*(x(1:end)-xb)-Depth);
        
    end
end
cfSAnow(bath<0)=1;
bath=bath.*cfSAnow;
Depthrunup=-bath;
% figure;
% plot(x,bath)

parslope(4)=indxb;
parslope(5)=closest(x,xshore);


