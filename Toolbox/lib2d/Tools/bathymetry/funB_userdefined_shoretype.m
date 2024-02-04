function Bathytype=funB_userdefined_shoretype(bathy,XX,YY,Hmin)
WaveChar=ones(size(bathy));
WaveChar(-bathy<Hmin)=zeros;


IndxZerofirstInX=funC_findIndexinColumn(WaveChar,'first');
IndxZerofirstInY=funC_findIndexinColumn(WaveChar.','first');
if all(IndxZerofirstInY)==0
    Bathytype='Userdefined (shore in x-axis)';
elseif all(IndxZerofirstInX)==0
    Bathytype='Userdefined (shore in y-axis)';
else
    Xdiff=max(XX(IndxZerofirstInX))-min(XX(IndxZerofirstInX));
    Ydiff=max(YY(IndxZerofirstInY))-min(YY(IndxZerofirstInY));
    if Xdiff < Ydiff
        Bathytype='Userdefined (shore in x-axis)';
    else
        Bathytype='Userdefined (shore in y-axis)';
    end
end
