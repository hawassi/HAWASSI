function Kwad=findG_kwadran_angle(theta1)
if theta1<=0 && theta1>-90
    Kwad=4;
elseif theta1<=-90 && theta1>-180
    Kwad=3;
elseif theta1<=-180 && theta1>-270
    Kwad=2;  
elseif theta1<=-270 && theta1>-360
    Kwad=2;    
elseif theta1>=0 && theta1<90
    Kwad=1;
elseif theta1>=90 && theta1<180
    Kwad=2;
elseif theta1>=180 && theta1<270
    Kwad=3;
elseif theta1>=270 && theta1<360
    Kwad=4;
elseif theta1>=360 && theta1<450
    Kwad=1;
elseif theta1>=450 && theta1<540
    Kwad=2;
end