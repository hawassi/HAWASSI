function [zeta0]=InitValEtaRunUp(bathy,H_min)
n           =length(bathy);
zeta0       =zeros(n,1);

H           =zeta0-bathy;
HeavH       =Heaviside(H-H_min);

indDI       =find(HeavH==0,1,'first');
indDF       =find(HeavH==0,1,'last');

zeta0(indDI:indDF)=bathy((indDI:indDF));

