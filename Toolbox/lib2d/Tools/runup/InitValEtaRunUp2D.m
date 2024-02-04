function [zeta0]=InitValEtaRunUp2D(bathy,H_min)
n           =size(bathy);
zeta0       =zeros(n);

H           =zeta0-bathy;
HeavH       =Heaviside(H-H_min);

indDI       =funC_findIndexinColumn(HeavH,'first');
indDF       =funC_findIndexinColumn(HeavH,'last');

zeta0(:,indDI:indDF)=bathy(:,(indDI:indDF));