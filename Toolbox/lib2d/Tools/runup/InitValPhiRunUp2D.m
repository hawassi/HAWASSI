function [zphi0]=InitValPhiRunUp2D(bathy,H_min,dx,dy)
n           =size(bathy);
zphi0       =zeros(n);

H           =zphi0-bathy;
HeavH       =Heaviside(H-H_min);

indDI       =funC_findIndexinColumn(HeavH,'first');
indDF       =funC_findIndexinColumn(HeavH,'last');

[dxbathy,~]=gradient(bathy,dx,dy);

zphi0(:,indDI:indDF)=-dxbathy(:,(indDI:indDF));
