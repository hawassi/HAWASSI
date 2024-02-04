function [zu0]=InitValURunUp(bathy,H_min,dx)
n           =length(bathy);
zu0         =zeros(n,1);

H           =zu0-bathy;
HeavH       =Heaviside(H-H_min);

indDI       =find(HeavH==0,1,'first');
indDF       =find(HeavH==0,1,'last');
dxbathy=gradient(bathy,dx);
zu0(indDI:indDF,1)=-dxbathy((indDI:indDF));
