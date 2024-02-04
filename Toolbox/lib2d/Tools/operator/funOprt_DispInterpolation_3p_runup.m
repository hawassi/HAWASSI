function [InterpDisp,H_minShore,H_min]= funOprt_DispInterpolation_3p_runup(g,influx,dom,Ug_fun,Up_fun,Om_fun,omAdd,nupeak,Depth,bath,...
                            Proj,dynmodel)

if strcmpi(dynmodel,'Hsdirect')
maxEtaInit= 1.5.*influx.par.Hs(1);
H_plus    = max(max(Depth))+maxEtaInit;
else
H_plus    = max(max(Depth));   
end
% nu        = mean(nupeak);

if strcmp(bath.type,'Shore (in x-axis)') || strcmp(bath.type,'Userdefined (shore in x-axis)')
    k_cut       = dom.kx(floor(length(dom.kx)/dom.cutfracx));  
%	k_hmin      = max(dom.kx)/dom.cutfracx ;  
    H_minShore  = bath.Hmin;%(3*nu/4/k_hmin)^2/g;
    H_min       = 5*H_minShore;
% 	k       = sort(reshape(sqrt(dom.Kx.^2+dom.Ky.^2),[],1),'ascend');
	k       	= dom.kx(1:floor(end/2));
  
elseif strcmp(bath.type,'Shore (in y-axis)') || strcmp(bath.type,'Userdefined (shore in y-axis)')
	k_cut       = dom.ky(floor(length(dom.ky)/dom.cutfracy));  
 %   k_hmin       = max(dom.ky)/dom.cutfracy ;  
    H_minShore  = bath.Hmin;%(3*nu/4/k_hmin)^2/g;
    H_min       = 5*H_minShore;
% 	k       = sort(reshape(sqrt(dom.Kx.^2+dom.Ky.^2),[],1),'ascend');
	k       	= dom.ky(1:floor(end/2));
end

if strcmpi(dynmodel,'Hs')
    Dperlambda=1/200; %% this is for very longwaves i.e T=100 s, if 1/20 (SWE) this leads to D~~200 (too much)
    KD=2*pi*Dperlambda;
    H_min=max(H_min,g*(KD*tanh(KD))/(influx.par.nu_p.^2)); % for windwave H_min=5*Hminshore it seems safe.
end

if ~isempty(bath.depthref)
    H_mid = bath.depthref;
else
    H_mid = (H_plus+H_min)/2;
end
Htot1  = linspace(H_min,H_mid,100);
Htot2  = linspace(H_mid,H_plus,100);
H_min1 = H_min;
H_plus1= H_mid;
H_mid1 = (H_min1+H_plus1)/2;
H_min2 = H_mid;
H_plus2= H_plus;
%H_mid2=H_min2+(H_plus2-H_min2)/6;
if strcmp(influx.par.category,'Deep water')
    Ddeep = invUpCgDeep(nupeak,H_min2,H_plus/300,0.51);
    if Ddeep==H_min2
    Ddeep=H_plus2;
    end
    H_mid2=H_min2+(Ddeep-H_min2)/6; %H_min2+(H_plus2-H_min2)/12;
   
else
    H_mid2=H_min2+(H_plus2-H_min2)/2;
end

wcut   = Om_fun(k_cut,H_plus,g,omAdd);

Ninterp     = length(Htot1);
Cp2_H1      = zeros(Ninterp,1);
Cp2_min1    = zeros(Ninterp,1);
Cp2_plus1   = zeros(Ninterp,1);
Cp2_mid1    = zeros(Ninterp,1);
Cp2kL_H1    = zeros(Ninterp,1);
Cp2kL_min1  = zeros(Ninterp,1);
Cp2kL_plus1 = zeros(Ninterp,1);
Cp2kL_mid1  = zeros(Ninterp,1);
Cp2IG_H1    = zeros(Ninterp,1);
Cp2IG_min1  = zeros(Ninterp,1);
Cp2IG_plus1 = zeros(Ninterp,1);
Cp2IG_mid1  = zeros(Ninterp,1);

Cp2_H2      = zeros(Ninterp,1);
Cp2_min2    = zeros(Ninterp,1);
Cp2_plus2   = zeros(Ninterp,1);
Cp2_mid2    = zeros(Ninterp,1);
Cp2kL_H2    = zeros(Ninterp,1);
Cp2kL_min2  = zeros(Ninterp,1);
Cp2kL_plus2 = zeros(Ninterp,1);
Cp2kL_mid2  = zeros(Ninterp,1);
Cp2IG_H2    = zeros(Ninterp,1);
Cp2IG_min2  = zeros(Ninterp,1);
Cp2IG_plus2 = zeros(Ninterp,1);
Cp2IG_mid2  = zeros(Ninterp,1);

knu1    = zeros(Ninterp,1);
knu2    = zeros(Ninterp,1);

for j= 1:Ninterp
    Hj1 = Htot1(j);
    Hj2 = Htot2(j);
    
    kappanu_H1   = funOprt_invOmExactVal(nupeak,Hj1,g);%interp1(Om_fun(k,Hj1,g,omAdd),k,nupeak,'spline');
    kappanuIG_H1 = funOprt_invOmExactVal(nupeak/10,Hj1,g);%interp1(Om_fun(k,Hj1,g,omAdd),k,nupeak/10,'spline');
    kappanuL_H1  = funOprt_invOmExactVal(wcut,Hj1,g);%interp1(Om_fun(k,Hj1,g,omAdd),k,wcut,'spline');
    kappanu_H2   = funOprt_invOmExactVal(nupeak,Hj2,g);%interp1(Om_fun(k,Hj2,g,omAdd),k,nupeak,'spline');
    kappanuIG_H2 = funOprt_invOmExactVal(nupeak/10,Hj2,g);%interp1(Om_fun(k,Hj2,g,omAdd),k,nupeak/10,'spline');
    kappanuL_H2  = funOprt_invOmExactVal(wcut,Hj2,g);%interp1(Om_fun(k,Hj2,g,omAdd),k,wcut,'spline');
    
    knu1(j) =  kappanu_H1;
    knu2(j) =  kappanu_H2;

    Cp2_min1(j) = Up_fun(kappanu_H1,H_min1,g,omAdd).^2;
    Cp2_plus1(j)= Up_fun(kappanu_H1,H_plus1,g,omAdd).^2;
    Cp2_mid1(j) = Up_fun(kappanu_H1,H_mid1,g,omAdd).^2;
    Cp2_H1(j)   = Up_fun(kappanu_H1,Hj1,g,omAdd).^2;

    Cp2kL_min1(j) = Up_fun(kappanuL_H1,H_min1,g,omAdd).^2;% sqrt(g*H_min) ;%
    Cp2kL_plus1(j)= Up_fun(kappanuL_H1,H_plus1,g,omAdd).^2;%%sqrt(g*H_plus) ;%
    Cp2kL_mid1(j) = Up_fun(kappanuL_H1,H_mid1,g,omAdd).^2;%sqrt(g*H_mid) ;%
    Cp2kL_H1(j)   = Up_fun(kappanuL_H1,Hj1,g,omAdd).^2;% sqrt(g*Hj) ;%    
    
    Cp2IG_min1(j) = Up_fun(kappanuIG_H1,H_min1,g,omAdd).^2;%
    Cp2IG_plus1(j)= Up_fun(kappanuIG_H1,H_plus1,g,omAdd).^2;%%
    Cp2IG_mid1(j) = Up_fun(kappanuIG_H1,H_mid1,g,omAdd).^2;%
    Cp2IG_H1(j)   = Up_fun(kappanuIG_H1,Hj1,g,omAdd).^2;%

    Cp2_min2(j) = Up_fun(kappanu_H2,H_min2,g,omAdd).^2;
    Cp2_plus2(j)= Up_fun(kappanu_H2,H_plus2,g,omAdd).^2;
    Cp2_mid2(j) = Up_fun(kappanu_H2,H_mid2,g,omAdd).^2;
    Cp2_H2(j)   = Up_fun(kappanu_H2,Hj2,g,omAdd).^2;

    Cp2kL_min2(j) = Up_fun(kappanuL_H2,H_min2,g,omAdd).^2;% sqrt(g*H_min) ;%
    Cp2kL_plus2(j)= Up_fun(kappanuL_H2,H_plus2,g,omAdd).^2;%%sqrt(g*H_plus) ;%
    Cp2kL_mid2(j) = Up_fun(kappanuL_H2,H_mid2,g,omAdd).^2;%sqrt(g*H_mid) ;%
    Cp2kL_H2(j)   = Up_fun(kappanuL_H2,Hj2,g,omAdd).^2;% sqrt(g*Hj) ;%

    Cp2IG_min2(j) = Up_fun(kappanuIG_H2,H_min2,g,omAdd).^2;%
    Cp2IG_plus2(j)= Up_fun(kappanuIG_H2,H_plus2,g,omAdd).^2;%%
    Cp2IG_mid2(j) = Up_fun(kappanuIG_H2,H_mid2,g,omAdd).^2;%
    Cp2IG_H2(j)   = Up_fun(kappanuIG_H2,Hj2,g,omAdd).^2;%
end

A1= (Cp2kL_mid1.*Cp2IG_plus1)-(Cp2kL_plus1.*Cp2IG_mid1);
B1=-(Cp2kL_min1.*Cp2IG_plus1)+(Cp2kL_plus1.*Cp2IG_min1);
C1= (Cp2kL_min1.*Cp2IG_mid1)-(Cp2kL_mid1.*Cp2IG_min1);
D1=-(Cp2_mid1.*Cp2IG_plus1)+(Cp2_plus1.*Cp2IG_mid1);
E1= (Cp2_min1.*Cp2IG_plus1)-(Cp2_plus1.*Cp2IG_min1);

F1=-(Cp2_min1.*Cp2IG_mid1)+(Cp2_mid1.*Cp2IG_min1);
G1= (Cp2_mid1.*Cp2kL_plus1)-(Cp2_plus1.*Cp2kL_mid1);
H1=-(Cp2_min1.*Cp2kL_plus1)+(Cp2_plus1.*Cp2kL_min1);
I1= (Cp2_min1.*Cp2kL_mid1)-(Cp2_mid1.*Cp2kL_min1);

det1=Cp2_min1.*A1+Cp2_mid1.*B1+Cp2_plus1.*C1;

gam_min1  = ( A1.*Cp2_H1+D1.*Cp2kL_H1+G1.*Cp2IG_H1 )./det1;
gam_mid1  = ( B1.*Cp2_H1+E1.*Cp2kL_H1+H1.*Cp2IG_H1 )./det1;
gam_plus1 = ( C1.*Cp2_H1+F1.*Cp2kL_H1+I1.*Cp2IG_H1 )./det1;

A2= (Cp2kL_mid2.*Cp2IG_plus2)-(Cp2kL_plus2.*Cp2IG_mid2);
B2=-(Cp2kL_min2.*Cp2IG_plus2)+(Cp2kL_plus2.*Cp2IG_min2);
C2= (Cp2kL_min2.*Cp2IG_mid2)-(Cp2kL_mid2.*Cp2IG_min2);
D2=-(Cp2_mid2.*Cp2IG_plus2)+(Cp2_plus2.*Cp2IG_mid2);
E2= (Cp2_min2.*Cp2IG_plus2)-(Cp2_plus2.*Cp2IG_min2);

F2=-(Cp2_min2.*Cp2IG_mid2)+(Cp2_mid2.*Cp2IG_min2);
G2= (Cp2_mid2.*Cp2kL_plus2)-(Cp2_plus2.*Cp2kL_mid2);
H2=-(Cp2_min2.*Cp2kL_plus2)+(Cp2_plus2.*Cp2kL_min2);
I2= (Cp2_min2.*Cp2kL_mid2)-(Cp2_mid2.*Cp2kL_min2);

det2=Cp2_min2.*A2+Cp2_mid2.*B2+Cp2_plus2.*C2;

gam_min2  = ( A2.*Cp2_H2+D2.*Cp2kL_H2+G2.*Cp2IG_H2 )./det2;
gam_mid2  = ( B2.*Cp2_H2+E2.*Cp2kL_H2+H2.*Cp2IG_H2 )./det2;
gam_plus2 = ( C2.*Cp2_H2+F2.*Cp2kL_H2+I2.*Cp2IG_H2 )./det2;

HtotComb    = [Htot1,Htot2(2:end)];  
indHmin1    = closest(HtotComb,H_min1);
indHplus1   = closest(HtotComb,H_plus1);
indHplus2   = closest(HtotComb,H_plus2);

gam_min1tot = zeros(size(HtotComb));
gam_mid1tot = zeros(size(HtotComb));
gam_plus1tot= zeros(size(HtotComb));
gam_min2tot = zeros(size(HtotComb));
gam_mid2tot = zeros(size(HtotComb));
gam_plus2tot= zeros(size(HtotComb));

gam_min1tot(indHmin1:indHplus1)     = gam_min1;
gam_mid1tot(indHmin1:indHplus1)     = gam_mid1;
gam_plus1tot(indHmin1:indHplus1)    = gam_plus1;
gam_min2tot(indHplus1+1:indHplus2)  = gam_min2(2:end).';
gam_mid2tot(indHplus1+1:indHplus2)  = gam_mid2(2:end).';
gam_plus2tot(indHplus1+1:indHplus2) = gam_plus2(2:end).';
if strcmpi(dynmodel,'Hsdirect')
gam_plus1   = spline(HtotComb,gam_plus1tot);
gam_min1    = spline(HtotComb,gam_min1tot);
gam_mid1    = spline(HtotComb,gam_mid1tot);
gam_plus2   = spline(HtotComb,gam_plus2tot);
gam_min2    = spline(HtotComb,gam_min2tot);
gam_mid2    = spline(HtotComb,gam_mid2tot);

InterpDisp.sp_gam_plus1 = gam_plus1;
InterpDisp.sp_gam_min1  = gam_min1;
InterpDisp.sp_gam_mid1  = gam_mid1;
InterpDisp.sp_gam_plus2 = gam_plus2;
InterpDisp.sp_gam_min2  = gam_min2;
InterpDisp.sp_gam_mid2  = gam_mid2;
else
InterpDisp.sp_gam_plus1 = interp1(HtotComb,gam_plus1tot,Depth,'spline');
InterpDisp.sp_gam_min1  = interp1(HtotComb,gam_min1tot,Depth,'spline');
InterpDisp.sp_gam_mid1  = interp1(HtotComb,gam_mid1tot,Depth,'spline');
InterpDisp.sp_gam_plus2 = interp1(HtotComb,gam_plus2tot,Depth,'spline');
InterpDisp.sp_gam_min2  = interp1(HtotComb,gam_min2tot,Depth,'spline');
InterpDisp.sp_gam_mid2  = interp1(HtotComb,gam_mid2tot,Depth,'spline');
InterpDisp.sp_gam_plus1(Depth<H_min)=0;
InterpDisp.sp_gam_plus2(Depth<H_min)=0;
InterpDisp.sp_gam_min1(Depth<H_min)=0;
InterpDisp.sp_gam_min2(Depth<H_min)=0;
InterpDisp.sp_gam_mid1(Depth<H_min)=0;
InterpDisp.sp_gam_mid2(Depth<H_min)=0;
% figure;
% subplot(4,1,1)
% surf(dom.XX,dom.YY,Depth,'edgecolor','none');view(2);colorbar;
% subplot(4,1,2)
% surf(dom.XX,dom.YY,InterpDisp.sp_gam_min2+InterpDisp.sp_gam_min1,'edgecolor','none');view(2);colorbar;
% subplot(4,1,3)
% surf(dom.XX,dom.YY,InterpDisp.sp_gam_mid2+InterpDisp.sp_gam_mid1,'edgecolor','none');view(2);colorbar;
% subplot(4,1,4)
% surf(dom.XX,dom.YY,InterpDisp.sp_gam_plus2+InterpDisp.sp_gam_plus1,'edgecolor','none');view(2);colorbar;
end
InterpDisp.H_min1  = H_min1;
InterpDisp.H_mid1  = H_mid1;
InterpDisp.H_plus1 = H_plus1;
InterpDisp.H_min2  = H_min2;
InterpDisp.H_mid2  = H_mid2;
InterpDisp.H_plus2 = H_plus2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%