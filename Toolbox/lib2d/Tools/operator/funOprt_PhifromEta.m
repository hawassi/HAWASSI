function phi=funOprt_PhifromEta(g,dom,model,meandepth,eta,propdir)

if model.current.check==0
Oprt_Om2d_fun= str2func(['funOprt_',model.dispersion]);
else
Oprt_Om2d_fun= str2func(['funOprt_',model.dispersion,'_current']);
end
if strcmpi(propdir,'South')
Mdir=[0 -1];    
elseif strcmpi(propdir,'North')
Mdir=[0 1];    
elseif strcmpi(propdir,'West')
Mdir=[-1 0];  
elseif strcmpi(propdir,'East')
Mdir=[1 0];
end
Om2d  =Oprt_Om2d_fun(dom.Kx,dom.Ky,meandepth,g,Mdir);
phihat=fft2(eta)./(1i.*Om2d);
phihat(abs(phihat)==Inf)=0;
phi=g*funC_ifft2(phihat);