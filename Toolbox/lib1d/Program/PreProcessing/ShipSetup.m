% !--------------------------------------------------------------------------------------
% !
% !    Copyright (C) 2024 - LabMath-Indonesia
% !
% !    This program is free software: you can redistribute it and/or modify
% !    it under the terms of the GNU General Public License as published by
% !    the Free Software Foundation, either version 3 of the License, or
% !    (at your option) any later version.
% !
% !    This program is distributed in the hope that it will be useful,
% !    but WITHOUT ANY WARRANTY; without even the implied warranty of
% !    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% !    GNU General Public License for more details.
% !
% !    You should have received a copy of the GNU General Public License
% !    along with this program.  If not, see <http://www.gnu.org/licenses/>.
% !
% !   Contributors list:
% !   - R. Kurnia
% !   - E. van Groesen
% !--------------------------------------------------------------------------------------

%%%%%%%Ship-configuration%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%Ship-Shape%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shippar.Nship=length(shippar.data(:,1));

FlagAllSFixed=1;
for i=1:shippar.Nship
    try
        SXc(i)=str2num(cell2mat(shippar.data(i,5)));
    catch
        SXc(i)=cell2mat(shippar.data(i,5));
    end
    try
        SZc(i)=str2num(cell2mat(shippar.data(i,6)));
    catch
        SZc(i)=cell2mat(shippar.data(i,6));
    end
    try
        SThetac(i)=ded2rad(str2num(cell2mat(shippar.data(i,7))));
    catch
        SThetac(i)=deg2rad(cell2mat(shippar.data(i,7)));
    end
    
    if  ~strcmpi(shippar.data(i,2),'Fixed')
    FlagAllSFixed=0;
    end
end

[shippar.form]= shapeformsetup(x,shippar,SXc,SZc,SThetac,bath.name);

nup=influx.nu_p;
gradDepth=gradient(-par.bathy,par.dx);
kappaS=zeros((shippar.Evmode+1),shippar.Nship*2);
Dwl=zeros(shippar.Nship,2);
dxDwl=zeros(shippar.Nship,2);
iterx=1;
shipform=shippar.form.shapeXcZc(end,:).';
Sdampcoef=zeros(length(par.x),1);
for jj=1:shippar.Nship
    indXc=closest(par.x,shippar.form.xShip(jj,2));
    SXc=par.x(indXc);
    indXL=closest(par.x,shippar.form.xShip(jj,1));
    indXR=closest(par.x,shippar.form.xShip(jj,3));
    DL=-par.bathy(indXL);DR=-par.bathy(indXR);
    dxDL=gradDepth(indXL);dxDR=gradDepth(indXR);
    kpL=invOmExact(nup,DL);
    kpR=invOmExact(nup,DR);
    sig0=nup.^2./g;
    if shippar.Evmode>=1
        kappaDL=invOmEvanescent(sig0*DL,shippar.Evmode)';
        kappaDR=invOmEvanescent(sig0*DR,shippar.Evmode)';
        kappaL=[-1i*kpL kappaDL./DL].';
        kappaR=[-1i*kpR kappaDR./DR].';
       
        kappaS(:,iterx:iterx+1)=[kappaL kappaR];
    else
        kappaS(1,iterx:iterx+1)=[-1i*kpL -1i*kpR];
    end
    Dwl(jj,:)=[DL DR];
    dxDwl(jj,:)=[dxDL dxDR];
    iterx=iterx+2;
    
    shipform(indXL)=shipform(indXL+1); %%% % due to the barge make sure the flux at ship bottom at water line is not zero.
    shipform(indXR)=shipform(indXR-1);
    
    DS=(DL+DR)/2;
    kpS=invOmExact(nup,DS);
    CpS=UpExact(kpS,DS);
    LS=shippar.form.Slength(jj);
%     if  2*pi/kpS/LS<10
     Sdampcoef(indXL:indXR)=7*CpS/(LS/3); 
%     else
%    Sdampcoef(indXL:indXR)=0.5*CpS/(LS/3);
%    end
end
shippar.form.kappa=kappaS;
shippar.form.Dwl=Dwl;
shippar.form.dxDwl=dxDwl;
shippar.form.Sdampcoef=Sdampcoef;



%%% Initial condition
if IVP.type~=1
    shippar.init.phi0=Ifft(IVP.zu0_hat./(1i.*k));
    shippar.init.zeta0       = Ifft(IVP.zeta0_hat);
else
    shippar.init.zeta0       =zeros(size(k));%
    shippar.init.phi0        =zeros(size(shippar.init.zeta0));
end


