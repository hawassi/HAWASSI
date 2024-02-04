function [cfSA,fbl]=funD_simulationAreaChar(dom,fbl)
dampdata=fbl.param;
cfSA=ones(size(dom.XX));
fbl.N=length(dampdata(:,1));
fblChar=zeros(size(dom.XX));
LL=7*ones(size(dom.XX));

for i=1:fbl.N
    if strcmp(dampdata(i,1),'Fourier Bdry.')
        fbl.l=cell2mat(dampdata(i,2));fbl.r=cell2mat(dampdata(i,3));
        fbl.b=cell2mat(dampdata(i,4));fbl.t=cell2mat(dampdata(i,5));
        cfSAii=funC_cfSA2d(dom.X,dom.Y,fbl); 
        
        indxl =funC_closest(dom.X,dom.X(1)+fbl.l); indxr   =funC_closest(dom.X,dom.X(end)-fbl.r);
        indyb =funC_closest(dom.Y,dom.Y(1)+fbl.b); indyt   =funC_closest(dom.Y,dom.Y(end)-fbl.t);
        
        if fbl.b==0&& fbl.t==0
            if fbl.l~=0
                LL(:,1:indxl)=fbl.l;
            else
                LL(:,1:indxl)=1;
            end
            if fbl.r~=0
                LL(:,indxr:end)=fbl.r;
            else
                LL(:,indxr:end)=1;
            end
        elseif fbl.l==0&& fbl.r==0
            if fbl.b~=0
                LL(1:indyb,:)=fbl.b;
            else
                LL(1:indyb,:)=1;
            end
            if fbl.t~=0
                LL(indyt:end,:)=fbl.t;
            else
                LL(indyt:end,:)=1;
            end
        else
            LL(:,1:indxl)=fbl.l; 
            LL(:,indxr:end)=fbl.r;
            LL(1:indyb,indxl+1:indxr-1)=fbl.b;
            LL(indyt:end,indxl+1:indxr-1)=fbl.t;
        end
        fblCharii=(1-cfSAii)./(0.9*LL); % factor 0.9 to damp faster before end of the domain
    else
        dampuser=fbl.userdata(i).bdry;
        smoothfact=cell2mat(dampdata(i,6));
        cfSAii  = funD_fourier_boundary_userdefined(dom,dampuser,smoothfact);
%       xx=dampuser(:,1);
%       yy=dampuser(:,2);
%       Lud=funC_minradius(xx,yy)
        %Lud=7;%
        xmin=min(dampuser(:,1));xmax=max(dampuser(:,1)); 
        ymin=min(dampuser(:,2));ymax=max(dampuser(:,2));
        if ymax-ymin<xmax-xmin
           dS=dom.dy; 
        else
           dS=dom.dx; 
        end
        Lud=smoothfact.*dS;
        fblCharii=(1-cfSAii)./(0.9*Lud); % factor 0.9 to damp faster before end of the domain
    end
    cfSA=cfSA.*cfSAii;
    fblChar=fblChar+fblCharii;
end
fbl.char=fblChar;


% figure(11)
% set(gcf,'Renderer','zbuffer'); %due to graphics driver
% surf(dom.XX,dom.YY,fblChar,'edgecolor','none');
% 
% figure(12)
% set(gcf,'Renderer','zbuffer'); %due to graphics driver
% surf(dom.XX,dom.YY,LL,'edgecolor','none');


