function par=funSA_preparation(wave,ship,bottom,calc)
g=9.811;

shipprof=funSA_shiplib(ship);

if ship.shape==4
    draft=abs(min(shipprof(:,2)));
    B=max(shipprof(:,1))-min(shipprof(:,1));
    par.S4draft=draft;
    par.S4B=B;    
else
    B=ship.width;
    draft=ship.draft;
end
Ndom=funC_Ndom(ship.Nx);

b=B./2;
NdxS=ship.Nx;

if bottom.slope>0
    IdSlope=1;
else
    IdSlope=0;
end

D=bottom.depth;


if wave.inputVar==1
wnorm=wave.input;          %wave.input=ww sqrt(B/2g)
w0=wnorm/sqrt(B/2/g);
k0=invOmExact(w0,D);
w0=w0.';
elseif wave.inputVar==2
wnorm=wave.input;          %wave.input=ww^2 B/2g
w0   =sqrt(2*g*wnorm/B);
k0=invOmExact(w0,D);
w0=w0.';
elseif wave.inputVar==3
lambda=wave.input.'*D;    %wave.input=lambda/D
k0=2*pi./lambda;   
w0=OmExact(k0,D);
elseif wave.inputVar==4
lambda=wave.input.'*B;    %wave.input=lambda/B
k0=2*pi./lambda;   
w0=OmExact(k0,D);
elseif wave.inputVar==5
k0=wave.input.'./D;      %wave.input=kD
w0=OmExact(k0,D);
elseif wave.inputVar==6
k0=wave.input.'./B;      %wave.input=kB
w0=OmExact(k0,D);
elseif wave.inputVar==7
k0=wave.input.'./b;      %wave.input=kb
w0=OmExact(k0,D);
elseif wave.inputVar==8
k0=wave.input.'./draft;      %wave.input=kT
w0=OmExact(k0,D);
end

Nk=length(w0);

IdDraftMax=closest(shipprof(:,2),min(shipprof(:,2)));
XshipL=sort(shipprof(1:IdDraftMax,1),'descend');
ZshipL=sort(shipprof(1:IdDraftMax,2));
XshipR=shipprof(IdDraftMax:end,1);
ZshipR=shipprof(IdDraftMax:end,2);

dxS=B/NdxS;
xx=(-20*B:B/NdxS:20*B).';
dxx=xx(2)-xx(1);
kk=fun_freqspace(xx.');
%     Cp=UpExact(kk,D);
%     OpF0=Cp.^2.*1i.*kk/g;
%     OpL0=-1i.*kk.*Cp.^2.*1i.*kk/g;
bf0=[B B];
cfSA    = cf(xx,xx(1),bf0(1))-cf(xx,xx(end)-bf0(2),bf0(2)); % cfSimulationArea

Dxx=ones(size(xx));
if IdSlope==1
    Lslope=B+2*dxS;
    Dslope=Lslope*bottom.slope;
    DR=D-Dslope/2; % water depth at Right bdy
    DL=D+Dslope/2; % water depth at Left bdy
    XslopeStart=-Lslope/2;% slope start position
    
    Dxx =-Bathy_slope(xx,[DL DR bottom.slope XslopeStart]).';
else
    Dxx=Dxx.*D;
end
dxDxx=gradient(Dxx,xx(2)-xx(1));%

% figure;
% plot(xx,Dxx,'r',xx,dxDxx,'b')

shipform=zeros(size(xx));

indx1=closest(xx,shipprof(1,1));
indx2=closest(xx,shipprof(end,1));

if xx(indx1)<shipprof(1,1)
    indx1=indx1+1;
end
%
if xx(indx2)>shipprof(end,1)
    indx2=indx2-1;
end

shipform(indx1:indx2)=interp1(shipprof(:,1),shipprof(:,2),xx(indx1:indx2));
% shipform(indx1)=0;%shipform(indx1+1);
% shipform(indx2)=0;%shipform(indx2-1);

if shipform(indx1)~=0
    indx1=indx1-1;
end
if shipform(indx2)~=0
    indx2=indx2+1;
end


indc=find(shipform==min(shipform));
Xc=xx(indc);
xShip=[xx(indx1) Xc xx(indx2)];
chiS=zeros(size(xx));
chiS(indx1+1:indx2-1)=1;
chiSwl=zeros(size(xx));
chiSwl(indx1:indx2)=1;


dzX_L=funOprt_FDGradient1d(XshipL,ZshipL(2)-ZshipL(1),1);
dzX_R=funOprt_FDGradient1d(XshipR,ZshipR(2)-ZshipR(1),1);


fun_Om=str2func('OmExact');
fun_Up=str2func('UpExact');

%         if IDTestcase~=13
%             DxL=ones(size(XshipL)).*D;
%             DxR=ones(size(XshipR)).*D;
%             dxDxL=gradient(DxL)./gradient(XshipL);
%             dxDxR=gradient(DxR)./gradient(XshipR);
%         end
dxZshipL=gradient(ZshipL)./gradient(XshipL);
dxZshipR=gradient(ZshipR)./gradient(XshipR);


indxBdy=zeros(1,Ndom-1);
indxBdy(1)=indx1;indxBdy(end)=indx2;
indxBdy(Ndom/2)=indc;
dIndL=floor((indc-indx1)/((Ndom-2)/2));
dIndR=floor((indx2-indc)/((Ndom-2)/2));

for ii=2:Ndom/2-1
    indxBdy(Ndom/2-(ii-1))=indc-dIndL*(ii-1);
    indxBdy(Ndom/2+(ii-1))=indc+dIndR*(ii-1);
end


%
% indxlc=indx1+floor((indc-indx1)/6);
% indxcr=closest(xx,-xx(indxlc));
%
% indxBdy=[indx1 indxlc indc indxcr indx2];

shipform_aprox=shipform;
for ii=1:Ndom/2-1
    shipform_aprox(indxBdy(ii)+1:indxBdy(ii+1))= shipform(indxBdy(ii+1));
    shipform_aprox(indxBdy(end-ii):indxBdy(end-ii+1)-1)= shipform(indxBdy(end-ii));
end


AreaO=-trapz(shipprof(:,1),shipprof(:,2));
Area=-trapz(xx(indx1:indx2),shipform_aprox(indx1:indx2));

Xbdy=xx(indxBdy);
dr=zeros(size(indxBdy));

for ii=1:(Ndom-2)/2
    dr(ii)=-shipform_aprox(indxBdy(ii+1));% shipform(indc) shipform(indc) shipform(indc) shipform(indxcr+1)];
    dr(end-ii+1)=-shipform_aprox(indxBdy(end-ii));
end
dr(Ndom/2)=-shipform_aprox(indxBdy(Ndom/2));



Dbdy=Dxx(indxBdy);%D*ones(size(indxBdy));%
dxDbdy=dxDxx(indxBdy);
dxDbdy(1)=0;
dxDbdy(end)=0;

if wave.amplVar==1   %wave.ampl=A
    kAinc=k0*wave.ampl;
    Ainc=wave.ampl*ones(size(kAinc));
else                 %wave.ampl=kA
    Ainc=wave.ampl./k0;
    kAinc=wave.ampl*ones(size(Ainc));
end

M=AreaO;
Zc=ship.cog(2);
% if ship.shape==1
%     I    = B*draft^3/3 +B*draft*(Zc)^2;
% elseif ship.shape==2
%     I    = pi.*draft.^4/8 +M*Zc^2;
% else
%     datS.x=shipprof(:,1);datS.y=shipprof(:,2);
%     I=MOMENT(datS,2,0)  ;
% end
% I
I=M*ship.gyradius^2;

M_mat=diag([M M I]);

par.M_mat=M_mat;
par.kAinc=kAinc;
par.Ainc=Ainc;
par.Nk=Nk;
par.bottom.IdSlope=IdSlope;
par.shipform=shipform;
par.shipform_aprox=shipform_aprox;
par.k0=k0;
par.w0=w0;
par.Dbdy=Dbdy;
par.g=g;
par.Dxx=Dxx;
par.dxx=dxx;
par.xx=xx;
par.indxBdy=indxBdy;
par.Ndom=Ndom;
par.dxDbdy=dxDbdy;
par.dr=dr;
par.Xbdy=Xbdy;
par.chiS=chiS;
par.XshipL=XshipL;
par.XshipR=XshipR;
par.dxZshipL=dxZshipL;
par.dxZshipR=dxZshipR;
par.ZshipL=ZshipL;
par.ZshipR=ZshipR;
par.Xc    =Xc;
par.AreaO=AreaO;
par.indx1=indx1;
par.indx2=indx2;
par.chiSwl=chiSwl;
par.dzX_L=dzX_L;
par.dzX_R=dzX_R;

%for progress bar
tic;
[statusbarObj]=JavaFrame_handling();
statusbarObj.setText('');
JavProgressBar;
[jbStop]=Java_stopbutton(statusbarObj);
set(jProgressBar,'Maximum',Nk, 'Value',0);
jProgressBar.setStringPainted( true );
ETA=remain_time(1,Nk);
par.jProgressBar=jProgressBar;
par.jbStop=jbStop;
par.statusbarObj=statusbarObj;
%
if calc.Id==1
    figure;
    plot(shipprof(:,1),shipprof(:,2),'b',xx,shipform,'r',xx(indxBdy),shipform(indxBdy),'or',xx,shipform_aprox,'--r',xx,-Dxx,'k');
    xlabel('x[m]'); ylabel('y[m]');
    
    xlim([-b-dxx*5 b+dxx*5]);
    %ylim([-draft 1]);
    % axis equal;
    title(['Area (b)=',num2str(AreaO),'Area (r)=',num2str(Area)])
    plot_properties;
end
end

