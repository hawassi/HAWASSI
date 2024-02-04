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

global Idstop
[jProgressBar,statusbarObj]=funGui_JavProgressBar(h);
jProgressBar.setStringPainted( true );
statusbarObj.setText('combining output files...');

output.time=timeSimul.interval';
output.dom =dom;
Nttot      =length(output.time);
Nx         =length(dom.X);
Ny         =length(dom.Y);

EtaAll     =zeros(Nttot,Ny,Nx);
if  model.phiForm==1
phiAll     =EtaAll;
else
uAll       =EtaAll;
vAll       =EtaAll;
end
timeAll    =zeros(Nttot,1);

if options.interior.check==1 
    Nt_IP=round((options.interior.time(2)-options.interior.time(1))/(options.interior.time(3).*timeSimul.dt))-1;
    dtetaAll=zeros(Nt_IP,Ny,Nx);
    dtphiAll=dtetaAll;
    etaIPAll=dtetaAll;
    phiIPAll  =dtetaAll;
    timeIPAll  =zeros(Nt_IP,1);
end
IndtI=1;IndtIntF=0;Indb=0; Nbt=0;Nbx=0; Indcb=0; Ncbt=0;Ncbx=0;
IndtbI=1;IndtbsI=1;

iter=1;
set(jProgressBar,'Maximum',Npartition, 'Value',iter);
datCBreak_temp=[];
datSBreak_temp=[];
if input.bdyassim.option==1
Npartition=ITERAssimSave-1;
end

for i=1:Npartition
    temp=load([Proj.workdir,Proj.savename,'_simul_part_',num2str(i),'.mat']);
    Ntnow=length(temp.output.time);
    IndtF=IndtI+(Ntnow-1);
    EtaAll(IndtI:IndtF,:,:)=temp.output.eta;
    if strcmpi(options.OnlyEta,'No')
        if  model.phiForm==1
            phiAll(IndtI:IndtF,:,:)  =temp.output.phi;
        else
            uAll(IndtI:IndtF,:,:)  =temp.output.u;
            vAll(IndtI:IndtF,:,:)  =temp.output.u;
        end
    end
    IndtI=IndtF;
    
    if ~isempty(temp.output.break_nodes)
        Ntbnow=length(temp.output.break_nodes(:,1));
        Nxbnow=length(temp.output.break_nodes(1,:));
        Indb=Indb+1;
        datBreak_temp(Indb).val=temp.output.break_nodes;
        datBreak_temp(Indb).Nt =Ntbnow;
        datBreak_temp(Indb).Nx =Nxbnow;
        Nbt=Nbt+Ntbnow;
        Nbx=max(Nbx,Nxbnow);
    end

    if ~isempty(temp.output.break_crest)
        Ntbnow=length(temp.output.break_crest);  
       
        for ii=1:Ntbnow
        datCBreak_temp(IndtbI).X=temp.output.break_crest(ii).X;
        datCBreak_temp(IndtbI).Y=temp.output.break_crest(ii).Y;
        datCBreak_temp(IndtbI).t=temp.output.break_crest(ii).t; 
        IndtbI=IndtbI+1;
        end
    end
    
    if ~isempty(temp.output.break_speed)
        Ntbsnow=length(temp.output.break_speed);  
        IndtbsF=IndtbsI+(Ntbsnow-1);
        datSBreak_temp(IndtbsI:IndtbsF,1)=temp.output.break_speed(1:end,1);
        datSBreak_temp(IndtbsI:IndtbsF,2)=temp.output.break_speed(1:end,2);
        datSBreak_temp(IndtbsI:IndtbsF,3)=temp.output.break_speed(1:end,3);
        IndtbsI=IndtbsF;     
    end

    if options.interior.check==1 
        if exist([Proj.workdir,Proj.savename,'_simul_Interior2D_part_',num2str(i),'.mat'],'file')
            tempInter=load([Proj.workdir,Proj.savename,'_simul_Interior2D_part_',num2str(i),'.mat']);

            NtIntnow=length(tempInter.outkinematic.time);
            IndtIntI=IndtIntF+1;
            IndtIntF=IndtIntI+(NtIntnow-1);
            dtetaAll(IndtIntI:IndtIntF,:,:)=tempInter.outkinematic.dteta;
            dtphiAll(IndtIntI:IndtIntF,:,:)=tempInter.outkinematic.dtphi;
            etaIPAll(IndtIntI:IndtIntF,:,:)=tempInter.outkinematic.eta;
            phiIPAll(IndtIntI:IndtIntF,:,:)=tempInter.outkinematic.phi;
            timeIPAll(IndtIntI:IndtIntF,1)=tempInter.outkinematic.time;
            clearvars tempInter;
        end
    end
    clearvars temp;


    set(jProgressBar,'Maximum',Npartition, 'Value',iter);
    statusbarObj.setText(['combining output files...']);
    iter=iter+1;

end
jProgressBar.setVisible(0);

if Indb~=0
    Ndat=Indb;
    dataBreakAll=zeros(Nbt,Nbx);
    Indtbf=0;
    for I=1:Ndat
        valbreak=datBreak_temp(I).val;
        Nbxn    =datBreak_temp(I).Nx;
        Nbtn    =datBreak_temp(I).Nt;
        IndtbI  =1+Indtbf;
        Indtbf  = IndtbI+(Nbtn-1);
        dataBreakAll(IndtbI:Indtbf,1:Nbxn)=valbreak;
    end
else
    dataBreakAll=[];
end


if options.interior.check==1 
    outkinematic.dteta      =dtetaAll;
    outkinematic.dtphi      =dtphiAll;
    outkinematic.eta        =etaIPAll;
    outkinematic.phi        =phiIPAll;
    outkinematic.time       =timeIPAll;
end

if  Idstop==1
    temp_t=output.time(1:IndtF,1);
    output.time=[];
    output.time=temp_t;
    temp_eta=EtaAll(1:IndtF,:,:);
    EtaAll=[];
    EtaAll=temp_eta;
    if strcmpi(options.OnlyEta,'No')
        if  model.phiForm==1
            temp_phi=phiAll(1:IndtF,:,:);
            phiAll=[];
            phiAll=temp_phi;
        else
            temp_u=uAll(1:IndtF,:,:);
            uAll=[];
            uAll=temp_u;
            temp_v=vAll(1:IndtF,:,:);
            vAll=[];
            vAll=temp_v;
            
        end
    end
end

output.eta   =EtaAll;
if strcmpi(options.OnlyEta,'No')
    if  model.phiForm==1
        output.phi   =phiAll;
    else
        output.u   =uAll;
        output.v   =vAll;
    end
end
output.break_nodes =dataBreakAll;
output.break_crest =datCBreak_temp;
output.break_speed =datSBreak_temp;

try
   statusbarObj.setText(['saving data...']);
    save ('-v7.3',[Proj.workdir,Proj.savename,'_simul','.mat'], ...
        'input','par','output','model','bath','dom','ivp','bdyassim','influx','Proj','timeSimul')
%   funC_savefast([Proj.workdir,Proj.savename,'_simul','.mat'], ...
%         'input','par','output','model','bath','input','dom','ivp','bdyassim','influx','Proj','timeSimul')

if options.interior.check==1
     save ('-v7.3',[Proj.workdir,Proj.savename,'_simul_Interior2D.mat'],'Proj','outkinematic','par','dom','Oprt','influx')        
end

Iddelete=1;
    statusbarObj.setText(['']);
catch
    statusbarObj.setText('Failed. Data cannot be combined (too large)');
    Iddelete=0;
end

% Iddelete=0;%never delete files;
if Iddelete==1
    for i=1:Npartition
        delete([Proj.workdir,Proj.savename,'_simul_part_',num2str(i),'.mat']);
        if options.interior.check==1
            if exist([Proj.workdir,Proj.savename,'_simul_Interior2D_part_',num2str(i),'.mat'],'file')
                delete([Proj.workdir,Proj.savename,'_simul_Interior2D_part_',num2str(i),'.mat']);
            end
        end
    end
    statusbarObj.setText('data has been combined');
end


