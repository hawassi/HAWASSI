statusbarObj.setText('combining output files...');
global Idstop

output.time=influx.gen.timesig';
output.x   =x;
Nttot      =length(output.time);
Nx         =length(x);

EtaAll     =zeros(Nttot,Nx);
uAll       =EtaAll;
timeAll    =zeros(Nttot,1);

if shippar.check== 1
    
    RBAll= zeros(Nttot,6*shippar.Nship);
    ExtAll=zeros(Nttot,6*shippar.Nship);
end

if options.interior.check==1
    Nt_IP=round((options.interior.time(2)-options.interior.time(1))/(options.interior.time(3).*par.dt))-1;
    dtetaAll=zeros(Nt_IP,length(x));
    dtphihatAll=dtetaAll;
    etaIPAll=dtetaAll;
    uIPAll  =dtetaAll;
    timeIP  =zeros(Nt_IP,1);
end
IndtI=1;IndtIntF=0;Indb=0; Nbt=0;Nbx=0; Indcb=0; Ncbt=0;Ncbx=0;


JavProgressBar;
jProgressBar.setStringPainted( true );
iter=1;
set(jProgressBar,'Maximum',Npartition, 'Value',iter);


for i=1:Npartition
    temp=load([Proj.Dir,Proj.savename,'_simul_part_',num2str(i),'.mat']);
    Ntnow=length(temp.output.time);
    IndtF=IndtI+(Ntnow-1);
    EtaAll(IndtI:IndtF,:)=temp.output.eta;
    uAll(IndtI:IndtF,:)  =temp.output.u;
    if shippar.check== 1
        RBAll(IndtI:IndtF,:)=temp.output.RB;
        ExtAll(IndtI:IndtF,:)=temp.output.Ext;
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
        Ntbnow=length(temp.output.break_crest(:,1));
        Nxbnow=length(temp.output.break_crest(1,:));
        Indcb=Indcb+1;
        datCBreak_temp(Indcb).val=temp.output.break_crest;
        datCBreak_temp(Indcb).Nt =Ntbnow;
        datCBreak_temp(Indcb).Nx =Nxbnow;
        Ncbt=Ncbt+Ntbnow;
        Ncbx=max(Ncbx,Nxbnow);
    end
    
    if options.interior.check==1
        if exist([Proj.Dir,Proj.savename,'_simul_Interior_part_',num2str(i),'.mat'],'file')
            tempInter=load([Proj.Dir,Proj.savename,'_simul_Interior_part_',num2str(i),'.mat']);
            
            NtIntnow=length(tempInter.InputInteriorCalc.time);
            IndtIntI=IndtIntF+1;
            IndtIntF=IndtIntI+(NtIntnow-1);
            dtetaAll(IndtIntI:IndtIntF,:)=tempInter.InputInteriorCalc.dteta;
            dtphihatAll(IndtIntI:IndtIntF,:)=tempInter.InputInteriorCalc.dtphihat;
            etaIPAll(IndtIntI:IndtIntF,:)=tempInter.InputInteriorCalc.eta;
            uIPAll(IndtIntI:IndtIntF,:)=tempInter.InputInteriorCalc.u;
            timeIP(IndtIntI:IndtIntF,1)=tempInter.InputInteriorCalc.time;
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

if Indcb~=0
    Ndat=Indcb;
    dataCBreakAll=zeros(Ncbt,Ncbx);
    Indtcbf=0;
    for I=1:Ndat
        valbreak=datCBreak_temp(I).val;
        Nbxn    =datCBreak_temp(I).Nx;
        Nbtn    =datCBreak_temp(I).Nt;
        IndtcbI  =1+Indtcbf;
        Indtcbf  = IndtcbI+(Nbtn-1);
        dataCBreakAll(IndtcbI:Indtcbf,1:Nbxn)=valbreak;
    end
else
    dataCBreakAll=[];
end


if options.interior.check==1
    dteta=dtetaAll;
    dtphihat=dtphihatAll;
end

if  Idstop==1
    temp_t=output.time(1:IndtF,1);
    output.time=[];
    output.time=temp_t;
    temp_eta=EtaAll(1:IndtF,:);
    EtaAll=[];
    EtaAll=temp_eta;
    temp_u=uAll(1:IndtF,:);
    uAll=[];
    uAll=temp_u;
    
     if shippar.check== 1
       temp_RBAll=RBAll(1:IndtF,:);  
       RBAll=temp_RBAll;
       temp_ExtAll=ExtAll(1:IndtF,:);  
       ExtAll=temp_ExtAll;
    end
end

output.eta =EtaAll;
output.u   =uAll;
output.break_nodes =dataBreakAll;
output.break_crest =dataCBreakAll;
if shippar.check== 1
output.RB=RBAll;
output.Ext=ExtAll;
end

try
    statusbarObj.setText(['saving data...']);
    save ('-v7.3',[Proj.Dir,Proj.savename,'_simul','.mat'], ...
        'output','model','bath','par','influx','Proj','shippar','Oprt')
    InputdataPostProcInternalFlowsCombines;
    
    if IndtIntF>0
        save ('-v7.3',[Proj.Dir,Proj.savename,'_simul_Interior.mat'],'InputInteriorCalc')
        checkInterior=1;
    end
    Iddelete=1;
    statusbarObj.setText(['']);
catch
    statusbarObj.setText('Combining file is fail. Data is too large');
    Iddelete=0;
end

if Iddelete==1
    for i=1:Npartition
        delete([Proj.Dir,Proj.savename,'_simul_part_',num2str(i),'.mat']);
        if options.interior.check==1
            if exist([Proj.Dir,Proj.savename,'_simul_Interior_part_',num2str(i),'.mat'],'file')
                delete([Proj.Dir,Proj.savename,'_simul_Interior_part_',num2str(i),'.mat']);
            end
        end
    end
    statusbarObj.setText('data has been combined');
end


