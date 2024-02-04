function  [timeSig_adj,Insig_adj]=funG_adjust_tsignal_to_tsimul(timesig,Insig,timeSimul,ramp2d)
dt=timesig(2)-timesig(1);
    
timeSig_adj=timesig;Insig_adj=Insig;
Insig=Insig.*ramp2d;
if timeSimul.t_init<timesig(1)
     Nt0=length(timesig);
    Nadd=round((timesig(1)-timeSimul.t_init)/dt)+1;
    timeSig_adj=zeros(1,Nt0+Nadd);
    timeSig_adj(1,Nadd+1:end)=timesig;
    timeSig_adj(1,1:Nadd)=timesig(1)-dt.*(Nadd:-1:1);
    Insig_adj=zeros(size(Insig,1)+Nadd,size(Insig,2));
    Insig_adj(Nadd+1:end,:)=Insig;
    Insig=Insig_adj;
    timesig=timeSig_adj;
end

if timeSimul.t_end>timesig(end)
    Nt0=length(timesig);
    Nadd=round((timeSimul.t_end-timesig(end))/dt)+1;
    timeSig_adj=zeros(1,Nt0+Nadd);
    timeSig_adj(1,1:Nt0)=timesig;
    timeSig_adj(1,Nt0+1:end)=timesig(end)+dt.*(1:Nadd);
    Insig_adj=zeros(size(Insig,1)+Nadd,size(Insig,2));
    Insig_adj(1:Nt0,:)=Insig;
end

end