[XX,YY]=meshgrid(output.X,output.Y);
Xvect=reshape(XX,[],1);
Yvect=reshape(YY,[],1);
indx=funC_closest(output.X,2100);% for runup removing land
EtaN=squeeze(output.eta(end,:,:));
EtaN(:,indx:end)=0;
Etavect=reshape(EtaN,[],1);
Uvect=reshape(squeeze(output.u(1,:,:)),[],1);
Vvect=reshape(squeeze(output.v(1,:,:)),[],1);
dat_initval=[Xvect Yvect Etavect Uvect Vvect];
save dat_initval dat_initval