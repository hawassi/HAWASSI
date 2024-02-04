function Indx=funW_indx_bdyCharInfl(MM,MMInfl)
[GMx,GMy]=gradient(MMInfl);
bb=(abs(GMx)+abs(GMy)).*(1-MMInfl);
MM(MM<1)=0; % make it zero for partial reflecting wall
[GMMx,GMMy]=gradient(MM);
cc=(abs(GMMx)+abs(GMMy)).*(1-MM);

bb(bb~=0 & cc==0)=0;
Indx=find(bb>0);


% figure(110);
% subplot(2,1,1)
% mesh(dom.XX,dom.YY,bb);
% xlabel('x');ylabel('y');
% view(2);
% plot_properties;
% 
% subplot(2,1,2)
% mesh(dom.XX,dom.YY,cc);
% xlabel('x');ylabel('y');
% view(2);
% plot_properties;

end