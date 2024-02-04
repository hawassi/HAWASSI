function Flag=funW_bdy_walldir_flag(dom)
Wall=dom.wall;
[GMx,GMy]=gradient(Wall.charInfl);
[GMxx,GMxy]=gradient(GMx);
[GMyx,GMyy]=gradient(GMy);
GMxbdyInfl=GMx(Wall.bdyInfl.index);GMybdyInfl=GMy(Wall.bdyInfl.index);
GMxxbdyInfl=GMxx(Wall.bdyInfl.index);GMxybdyInfl=GMxy(Wall.bdyInfl.index);
GMyxbdyInfl=GMyx(Wall.bdyInfl.index);GMyybdyInfl=GMyy(Wall.bdyInfl.index);
% xc=mean(dom.XX(Wall.bdyInfl.index));
% yc=mean(dom.YY(Wall.bdyInfl.index));
if dom.Nx>=dom.Ny
Const=1;
else
Const=2;    
end
Flag=Const*ones(size(Wall.bdyInfl.index));

% % Flag(GMxbdyInfl>GMybdyInfl)=2;Flag(GMxbdyInfl<GMybdyInfl)=1;
% Flag(GMybdyInfl==0)=1;Flag(GMxbdyInfl==0)=2;
% %Flag(GMxybdyInfl~=0)=2;%Flag(GMyxbdyInfl~=0)=2;
%Flag(abs(GMxbdyInfl)+abs(GMybdyInfl)==1)=Const;
% 

Flag(GMybdyInfl==0)=1;%*sign(GMxbdyInfl(GMybdyInfl==0));
Flag(GMxbdyInfl==0)=2;%*sign(GMybdyInfl(GMxbdyInfl==0));
Flag(GMxybdyInfl~=0)=Const*sign(GMxbdyInfl(GMxybdyInfl~=0));
Flag(Flag==0)=Const*sign(GMyybdyInfl(Flag==0));

indFlag= Flag==2 ;
indSign= sign(GMybdyInfl)==-1;
Flag(indFlag.*indSign==1)=-2;

indFlag= Flag==1 ;
indSign= sign(GMxbdyInfl)==-1;
Flag(indFlag.*indSign==1)=-1;


% NN=1:1:length(Wall.bdyInfl.index);
% WallC=zeros(size(dom.XX));
% WallC(Wall.bdyInfl.index)=1;
% figure(33);
% subplot(3,1,1)
% plot(NN,dom.XX(Wall.bdyInfl.index),'r',NN,dom.YY(Wall.bdyInfl.index),'--b');
% subplot(3,1,2)
% plot(NN,GMxbdyInfl,'r',NN,GMybdyInfl,'--b',NN,Flag,'k');
% subplot(3,1,3)
% surf(dom.XX,dom.YY,WallC,'edgecolor','none');view(2);
% xlim([min(min(dom.XX)) max(max(dom.XX))]);
% ylim([min(min(dom.YY)) max(max(dom.YY))]);
% 
% figure(34)
% subplot(3,1,1)
% surf(dom.XX,dom.YY,GMx,'edgecolor','none');view(2);
% xlim([min(min(dom.XX)) max(max(dom.XX))]);axis equal;
% ylim([min(min(dom.YY)) max(max(dom.YY))]);
% subplot(3,1,2)
% surf(dom.XX,dom.YY,GMy,'edgecolor','none');view(2);
% xlim([min(min(dom.XX)) max(max(dom.XX))]);axis equal;
% ylim([min(min(dom.YY)) max(max(dom.YY))]);
% subplot(3,1,3)
% surf(dom.XX,dom.YY,abs(GMx)+abs(GMy),'edgecolor','none');view(2);
% xlim([min(min(dom.XX)) max(max(dom.XX))]);axis equal;
% ylim([min(min(dom.YY)) max(max(dom.YY))]);
% figure(35)
% subplot(2,1,1)
% surf(dom.XX,dom.YY,GMxx,'edgecolor','none');view(2);
% xlim([min(min(dom.XX)) max(max(dom.XX))]);
% ylim([min(min(dom.YY)) max(max(dom.YY))]);
% subplot(2,1,2)
% surf(dom.XX,dom.YY,GMxy,'edgecolor','none');view(2);
% xlim([min(min(dom.XX)) max(max(dom.XX))]);
% ylim([min(min(dom.YY)) max(max(dom.YY))]);
% figure(36)
% subplot(2,1,1)
% surf(dom.XX,dom.YY,GMyx,'edgecolor','none');view(2);
% xlim([min(min(dom.XX)) max(max(dom.XX))]);
% ylim([min(min(dom.YY)) max(max(dom.YY))]);
% subplot(2,1,2)
% surf(dom.XX,dom.YY,GMyy,'edgecolor','none');view(2);
% xlim([min(min(dom.XX)) max(max(dom.XX))]);
% ylim([min(min(dom.YY)) max(max(dom.YY))]);
