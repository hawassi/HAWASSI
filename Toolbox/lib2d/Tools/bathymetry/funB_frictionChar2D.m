function frictionChar=funB_frictionChar2D(dom,friction)
fricdata=friction.param;
Ni=length(fricdata(:,1));
frictionChar=0;g=9.81;
for i=1:Ni
    if strcmp(fricdata(i,1),'Rectangle')
         bb=cell2mat(fricdata(i,4))-dom.Y(1);tt=dom.Y(end)-cell2mat(fricdata(i,5));
         ll=cell2mat(fricdata(i,2))-dom.X(1);rr=dom.X(end)-cell2mat(fricdata(i,3));
         if bb==0 &&tt==0 && ll==0 && rr==0
          chi=ones(size(dom.XX)).*cell2mat(fricdata(i,6)).^2.*g;     
         else
         chi=funC_HeavCuboid2D(dom.X,dom.Y,ll,rr,bb,tt).*cell2mat(fricdata(i,6)).^2.*g;
         end
    else
        fricuser=friction.userdata(i).bdry;
        chi=funB_bottom_friction_userdefined(dom,fricuser).*cell2mat(fricdata(i,6)).^2.*g;  
    end
        frictionChar=frictionChar+chi;
end
% figure;
% surf(dom.XX,dom.YY,frictionChar,'edgecolor','none');