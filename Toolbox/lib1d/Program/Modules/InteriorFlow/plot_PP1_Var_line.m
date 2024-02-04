function plot_PP1_Var_line(X,Var,label_y,label_x,SetAx,~)

if SetAx.Xlim.Id==1
indx1=closest(X,SetAx.Xlim.val(1));
indx2=closest(X,SetAx.Xlim.val(2));    
plot(X(indx1:indx2),Var(indx1:indx2));
xlim([X(indx1) X(indx2)]);
else
plot(X,Var);
xlim([min(X(1),X(end)) max(X(1),X(end))]);
end

if SetAx.Ylim.Id==1
indz1=closest(Var,SetAx.Ylim.val(1));
indz2=closest(Var,SetAx.Ylim.val(2));    
ylim([Var(indz1) Var(indz2)]);
end

xlabel(label_x);
ylabel(label_y);
