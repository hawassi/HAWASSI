function om2k=invom_interp(k,ww,depthinflux)
Om_k       = OmExact(k,depthinflux);
if max(ww)<=max(Om_k)
    om2k            = interp1(OmExact(k,depthinflux),k,ww,'spline');
else
    om2k        =zeros(length(ww),1);
    Ind_half    = floor(length(k)/2);
    k_half1      =k(1:Ind_half);
    k_half2      =k(Ind_half+1:end);
    
    Omega1       =OmExact(k_half1,depthinflux);
    [~,Idnu1]    =max(Omega1);
    
    Omega2       =OmExact(k_half2,depthinflux);
    [~,Idnu2]    =min(Omega2);
    
    Idwwmax      =closest(ww,max(Omega1));
    om2k(1:Idwwmax)=interp1(Omega1(1:Idnu1),k_half1(1:Idnu1),ww(1:Idwwmax));% invOmExact(2*pi/nupeak,depth);
    Idwwmin      =closest(ww,min(Omega2));
    om2k(Idwwmin:end)= interp1(Omega2(Idnu2:end),k_half2(Idnu2:end),ww(Idwwmin:end));% invOmExact(2*pi/nupeak,depth);
    om2k(isnan(om2k))=0;
    %         figure
    %         subplot(3,1,1)
    %         plot(k,Om(k,depthinflux,omAdd));
    %         subplot(3,1,2)
    %         plot(ww)
    %         subplot(3,1,3)
    %         plot(ww,om2k)
end
end