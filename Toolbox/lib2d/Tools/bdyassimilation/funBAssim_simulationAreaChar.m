function [cfSAbdyassim]=funBAssim_simulationAreaChar(dom,fbl,lambda_p)
dampdata=fbl.param;

if strcmpi(fbl.option,'Yes')
    fbl.l=cell2mat(dampdata(1,2));
    fbl.r=cell2mat(dampdata(1,3));
    fbl.b=cell2mat(dampdata(1,4));
    fbl.t=cell2mat(dampdata(1,5));
    if fbl.l==0
        fbl.l=lambda_p/3;
    end
    if fbl.r==0
        fbl.r=lambda_p/3;
    end
    if fbl.b==0
        fbl.b=lambda_p/3;
    end
    if fbl.t==0
        fbl.t=lambda_p/3;
    end
else
    fbl.l=lambda_p/3;
    fbl.r=lambda_p/3;
    fbl.b=lambda_p/3;
    fbl.t=lambda_p/3;
end
cfSAbdyassim=funC_cfSA2d(dom.X,dom.Y,fbl);

% ff=figure;
% set(ff,'Renderer','zbuffer')
% surf(dom.XX,dom.YY,cfSAbdyassim,'edgecolor','none');

end
