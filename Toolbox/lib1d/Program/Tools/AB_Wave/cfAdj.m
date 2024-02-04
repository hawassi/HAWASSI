%Characteristic function for Nonlinear influxing adjustment 
%ChiAdj.R : cf for propagation to the right
%ChiAdj.L : cf for propagation to the left
%ChiAdj.RL : cf for bidirection propagation
function ChiAdj=cfAdj(x,Xinflux,lambda_p,adjcoef,Indirection)
if strcmp(Indirection,'Uni+')
ChiAdj    = cf(x,Xinflux,adjcoef*lambda_p)';
elseif strcmp(Indirection,'Uni-')
ChiAdj    = 1-cf(x,Xinflux-adjcoef*lambda_p,adjcoef*lambda_p)';
elseif strcmp(Indirection,'Bi')
ChiAdjL    = cf(x,Xinflux,adjcoef*lambda_p)';
ChiAdjR    = 1-cf(x,Xinflux-adjcoef*lambda_p,adjcoef*lambda_p)';
ChiAdj    = ChiAdjR+ChiAdjL;
end     