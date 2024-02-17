%%% 9 Feb 2012 infinitely smooth char function (not inf at x=w)

function f=SmChar(x,x0,w)
hw=w/2;
f= Heaviside(x-x0).*min(((x-x0).*(1+tanh(x-x0)))/w/(1+tanh(x0+w)),Heaviside(x0+w-x)+1);
% exp(-(w./(100.*(x-x0))).^2+1/100^2)