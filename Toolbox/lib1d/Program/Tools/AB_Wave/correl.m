function corr=correl(s,m) %s=simulation, m = measurement
slen    = length(s);
mlen    = length(m);
coarse  = round(mlen/slen);%nearest(mlen/slen);
inp=0;
    for j=1:coarse*floor(slen/coarse)
        inp=inp+s(j)*m(1+coarse*(j-1));
    end;
    ss=s'*s;
    mm= m'*m/coarse;
corr=inp/sqrt(ss*mm);
% corr2=inp*dt/sqrt(Emeas*Esimul)/sqrt(coarse)

