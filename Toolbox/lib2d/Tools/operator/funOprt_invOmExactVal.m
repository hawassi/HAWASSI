% invOm, version March 2012
% % % % Software developed by 
% % % % LabMath-Indonesia & University of Twente
function k = funOprt_invOmExactVal(v,d,g)
M=length(v);
J=length(d);
k=zeros(M,J);
for j=1:J
    k(1,j)= fzero(@(k) v(1)-OmExact(k,d(j)),v(1)/sqrt(g*d(j)));
    for m=2:M
        k(m,j) = fzero(@(k) v(m)-OmExact(k,d(j)),v(m)/sqrt(9.81*d(j)),optimset('Display','off'));
    end
end