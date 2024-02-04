% invOm, version March 2012
% % % % Software developed by 
% % % % LabMath-Indonesia & University of Twente
function k = invOm(v,d,Om,omAdd)
M=length(v);
J=length(d);
k=zeros(M,J);
for j=1:J
    k(1,j)= fzero(@(k) v(1)-Om(k,d(j),omAdd),v(1)/sqrt(9.81*d(j)));
    for m=2:M
        k(m,j) = fzero(@(k) v(m)-Om(k,d(j),omAdd),k(m-1,j));%v(m)/sqrt(9.81*d(j)),optimset('Display','off'));
    end
end