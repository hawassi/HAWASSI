function w = OmExact(k,d,~) 
w = zeros(length(k),length(d));
for j=1:length(d)
    w(:,j) = sign(k).*sqrt(9.81*k.*tanh(d(j)*k)); 
end
end

% function w = OmExact(k,d,~) 
% Nk=length(k);Nd=length(d);
% d =repmat(d.',Nk,1);
% k =repmat(k,1,Nd);
% w=sign(k).*sqrt(9.81*k.*tanh(d.*k));
% end