function idx=funC_findIndexinColumn(input,arg)
nr  = size(input,1);
idx = zeros(nr,1);
for j=1:nr
    try
    idx(j) = find(input(j,:) == 0, 1,arg);  
    catch
       continue; 
    end
end

