function path = decode(x)
% 染色体解码为方案

x = [0;x;0];

ind = find(x==0);
m = numel(ind)-1;

for i = 1 : m
    if ind(i+1)-ind(i)==1
        path{i} = [0;0];
    else
        path{i} = [0;x(ind(i)+1:ind(i+1)-1);0];
    end
end


end