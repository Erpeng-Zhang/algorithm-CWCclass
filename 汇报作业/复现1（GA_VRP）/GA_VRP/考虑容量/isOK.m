function output = isOK(x,m,n,demands,volume)
%判断解是否可行
%   


output = false;



path = decode(x);

for i = 1 : m
    p = path{i};
    p(p==0) = []; 
    if numel(p)==0 || sum(demands(p)) > volume 
        return;
    end
end

%判断是否经过所有客户
if sum(x) == sum(1:1:n)
    output = true;
else
    return;
end



end