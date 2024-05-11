function [ p ] = mutatePop( x, mu ,m,n,demands,volume)
%mutatePop 变异
%   mu:变异概率
    condition = true;
    while condition
        if rand <= mu
            nVar = numel(x);
            % 交换两个基因位置
            index = randperm(nVar,2);
            v = x(index(1));
            x(index(1)) = x(index(2));
            x(index(2)) = v;
        
        end
        p=x;

        condition = ~isOK(p,m,n,demands,volume);
    end
    

end

