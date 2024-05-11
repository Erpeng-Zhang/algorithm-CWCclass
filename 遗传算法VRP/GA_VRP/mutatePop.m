function [ p ] = mutatePop( x, mu )
%mutatePop 变异
%   mu:变异概率

    if rand <= mu
        nVar = numel(x);
        % 交换两个基因位置
        index = randperm(nVar,2);
        v = x(index(1));
        x(index(1)) = x(index(2));
        x(index(2)) = v;
        
    end
    p=x;

end

