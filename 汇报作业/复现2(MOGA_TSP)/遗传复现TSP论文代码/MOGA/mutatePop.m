function [ p ] = mutatePop( x, mu )
    %mutatePop 变异
    %   mu:变异概率
    
        if rand <= mu
            nVar = numel(x);
            % 随机生成变异片段
            dots = 1+randperm(nVar-1,2);
            dots = sort(dots);
            
            value = x(dots(1):dots(2));
            value = flip(value);
            x(dots(1):dots(2)) = value;
        end
        p=x;
    end
    
    

