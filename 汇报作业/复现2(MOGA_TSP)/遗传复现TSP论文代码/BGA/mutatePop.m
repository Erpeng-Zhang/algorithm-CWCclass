function [ p ] = mutatePop( x, mu )
%mutatePop ����
%   mu:�������

    if rand <= mu
        nVar = numel(x);
        % ������ɱ���Ƭ��
        dots = 1+randperm(nVar-1,2);
        dots = sort(dots);
        
        value = x(dots(1):dots(2));
        value = flip(value);
        x(dots(1):dots(2)) = value;
    end
    p=x;
end

