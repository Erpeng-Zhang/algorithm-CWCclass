function [ p ] = mutatePop( x, mu )
%mutatePop ����
%   mu:�������

    if rand <= mu
        nVar = numel(x);
        % ������������λ��
        index = randperm(nVar,2);
        v = x(index(1));
        x(index(1)) = x(index(2));
        x(index(2)) = v;
        
    end
    p=x;

end

