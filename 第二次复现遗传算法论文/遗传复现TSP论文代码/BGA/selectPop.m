function [ id,p ] = selectPop( Parent, It )
%UNTITLED3 ˫�����̶�ѡ��
%   �˴���ʾ��ϸ˵��
    p = Parent(1);

    n = numel(Parent);

    %��������
    wheel = zeros(n,1);
    if It/2 ~= 0
        % ������������ʹ��F1
        beta = rand * 0.3;
        for i = 1 : n
            wheel(i) = beta*(1-beta)^(i-1);
        end
        wheel(:) = wheel(:) / sum(wheel(:));
    else
        % ��������ż��ʹ��F2
        for i = 1 : n
            wheel(i) = (n-i+1)/n;
        end
        wheel(:) = wheel(:) / sum(wheel(:));
    end



    %�������һ��ֵ
    index = rand;
    for i = 2 : n
        wheel(i) = wheel(i)+ wheel(i-1);
    end

    id = find(wheel<=index, 1, 'last' );
      
    if id ~= 0
        id = id + 1;
        p = Parent(id);
    end

        

end

