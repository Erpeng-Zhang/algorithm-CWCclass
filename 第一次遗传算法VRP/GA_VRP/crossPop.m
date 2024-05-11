function [ y1, y2 ] = crossPop( x1, x2 )
%crossPop ���㽻��
%   �˴���ʾ��ϸ˵��

    nVar = numel(x1);
    
    %������ɽ����
    index = randi([1 nVar]);
    
    % ����Ⱦɫ��
    v1 = x1;v2=x2;

    %��������
    v1(1:index) = x2(1:index);
    v2(1:index) = x1(1:index);

    %��һ������ʹȾɫ�����Ľ����
    ind_v1 = [];
    ind_v2 = [];
    for i = 1:index
        % ��¼v1�ظ�Ԫ��λ��
        if v1(i) ~= 0
            ind_v1 = [ind_v1,find(v1(index+1:end)==v1(i))];
        end

        % ��¼v2�ظ�Ԫ��
        if v2(i) ~= 0
            ind_v2 = [ind_v2,find(v2(index+1:end)==v2(i))];
        end
    end

    j=1;
    while numel(ind_v1) < numel(ind_v2)
        
        i = find(v1(index+1:end)==0);
        ind_v1 = [ind_v1,i(j)];
        j = j+1;
    end
    
    j=1;
    while numel(ind_v1) > numel(ind_v2)        
        i = find(v2(index+1:end)==0);
        ind_v2 = [ind_v2,i(j)];
        j = j+1;
    end

    ind_v1 = index + ind_v1; ind_v2 = index + ind_v2;

    y1 = v1;
    y2 = v2;

    % ���н���
    for i = 1 : numel(ind_v1)
        y1(ind_v1(i)) = v2(ind_v2(i));
        y2(ind_v2(i)) = v1(ind_v1(i));
    end



end

