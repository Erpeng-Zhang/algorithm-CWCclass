function [C1,C2] = crossPop( x1, x2 )
%crossPop ���㽻��
%   ������������4���Ӵ�

    nVar = numel(x1);
    
    A = x1;
    B = x2;
    %������ɽ���� �����ǵ�һ�������һ����
    dots = 1+randperm(nVar-1,2);
    dots = sort(dots);

    
    %% ����ǰ�����Ӵ�
    % ����Ⱦɫ��
    C1 = zeros(nVar,1);C2 = zeros(nVar,1);
    %�����뻻λ
    C1(dots(1):dots(2)) = B(dots(1):dots(2));
    C2(dots(1):dots(2)) = A(dots(1):dots(2));
    A = [A(dots(2)+1:end);A(1:dots(1)-1);A(dots(1):dots(2))];
    B = [B(dots(2)+1:end);B(1:dots(1)-1);B(dots(1):dots(2))];
    %ȥ��
    [~, ia] = setdiff(A, C1);
    A = A(sort(ia));
    %��δ��ֵ�������ݲ�ȫ
    C1(dots(2)+1:end) = A(1:nVar-dots(2));
    C1(1:dots(1)-1) = A(nVar-dots(2)+1:end);  

    %ȥ��
    [~, ia] = setdiff(B, C2);
    B = B(sort(ia));
    %��δ��ֵ�������ݲ�ȫ
    C2(dots(2)+1:end) = B(1:nVar-dots(2));
    C2(1:dots(1)-1) = B(nVar-dots(2)+1:end);

     


end

