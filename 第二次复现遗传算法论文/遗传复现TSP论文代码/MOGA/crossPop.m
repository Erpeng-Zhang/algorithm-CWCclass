function [C1,C2,C3,C4 ] = crossPop( x1, x2 )
    %crossPop 两点交叉
    %   两个父代生成4个子代
    
        nVar = numel(x1);
        
        A = x1;
        B = x2;
        %随机生成交叉点 不能是第一个和最后一个点
        dots = 1+randperm(nVar-1,2);
        dots = sort(dots);
    
        
        %% 生成前两个子代
        % 复制染色体
        C1 = zeros(nVar,1);C2 = zeros(nVar,1);
        %交叉与换位
        C1(dots(1):dots(2)) = B(dots(1):dots(2));
        C2(dots(1):dots(2)) = A(dots(1):dots(2));
        A = [A(dots(2)+1:end);A(1:dots(1)-1);A(dots(1):dots(2))];
        B = [B(dots(2)+1:end);B(1:dots(1)-1);B(dots(1):dots(2))];
        %去重
        [~, ia] = setdiff(A, C1);
        A = A(sort(ia));
        %将未赋值基因数据补全
        C1(dots(2)+1:end) = A(1:nVar-dots(2));
        C1(1:dots(1)-1) = A(nVar-dots(2)+1:end);  
    
        %去重
        [~, ia] = setdiff(B, C2);
        B = B(sort(ia));
        %将未赋值基因数据补全
        C2(dots(2)+1:end) = B(1:nVar-dots(2));
        C2(1:dots(1)-1) = B(nVar-dots(2)+1:end);
    
        %% 生成后两个子代
        % 复制染色体
        C3 = zeros(nVar,1);C4 = zeros(nVar,1);
        AA = x1;
        BB = x2;
        C3(dots(2)+1:end) = BB(dots(2)+1:end);
        C4(dots(2)+1:end) = AA(dots(2)+1:end);
        % 交换前两部分基因位置
        AA = [AA(dots(1):dots(2));AA(1:dots(1)-1);AA(dots(2)+1:end)];
        BB = [BB(dots(1):dots(2));BB(1:dots(1)-1);BB(dots(2)+1:end)];
    
        %找出AA中BB的第三部分元素删除
        [~, ia] = setdiff(AA, BB(dots(2)+1:end));
    
        C3(1:dots(2)) = AA(sort(ia));
        %找出BB中AA的第三部分元素删除
        [~, ia] = setdiff(BB, AA(dots(2)+1:end));
        C4(1:dots(2)) = BB(sort(ia));
    
    
       
    
    
    end
    
    