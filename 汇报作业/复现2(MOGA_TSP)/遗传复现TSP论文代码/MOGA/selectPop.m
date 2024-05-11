function [ id,p ] = selectPop( Parent, It )
    %UNTITLED3 双重轮盘赌选择法
    %   此处显示详细说明
        p = Parent(1);
    
        n = numel(Parent);
    
        %生成轮盘
        wheel = zeros(n,1);
        if It/2 ~= 0
            % 迭代次数奇数使用F1
            beta = rand * 0.3;
            for i = 1 : n
                wheel(i) = beta*(1-beta)^(i-1);
            end
            wheel(:) = wheel(:) / sum(wheel(:));
        else
            % 迭代次数偶数使用F2
            for i = 1 : n
                wheel(i) = (n-i+1)/n;
            end
            wheel(:) = wheel(:) / sum(wheel(:));
        end
    
    
    
        %随机生成一个值
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
    
    