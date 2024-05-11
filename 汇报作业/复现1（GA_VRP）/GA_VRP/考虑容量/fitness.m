function y = fitness(x,D,n)
%fitness 计算适应度
    x(x ==0) = n+1;
    y = 0;


    s = [n+1;x];
    t = [x;n+1];
    for i = 1 : numel(s)
        y = y + D(s(i),t(i));
    end

    
end