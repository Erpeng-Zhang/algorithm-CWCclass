function y = fitness(x,D)
%fitness 计算适应度
    y = 0;
    
    % 路线
    s = x;
    t = [x(2:end);x(1)];
    for i = 1 : numel(s)
        y = y + D(s(i),t(i));
    end    
end