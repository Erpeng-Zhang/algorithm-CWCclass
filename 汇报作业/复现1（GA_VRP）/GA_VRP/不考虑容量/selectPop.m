function [ p ] = selectPop( Parent )
%UNTITLED3 锦标赛选择法
%   此处显示详细说明

    n = numel(Parent);
    index = randperm(n);
    p1 = Parent(index(1));
    p2 = Parent(index(2));
    
    if p1.y <= p2.y
        p = p1;
    else
        p = p2;
    end
        

end

