function [result,time,It] = MOGA(data,opt_value)

% 开始计时
tic;

%城市数
n = size(data,1);

%位置坐标
X = (data(:,2:3));

D=Distanse(X);  %生成距离矩阵

nVar = n;%染色体数

nPop = 100;  %种群规模
maxIt = 1000;    %最大迭代次数
nPc = 1;      %子代比例
nC = round(nPop* nPc/4)*4; %子代种群规模
mu = 0.3;  %变异概率

template.x = [];% 位置
template.y = [];% 适应度函数值

%生成初始种群
Parent = repmat(template, nPop, 1);

%适应度变化曲线数据
fig_data = [];

%初始化种群
for i = 1:nPop
    % 生成染色体
    Parent(i).x = randperm(nVar)';
    %计算适应度
    Parent(i).y = fitness(Parent(i).x,D);    
end

% 按适应度升序排序
[~, so] = sort([Parent.y], 'ascend');
Parent = Parent(so);
opt_old = Parent(1).y;%旧值
shoulian_time = 0;

% 主循环
for It = 1:maxIt
    %生成子代种群
    Offspring =repmat(template, nC/2, 4);
    
    selected_id = [];%被选中精英ID
    for j = 1: nC / 2 
        % 选择
        [id,p1] = selectPop(Parent,It); 
        selected_id = [selected_id,id];
        
        [id,p2] = selectPop(Parent,It); 
        selected_id = [selected_id,id];
        
        % 交叉
        [Offspring(j,1).x, Offspring(j,2).x,Offspring(j,3).x, Offspring(j,4).x] = crossPop(p1.x, p2.x);        
    end

    %选择q个精英
    selected_id = unique(selected_id);%去重
    elite = Parent(selected_id);

    % 两列变一列
    Offspring = Offspring(:);
    
    for k = 1 : numel(Offspring)    
        %变异
        Offspring(k).x = mutatePop( Offspring(k).x, mu);   
        %计算
        Offspring(k).y = fitness(Offspring(k).x,D);        
    end
    
    % 新种群=精英+子代
    newPop = [elite; Offspring];
    % 按适应度进行排序
    [~, so] = sort([newPop.y], 'ascend');
    newPop = newPop(so);
    %末位淘汰 保持种群规模
    Parent = newPop(1: nPop);

    % 中止判定
    %连续20代未改变跳出
    opt = Parent(1).y;
    if opt == opt_old
        shoulian_time = shoulian_time+1;
    else
        shoulian_time = 1;
    end
    
    if shoulian_time >= 20
        break;
    end
    %达到最优值退出循环
    if round(opt,4)<=opt_value
        break;
    end
      
end

time = toc;%运行时间

result = Parent(1).y;
end

%% 计算两两城市之间的距离（x1-x2）^2+(y1-y2)^2
%输入 a  各城市的位置坐标
%输出 D  两两城市之间的距离
function D=Distanse(a)
    row=size(a,1);
    D=zeros(row);
    for i=1:row
        for j=i+1:row
            D(i,j)=((a(i,1)-a(j,1))^2+(a(i,2)-a(j,2))^2)^0.5;
            D(j,i)=D(i,j);
        end
    end
end


function [ id,p ] = selectPop( Parent, It )
%selectPop 双重轮盘赌选择法
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

function [ p ] = mutatePop( x, mu )
%mutatePop 变异
%   mu:变异概率
    if rand <= mu
        nVar = numel(x);
        % 随机生成变异片段
        dots = 1+randperm(nVar-1,2);
        dots = sort(dots);
        
        value = x(dots(1):dots(2));
        value = flip(value);
        x(dots(1):dots(2)) = value;
    end
    p=x;
end

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






