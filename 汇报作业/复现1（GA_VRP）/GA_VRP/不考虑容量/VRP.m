clc;clear;

%% TSP拿来的数据
% m = 5; %汽车数
% n = 30; %城市数
% 
% %导入城市坐标位置
% load('CityPosition3.mat');

% 所罗门数据集的一个数据
m = 5;%汽车数

%CUST NO.   XCOORD.   YCOORD.   DEMAND    READY TIME   DUE DATE   SERVICE TIME
data = [1      45         68         10        885        994         90   
    2      45         70         30        802        893         90   
    3      42         66         10         25        186         90   
    4      42         68         10        699        810         90   
    5      42         65         10         15        120         90   
    6      40         69         20        580        743         90   
    7      40         66         20        142        253         90   
    8      38         68         20        220        359         90   
    9      38         70         10        499        640         90   
   10      35         66         10        331        436         90   
   11      35         69         10        420        533         90   
   12      25         85         20        617        756         90   
   13      22         75         30         30        155         90   
   14      22         85         10        541        646         90   
   15      20         80         40        362        451         90   
   16      20         85         40        448        555         90   
   17      18         75         20         75        172         90   
   18      15         75         20        142        291         90   
   19      15         80         10        244        379         90   
   20      30         50         10         10        137         90   
   21      30         52         20        888        991         90   
   22      28         52         20        776        919         90   
   23      28         55         10        709        800         90   
   24      25         50         10         25        184         90   
   25      25         52         40        142        251         90   
   0      40         50          0          0       1236          0];

%客户数
n = size(data,1)-1;

%位置坐标
X = (data(:,2:3));


D=Distanse(X);  %生成距离矩阵

% 在二维图上画出所有坐标点
figure(1)
plot(X(:,1), X(:,2),'red*');hold on; % 保持图形窗口开启
% 在点旁边标上序号
for i = 1:numel(X(:,1))-1
    text(X(i,1), X(i,2), num2str(i));
end
text(X(end,1), X(end,2), '仓库');
hold on;



nVar = m+n-1;%染色体数

nPop = 30;  %种群规模
maxIt = 1000;    %最大迭代次数
nPc = 0.8;      %子代比例
nC = round(nPop* nPc/2)*2; %子代种群规模
mu = 0.1;  %变异概率

template.x = [];
template.y = [];

%生成初始种群
Parent = repmat(template, nPop, 1);


%适应度变化曲线数据
fig_data = [];

%初始化种群
for i = 1:nPop
    % 生成染色体
    Parent(i).x = ini_x(nVar,m,n);
    %计算适应度
    Parent(i).y = fitness(Parent(i).x,D,n);
    
end

% 主循环
for It = 1:maxIt
    %生成子代种群
    Offspring =repmat(template, nC/2, 2);
    
    for j = 1: nC / 2 
        % 选择
        p1 = selectPop(Parent);
        p2 = selectPop(Parent);
        % 交叉
        [Offspring(j,1).x, Offspring(j,2).x] = crossPop(p1.x, p2.x);        
    end
    
    % 两列变一列
    Offspring = Offspring(:);
    
    for k = 1 : nC     
        %变异
        Offspring(k).x = mutatePop( Offspring(k).x, mu);   
        %计算
        Offspring(k).y = fitness(Offspring(k).x,D,n);        
    end
    
    % 新种群=父代+子代
    newPop = [Parent; Offspring];
    % 按适应度进行排序
    [~, so] = sort([newPop.y], 'ascend');
    newPop = newPop(so);
    %末位淘汰 保持种群规模
    Parent = newPop(1: nPop);
    
    disp(['迭代次数:', num2str(It), ',最小值为：', num2str(Parent(1).y)])
    
    fig_data(It,1) = It;fig_data(It,2) = Parent(1).y;
    
end


%% 最优解的路线图
x = Parent(1).x;x = [0;x;0];
ind = find(x==0);lable = [];
for i = 1 : m
    if ind(i+1)-ind(i)==1
        path{i} = [n+1;n+1];
    else
        path{i} = [n+1;x(ind(i)+1:ind(i+1)-1);n+1];
    end
    plot(X(path{i},1),X(path{i},2),'-',LineWidth=i*.5);hold on;
end
title('路径图');
hold off;




% 画出适应度迭代曲线
figure(2)
plot(fig_data(:,1),fig_data(:,2));
xlabel('迭代次数');
ylabel('适应度')
title('适应度迭代曲线');


%写出路线规划图
for i = 1 : m
    fprintf('汽车%d路径:仓库-->',i);
    for j = 2 : numel(path{i})-1
        fprintf('客户%d-->',path{i}(j));
    end
    fprintf('仓库\n');
end
fprintf('最短路径长度为%d',round(Parent(1).y,4));






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

function x = ini_x(nVar,m,n)
    x = zeros(nVar,1);
    value = randperm(nVar,n);
    for i = 1 : n
        x(value(i)) = i;
    end
end



