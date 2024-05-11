clc;clear;

%写入数据
[m,n,X,demands,volume] = data();


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
maxIt = 2000;    %最大迭代次数
nPc = 0.8;      %子代比例
nC = round(nPop* nPc/2)*2; %子代种群规模
mu = 0.01;  %变异概率

template.x = [];
template.y = [];

%生成初始种群
Parent = repmat(template, nPop, 1);


%适应度变化曲线数据
fig_data = [];

%初始化种群
for i = 1:nPop
    %生成变量
    Parent(i).x = ini_x(nVar,m,n,demands,volume);
    %计算适应度
    Parent(i).y = fitness(Parent(i).x,D,n);
    
end

for It = 1:maxIt
    %生成子代种群
    Offspring =repmat(template, nC/2, 2);
    
    for j = 1: nC / 2
        %选择
        p1 = selectPop(Parent);
        p2 = selectPop(Parent);
        %交叉
        [Offspring(j,1).x, Offspring(j,2).x] = crossPop(p1.x, p2.x,m,n,demands,volume);      
    end

    %两列变一列
    Offspring = Offspring(:);
    
    for k = 1 : nC
        %变异
        Offspring(k).x = mutatePop( Offspring(k).x, mu,m,n,demands,volume);  
        Offspring(k).y = fitness(Offspring(k).x,D,n);        
    end
    
    %新种群
    newPop = [Parent; Offspring];
    %按适应度排序
    [~, so] = sort([newPop.y], 'ascend');
    newPop = newPop(so);
    %末位淘汰；使种群规模稳定
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


figure(2)
plot(fig_data(:,1),fig_data(:,2));
xlabel('迭代次数');
ylabel('适应度')
title('适应度迭代曲线');


%写出路线图
for i = 1 : m
    fprintf('汽车%d路径:仓库-->',i);
    for j = 2 : numel(path{i})-1
        fprintf('客户%d-->',path{i}(j));
    end
    fprintf('仓库;载重%d\n',sum(demands(path{i}(2:end-1))));
end

fprintf('最短路径长度为%d',round(Parent(1).y,4));






%% 计算两两城市之间的距离
%输入 a  各城市的位置坐标
%输出 D  两两城市之间的距离
function D=Distanse(a)
    row=size(a,1);
    D=zeros(row,row);
    for i=1:row
        for j=i+1:row
            D(i,j)=((a(i,1)-a(j,1))^2+(a(i,2)-a(j,2))^2)^0.5;
            D(j,i)=D(i,j);
        end
    end
end

% 生成染色体
function x = ini_x(nVar,m,n,demands,volume)  

    condition = true;
    while condition 
        x = zeros(nVar,1);
        value = randperm(nVar,n);
        for i = 1 : n
            x(value(i)) = i;
        end
        condition  = ~isOK(x,m,n,demands,volume);
    end
     
end



