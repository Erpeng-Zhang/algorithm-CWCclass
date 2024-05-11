clc;clear;

% 开始计时
tic;

%% 导入数据
% burma14数据
data = [ 1  16.47       96.10
   2  16.47       94.44
   3  20.09       92.54
   4  22.39       93.37
   5  25.23       97.24
   6  22.00       96.05
   7  20.47       97.02
   8  17.20       96.29
   9  16.30       97.38
  10  14.05       98.12
  11  16.53       97.38
  12  21.52       95.59
  13  19.41       97.13
  14  20.09       94.55];opt_value = 30.8785;

% data_name = 'eil51.tsp';opt_value = 442.0471;
% data_name = 'kroB100.tsp';opt_value = 22141;
% data_name = 'pr1002.tsp';opt_value = 259045;

% data_name = 'kroC100.tsp';
% 
% data = dlmread(data_name, ' ', 6, 0);

%城市数
n = size(data,1);

%位置坐标
X = (data(:,2:3));


D=Distanse(X);  %生成距离矩阵

% 在二维图上画出所有坐标点
figure(1)
plot(X(:,1), X(:,2),'red*');hold on; % 保持图形窗口开启
% 在点旁边标上序号
for i = 1:numel(X(:,1))
    text(X(i,1), X(i,2), num2str(i));
end
hold on;



nVar = n;%染色体数

nPop = 100;  %种群规模
maxIt = 1000;    %最大迭代次数
nPc = 1;      %子代比例
nC = round(nPop* nPc/2)*4; %子代种群规模
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


% 主循环
for It = 1:maxIt
    %生成子代种群
    Offspring =repmat(template, nC/4, 4);
    
    selected_id = [];%被选中精英ID
    for j = 1: size(Offspring,1)
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

    % 四列变一列
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
   
    disp(['迭代次数:', num2str(It), ',最小值为：', num2str(Parent(1).y)])
    
    fig_data(It,1) = It;fig_data(It,2) = Parent(1).y;

    if round(Parent(1).y,4)<=opt_value
        break;
    end
    
end

toc;


%% 最优解的路线图
x = Parent(1).x;
path = [x;x(1)];
plot(X(path,1),X(path,2),'-',LineWidth=.5);hold on;
title('路径图');
xlabel('x');ylabel('y');
hold off;

% 画出适应度迭代曲线
figure(2)
plot(fig_data(:,1),fig_data(:,2));
xlabel('迭代次数');
ylabel('适应度')
title('适应度迭代曲线');


%写出路线规划图
fprintf('最短路径为：%d',path(1));
for i = 2 : numel(path)
    fprintf('-->%d',path(i));
end
fprintf('\n最短路径长度为%d',round(Parent(1).y,4));


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





