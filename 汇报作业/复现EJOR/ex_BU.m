%%
clc;clear;


results = zeros(10*4,4);
V = 2 ;
dataset = data(3,3,3);

datas = dataset(3).data;
times = dataset(3).time;
Optimals = dataset(3).Optimal;



for i = 1 : 10
    for j = 1 : 4
        output = PSOSA(datas{i},times{j},V,Optimals(i,j));
        fprintf('工作%d，布局%d，已知最优值%d---算法最优值%d,耗时：%d秒,GAP:%d,LB:%d\n',i,j,Optimals(i,j),output(1),round(output(2),2),output(3),output(4));
        results((i-1)*4 + j,:) = output;
    end
end

% results输出为csv
csvwrite('results.csv',results);
fprintf('done1\n');

%% 重复150次，取平均值
results_150 = zeros(10*4,4);
for i = 1 : 10
    for j = 1 : 4
        outputs = [];
        for k = 1 : 150            
            output = PSOSA(datas{i},times{j},V,Optimals(i,j));
            outputs = [outputs;otuput];
        end
        output = mean(outputs);
        fprintf('150次平均：工作%d，布局%d，已知最优值%d---算法最优值%d,耗时：%d,GAP:%d,LB:\n',i,j,Optimals(i,j),output(1),output(2),output(3),output(4));
        results_150((i-1)*4 + j,:) = mean(outputs);
    end
end

% results输出为csv
csvwrite('results_150.csv',results_150);
fprintf('done2\n');


%% 函数化的PSOSA
function output = PSOSA(data,time,V,Optimal)
    N = size(data,2); % 任务数
    LB = cal_LB(data, time); %可行解下界
    
    tic
    %% 算法参数
    w = 1 / (2 * log(2)); % 惯性权重
    c1 = 0.5 + log(2);    % 个体学习因子
    c2 = 0.5 + log(2);    % 群体学习因子
    V_min = -2;           % 速度最小值
    V_max = 2;            % 速度最大值
    S_max = 40;           % 粒子群规模
    
    TI = 100;             % 初始温度
    TF = 0.01;            % 最小温度
    R = 50;               % 模拟退火搜索次数
    alpha = 0.9999;       % 退火系数
    
    % 粒子结构体
    model.x = [];
    model.y = [];
    model.v = [];
    
    
    % 初始化粒子群
    particles = repmat(model, S_max, 1);
    
    for i = 1:S_max
        particles(i).x = rand(1,2*N);
        particles(i).y = fitness(particles(i).x,V, data, time);
        particles(i).v = rand(1,2*N);
    end
    
    % 初始化个体最优
    Pbest = particles;
    
    % 找出初始化全局最优下标a
    [~, a] = min([particles.y]);
    
    % 初始化全局最优
    Gbest = particles(a);
    
    % 模拟退火解初始化
    X_SA = particles(a).x;
    F_SA = particles(a).y;
    
    % 初始化温度和迭代次数
    T = TI; g = 0;
    %% 主循环
    while T >= TF && Gbest.y > Optimal
        for r = 1 : R
            % 产生邻域解
            Y = cal_larboslusion(X_SA,N);
            F_Y = fitness(Y, V, data, time);
            % 计算适应度变化
            delta_F = F_SA - F_Y;
            if delta_F > 0 %邻域解更好
                X_SA = Y;
                F_SA = F_Y;
            else
                if rand < exp(delta_F / T)
                    X_SA = Y;
                    F_SA = F_Y;
                    % 更新全局最优
                    if F_SA < Gbest.y
                        Gbest.x = X_SA;
                        Gbest.y = F_SA;
                    end
                else 
                    % 计算新解接受度
                    Pr = exp(delta_F / T);
                    if Pr > rand
                        X_SA = Y;
                        F_SA = F_Y;
                    end
                end
            end
        end
    
        % 更新粒子
        for s = 1 : S_max
            % 更新粒子速度
            particles(s).v = update_v(w, c1, c2, particles(s).v, particles(s).x, Pbest(s).x, Gbest.x, V_min, V_max);
            % TODO:更新粒子位置
            particles(s).x = particles(s).x + particles(s).v;
            % 更新粒子适应度
            particles(s).y = fitness(particles(s).x, V, data, time);
            % 更新个体最优
            if particles(s).y < Pbest(s).y
                Pbest(s).x = particles(s).x;
                Pbest(s).y = particles(s).y;
            end
        end
    
        % 更新全局最优
        [~, a] = min([particles.y]);
        if particles(a).y < Gbest.y
            Gbest.x = particles(a).x;
            Gbest.y = particles(a).y;
        end
    
        % 更新X_SA和F_SA
        X_SA = Gbest.x;
        F_SA = Gbest.y;
    
        % 更新温度和迭代次数
        T = alpha * T;
        g = g + 1;
    end
    t = toc;

    GAP = (Gbest.y-Optimal) / Optimal;
    output = [Gbest.y, t, GAP, LB];
end

%% 计算可行解下界
function lb = cal_LB(data, time)
    J = data(1,:); % 工件
    n_J = max(J); %工件数
    exit_time = zeros(1,n_J);
    for i = 1 : n_J
        [~,ind] = find(J == i);
        value = data(:,ind);
        exit_time(i) = sum(value(4,:));
        ms = value(3,:)+1;
        for k = 1 : numel(ms)-1
            exit_time(i) = exit_time(i) + time(ms(k),ms(k+1));
        end
    end
    lb = max(exit_time);
end

%% 产生邻域解
function y = cal_larboslusion(x,N)
    x1 = x(1 : N);x2=x(N+1:end);
    if rand < 0.5
        ind = randperm(N,2);
        dots = x1(ind);
        dots = flip(dots);
        x1(ind) = dots;
    else
        ind = randperm(N,2);
        dots = x2(ind);
        dots = flip(dots);
        x2(ind) = dots;
    end
    y = [x1 x2];
end

%% 更新速度
function v = update_v(w, c1, c2, v, x, pbest, gbest, V_min, V_max)
    p_sg = x + c1 .* rand(size(x)) .* (pbest - x); % 个体最优方向上一点
    l_sg = x + c2 .* rand(size(x)) .* (gbest - x); % 全局最优方向上一点
    G_sg = (x + p_sg + l_sg) / 3; % 超球中心点
    r = (sum((x-G_sg).^2)).^0.5; % 超球半径; r = ||x-G_sg||_1;范数，记p=2
    X = G_sg + r;         % 超球表面点
    v = w * v + X - x; % 速度更新
    v(v < V_min) = V_min;
    v(v > V_max) = V_max;
end