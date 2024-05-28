%% 对《A hybrid particle swarm optimization and simulated annealing algorithm for the job shop scheduling problem with transport resources》
%  中的算法进行复现

clc;clear;


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


model.x = [];
model.y = [];
model.v = [];

% 问题数据
N = 12; % 任务数
M = 3;  % 机器数
V = 2;  % 小车数
J1 = [1,3,2,5;
    10,16,18,0]; % 工件1
J2 = [3,1,4,5;
    15,22,18,0];  % 工件2
J3 = [2,4,1,5;
    18,16,22,0];  % 工件3
data ={J1,J2,J3} ; % 任务数据——任务号，机器号，处理时间

LB = cal_LB(); %可行解下界

% 小车转运时间数据
time = [0 6 8 10 12
        12 0 6 8 10
        10 6 0 6 8
        8 8 6 0 6
        6 10 8 6 0];


% 初始化粒子群
particles = repmat(model, S_max, 1);

for i = 1:S_max
    particles(i).x = rand(1,2*N);
    particles(i).y = fitness(particles(i).x, data, time);
    particles(i).v = rand(1,2*N);
end

% 初始化个体最优
Pbest = particles;

% 找出初始化全局最优下标a
[~, a] = min([particles.y]);

% 初始化全局最优
Gbest = particles(a);

X_SA = particles(a).x;
F_SA = particles(a).y;

% 初始化温度和迭代次数
T = TI; g = 0;
d_f_p = [];%存储迭代过程中的最优解  用于绘图
%% 主循环
while (T > TF && g<100) %|| Gbest.y > LB
    for r = 1 : R
        % TODO:产生邻域解
        Y = cal_larboslusion(X_SA);
        F_Y = fitness(Y, data, V, time);
        % 计算适应度变化
        delta_F = F_SA - F_Y;
        if delta_F > 0
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

    % 每次迭代输出最优解
    disp(['迭代次数：', num2str(g), '，当前温度：', num2str(T), '，最优解：', num2str(Gbest.y)]);
    % 存储绘图
    d_f_p = [d_f_p; g,Gbest.y];
end

figure(1)
plot(d_f_p(:,1),d_f_p(:,2),'-');
xlabel('迭代次数');ylabel('适应度值');
title('适应度迭代曲线');




%% TODO:计算可行解下界
function lb = cal_LB()
    lb = 0;
end

%% 产生邻域解
function y = cal_larboslusion(x)
    value = x(1 : numel(x)/2);
    [max_v,max_i] = max(value);
    [min_v,min_i] = min(value);
    x(max_i) = min_v;
    x(min_i) = max_v;
    y = x;
end

%% 更新速度
function v = update_v(w, c1, c2, v, x, pbest, gbest, V_min, V_max)
    p_sg = x + c1 .* rand(size(x)) .* (pbest - x); % 个体最优方向上一点
    l_sg = x + c2 .* rand(size(x)) .* (gbest - x); % 全局最优方向上一点
    G_sg = (x + p_sg + l_sg) / 3; % 超球中心点
    r = sum(abs(x-G_sg)); % 超球半径; r = ||x-G_sg||_1;范数，记p=1
    X = G_sg + r;         % 超球表面点
    v = w * v + X - x; % 速度更新
    v(v < V_min) = V_min;
    v(v > V_max) = V_max;
end