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
R = 50;               % 
alpha = 0.9999;       % 退火系数


model.x = [];
model.y = [];
model.v = [];

% 问题数据
N = 12; % 任务数
M = 3;  % 机器数
V = 2;  % 小车数
data = []; % 任务数据——任务号，机器号，处理时间

LB = cal_LB(); %可行解下界

% 小车转运时间数据
time = [];


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
d_f_p = [];
%% 主循环
while (T > TF && g<100) || Gbest.y > LB
    for r = 1 : R
        % TODO:产生邻域解
        Y = cal_larboslusion(X_SA);
        F_Y = fitness(Y, data, time);
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
        % TODO:更新粒子速度
        particles(s).v = update_v(w, c1, c2, particles(s).v, particles(s).x, Pbest(s).x, Gbest.x, V_min, V_max);
        % TODO:更新粒子位置
        particles(s).x = particles(s).x + particles(s).v;
        % 更新粒子适应度
        particles(s).y = fitness(particles(s).x, data, time);
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


%% TODO:适应度函数
function y = fitness(x, data, time)
    y = sum(x);
end

%% TODO:计算可行解下界
function lb = cal_LB()
    lb = 0;
end

%% TODO:产生邻域解
function y = cal_larboslusion(x)
    y = x;
end

%% TODO:更新速度
function v = update_v(w, c1, c2, v, x, pbest, gbest, V_min, V_max)
    v = w * v + c1 * rand * (pbest - x) + c2 * rand * (gbest - x);
    v(v < V_min) = V_min;
    v(v > V_max) = V_max;
end