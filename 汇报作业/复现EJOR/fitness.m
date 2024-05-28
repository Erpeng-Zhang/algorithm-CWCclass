% function y = fitness(x,V, data, time)
    x = rand(1,24);
    V = 2;
    data = [1 1 1 1 2 2 2 2 3 3 3 3  
            1 2 3 4 5 6 7 8 9 10 11 12
            1 3 2 5 3 1 4 5 2 4 1 5
            10 16 18 0 15 22 18 0 18 16 22 0];
    time = [0 6 8 10 12
            12 0 6 8 10
            10 6 0 6 8
            8 8 6 0 6
            6 10 8 6 0];


    %% TODO:适应度函数
    J = data(1,:); % 工件
    n_J = max(J); %工件数
    j = data(2,:); % 工序
    M = data(3,:); % 工序对应机器
    process_time = data(4,:); %操作时间
    N = numel(x)/2;
    
    % 机器工序方案
    x1 = x(1 : N);     
    [~,index] = sort(x1,'ascend');
    p = J(index);
    % 方案:每行含义：工件；工序；机器；操作时间；小车
    project = zeros(5,N);
    project(1,:) = p;
    for i = 1 : n_J
           [~,ind] = find(p == i);
           project(2,ind) = j(J == i);
           project(3,ind) = M(J == i);
           project(4,ind) = process_time(J == i);
    end

    
    % 小车方案
    x2 = x(N+1 : end); 
    n = floor(numel(x2) / V); %小车负责到的数量
    [~,index2] = sort(x2,'ascend');
    p_v = zeros(1,N);
    for i = 1 : V-1
        p_v(index2( n*(i-1)+1 : n*i)) = i;
    end
    p_v(index2( n*i+1 : N) ) = V;

    project(5,:) = p_v;

    
    % 各工件退出时间
    Exit_time = zeros(1,n_J);
    machine_time = zeros(1,5);

    for i = 1 : N
        job = project(1,i);%工件
        sub_job = project(2,i);%工序

    end
    

    time = zeros(1,J);
    y = sum(x);
