function y = fitness(x,V, data, time)
    % x = rand(1,24);
    % V = 2;
    % data = [1 1 1 1 2 2 2 2 3 3 3 3  
    %         1 2 3 4 5 6 7 8 9 10 11 12
    %         1 3 2 0 3 1 4 0 2 4 1 0
    %         10 16 18 0 15 22 18 0 18 16 22 0];
    % time = [0 6 8 10 12
    %         12 0 6 8 10
    %         10 6 0 6 8
    %         8 8 6 0 6
    %         6 10 8 6 0];


    %% TODO:适应度函数
    J = data(1,:); % 工件
    n_J = max(J); %工件数
    j = data(2,:); % 工序

    M = data(3,:)+1; % 工序对应机器
    n_M = max(M); %机器数
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
    % n = floor(numel(x2) / V); %小车负责到的数量
    % [~,index2] = sort(x2,'ascend');
    % p_v = zeros(1,N);
    % for i = 1 : V-1
    %     p_v(index2( n*(i-1)+1 : n*i)) = i;
    % end
    % p_v(index2( n*i+1 : N) ) = V;

    p_v = ceil(x2);
    project(5,:) = p_v;

    if max(p_v) > 2
        a = 1;
    end

    % 各工件退出时间
    machine_time = zeros(1,n_M);
    car_time = zeros(1,V);
    car_location = ones(1,V);%开始时小车都在转运中心
    location = ones(n_J);

    for i = 1 : N
        job = project(1,i);%工件
        sub_job = project(2,i);%工序
        machine = project(3,i);%机器
        car = project(5,i);%小车

        % 小车转运
        % 转运时间 = 取+送
        get = max(car_time(car) + time(car_location(car),location(job)) , machine_time(location(job)));
        % 送达
        car_time(car) = get + time(car_location(car),machine);
        % 机器生产 + max(到达时间,机器空闲时刻)
        machine_time(machine) = max(car_time(car),machine_time(machine)) + project(4,i);
        % 小车位置更新
        car_location(car) = machine;
        % 工件位置更新
        location(job) = machine;
    end
    
    Exit_time = max([machine_time car_time]);
    y = Exit_time;
