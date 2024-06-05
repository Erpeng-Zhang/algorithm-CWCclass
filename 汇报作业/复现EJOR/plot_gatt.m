function y = plot_gatt(x,V, data, time)


    %% 画出甘特图
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
    n = floor(numel(x2) / V); %小车负责到的数量
    [~,index2] = sort(x2,'ascend');
    p_v = zeros(1,N);
    for i = 1 : V-1
        p_v(index2( n*(i-1)+1 : n*i)) = i;
    end
    p_v(index2( n*i+1 : N) ) = V;

    project(5,:) = p_v;

    %% 
    machine_gatt = zeros(n_M,N*2);
    car_gatt = zeros(V,N*4);
    
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
        
        car_gatt(car,i*4-3) = car_time(car);%小车可以开始下一轮任务时间
        car_gatt(car,i*4-2) = car_time(car) + time(car_location(car),location(job));%小车达到取货机器时间
        
        get = max(car_time(car) + time(car_location(car),location(job)) , machine_time(location(job)));

        car_gatt(car,i*4-1) = get;%小车取到货物可以开始运输时间

        % 送达
        car_time(car) = get + time(car_location(car),machine);
        car_gatt(car,i*4) = car_time(car);%小车送达时间
        % 机器生产 + max(到达时间,机器空闲时刻)

        machine_gatt(machine,i*2-1) = max(car_time(car),machine_time(machine));%机器开始工作时间
        machine_time(machine) = max(car_time(car),machine_time(machine)) + project(4,i);
        machine_gatt(machine,i*2) = machine_time(machine);%机器结束工作时间

        % 小车位置更新
        car_location(car) = machine;
        % 工件位置更新
        location(job) = machine;
    end
    
    Exit_time = max([machine_time car_time]);
    y = Exit_time;

    % 画出甘特图
    colormap(jet); % 使用jet颜色映射
    color = jet(N); % 从jet颜色映射中获取N个颜色
    figure(2)
    
    M_Car_lable = cell(1,n_M+V);
    for i = 2 : n_M
        M_Car_lable{i} = strcat('M',num2str(i-1));
        for j = 1 : 2 : (N-1)*2

            if machine_gatt(i,j) ~= machine_gatt(i,j+1)
                plot([machine_gatt(i,j) machine_gatt(i,j+1)],[i i],'Color',color(ceil(j/2),:),'LineWidth',20);hold on;
            end
        end
    end
    for i = 1 : V 
        M_Car_lable{n_M+i} = strcat('V',num2str(i));
        for j = 1 : 2 : (N-1)*4
            if car_gatt(i,j) ~= car_gatt(i,j+1)
                plot([car_gatt(i,j) car_gatt(i,j+1)],[i+n_M i+n_M],'Color',color(ceil(j/4),:),'LineWidth',20);hold on;
            end
            if car_gatt(i,j+2) ~= car_gatt(i,j+3)
                plot([car_gatt(i,j) car_gatt(i,j+1)],[i+n_M i+n_M],'Color',color(ceil(j/4),:),'LineWidth',20);hold on;
            end
        end
    end
    
    % set(gca,'yticklabel',{'','M1','M2','M3','M4','V1','V2'});
    set(gca,'yticklabel',M_Car_lable);
    xlabel('时间');ylabel('机器或小车');title('甘特图');
    ylim([1, n_M+V+1]);
    hold off;

    
