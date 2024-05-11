clc;clear;

%% 导入数据
%  burma14数据
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
  14  20.09       94.55];
opt_value = 30.8785;
%运行100次记录结果
ex_t = 100;
result = zeros(ex_t,3);%记录实验结果，第一列最优值，第二列时间,3列迭代次数
parfor i = 1 : ex_t
    [best_value,time,it] = MOGA(data,opt_value);
    result(i,:) = [best_value,time,it];
end
result_burma14.optvalue = opt_value;
result_burma14.data_moga = result;
result_burma14.mean_value_moga = mean(result);

result = zeros(ex_t,3);%记录实验结果，第一列最优值，第二列时间,3列迭代次数
parfor i = 1 : ex_t
    [best_value,time,it] = BGA(data,opt_value);
    result(i,:) = [best_value,time,it];
end
result_burma14.data_bga = result;
result_burma14.mean_value_bga = mean(result);


%% eil51实验
data_name = 'eil51.tsp';
data = dlmread(data_name, ' ', 6, 0);

opt_value = 442.0471;
%运行100次记录结果
ex_t = 100;
result = zeros(ex_t,3);%记录实验结果，第一列最优值，第二列时间,3列迭代次数
parfor i = 1 : ex_t
    [best_value,time,it] = MOGA(data,opt_value);
    result(i,:) = [best_value,time,it];
end
result_eil51.optvalue = opt_value;
result_eil51.data_moga = result;
result_eil51.mean_value_moga = mean(result);

result = zeros(ex_t,3);%记录实验结果，第一列最优值，第二列时间,3列迭代次数
parfor i = 1 : ex_t
    [best_value,time,it] = BGA(data,opt_value);
    result(i,:) = [best_value,time,it];
end
result_eil51.data_bga = result;
result_eil51.mean_value_bga = mean(result);



%% kroB100实验
data_name = 'kroB100.tsp';
data = dlmread(data_name, ' ', 6, 0);

opt_value = 22141;
%运行100次记录结果
ex_t = 100;
result = zeros(ex_t,3);%记录实验结果，第一列最优值，第二列时间,3列迭代次数
parfor i = 1 : ex_t
    [best_value,time,it] = MOGA(data,opt_value);
    result(i,:) = [best_value,time,it];
end
result_kroB100.optvalue = opt_value;
result_kroB100.data_moga = result;
result_kroB100.mean_value_moga = mean(result);

result = zeros(ex_t,3);%记录实验结果，第一列最优值，第二列时间,3列迭代次数
parfor i = 1 : ex_t
    [best_value,time,it] = BGA(data,opt_value);
    result(i,:) = [best_value,time,it];
end
result_kroB100.data_bga = result;
result_kroB100.mean_value_bga = mean(result);

% 保存结果
save('result.mat','result_burma14','result_eil51','result_kroB100')


%导出
value1 = [result_burma14.optvalue 0 0 0 0 0
    result_burma14.mean_value_moga result_burma14.mean_value_bga
        result_burma14.data_moga result_burma14.data_bga];

value2 = [result_eil51.optvalue 0 0 0 0 0
    result_eil51.mean_value_moga result_eil51.mean_value_bga
        result_eil51.data_moga result_eil51.data_bga];

value3 = [result_kroB100.optvalue 0 0 0 0 0
    result_kroB100.mean_value_moga result_kroB100.mean_value_bga
        result_kroB100.data_moga result_kroB100.data_bga];


%第一行已知最优值，第二行平均值，后面每一次数据
writematrix([value1 value2 value3],'结果.csv');

