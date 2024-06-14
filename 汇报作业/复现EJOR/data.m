function dataset = data(M,N,J)
% 数据集：包括案例数据和随机生成数据和BU数据集数据

% 参数为随机生成数据集所需参数
% M: 机器数
% N：工件数
% J：各工件工序数

dataset.data = [];
dataset.time = [];
dataset.Optimal = [];
%% 论文案例数据
%% 任务数据每行含义——工件号，任务号，机器号，处理时间
dataset(1).data = [1 1 1 1 2 2 2 2 3 3 3 3
                    1 2 3 4 5 6 7 8 9 10 11 12
                    1 3 2 0 3 1 4 0 2 4 1 0
                    10 16 18 0 15 22 18 0 18 16 22 0];

dataset(1).time = [0 6 8 10 12
                    12 0 6 8 10
                    10 6 0 6 8
                    8 8 6 0 6
                    6 10 8 6 0];
dataset(1).Optimal = 108;
%% 随机生成数据
data = zeros(4,N*J);

for i = 1 : N
    data(1,(i-1)*(J+1)+1:i*(J+1)) = i;
    data(2,(i-1)*(J+1)+1:i*(J+1)) = [(i-1)*J+1:i*J,0];
    for j = (i-1)*(J+1)+1:i*(J+1)-1
        data(3,j) = randi(M);
        data(4,j) = randi(100);
    end
    data(3,j+1) = 0;
    data(4,j+1) = 0;
end

time = rand(M+1);
time = round(time.*30);
time(time==0) = 1;
for j = 1 : M+1
    time(j,j) = 0;
end

dataset(2).data = data;
dataset(2).time = time;
dataset(2).Optimal = [];


%% BU数据集
job_set1=[1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,5,5,5;
                   1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18;
                   1,2,4,0,1,3,2,0,3,4,1,0,4,2,0,3,1,0;
                   8,16,12,0,20,10,18,0,12,8,15,0,14,18,0,10,15,0];
job_set2=[1,1,1,2,2,2,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6;
                  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21;
                  1,4,0,2,4,0,1,3,0,2,3,4,0,1,2,4,0,1,2,3,0;
                  10,18,0,10,18,0,10,20,0,10,15,12,0,10,15,12,0,10,15,12,0];
job_set3=[1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,5,5,6,6,6,6,6;
                  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22;
                  1,3,0,2,4,0,1,2,0,3,4,0,1,2,3,4,0,2,3,4,1,0;
                  16,15,0,18,15,0,20,10,0,15,10,0,8,10,15,17,0,10,15,8,15,0];
job_set4=[1,1,1,1,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5,5;
                  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24;
                  4,1,2,0,3,2,4,0,2,3,1,3,0,2,4,1,2,0,1,2,4,2,3,0;
                  11,10,7,0,12,10,8,0,7,10,9,8,0,7,8,12,6,0,9,7,8,10,8,0];
job_set5=[1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,5,5,5;
                  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18;
                  1,2,4,0,1,3,2,0,3,4,1,0,4,2,0,3,1,0;
                  6,12,9,0,18,6,15,0,9,3,12,0,6,15,0,3,9,0];
job_set6=[1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6;
                  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24;
                  1,2,4,0,1,2,4,0,2,3,4,0,2,3,4,0,1,3,4,0,1,3,4,0;
                  9,11,7,0,19,20,13,0,14,20,9,0,14,20,9,0,11,16,8,0,10,12,10,0];
job_set7=[1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,6,7,7,7,7,8,8,8,8;
                  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27;
                  1,4,0,2,4,0,2,4,0,3,4,0,1,3,0,2,3,4,0,1,2,3,0,1,2,4,0;
                  6,6,0,11,9,0,9,7,0,16,7,0,9,18,0,13,19,6,0,10,9,13,0,11,9,8,0];
job_set8=[1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,5,6,6,6,6,6;
                  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26;
                  2,3,4,0,2,3,4,0,2,3,4,0,2,3,4,0,1,2,3,4,0,1,2,3,4,0;
                  12,21,11,0,12,21,11,0,12,21,11,0,12,21,11,0,10,14,18,9,0,10,14,18,9,0];
job_set9=[1,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,5;
                  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,19,17,18,19,20,21,22;
                  3,1,2,4,0,3,2,4,0,1,2,4,0,2,3,4,0,3,1,3,4,0;
                  9,12,9,6,0,16,11,9,0,21,18,7,0,20,22,11,0,14,16,13,9,0];
job_set10=[1,1,1,1,1,2,2,2,2,3,3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6,6;
                  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27;
                  1,3,2,4,0,2,3,4,0,3,2,1,4,0,2,3,4,0,1,3,4,0,2,1,3,4,0;
                  11,19,16,13,0,21,16,14,0,8,10,14,9,0,13,20,10,0,9,16,18,0,19,21,11,15,0];
layout_1=[0,6,8,10,12;
          12,0,6,8,10;
          10,6,0,6,8;
          8,8,6,0,6;
          6,10,8,6,0];
layout_2=[0,4,6,8,6;
          6,0,2,4,2;
          8,12,0,2,4;
          6,10,12,0,2;
          4,8,10,12,0];
layout_3=[0,2,4,10,12;
          12,0,2,8,10;
          10,12,0,6,8;
          4,6,8,0,2;
          2,4,6,12,0];
layout_4=[0,4,8,10,14;
          18,0,4,6,10;
          20,14,0,8,6;
          12,8,6,0,6;
          14,14,12,6,0];

Optimals =[114,90,98,140;
          116,82,89,134;
          121,89,96,148;
          136,100,102,163;
          110,81,89,134;
          129,102,105,151;
          146,86,93,169;
          167,155,155,178;
          127,106,107,149;
          153,139,139,183];

datas = {job_set1,job_set2,job_set3,job_set4,job_set5,job_set6,job_set7,job_set8,job_set9,job_set10};
layouts = {layout_1,layout_2,layout_3,layout_4};

dataset(3).data = datas;
dataset(3).time = layouts;
dataset(3).Optimal = Optimals;

end



