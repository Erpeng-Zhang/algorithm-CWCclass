function dataset = data(M,N,J)
% 数据集
% M: 机器数
% N：工件数
% J：工序数

dataset.data = [];
dataset.time = [];


dataset(1).data = [1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,5,5,5;
                   1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18;
                   1,2,4,0,1,3,2,0,3,4,1,0,4,2,0,3,1,0;
                   8,16,12,0,20,10,18,0,12,8,15,0,14,18,0,10,15,0];

dataset(1).time = [0,22,18,6,8,20,7,25,31,12,27;
                    22,0,11,27,30,32,20,10,28,21,7;
                    18,11,0,22,26,22,21,7,18,12,19;
                    6,27,22,0,14,8,12,26,10,16,21;
                    8,30,26,14,0,18,24,32,22,28,16;
                    20,32,22,8,18,0,14,26,12,10,22;
                    7,20,21,12,24,14,0,18,22,16,10;
                    25,10,7,26,32,26,18,0,21,14,22;
                    31,28,18,10,22,12,22,21,0,8,14;
                    12,21,12,16,28,10,16,14,8,0,18;
                    27,7,19,21,16,22,10,22,14,18,0];
J = J;
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
