
%% 该代码为基于遗传算法神经网络的预测代码
% 清空环境变量
clc
clear
% 

%%%%%问题1.选择
%%%%%问题2.同一个体是否可以变异两次
%% 网络结构建立
%读取数据
load SCRdata input output

%节点个数
inputnum=4;    %4;                                            %%%%%%%%%%%%%%%%%%%%%%%可改
hiddennum= 7;      %=7;                                          %%%%%%%%%%%%%%%%%%%%%%%可改
outputnum =1;        %=1;                                           %%%%%%%%%%%%%%%%%%%%%%%可改

%训练数据和预测数据
output=output';%%%%%%%%%%%%%%%%将输出改为列向量；输入矩阵必须也是  输入节点个数*数据组数
input_train=input(1:190,:)';                                          %%%%%%%%%%%%%%%%%%%%%%%可改
input_test=input(191:200,:)';                                          %%%%%%%%%%%%%%%%%%%%%%%可改
output_train=output(1:190)';                                           %%%%%%%%%%%%%%%%%%%%%%%可改
output_test=output(191:200)';                                           %%%%%%%%%%%%%%%%%%%%%%%可改

%选连样本输入输出数据归一化
[inputn,inputps]=mapminmax(input_train);%%%（ai-a_min）/（a_max-a_min）
[outputn,outputps]=mapminmax(output_train);

%构建网络
net=newff(inputn,outputn,hiddennum);

%% 遗传算法参数初始化
maxgen=100;                         %进化代数，即迭代次数  100          %%可改
sizepop=10;                        %种群规模     10                 %%可改
pcross=0.8;                       %交叉概率选择，0和1之间           %%可改
pmutation=0.1;                    %变异概率选择，0和1之间           %%可改

%节点总数
numsum=inputnum*hiddennum+hiddennum+hiddennum*outputnum+outputnum;%%个体中所含实数个数，即权值与阈值数量之和

lenchrom=ones(1,numsum);%个体长度        
bound=[-3*ones(numsum,1) 3*ones(numsum,1)];    %数据范围        %%可改

%------------------------------------------------------种群初始化--------------------------------------------------------
    individuals=struct('fitness',zeros(1,sizepop), 'chrom',[]); 
    %将种群信息定义为一个结构体%%fitness适应度设为0，chrom为种群个体
avgfitness=[];                     
%每一代种群的平均适应度
bestfitness=[];                    
%每一代种群的最佳适应度
    bestchrom=[];     %适应度最好的染色体
    
%初始化种群
for i=1:sizepop
    %随机产生一个种群
    individuals.chrom(i,:)=Code(lenchrom,bound);    %编码（binary和grey的编码结果为一个实数，float的编码结果为一个实数向量）
    x=individuals.chrom(i,:);%%便于输入fun函数
    %计算适应度
    individuals.fitness(i)=fun(x,inputnum,hiddennum,outputnum,net,inputn,outputn);   %染色体的适应度
end
FitRecord=[];
%找最好的染色体（个体）
[bestfitness,bestindex]=min(individuals.fitness);%%左边为最小的适应度值，右边index为对应的位置（染色体位置）
bestchrom=individuals.chrom(bestindex,:);  %最好的染色体
avgfitness=sum(individuals.fitness)/sizepop; %染色体的平均适应度
% 记录每一代进化中最好的适应度和平均适应度
trace=[avgfitness bestfitness]; 
 
%% 迭代求解最佳初始阀值和权值
% 进化开始
for i=1:maxgen
    
    %选择
    individuals=Select(individuals,sizepop); 
%     avgfitness=sum(individuals.fitness)/sizepop;
    %交叉
    individuals.chrom=Cross(pcross,lenchrom,individuals.chrom,sizepop,bound);
    % 变异
    individuals.chrom=Mutation(pmutation,lenchrom,individuals.chrom,sizepop,i,maxgen,bound);
    
    % 计算适应度 
    for j=1:sizepop
        x=individuals.chrom(j,:); %解码
        individuals.fitness(j)=fun(x,inputnum,hiddennum,outputnum,net,inputn,outputn);   
    end
    
  %找到最小和最大适应度的染色体及它们在种群中的位置
    [newbestfitness,newbestindex]=min(individuals.fitness);
    [worestfitness,worestindex]=max(individuals.fitness);
    % 代替上一次进化中最好的染色体
    if bestfitness>newbestfitness%%若本文适应度越小越好
        bestfitness=newbestfitness;
        bestchrom=individuals.chrom(newbestindex,:);
    end
    individuals.chrom(worestindex,:)=bestchrom;%%最差的染色体用最好的代替
    individuals.fitness(worestindex)=bestfitness;
    
    avgfitness=sum(individuals.fitness)/sizepop;
    
    trace=[trace;avgfitness bestfitness]; %记录每一代进化中最好的适应度和平均适应度
    FitRecord=[FitRecord;individuals.fitness];%每一代的适应度
end

%% 遗传算法结果分析 
figure(1)
[r,c]=size(trace);
plot([1:r]',trace(:,2),'b--');
title(['适应度曲线  ' '终止代数＝' num2str(maxgen)]);
xlabel('进化代数');ylabel('适应度');
legend('平均适应度','最佳适应度');
disp('适应度                   变量');

%% 把最优初始阀值权值赋予网络预测
% %用遗传算法优化的BP网络进行值预测
x=bestchrom;%%%%确定最优个体（染色体）
w1=x(1:inputnum*hiddennum);
B1=x(inputnum*hiddennum+1:inputnum*hiddennum+hiddennum);
w2=x(inputnum*hiddennum+hiddennum+1:inputnum*hiddennum+hiddennum+hiddennum*outputnum);
B2=x(inputnum*hiddennum+hiddennum+hiddennum*outputnum+1:inputnum*hiddennum+hiddennum+hiddennum*outputnum+outputnum);

net.iw{1,1}=reshape(w1,hiddennum,inputnum);
net.lw{2,1}=reshape(w2,outputnum,hiddennum);
net.b{1}=reshape(B1,hiddennum,1);
net.b{2}=B2;

%% BP网络训练
%网络进化参数
net.trainParam.epochs = 500;      %=100;
net.trainParam.lr=0.1;
net.trainParam.goal = 0.00001;    %=0.00001;

%网络训练
[net,per2]=train(net,inputn,outputn);%%%per2为训练过程记录

%% BP网络预测
%数据归一化
inputn_test=mapminmax('apply',input_test,inputps);
an=sim(net,inputn_test);
test_simu=mapminmax('reverse',an,outputps);%%输出结果反归一化
error=sum(abs(test_simu-output_test))
