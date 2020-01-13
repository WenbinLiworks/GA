function ret=select(individuals,sizepop)
% 本函数对每一代种群中的染色体进行选择，以进行后面的交叉和变异
% individuals input  : 种群信息
% sizepop     input  : 种群规模
% ret         output : 经过选择后的种群

%根据个体适应度值进行排序
fitness1=10./individuals.fitness;

sumfitness=sum(fitness1);%%%个体的适应度倒数之和
sumf=fitness1./sumfitness;%%各个个体的概率
index=[];%%空矩阵
for i=1:sizepop   %转sizepop次轮盘
    pick=rand;%%rand为生成一个[0，1]的随机数，相当于rand（）
    while pick==0 %%生成0时重新随机   
        pick=rand;        
    end
    for j=1:sizepop    
        pick=pick-sumf(j); %%随机数与个体概率之差       
        if pick<0%%%随机数比个体的概率小，则将该个体加入新种群。（概率越大的个体被选中的几率越大）        
            index=[index j];%[index j]为将j个体加入新种群            
            break;  %寻找落入的区间，此次转轮盘选中了染色体i，注意：在转sizepop次轮盘的过程中，有可能会重复选择某些染色体
        end
    end
end
individuals.chrom=individuals.chrom(index,:);%individuals.chrom为种群中个体  （index，：）中index为一维矩阵，表示里面数字的行数，种群的个体变为原来种群的某些个体
individuals.fitness=individuals.fitness(index);%选择出的新个体的适应度为新种群适应度
ret=individuals;
 