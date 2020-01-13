function ret=Mutation(pmutation,lenchrom,chrom,sizepop,num,maxgen,bound)
% 本函数完成变异操作
% pcorss                input  : 变异概率
% lenchrom              input  : 染色体长度
% chrom     input  : 染色体群
% sizepop               input  : 种群规模
% opts                  input  : 变异方法的选择
% pop                   input  : 当前种群的进化代数和最大的进化代数信息
% bound                 input  : 每个个体的上届和下届
% maxgen                input  ：最大迭代次数
% num                   input  : 当前迭代次数
% ret                   output : 变异后的染色体

for i=1:sizepop   %每一轮for循环中，可能会进行一次变异操作，染色体是随机选择的，变异位置也是随机选择的，
    %但该轮for循环中是否进行变异操作则由变异概率决定（continue控制）
    % 随机选择一个染色体进行变异
    pick=rand;%%rand为生成一个[0，1]的随机数，相当于rand（）
    while pick==0%%生成0时重新随机   
        pick=rand;
    end
    index=ceil(pick*sizepop);%%确定第几个染色体%%%一个个体是否可以变异多次
    % 变异概率决定该轮循环是否进行变异
    pick=rand;
    if pick>pmutation  %%若生成的随机数大于变异概率，则此轮循环不进行变异
        continue;
    end
    flag=0;
    while flag==0
        % 变异位置
        pick=rand;
        while pick==0      
            pick=rand;
        end
        pos=ceil(pick*sum(lenchrom));  %随机选择了染色体变异的位置，即选择了第pos个变量进行变异
        %%变异操作
        pick=rand; %变异开始     
        fg=(rand*(1-num/maxgen))^2;%%变异公式maxgen为最大迭代次数（最大变异次数）
        if pick>0.5
            chrom(index,pos)=chrom(index,pos)+(chrom(index,pos)-bound(pos,2))*fg;%bound(pos,2)为该个体（染色体）基因的上届%%%修改过index，因为选择哪个个体变异
        else
            chrom(index,pos)=chrom(index,pos)-(chrom(index,pos)-bound(pos,1))*fg;%bound(pos,1)为该个体（染色体）基因的下届
        end   %变异结束
        flag=test(lenchrom,bound,chrom(i,:));     %检验染色体的可行性
    end
end
ret=chrom;