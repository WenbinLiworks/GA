function ret=select(individuals,sizepop)
% ��������ÿһ����Ⱥ�е�Ⱦɫ�����ѡ���Խ��к���Ľ���ͱ���
% individuals input  : ��Ⱥ��Ϣ
% sizepop     input  : ��Ⱥ��ģ
% ret         output : ����ѡ������Ⱥ

%���ݸ�����Ӧ��ֵ��������
fitness1=10./individuals.fitness;

sumfitness=sum(fitness1);%%%�������Ӧ�ȵ���֮��
sumf=fitness1./sumfitness;%%��������ĸ���
index=[];%%�վ���
for i=1:sizepop   %תsizepop������
    pick=rand;%%randΪ����һ��[0��1]����������൱��rand����
    while pick==0 %%����0ʱ�������   
        pick=rand;        
    end
    for j=1:sizepop    
        pick=pick-sumf(j); %%�������������֮��       
        if pick<0%%%������ȸ���ĸ���С���򽫸ø����������Ⱥ��������Խ��ĸ��屻ѡ�еļ���Խ��        
            index=[index j];%[index j]Ϊ��j�����������Ⱥ            
            break;  %Ѱ����������䣬�˴�ת����ѡ����Ⱦɫ��i��ע�⣺��תsizepop�����̵Ĺ����У��п��ܻ��ظ�ѡ��ĳЩȾɫ��
        end
    end
end
individuals.chrom=individuals.chrom(index,:);%individuals.chromΪ��Ⱥ�и���  ��index��������indexΪһά���󣬱�ʾ�������ֵ���������Ⱥ�ĸ����Ϊԭ����Ⱥ��ĳЩ����
individuals.fitness=individuals.fitness(index);%ѡ������¸������Ӧ��Ϊ����Ⱥ��Ӧ��
ret=individuals;
 