function ret=Mutation(pmutation,lenchrom,chrom,sizepop,num,maxgen,bound)
% ��������ɱ������
% pcorss                input  : �������
% lenchrom              input  : Ⱦɫ�峤��
% chrom     input  : Ⱦɫ��Ⱥ
% sizepop               input  : ��Ⱥ��ģ
% opts                  input  : ���췽����ѡ��
% pop                   input  : ��ǰ��Ⱥ�Ľ������������Ľ���������Ϣ
% bound                 input  : ÿ��������Ͻ���½�
% maxgen                input  ������������
% num                   input  : ��ǰ��������
% ret                   output : ������Ⱦɫ��

for i=1:sizepop   %ÿһ��forѭ���У����ܻ����һ�α��������Ⱦɫ�������ѡ��ģ�����λ��Ҳ�����ѡ��ģ�
    %������forѭ�����Ƿ���б���������ɱ�����ʾ�����continue���ƣ�
    % ���ѡ��һ��Ⱦɫ����б���
    pick=rand;%%randΪ����һ��[0��1]����������൱��rand����
    while pick==0%%����0ʱ�������   
        pick=rand;
    end
    index=ceil(pick*sizepop);%%ȷ���ڼ���Ⱦɫ��%%%һ�������Ƿ���Ա�����
    % ������ʾ�������ѭ���Ƿ���б���
    pick=rand;
    if pick>pmutation  %%�����ɵ���������ڱ�����ʣ������ѭ�������б���
        continue;
    end
    flag=0;
    while flag==0
        % ����λ��
        pick=rand;
        while pick==0      
            pick=rand;
        end
        pos=ceil(pick*sum(lenchrom));  %���ѡ����Ⱦɫ������λ�ã���ѡ���˵�pos���������б���
        %%�������
        pick=rand; %���쿪ʼ     
        fg=(rand*(1-num/maxgen))^2;%%���칫ʽmaxgenΪ�����������������������
        if pick>0.5
            chrom(index,pos)=chrom(index,pos)+(chrom(index,pos)-bound(pos,2))*fg;%bound(pos,2)Ϊ�ø��壨Ⱦɫ�壩������Ͻ�%%%�޸Ĺ�index����Ϊѡ���ĸ��������
        else
            chrom(index,pos)=chrom(index,pos)-(chrom(index,pos)-bound(pos,1))*fg;%bound(pos,1)Ϊ�ø��壨Ⱦɫ�壩������½�
        end   %�������
        flag=test(lenchrom,bound,chrom(i,:));     %����Ⱦɫ��Ŀ�����
    end
end
ret=chrom;