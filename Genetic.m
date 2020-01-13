
%% �ô���Ϊ�����Ŵ��㷨�������Ԥ�����
% ��ջ�������
clc
clear
% 

%%%%%����1.ѡ��
%%%%%����2.ͬһ�����Ƿ���Ա�������
%% ����ṹ����
%��ȡ����
load SCRdata input output

%�ڵ����
inputnum=4;    %4;                                            %%%%%%%%%%%%%%%%%%%%%%%�ɸ�
hiddennum= 7;      %=7;                                          %%%%%%%%%%%%%%%%%%%%%%%�ɸ�
outputnum =1;        %=1;                                           %%%%%%%%%%%%%%%%%%%%%%%�ɸ�

%ѵ�����ݺ�Ԥ������
output=output';%%%%%%%%%%%%%%%%�������Ϊ������������������Ҳ��  ����ڵ����*��������
input_train=input(1:190,:)';                                          %%%%%%%%%%%%%%%%%%%%%%%�ɸ�
input_test=input(191:200,:)';                                          %%%%%%%%%%%%%%%%%%%%%%%�ɸ�
output_train=output(1:190)';                                           %%%%%%%%%%%%%%%%%%%%%%%�ɸ�
output_test=output(191:200)';                                           %%%%%%%%%%%%%%%%%%%%%%%�ɸ�

%ѡ����������������ݹ�һ��
[inputn,inputps]=mapminmax(input_train);%%%��ai-a_min��/��a_max-a_min��
[outputn,outputps]=mapminmax(output_train);

%��������
net=newff(inputn,outputn,hiddennum);

%% �Ŵ��㷨������ʼ��
maxgen=100;                         %��������������������  100          %%�ɸ�
sizepop=10;                        %��Ⱥ��ģ     10                 %%�ɸ�
pcross=0.8;                       %�������ѡ��0��1֮��           %%�ɸ�
pmutation=0.1;                    %�������ѡ��0��1֮��           %%�ɸ�

%�ڵ�����
numsum=inputnum*hiddennum+hiddennum+hiddennum*outputnum+outputnum;%%����������ʵ����������Ȩֵ����ֵ����֮��

lenchrom=ones(1,numsum);%���峤��        
bound=[-3*ones(numsum,1) 3*ones(numsum,1)];    %���ݷ�Χ        %%�ɸ�

%------------------------------------------------------��Ⱥ��ʼ��--------------------------------------------------------
    individuals=struct('fitness',zeros(1,sizepop), 'chrom',[]); 
    %����Ⱥ��Ϣ����Ϊһ���ṹ��%%fitness��Ӧ����Ϊ0��chromΪ��Ⱥ����
avgfitness=[];                     
%ÿһ����Ⱥ��ƽ����Ӧ��
bestfitness=[];                    
%ÿһ����Ⱥ�������Ӧ��
    bestchrom=[];     %��Ӧ����õ�Ⱦɫ��
    
%��ʼ����Ⱥ
for i=1:sizepop
    %�������һ����Ⱥ
    individuals.chrom(i,:)=Code(lenchrom,bound);    %���루binary��grey�ı�����Ϊһ��ʵ����float�ı�����Ϊһ��ʵ��������
    x=individuals.chrom(i,:);%%��������fun����
    %������Ӧ��
    individuals.fitness(i)=fun(x,inputnum,hiddennum,outputnum,net,inputn,outputn);   %Ⱦɫ�����Ӧ��
end
FitRecord=[];
%����õ�Ⱦɫ�壨���壩
[bestfitness,bestindex]=min(individuals.fitness);%%���Ϊ��С����Ӧ��ֵ���ұ�indexΪ��Ӧ��λ�ã�Ⱦɫ��λ�ã�
bestchrom=individuals.chrom(bestindex,:);  %��õ�Ⱦɫ��
avgfitness=sum(individuals.fitness)/sizepop; %Ⱦɫ���ƽ����Ӧ��
% ��¼ÿһ����������õ���Ӧ�Ⱥ�ƽ����Ӧ��
trace=[avgfitness bestfitness]; 
 
%% ���������ѳ�ʼ��ֵ��Ȩֵ
% ������ʼ
for i=1:maxgen
    
    %ѡ��
    individuals=Select(individuals,sizepop); 
%     avgfitness=sum(individuals.fitness)/sizepop;
    %����
    individuals.chrom=Cross(pcross,lenchrom,individuals.chrom,sizepop,bound);
    % ����
    individuals.chrom=Mutation(pmutation,lenchrom,individuals.chrom,sizepop,i,maxgen,bound);
    
    % ������Ӧ�� 
    for j=1:sizepop
        x=individuals.chrom(j,:); %����
        individuals.fitness(j)=fun(x,inputnum,hiddennum,outputnum,net,inputn,outputn);   
    end
    
  %�ҵ���С�������Ӧ�ȵ�Ⱦɫ�弰��������Ⱥ�е�λ��
    [newbestfitness,newbestindex]=min(individuals.fitness);
    [worestfitness,worestindex]=max(individuals.fitness);
    % ������һ�ν�������õ�Ⱦɫ��
    if bestfitness>newbestfitness%%��������Ӧ��ԽСԽ��
        bestfitness=newbestfitness;
        bestchrom=individuals.chrom(newbestindex,:);
    end
    individuals.chrom(worestindex,:)=bestchrom;%%����Ⱦɫ������õĴ���
    individuals.fitness(worestindex)=bestfitness;
    
    avgfitness=sum(individuals.fitness)/sizepop;
    
    trace=[trace;avgfitness bestfitness]; %��¼ÿһ����������õ���Ӧ�Ⱥ�ƽ����Ӧ��
    FitRecord=[FitRecord;individuals.fitness];%ÿһ������Ӧ��
end

%% �Ŵ��㷨������� 
figure(1)
[r,c]=size(trace);
plot([1:r]',trace(:,2),'b--');
title(['��Ӧ������  ' '��ֹ������' num2str(maxgen)]);
xlabel('��������');ylabel('��Ӧ��');
legend('ƽ����Ӧ��','�����Ӧ��');
disp('��Ӧ��                   ����');

%% �����ų�ʼ��ֵȨֵ��������Ԥ��
% %���Ŵ��㷨�Ż���BP�������ֵԤ��
x=bestchrom;%%%%ȷ�����Ÿ��壨Ⱦɫ�壩
w1=x(1:inputnum*hiddennum);
B1=x(inputnum*hiddennum+1:inputnum*hiddennum+hiddennum);
w2=x(inputnum*hiddennum+hiddennum+1:inputnum*hiddennum+hiddennum+hiddennum*outputnum);
B2=x(inputnum*hiddennum+hiddennum+hiddennum*outputnum+1:inputnum*hiddennum+hiddennum+hiddennum*outputnum+outputnum);

net.iw{1,1}=reshape(w1,hiddennum,inputnum);
net.lw{2,1}=reshape(w2,outputnum,hiddennum);
net.b{1}=reshape(B1,hiddennum,1);
net.b{2}=B2;

%% BP����ѵ��
%�����������
net.trainParam.epochs = 500;      %=100;
net.trainParam.lr=0.1;
net.trainParam.goal = 0.00001;    %=0.00001;

%����ѵ��
[net,per2]=train(net,inputn,outputn);%%%per2Ϊѵ�����̼�¼

%% BP����Ԥ��
%���ݹ�һ��
inputn_test=mapminmax('apply',input_test,inputps);
an=sim(net,inputn_test);
test_simu=mapminmax('reverse',an,outputps);%%����������һ��
error=sum(abs(test_simu-output_test))
