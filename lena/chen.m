function [ sample ] = chen( N, W,K,  M )
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here
% NΪ�ڵ���
%WΪ�ڽӾ���
%KΪ����
%MΪ��������
sample=[];%��ʼ������
D=zeros(N,N);%�Խ����ʼ��
D=diag(sum(W));%�Ⱦ���
L=D-W;%������˹����
[VL, DL] = eig(L); %ͼ����Ҷ�任
VLR=VL(:,[1:K]);%ȡǰK��������ɾ���
for i=1:M
    for j=1:N
        if ismember(j,sample)==1
            sigema(j)=-10000;
        else
            sample(i)=j;
            FUK=VLR(sample,:);
            AAH=FUK*FUK';
            [VAAH, DAAH] = eig(AAH); 
            for m=1:length(sample)
                DAv(m)=DAAH(m,m);
            end
            [DAvm,DAvn]=min(DAv);
            sigema(j)=DAvm;
        end
            [fum,fun]=max(sigema);
    end
    sample(i)=fun;
    i=i+1;   
end
end

