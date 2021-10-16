function [ sample ] = chen( N, W,K,  M )
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here
% N为节点数
%W为邻接矩阵
%K为带宽
%M为采样点数
sample=[];%初始采样集
D=zeros(N,N);%对角阵初始化
D=diag(sum(W));%度矩阵
L=D-W;%拉普拉斯矩阵；
[VL, DL] = eig(L); %图傅里叶变换
VLR=VL(:,[1:K]);%取前K个向量组成矩阵
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

