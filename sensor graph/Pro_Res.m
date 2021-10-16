function [xres] = Pro_Res(W, x, Rvector)
%UNTITLED5 此处显示有关此函数的摘要
%   此处显示详细说明 待修改
N = length(W);   %节点数
K = 20;    %频域非零值个数
M = K;
x_f = zeros(M,1);    %信号的频域分量初始化
%load W;%导入图的邻接矩阵
D=zeros(N,N);%对角阵
D=diag(sum(W));%度矩阵
L=D-W;%拉普拉斯矩阵；
[VL, DL] = eig(L); %图傅里叶变换
UVR=VL(:,[1:K]);%前K列傅里叶变换矩阵备用
sigma0 = 0.1;
hatLambda=zeros(M,M);
I=zeros(N,N);%单位矩阵初始化
for i=1:N
    I(i,i)=1;%N维单位矩阵定义
end
for i=2:M
    hatLambda(i,i)=1/DL(i,i);
end
hatLambda(1,1)=35;%频域分量的假定
E=UVR*hatLambda*UVR'+sigma0*I;%图信号的协方差矩阵
GammaX=hatLambda*UVR'*E^(-1);%最佳线性恢复
[U A V]=svd(GammaX*E^(0.5));%对tilde Gamma矩阵的svd分解
Mi1=zeros(M,1);%初始化比特分配
Lda=zeros(M,1);%初始化tilde Gamma矩阵的特征值组成的对角阵
for i=1:M
    Lda(i)=A(i,i);%tilde Gamma矩阵的特征值组成的对角阵
end
R=Rvector;
for i=1:M
    Mi2(i)=2^floor(R/M);%平均分配比特
end
eta=3;%均匀量化下高斯信号的范围参数
Mi1=BitAllocation(Lda,K,2^R,eta);%调用比特分配函数
p0=length(Mi1(Mi1>1));%获取预处理矩阵的维度
G1=zeros(p0,p0);%初始化量化误差的协方差矩阵
AI=zeros(p0,N);%随便定义一个对角非方阵且对角元素不为0
for i=1:p0
    AI(i,i)=1;
end
PsiX=AI*real(V'*E^(0.5));%获取最优的预处理矩阵
MX=PsiX*E*PsiX';%预处理后图信号的协方差矩阵
LamndaF=zeros(1,K);
for i=1:K
    LamndaF(i)=hatLambda(i,i)+sigma0;
    %LamndaF(i)=1;
end
avector = zeros(1,K);
Msum=2^R;
Mi1=zeros(1,K);
for i = 1:K
    avector(i) = 2*eta^2*MX(i,i)/3;
end
alphaM=prod(avector)^(1/K)/Msum^(2/K);
for i =1:K
    Mi1(i)=sqrt(avector(i)/alphaM);
end
for i=1:p0
    G1(i,i)=(4*eta^2*MX(i,i))/Mi1(i)^2/12;%量化误差的协方差矩阵定义(均匀量化+最优分配比特)
end
for i=1:p0
    G2(i,i)=(pi*sqrt(3)*MX(i,i))/Mi1(i)^2/2;%量化误差的协方差矩阵定义（非均匀量化）
end
PhiX1=GammaX*E*PsiX'*(PsiX*E*PsiX'+G1)^(-1);%均匀量化下的恢复矩阵
%PhiX2=GammaX*E*PsiX'*(PsiX*E*PsiX'+G2)^(-1);%非均匀量化下的恢复矩阵
e0=normrnd(0,sigma0,[N 1]);%加性白噪声向量
xnoise = x + e0;
PX=PsiX*xnoise;
PXQ=zeros(p0,1);%预处理后的信号进行量化
PXQ2=zeros(p0,1);
for xi=1:p0
%    [C,PXQ(xi)]=Non_uniformQ(PX(xi),sqrt(MX(xi,xi)),Mi1(xi));%非均匀量化
    [C2, PXQ2(xi)]=UniformQ(PX(xi),sqrt(MX(xi,xi)),Mi1(xi),3);%均匀量化
end
xres = PhiX1*PXQ2;
end

