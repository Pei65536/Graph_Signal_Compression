function [xres] = HDL_Res(W,x,Rvector)
%rand('state', 1);   % seed first  随机数生成模式
%randn('state', 2);   % seed second(normal distribution)
%limited graph signal 
N = 100;   %节点数
K = 20;    %频域非零值个数
M=K;
x_f = zeros(M,1);    %信号的频域分量初始化
%load W;%导入图的邻接矩阵
D=zeros(N,N);%对角阵
D=diag(sum(W));%度矩阵
L=D-W;%拉普拉斯矩阵；
[VL, DL] = eig(L); %图傅里叶变换
UVR=VL(:,[1:K]);%前K列傅里叶变换矩阵备用
sigma0=0.1;%高斯白噪声的方差
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
R=Rvector;
Mi1=ones(M,1)*2^floor(R/M); %
[PsiX, s_fDynRange, m_fB] = v_fGetHLQuant(E, GammaX, 2^R);
MX=PsiX*E*PsiX';%预处理后图信号的协方差矩阵
eta=3;%均匀量化下高斯信号的范围参数
for i=1:M
    G1(i,i)=(4*eta^2*MX(i,i))/Mi1(i)^2/12;%量化误差的协方差矩阵定义(均匀量化+最优分配比特)
end
PhiX1=GammaX*E*PsiX'*(PsiX*E*PsiX'+G1)^(-1);%均匀量化下的恢复矩阵
e0=normrnd(0,sigma0,[N 1]);%加性白噪声向量
xnoise = x + e0;
PX=PsiX*xnoise;
PXQ2=zeros(M,1);
for xi=1:M
%    [C,PXQ(xi)]=Non_uniformQ(PX(xi),sqrt(MX(xi,xi)),Mi1(xi));%非均匀量化
    [C2, PXQ2(xi)]=UniformQ(PX(xi),sqrt(MX(xi,xi)),Mi1(xi),3);%均匀量化
end
xres = PhiX1*PXQ2;
end


