function [xfres] = Sam_Res(W,x,Rvector)
N = 100;   %节点数
K = 20;    %频域非零值个数
M=K;
x_f = zeros(M,1);    %信号的频域分量初始化
SamChen=chen( N, W,K,  M );
SamChen=sort(SamChen);
D=zeros(N,N);%对角阵
D=diag(sum(W));%度矩阵
L=D-W;%拉普拉斯矩阵；
[VL, DL] = eig(L); %图傅里叶变换
UVR=VL(:,[1:K]);
sigma0=0.1;
hatLambda=zeros(M,M);
I=zeros(N,N);
for i=1:N
    I(i,i)=1;
end
for i=2:M
    hatLambda(i,i)=1/DL(i,i);
end
hatLambda(1,1)=35;
USVRI=UVR(SamChen,:)^(-1);
colA=zeros(M,1);
for i=1:M
    colA(i)=norm(USVRI(:,i));
end
SamMtrix=zeros(M,N);
for i=1:length(SamChen)
    SamMtrix(i,SamChen(i))=1;
end
E=UVR*hatLambda*UVR';
Imax=zeros(K,1);
eta=3;
for i=1:K
    Imax(i)=2*eta*sqrt(E(SamChen(i),SamChen(i)));
end
%Rvector=[20 40 60 80 100 120 140];
RT=Rvector;
Delt_S=DeltS_Compute( N, W, Imax, M ,RT,SamChen);
Delt_S=Delt_S';
Ri=zeros(M,1);
Mi=zeros(M,1);
for i=1:M
    Mi(i)=Imax(i)/Delt_S(i);
    Ri(i)=log2(Mi(i));
end
Rie=round(Ri);
Rie=ones(M,1)*RT/M;
for p=1:M
    x_f(p)=normrnd(0,sqrt(hatLambda(p,p)),[1 1]);
end
e0=normrnd(0,sigma0,[N 1]);%加性白噪声向量
xnoise=x+e0;
xs=SamMtrix*xnoise;
xsQ=zeros(M,1);
for xi=1:M
    [C,xsQ(xi)]=UniformQ(xs(xi),sqrt(E(SamChen(xi),SamChen(xi))),2^Rie(xi),eta);
end
xres=UVR*USVRI*xsQ;
xfres=UVR'*xres;
end