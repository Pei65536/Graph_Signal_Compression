function [xres] = Pro_Res(W, x, Rvector)
%UNTITLED5 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵�� ���޸�
N = length(W);   %�ڵ���
K = 20;    %Ƶ�����ֵ����
M = K;
x_f = zeros(M,1);    %�źŵ�Ƶ�������ʼ��
%load W;%����ͼ���ڽӾ���
D=zeros(N,N);%�Խ���
D=diag(sum(W));%�Ⱦ���
L=D-W;%������˹����
[VL, DL] = eig(L); %ͼ����Ҷ�任
UVR=VL(:,[1:K]);%ǰK�и���Ҷ�任������
sigma0 = 0.1;
hatLambda=zeros(M,M);
I=zeros(N,N);%��λ�����ʼ��
for i=1:N
    I(i,i)=1;%Nά��λ������
end
for i=2:M
    hatLambda(i,i)=1/DL(i,i);
end
hatLambda(1,1)=35;%Ƶ������ļٶ�
E=UVR*hatLambda*UVR'+sigma0*I;%ͼ�źŵ�Э�������
GammaX=hatLambda*UVR'*E^(-1);%������Իָ�
[U A V]=svd(GammaX*E^(0.5));%��tilde Gamma�����svd�ֽ�
Mi1=zeros(M,1);%��ʼ�����ط���
Lda=zeros(M,1);%��ʼ��tilde Gamma���������ֵ��ɵĶԽ���
for i=1:M
    Lda(i)=A(i,i);%tilde Gamma���������ֵ��ɵĶԽ���
end
R=Rvector;
for i=1:M
    Mi2(i)=2^floor(R/M);%ƽ���������
end
eta=3;%���������¸�˹�źŵķ�Χ����
Mi1=BitAllocation(Lda,K,2^R,eta);%���ñ��ط��亯��
p0=length(Mi1(Mi1>1));%��ȡԤ��������ά��
G1=zeros(p0,p0);%��ʼ����������Э�������
AI=zeros(p0,N);%��㶨��һ���ԽǷǷ����ҶԽ�Ԫ�ز�Ϊ0
for i=1:p0
    AI(i,i)=1;
end
PsiX=AI*real(V'*E^(0.5));%��ȡ���ŵ�Ԥ�������
MX=PsiX*E*PsiX';%Ԥ�����ͼ�źŵ�Э�������
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
    G1(i,i)=(4*eta^2*MX(i,i))/Mi1(i)^2/12;%��������Э���������(��������+���ŷ������)
end
for i=1:p0
    G2(i,i)=(pi*sqrt(3)*MX(i,i))/Mi1(i)^2/2;%��������Э��������壨�Ǿ���������
end
PhiX1=GammaX*E*PsiX'*(PsiX*E*PsiX'+G1)^(-1);%���������µĻָ�����
%PhiX2=GammaX*E*PsiX'*(PsiX*E*PsiX'+G2)^(-1);%�Ǿ��������µĻָ�����
e0=normrnd(0,sigma0,[N 1]);%���԰���������
xnoise = x + e0;
PX=PsiX*xnoise;
PXQ=zeros(p0,1);%Ԥ�������źŽ�������
PXQ2=zeros(p0,1);
for xi=1:p0
%    [C,PXQ(xi)]=Non_uniformQ(PX(xi),sqrt(MX(xi,xi)),Mi1(xi));%�Ǿ�������
    [C2, PXQ2(xi)]=UniformQ(PX(xi),sqrt(MX(xi,xi)),Mi1(xi),3);%��������
end
xres = PhiX1*PXQ2;
end

