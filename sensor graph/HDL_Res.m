function [xres] = HDL_Res(W,x,Rvector)
%rand('state', 1);   % seed first  ���������ģʽ
%randn('state', 2);   % seed second(normal distribution)
%limited graph signal 
N = 100;   %�ڵ���
K = 20;    %Ƶ�����ֵ����
M=K;
x_f = zeros(M,1);    %�źŵ�Ƶ�������ʼ��
%load W;%����ͼ���ڽӾ���
D=zeros(N,N);%�Խ���
D=diag(sum(W));%�Ⱦ���
L=D-W;%������˹����
[VL, DL] = eig(L); %ͼ����Ҷ�任
UVR=VL(:,[1:K]);%ǰK�и���Ҷ�任������
sigma0=0.1;%��˹�������ķ���
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
R=Rvector;
Mi1=ones(M,1)*2^floor(R/M); %
[PsiX, s_fDynRange, m_fB] = v_fGetHLQuant(E, GammaX, 2^R);
MX=PsiX*E*PsiX';%Ԥ�����ͼ�źŵ�Э�������
eta=3;%���������¸�˹�źŵķ�Χ����
for i=1:M
    G1(i,i)=(4*eta^2*MX(i,i))/Mi1(i)^2/12;%��������Э���������(��������+���ŷ������)
end
PhiX1=GammaX*E*PsiX'*(PsiX*E*PsiX'+G1)^(-1);%���������µĻָ�����
e0=normrnd(0,sigma0,[N 1]);%���԰���������
xnoise = x + e0;
PX=PsiX*xnoise;
PXQ2=zeros(M,1);
for xi=1:M
%    [C,PXQ(xi)]=Non_uniformQ(PX(xi),sqrt(MX(xi,xi)),Mi1(xi));%�Ǿ�������
    [C2, PXQ2(xi)]=UniformQ(PX(xi),sqrt(MX(xi,xi)),Mi1(xi),3);%��������
end
xres = PhiX1*PXQ2;
end


