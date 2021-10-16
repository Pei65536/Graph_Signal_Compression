clc
clear

N = 4;
W = CreateW(N);
D=zeros(N,N);%�Խ���
D=diag(sum(W));%�Ⱦ���
L=D-W;%������˹����
[VL, DL] = eig(L); %ͼ����Ҷ�任

A=imread('D:\fig\lena.jpg');
I=rgb2gray(A);
I=im2double(I);%��ͼ��ת��Ϊ˫���ȸ�ʽ
%T=dctmtx(8);%����һ��8*8��DCT�任����

B=blkproc(I,[4 4],'ImageFT',VL');

Xin = B;
[AverageX, VarX] = AverVarComputation(Xin, 4);
eta = 3;

% [bout] = MatrixToVector(W);
% 
% [Bout] = VectorToMatrix(bout);
UF = VL;
N = 16;
Lda = zeros(N,1);
sigm0 = 0;
[GammaRec, Cx] = OperatorComputation(UF, VarX, sigm0);

HatGammaRec = GammaRec * Cx^(-1/2);

[s,v,d] = svd(HatGammaRec);

for i=1:N
    Lda(i) = v(i,i);  %ͼ������˹�����Ӧ������ֵ����
end
M = 2^64;
MV = BitAllocation(Lda,N,M,eta); %���ط�����
NumMV = sum(MV > 1);

Lnaix = 3*MV.^2/(2*eta^2);

lambdai = Lda(1:NumMV);

Lnai = Lnaix(1:NumMV);

[ Aplhai ] = AlhpaSolve(lambdai, Lnai);

[Usum,Tsum] = ComputeTk(Lnai,Aplhai);

AplhaiMatrix=zeros(NumMV, N);

for i = 1:NumMV
    AplhaiMatrix(i,i) = sqrt(Aplhai(i));
end


MSAM = Usum * AplhaiMatrix * d * Cx^(-1/2);
%MSAM = AplhaiMatrix * d * Cx^(-1/2);
g = zeros(NumMV,1);
EM = MSAM  * Cx  * MSAM';
for i =1:NumMV
    g(i) = 2*eta^2*EM(i,i)/(3*MV(i)^2);
end

G = diag(g);
%G = zeros(NumMV, NumMV);
% 
MREC = GammaRec * Cx * MSAM' * (MSAM  * Cx  * MSAM' + G)^(-1);

Idi = ones(512/4, 512/4);

Di = kron(Idi,AverageX);
IDi2=blkproc(Di,[4 4],'ImageInvFT',VL);


Itry = I - IDi2;

B1=blkproc(Itry,[4 4],'Sampling_Pei',MSAM);

QB = VectorToMatrix(MV);
B2=blkproc(B1,[4 4],'ImageQuantization_Pei',QB, eta, AverageX, VarX);
%��������




B3=blkproc(B1,[4 4],'Recovery_Pei',MREC);

I2=blkproc(B3,[4 4],'ImageInvFT',VL);

Iend = I2 + IDi2;
%���з��任���õ�ѹ�����ͼ��
imshow(I);  
title('original image')
figure;                            
imshow(Iend); 
title('compressed image')


%Aplhai  = AlhpaSolve(lambdai, Lnai); 