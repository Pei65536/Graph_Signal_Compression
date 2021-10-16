clc
clear

N = 8;
W = CreateW(N);
D=zeros(N,N);%�Խ���
D=diag(sum(W));%�Ⱦ���
L=D-W;%������˹����
[VL, DL] = eig(L); %ͼ����Ҷ�任

A=imread('D:\fig\lena.jpg');
I=rgb2gray(A);
I=im2double(I);%��ͼ��ת��Ϊ˫���ȸ�ʽ
%T=dctmtx(8);%����һ��8*8��DCT�任����

B=blkproc(I,[N N],'ImageFT',VL');

Xin = B;
[AverageX, VarX] = AverVarComputation(Xin, N);
eta = 3;
QB = ones(N,N) * 16;

B2=blkproc(B,[N N],'ImageQuantization_Pei',QB, eta, AverageX, VarX);
%��������


I2=blkproc(B2,[N N],'ImageInvFT',VL);
%���з��任���õ�ѹ�����ͼ��
imshow(I);   
figure;                            
imshow(I2); 

PSN = norm(I-I2);
