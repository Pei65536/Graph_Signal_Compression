clc
clear
A=imread('D:\fig\lena.jpg');
N = 4;
I=rgb2gray(A);
I=im2double(I);%��ͼ��ת��Ϊ˫���ȸ�ʽ
T=dctmtx(N);%����һ��8*8��DCT�任����
B=blkproc(I,[N N],'P1*x*P2',T,T');%��ԭͼ�����DCT�任
I8 = ones(N);
Xin = B;
[AverageX, VarX] = AverVarComputation(Xin, N);
eta = 3;
QB = ones(N,N) * 16;

B2=blkproc(B,[N N],'ImageQuantization_Pei',QB, eta, AverageX, VarX);
I2=blkproc(B2,[N N],'P1*x*P2',T',T);
%����DCT���任���õ�ѹ�����ͼ��
imshow(I);   
figure;                            
imshow(I2); 







