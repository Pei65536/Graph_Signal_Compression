clc
clear
A=imread('D:\fig\lena.jpg');
N = 4;
I=rgb2gray(A);
I=im2double(I);%将图像转换为双精度格式
T=dctmtx(N);%返回一个8*8的DCT变换矩阵
B=blkproc(I,[N N],'P1*x*P2',T,T');%对原图像进行DCT变换
I8 = ones(N);
Xin = B;
[AverageX, VarX] = AverVarComputation(Xin, N);
eta = 3;
QB = ones(N,N) * 16;

B2=blkproc(B,[N N],'ImageQuantization_Pei',QB, eta, AverageX, VarX);
I2=blkproc(B2,[N N],'P1*x*P2',T',T);
%进行DCT反变换，得到压缩后的图像
imshow(I);   
figure;                            
imshow(I2); 







