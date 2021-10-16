clc
clear

N = 8;
W = CreateW(N);
D=zeros(N,N);%对角阵
D=diag(sum(W));%度矩阵
L=D-W;%拉普拉斯矩阵；
[VL, DL] = eig(L); %图傅里叶变换

A=imread('D:\fig\lena.jpg');
I=rgb2gray(A);
I=im2double(I);%将图像转换为双精度格式
%T=dctmtx(8);%返回一个8*8的DCT变换矩阵

B=blkproc(I,[N N],'ImageFT',VL');

Xin = B;
[AverageX, VarX] = AverVarComputation(Xin, N);
eta = 3;
QB = ones(N,N) * 16;

B2=blkproc(B,[N N],'ImageQuantization_Pei',QB, eta, AverageX, VarX);
%数据量化


I2=blkproc(B2,[N N],'ImageInvFT',VL);
%进行反变换，得到压缩后的图像
imshow(I);   
figure;                            
imshow(I2); 

PSN = norm(I-I2);
