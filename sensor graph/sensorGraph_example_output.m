addpath(genpath(pwd));%引入图函数
clear all
close all 
clc

N = 100;    %numbers of nodes
K = 20;    %numbers of no-zero of F field
M = 20;
sigma0 = 0.1;
x_f = zeros(N,1);    % N*1 dimensional zero-vector
q = randperm(N);    % disrupt the order
x_f(q(1:K)) = sign(randn(K,1));    % the original signal(T non-zeros, + - 1)
param.distribute = 1;
G = gsp_sensor(N, param);
W = full(G.W);    %sparse->full for matrix
D=zeros(N,N);%对角阵
D=diag(sum(W));%度矩阵
L=D-W;%拉普拉斯矩阵；
[VL, DL] = eig(L); %图傅里叶变换
UVR=VL(:,[1:K]);%前K列傅里叶变换矩阵备用
Rvector = 60;
hatLambda=zeros(M,M);
I=zeros(N,N);%单位矩阵初始化
for i=1:N
    I(i,i)=1;%N维单位矩阵定义
end
for i=2:M
    hatLambda(i,i)=1/DL(i,i);
end
hatLambda(1,1)=35;%频域分量的假定
x_fp = zeros(M,1);
    for p=1:M
        x_fp(p)=normrnd(0,sqrt(hatLambda(p,p)),[1 1]);%定义信号的频域参数
    end
x_try = UVR*x_fp;

[xres_pro] = Pro_Res(W, x_try, Rvector);
[xres_hdl] = HDL_Res(W, x_try, Rvector);
[xres_sam] = Sam_Res(W, x_try, Rvector);

figure(1);
gsp_plot_signal(G, UVR*x_fp);
title('Original Graph Signal');

figure(2);
gsp_plot_signal(G, UVR*xres_pro);
title('Joint sampling and quantization, non-identical quantizers (Algorithm 1)');

figure(3);
gsp_plot_signal(G, UVR*xres_hdl);
title('Joint sampling and quantization, identical quantizers');

figure(4);
gsp_plot_signal(G, UVR*xres_sam);
title('Separate sampling and quantization, non-identical quantizers');


