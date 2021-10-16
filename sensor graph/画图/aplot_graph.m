clear all
rand('state', 1);   % seed first
randn('state', 2);   % seed second(normal distribution)

%limited graph signal 
N = 512;    %numbers of nodes
K = 20;    %numbers of no-zero of F field
x_f = zeros(N,1);    % N*1 dimensional zero-vector
q = randperm(N);    % disrupt the order
x_f(q(1:K)) = sign(randn(K,1));    % the original signal(T non-zeros, + - 1)
%x_f(1:20)=[1,2,3,4,5,-5,-4,-3,-2,-1,6,7,8,9,10,-6,-7,-8,-9,-10];
%x_f(1:20)=[1,2,3,4,5,5,4,3,2,1,1,2,3,4,5,5,4,3,2,1];
param.distribute = 1;
G1 = gsp_sensor(N, param);
W1 = full(G1.W);    %sparse->full for matrix
[V1, D1] = eig(W1); 
x1 = V1*x_f;    %inverse graph F

[G2, tt]=gsp_swiss_roll(N)
W2 = full(G2.W);    %sparse->full for matrix
[V2, D2] = eig(W2); 
x2 = V2*x_f;    %inverse graph F

G3 = gsp_community(N);
W3 = full(G3.W);    %sparse->full for matrix
[V3, D3] = eig(W3); 
x3 = V3*x_f;    %inverse graph F
% Draw a picture corresponding to figure 2 in the paper
figure(1);
gsp_plot_signal(G1, x1);
%title('the original signal')

figure(2);
gsp_plot_signal(G2, x2);
title('the recovery signal using BCS')

figure(3);
gsp_plot_signal(G3, x3);
title('the recovery signal using BP')
% 
% figure(4);
% gsp_plot_signal(G, V*x_f_OMP);
% title('the recovery signal using OMP')
% 
% figure(5);
% gsp_plot_signal(G, V*x_f_CoSaMP);
% title('the recovery signal using CoSaMP')
% 
% figure(6);
% gsp_plot_signal(G, V*x_f_AMP);
% title('the recovery signal using AMP')
% 
% figure(7);
% gsp_plot_signal(G, V*x_f_MSBL);
% title('the recovery signal using MSBL')