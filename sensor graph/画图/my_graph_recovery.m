%------------------------------------------------------
% This code generates Figure 2 of the following paper: 
% "Bayesian Compressive Sensing" (Preprint, 2007)
% This example is modified from l1qc_example.m, an example 
% from l1magic.
% Coded by: Shihao Ji, ECE, Duke University
% last change: Jan. 2, 2007
%------------------------------------------------------
addpath(genpath(pwd));


clear all
close all 
clc

clear all
rand('state', 1);   % seed first
randn('state', 2);   % seed second(normal distribution)

%limited graph signal 
N = 512;    %numbers of nodes
K = 20;    %numbers of no-zero of F field
M = 150;
x_f = zeros(N,1);    % N*1 dimensional zero-vector
q = randperm(N);    % disrupt the order
x_f(q(1:K)) = sign(randn(K,1));    % the original signal(T non-zeros, + - 1)

param.distribute = 1;
G = gsp_sensor(N, param);
W = full(G.W);    %sparse->full for matrix
[V, D] = eig(W); 
x = V*x_f;    %inverse graph F

% Generate dictionary matrix
A = randn(M,N);    % initialize it to a normal distribution
A = A./repmat(sqrt(sum(A.^2,2)),[1,N]);	    % normalization
% Add noise
sigma = 0.005;
e = sigma*randn(M,1);
y = A*x_f + e;    % noise

signal = A*x_f; noise = e;
SNR_rec = 20*log10(norm(signal,'fro')/norm(noise,'fro'));


% BCS algorithm 
initsigma2 = std(y)^2/1e2;    % the initial variance of noise
%initsigma2 = std(y)^2/1e6;    % the initial variance of non-noise
tic;
[weights,used,sigma2,errbars] = BCS_fast_rvm(A,y,initsigma2,1e-8);
t_BCS = toc;    % time starts
x_f_BCS = zeros(N,1); err = zeros(N,1);
x_f_BCS(used) = weights; err(used) = errbars;    % assigned weights and error bound to corresponding positions
fprintf(1,'BCS number of nonzero weights: %d\n',length(used));    % print the number of non-zero weights

% BP algorithm 
x0 = A'*inv(A*A')*y;
% take epsilon a little bigger than sigma*sqrt(K)
epsilon =  sigma*sqrt(M)*sqrt(1 + 2*sqrt(2)/sqrt(M));                                                                                                              
tic;
x_f_BP = l1qc_logbarrier(x0, A, [], y, epsilon, 1e-3);
t_BP = toc;
fprintf(1,'BP number of nonzero weights: %d\n',sum(x_f_BP~=0));

% % OMP algorithm
% tic;
% PHI = A*V;
% %e1 = sigma*randn(N,1);
% y_omp = PHI*x_f+e;
% x_f_OMP = cs_omp(y_omp, PHI, K);
% t_OMP = toc;
% fprintf(1,'OMP number of nonzero weights: %d\n',sum(x_f_OMP~=0));

%CoSaMP algorithm
tic;
[s_hat, Supp_A] = CoSaMP(A, y, K, 300);
x_f_CoSaMP = s_hat;
t_CoSaMP = toc;
fprintf(1,'CoSaMP number of nonzero weights: %d\n',sum(s_hat~=0));

%AMP algorithm
tic;
T = 10000; tol = 0.00001;    %iterations and threshold
x_f_AMP = reconstructAmp(A, y, T, tol);
t_AMP = toc;
fprintf(1,'AMP number of nonzero weights: %d\n',sum(x_f_AMP~=0));

%MSBL algorithm
tic;
[x_f_MSBL] = MSBL(A, y, sigma, 0);    
t_MSBL = toc;
fprintf(1,'MSBL number of nonzero weights: %d\n',sum(x_f_MSBL~=0));

% Error ratio
E_BP = norm(x-V*x_f_BP)/norm(x);
E_BCS = norm(x-V*x_f_BCS)/norm(x);
% E_OMP = norm(x-V*x_f_OMP)/norm(x);
E_CoSaMP = norm(x-V*x_f_CoSaMP)/norm(x);
E_AMP = norm(x-V*x_f_AMP)/norm(x);
E_MSBL = norm(x-V*x_f_MSBL)/norm(x);

% Draw a picture corresponding to figure 2 in the paper
figure(1);
gsp_plot_signal(G, x);
title('the original signal')

figure(2);
gsp_plot_signal(G, V*x_f_BCS);
title('the recovery signal using BCS')

figure(3);
gsp_plot_signal(G, V*x_f_BP);
title('the recovery signal using BP')

% figure(4);
% gsp_plot_signal(G, V*x_f_OMP);
% title('the recovery signal using OMP')

figure(5);
gsp_plot_signal(G, V*x_f_CoSaMP);
title('the recovery signal using CoSaMP')

figure(6);
gsp_plot_signal(G, V*x_f_AMP);
title('the recovery signal using AMP')

figure(7);
gsp_plot_signal(G, V*x_f_MSBL);
title('the recovery signal using MSBL')

fprintf('\nSNR = %g dB  \n',SNR_rec);
disp(['BCS: ||I_hat-I||/||I|| = ' num2str(E_BCS) ', time = ' num2str(t_BCS) ' secs']); 
disp(['BP: ||I_hat-I||/||I|| = ' num2str(E_BP) ', time = ' num2str(t_BP) ' secs']);
% disp(['OMP: ||I_hat-I||/||I|| = ' num2str(E_OMP) ', time = ' num2str(t_OMP) ' secs']);
disp(['CoSaMP: ||I_hat-I||/||I|| = ' num2str(E_CoSaMP) ', time = ' num2str(t_CoSaMP) ' secs']);
disp(['AMP: ||I_hat-I||/||I|| = ' num2str(E_AMP) ', time = ' num2str(t_AMP) ' secs']);
disp(['MSBL: ||I_hat-I||/||I|| = ' num2str(E_MSBL) ', time = ' num2str(t_MSBL) ' secs']);
