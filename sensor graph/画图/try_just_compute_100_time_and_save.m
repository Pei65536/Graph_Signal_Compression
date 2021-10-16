%------------------------------------------------------
% This code generates Figure 2 of the following paper: 
% "Bayesian Compressive Sensing" (Preprint, 2007)
% This example is modified from l1qc_example.m, an example 
% from l1magic.
% Coded by: Shihao Ji, ECE, Duke University
% last change: Jan. 2, 2007
%------------------------------------------------------
clc;
clear all
total_count = 100;
samples = [71 : 1 : 150];
L = length(samples);
for count = 1:total_count
    count
    rand('state', count);
    randn('state', 2*count);
	for i = 1 : L
		%limited graph signal 
		N = 512;    %numbers of nodes
		K = 20;    %numbers of no-zero of F field
	%	M = 150;
		M = samples(i);
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
		[weights,used,sigma2,errbars] = BCS_fast_rvm(A,y,initsigma2,1e-8);
		x_f_BCS = zeros(N,1); err = zeros(N,1);
		x_f_BCS(used) = weights; err(used) = errbars;    % assigned weights and error bound to corresponding positions

		% BP algorithm 
		x0 = A'*inv(A*A')*y;
		% take epsilon a little bigger than sigma*sqrt(K)
		epsilon =  sigma*sqrt(M)*sqrt(1 + 2*sqrt(2)/sqrt(M));                                                                                                              
		x_f_BP = l1qc_logbarrier(x0, A, [], y, epsilon, 1e-3);

		% OMP algorithm
		PHI = A*V;
		y_omp = PHI*x_f+e;
		x_f_OMP = cs_omp(y_omp, PHI, K);
		
		%CoSaMP algorithm
		[s_hat, Supp_A] = CoSaMP(A, y, K, 300);
		x_f_CoSaMP = s_hat;
		
		%AMP algorithm
		T = 1000; tol = 0.001;    %iterations and threshold
		x_f_AMP = reconstructAmp(A, y, T, tol);

		%MSBL algorithm
		[x_f_MSBL] = MSBL(A, y, sigma, 0);    

		% Error ratio
		E_BP(count,i) = norm(x-V*x_f_BP)/norm(x);
		E_BCS(count,i) = norm(x-V*x_f_BCS)/norm(x);
		E_OMP(count,i) = norm(x-V*x_f_OMP)/norm(x);
		E_CoSaMP(count,i) = norm(x-V*x_f_CoSaMP)/norm(x);
		E_AMP(count,i) = norm(x-V*x_f_AMP)/norm(x);
		E_MSBL(count,i) = norm(x-V*x_f_MSBL)/norm(x);
	end
end
save BCS_results.mat E_BCS;
save BP_results.mat E_BP;
save CoSaMP_results.mat E_CoSaMP;
save AMP_results.mat E_AMP;
save MSBL_results.mat E_MSBL;
save OMP_results.mat E_OMP;