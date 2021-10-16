%------------------------------------------------------
% This code generates Figure 2 of the following paper: 
% "Bayesian Compressive Sensing" (Preprint, 2007)
% This example is modified from l1qc_example.m, an example 
% from l1magic.
% Coded by: Shihao Ji, ECE, Duke University
% last change: Jan. 2, 2007
%------------------------------------------------------
clc;
samples = 50 : 25 : 250;
total_count = 10;
L = length(samples);
t_BCS = zeros(total_count,L);
t_BP = zeros(total_count,L);
t_AMP = zeros(total_count,L);
t_CoSaMP = zeros(total_count,L);
t_MSBL = zeros(total_count,L);
for count = 1:total_count
    count
    rand('state', count);
    randn('state', 2*count);
	for i = 1 : L
		%limited graph signal 
		N = 512;    %numbers of nodes
		K = 20;    %numbers of no-zero of F field
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

% 		signal = A*x_f; noise = e;
% 		SNR_rec = 20*log10(norm(signal,'fro')/norm(noise,'fro'));
% 		% BCS algorithm 
% 		initsigma2 = std(y)^2/1e2;    % the initial variance of noise
% 		%initsigma2 = std(y)^2/1e6;    % the initial variance of non-noise
% 		tic;
% 		[weights,used,sigma2,errbars] = BCS_fast_rvm(A,y,initsigma2,1e-8);
% 		t_BCS(count,i) = toc;    % time starts
% 		fprintf(1,'BCS number of nonzero weights: %d\n',length(used));    % print the number of non-zero weights

% 		% BP algorithm 
% 		x0 = A'/(A*A')*y;
% 		% take epsilon a little bigger than sigma*sqrt(K)
% 		epsilon =  sigma*sqrt(M)*sqrt(1 + 2*sqrt(2)/sqrt(M));                                                                                                              
% 		tic;
% 		x_f_BP = l1qc_logbarrier(x0, A, [], y, epsilon, 1e-3);
% 		t_BP(count,i) = toc;
% 		fprintf(1,'BP number of nonzero weights: %d\n',sum(x_f_BP~=0));

		%CoSaMP algorithm
		tic;
		[s_hat, Supp_A] = CoSaMP(A, y, K, 300);
		x_f_CoSaMP = s_hat;
		t_CoSaMP(count,i) = toc;
		fprintf(1,'CoSaMP number of nonzero weights: %d\n',sum(s_hat~=0));

% 		%AMP algorithm
% 		tic;
% 		T = 1000; tol = 0.001;    %iterations and threshold
% 		x_f_AMP = reconstructAmp(A, y, T, tol);
% 		t_AMP(count,i) = toc;
% 		fprintf(1,'AMP number of nonzero weights: %d\n',sum(x_f_AMP~=0));

% 		%MSBL algorithm
% 		tic;
% 		[x_f_MSBL] = MSBL(A, y, sigma, 0);    
% 		t_MSBL(count,i) = toc;
% 		fprintf(1,'MSBL number of nonzero weights: %d\n',sum(x_f_MSBL~=0));
	end
end
% t_BCS1 = mean(t_BCS);
% t_AMP1 = mean(t_AMP);
% t_BP1 = mean(t_BP);
 t_CoSaMP1 = mean(t_CoSaMP);
% t_MSBL1 = mean(t_MSBL);