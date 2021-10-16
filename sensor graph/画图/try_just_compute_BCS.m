% just compute BCS
clc;
clear all
total_count = 100;
samples = [71 : 1 : 150];
L = length(samples);
E_BCS = zeros(total_count,L);
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

		% BCS algorithm 
		initsigma2 = std(y)^2/1e2;    % the initial variance of noise
		[weights,used,sigma2,errbars] = BCS_fast_rvm(A,y,initsigma2,1e-8);
		x_f_BCS = zeros(N,1); err = zeros(N,1);
		x_f_BCS(used) = weights; err(used) = errbars;    % assigned weights and error bound to corresponding positions

		% Error ratio
		E_BCS(count,i) = norm(x-V*x_f_BCS)/norm(x);
	end
end
save BCS_results.mat E_BCS;
disp('Done!');