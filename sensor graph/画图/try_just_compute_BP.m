% just compute BP
clc;
clear all
total_count = 32;
samples = [41 : 1 : 150];
L = length(samples);
E_BP = zeros(total_count,L);
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

		% BP algorithm 
		x0 = A'/(A*A')*y;
		% take epsilon a little bigger than sigma*sqrt(K)
		epsilon =  sigma*sqrt(M)*sqrt(1 + 2*sqrt(2)/sqrt(M));                                                                                                              
		x_f_BP = l1qc_logbarrier(x0, A, [], y, epsilon, 1e-3);

		% Error ratio
		E_BP(count,i) = norm(x-V*x_f_BP)/norm(x);
	end
end
save BP_results.mat E_BP;
disp('Done!');
