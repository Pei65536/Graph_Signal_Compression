% just compute MSBL
clc;
total_count = 32;
samples = 41 : 1 : 140;
L = length(samples);
E_MSBL = zeros(total_count,L);
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

		%MSBL algorithm
		[x_f_MSBL] = MSBL(A, y, sigma, 0);    

		% Error ratio
		E_MSBL(count,i) = norm(x-V*x_f_MSBL)/norm(x);
	end
end
save MSBL_results.mat E_MSBL;
disp('Done!');