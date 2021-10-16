function [GammaRec, Cx] = OperatorComputation(UF, VarX, sigm0)
%对参数矩阵的计算
%UF - 图傅里叶基，N^2*N^2  VarX - 协方差矩阵-实际上是向量的协方差进行矩阵化的结果   
%sigm0 - 加性高斯噪声的大小
% GammaRec - MMSE估计矩阵
% Cx - 实际使用的协方差矩阵， N^2*N^2

N2  = length(UF); %取矩阵的维度

N = sqrt(N2);

I  = eye(N2);

varx = MatrixToVector(VarX);

%varxcopy = sort(varx);

Cx = UF * diag(varx) * UF' + sigm0 * I;

GammaRec = diag(varx) * UF' * Cx^(-1);




end