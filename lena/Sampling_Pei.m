function [SamedB] = Sampling_Pei(Bin, MSAM)
%对图信号进行压缩采样
%U-正交矩阵 V-对角矩阵 M-量化器个数 K-带宽 eta-参数 sigm0 - 高斯噪声大小
%SamOperator-采样矩阵
%
bin =  MatrixToVector(Bin);

M = size(MSAM,1);
N = size(MSAM,2);

Id = zeros(N-M, 1);

Samedb = MSAM * bin;

Samedb = [Samedb', Id']';

SamedB = VectorToMatrix(Samedb);

end