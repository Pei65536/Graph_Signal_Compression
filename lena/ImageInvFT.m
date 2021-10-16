function [Bsamll] = ImageInvFT(BFsamll,InvIFT)
%对图像的一个小块进行傅里叶变换，输出矩阵形式的频域分量。
%Bsamll - 输入的区块矩阵，IFT  - 图傅里叶变换矩阵
%BFsamll-输出的矩阵形式的频域分量

bFsmall = MatrixToVector(BFsamll);

bsmall = InvIFT * bFsmall;

Bsamll = VectorToMatrix(bsmall);

end