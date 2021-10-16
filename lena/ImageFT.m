function [BFsamll] = ImageFT(Bsamll,IFT)
%对图像的一个小块进行傅里叶变换，输出矩阵形式的频域分量。
%Bsamll - 输入的区块矩阵，IFT  - 图傅里叶变换矩阵
%BFsamll-输出的矩阵形式的频域分量

bsmall = MatrixToVector(Bsamll);

bFsmall = IFT * bsmall;

BFsamll = VectorToMatrix(bFsmall);

end