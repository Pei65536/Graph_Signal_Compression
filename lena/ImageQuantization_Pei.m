function [QuantizationB] = ImageQuantization_Pei(Bin, QB, eta, Averagex, Varx)
%对图信号进行压缩采样
%Bin - 输入的图像频域矩阵 QB - 比特分配矩阵 eta - 参数 Averagex - 频域的均值 Varx - 频域的方差
%QuantizationB-输出的频域结果矩阵
N1 = size(Bin, 1);
N2 = size(Bin, 2);

for i = 1:N1
    for j = 1:N2
        QuantizationB(i,j) = UniformQuantization(Bin(i,j),QB(i,j),eta,Averagex(i,j),Varx(i,j));
    end
end
end