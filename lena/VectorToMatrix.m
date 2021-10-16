function [Bout] = VectorToMatrix(bin)
%将向量形式的频域信息转换为矩阵形式
%Bin - 输入的图像频域矩阵
%bout-输出的频域向量形式

Nx = length(bin);
Bout = zeros(sqrt(Nx),sqrt(Nx));
for i = 1:sqrt(Nx)
    for j =1:sqrt(Nx)
        if mod(i,2) == 1
            Bout(i,j) = bin((i-1)*sqrt(Nx)+j);
        else
            Bout(i,j) = bin((i-1)*sqrt(Nx)+sqrt(Nx)-j+1);
        end
    end
end

end