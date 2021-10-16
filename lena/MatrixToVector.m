function [bout] = MatrixToVector(Bin)
%将矩阵形式的频域信息转换为向量形式
%Bin - 输入的图像频域矩阵
%bout-输出的频域向量形式

N = length(Bin);
bout = zeros(N^2,1);
for i = 1:N
    for j =1:N
        if mod(i,2) == 1
            bout((i-1)*N+j) = Bin(i,j);
        else
            bout((i-1)*N+N-j+1) = Bin(i,j);
        end
    end
end

end