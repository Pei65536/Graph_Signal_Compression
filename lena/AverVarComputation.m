function [AverageX, VarX] = AverVarComputation(Xin, NumK)
%对输入的图像信号计算均值与方差
%Xin输入的图像信号（频域）
%AverageX - 均值  VarX - 方差
%

%NumK = 8;

N = length(Xin);
M = N/NumK;
XQ1 = zeros(NumK,NumK,M,M);

for i=1:M
    for j =1:M
        XQ1(:,:,i,j) = Xin(NumK*(i-1)+1:NumK*i,NumK*(j-1)+1:NumK*j);  %将信号以4-D矩阵的形式留存
    end
end

AverageX = zeros(NumK,NumK);
VarX = zeros(NumK,NumK);
midX = zeros(M,M);

for i=1:NumK
    for j=1:NumK
        for m =1:M
            for n = 1:M
                midX(m,n) = XQ1(i,j,m,n);
            end
        end
        AverageX(i,j) = mean(mean(midX));
        minX2 = reshape(midX,1,M^2);
        VarX(i,j) = var(minX2);
    end
end



end