function [W] = CreateW(N)
%邻接矩阵的生成
%N - 输入的节点个数，一般是N^2
%W-输出的矩阵，为N^2*N^2大小
A = zeros(N^2,N^2);
for i = 1:N^2
    for j = 1:N^2
        if abs(j - i) == 1
            A(i,j) = 1;
        end
        if i+j == ceil(i/N)*2*N+1
            A(i,j) = 1 ;
        end
    end
end

for i = 1:N^2
    for j = 1:i
        A(i,j) = A(j,i);
    end
end
W = A;

end