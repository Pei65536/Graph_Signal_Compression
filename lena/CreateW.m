function [W] = CreateW(N)
%�ڽӾ��������
%N - ����Ľڵ������һ����N^2
%W-����ľ���ΪN^2*N^2��С
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