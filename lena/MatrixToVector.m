function [bout] = MatrixToVector(Bin)
%��������ʽ��Ƶ����Ϣת��Ϊ������ʽ
%Bin - �����ͼ��Ƶ�����
%bout-�����Ƶ��������ʽ

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