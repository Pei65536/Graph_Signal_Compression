function [Bout] = VectorToMatrix(bin)
%��������ʽ��Ƶ����Ϣת��Ϊ������ʽ
%Bin - �����ͼ��Ƶ�����
%bout-�����Ƶ��������ʽ

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