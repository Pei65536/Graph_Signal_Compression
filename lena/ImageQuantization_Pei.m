function [QuantizationB] = ImageQuantization_Pei(Bin, QB, eta, Averagex, Varx)
%��ͼ�źŽ���ѹ������
%Bin - �����ͼ��Ƶ����� QB - ���ط������ eta - ���� Averagex - Ƶ��ľ�ֵ Varx - Ƶ��ķ���
%QuantizationB-�����Ƶ��������
N1 = size(Bin, 1);
N2 = size(Bin, 2);

for i = 1:N1
    for j = 1:N2
        QuantizationB(i,j) = UniformQuantization(Bin(i,j),QB(i,j),eta,Averagex(i,j),Varx(i,j));
    end
end
end