function [RecB] = Recovery_Pei(Bin, MREC)
%��ͼ�źŽ���ѹ������
%U-�������� V-�ԽǾ��� M-���������� K-���� eta-���� sigm0 - ��˹������С
%SamOperator-��������
%
bin =  MatrixToVector(Bin);

M = size(MREC,2);

bin = bin(1:M);


Recb = MREC * bin;

RecB = VectorToMatrix(Recb);

end