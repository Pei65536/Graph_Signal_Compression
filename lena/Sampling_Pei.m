function [SamedB] = Sampling_Pei(Bin, MSAM)
%��ͼ�źŽ���ѹ������
%U-�������� V-�ԽǾ��� M-���������� K-���� eta-���� sigm0 - ��˹������С
%SamOperator-��������
%
bin =  MatrixToVector(Bin);

M = size(MSAM,1);
N = size(MSAM,2);

Id = zeros(N-M, 1);

Samedb = MSAM * bin;

Samedb = [Samedb', Id']';

SamedB = VectorToMatrix(Samedb);

end