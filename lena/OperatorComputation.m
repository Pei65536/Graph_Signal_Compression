function [GammaRec, Cx] = OperatorComputation(UF, VarX, sigm0)
%�Բ�������ļ���
%UF - ͼ����Ҷ����N^2*N^2  VarX - Э�������-ʵ������������Э������о��󻯵Ľ��   
%sigm0 - ���Ը�˹�����Ĵ�С
% GammaRec - MMSE���ƾ���
% Cx - ʵ��ʹ�õ�Э������� N^2*N^2

N2  = length(UF); %ȡ�����ά��

N = sqrt(N2);

I  = eye(N2);

varx = MatrixToVector(VarX);

%varxcopy = sort(varx);

Cx = UF * diag(varx) * UF' + sigm0 * I;

GammaRec = diag(varx) * UF' * Cx^(-1);




end