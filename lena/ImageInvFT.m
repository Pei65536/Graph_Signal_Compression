function [Bsamll] = ImageInvFT(BFsamll,InvIFT)
%��ͼ���һ��С����и���Ҷ�任�����������ʽ��Ƶ�������
%Bsamll - ������������IFT  - ͼ����Ҷ�任����
%BFsamll-����ľ�����ʽ��Ƶ�����

bFsmall = MatrixToVector(BFsamll);

bsmall = InvIFT * bFsmall;

Bsamll = VectorToMatrix(bsmall);

end