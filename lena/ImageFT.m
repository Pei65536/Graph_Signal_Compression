function [BFsamll] = ImageFT(Bsamll,IFT)
%��ͼ���һ��С����и���Ҷ�任�����������ʽ��Ƶ�������
%Bsamll - ������������IFT  - ͼ����Ҷ�任����
%BFsamll-����ľ�����ʽ��Ƶ�����

bsmall = MatrixToVector(Bsamll);

bFsmall = IFT * bsmall;

BFsamll = VectorToMatrix(bFsmall);

end