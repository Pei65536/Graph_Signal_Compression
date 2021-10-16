function [xQ] = UniformQuantization(x,M,eta,Averagex,Varx)
%�������źŽ��о�����������
%x-�����ź� Varx-�źŵķ��� Averagex-�źŵľ�ֵ M-����״̬��Ŀ eta-����
%xQ-��������ź�
%
Range_Low = Averagex - eta*Varx;
%Range_Upper = Averagex + eta*Varx;

%GammaMax=eta*sign;%���ֵ��Χ
enn=[];
%enn(1)=-1*GammaMax;
for i=1:M+1
    enn(i)=Range_Low+(i-1)*2*eta*Varx/M;
end
C = ones(M,1);
for i=1:M
    C(i)=(enn(i)+enn(i+1))/2;
end
xres=C-x;
[min_ax,indexx]=min(abs(xres));
xQ=C(indexx);


end

