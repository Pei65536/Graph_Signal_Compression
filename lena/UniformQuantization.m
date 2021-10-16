function [xQ] = UniformQuantization(x,M,eta,Averagex,Varx)
%对输入信号进行均匀量化处理
%x-输入信号 Varx-信号的方差 Averagex-信号的均值 M-量化状态数目 eta-参数
%xQ-量化后的信号
%
Range_Low = Averagex - eta*Varx;
%Range_Upper = Averagex + eta*Varx;

%GammaMax=eta*sign;%最大值范围
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

