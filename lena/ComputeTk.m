function [Usum,Tsum] = ComputeTk(x,y)
%已完成  输入：x (<m)  y  输出  Usum为正交矩阵满足 Usum'*diag(y)*Usum 结果的对角元素为x
%Tsum满足Tsum'*y=x
k = 1;
N=length(y);
I = eye(N);
Tsum = I;
Usum = I;
while(k<1000)
Mmatrix = eye(N,N);
k = 1;
U = zeros(N,N);
largey = [];
largex = [];
ix1 = 0;
iy1 = 0;
for i1 = 1:N
    if y(i1)>x(i1)
        iy1 = iy1 + 1;
        largey(iy1) = i1;       
    elseif y(i1)<x(i1)
        ix1 = ix1 + 1;
        largex(ix1) = i1;   
    end
end
if iy1 == 0 
        break;
    elseif ix1 == 0
        break;
end
largeynum = [];
largexnum = [];
for i2 = 1:iy1
    largeynum(i2) = y(largey(i2));
end

for i3 = 1:ix1
    largexnum(i3) = y(largex(i3));
end
[maxVal, mi] = max(largeynum);
[minVal, mj] = min(largexnum);
i = largey(mi);
j = largex(mj);
Delt = min(x(j) - y(j), y(i) - x(i));
Alph = 1-Delt/(y(i) - y(j));
Mmatrix(i,i) = 0;
Mmatrix(j,j) = 0;
Mmatrix(i,j) = 1;
Mmatrix(j,i) = 1;
T = Alph * I + (1-Alph) * Mmatrix;
for m1 = 1:N
    for m2 = 1:N
        if m1<m2
            U(m1,m2) = sqrt(T(m1,m2));
        else
            U(m1,m2) = -1 * sqrt(T(m1,m2));
        end
    end
end
y = T * y;
Usum =  Usum * U;
Tsum =Tsum*T; 
if norm(y - x)>= 0.001
    k = k+1;
else
    break;
end
end

end