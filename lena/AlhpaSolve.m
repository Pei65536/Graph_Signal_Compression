function [ Aplhai ] = AlhpaSolve(lambdai, Lnai)
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here
%lambdaiΪ����Ҫ������ֵ����
%LnaiΪ���������ص�����

N = length(lambdai);
cvx_begin
    variable x(N)
    for i=1:N
        y(i) = lambdai(i)^2*inv_pos(x(i)+1);
        s(i) = sum(x(1:i));
    end
    minimize( sum(y) )	%Ŀ�꺯��
    subject to
        for j = 1:N
            s(j)>=sum(Lnai(1:j));
        end 
        sum(x) == sum(Lnai);
cvx_end

Aplhai = x;
end

