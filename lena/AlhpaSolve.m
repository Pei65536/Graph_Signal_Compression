function [ Aplhai ] = AlhpaSolve(lambdai, Lnai)
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here
%lambdai为所需要的奇异值序列
%Lnai为与比特数相关的序列

N = length(lambdai);
cvx_begin
    variable x(N)
    for i=1:N
        y(i) = lambdai(i)^2*inv_pos(x(i)+1);
        s(i) = sum(x(1:i));
    end
    minimize( sum(y) )	%目标函数
    subject to
        for j = 1:N
            s(j)>=sum(Lnai(1:j));
        end 
        sum(x) == sum(Lnai);
cvx_end

Aplhai = x;
end

