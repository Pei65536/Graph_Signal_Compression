function [MV] = BitAllocation(Lda,K,M,eta)

%输入Lda-特征值降序排列 K-量化器个数 M-总状态数 eta-参数


MV=ones(K,1);
%eta=3;
for i=1:K
    DeltaG(i)=6*((2*MV(i)+1)*eta^2*Lda(i)^2)/(3*(MV(i)+1)^2+2*eta^2)/(3*MV(i)^2+2*eta^2)/log2((MV(i)+1)/MV(i));
end

while sum(log2(MV))<log2(M)
    [p,q]=max(DeltaG);
    MV(q)=MV(q)+1;
    DeltaG(q)=6*((2*MV(q)+1)*eta^2*Lda(q)^2)/(3*(MV(q)+1)^2+2*eta^2)/(3*MV(q)^2+2*eta^2)/log2((MV(q)+1)/MV(q));
end

end

