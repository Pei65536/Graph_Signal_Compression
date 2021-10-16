function [s_hat,Supp_A]=CoSaMP(H,y,T,IterNum)
% CoSaMP is designed by D. Needell and J. A. Tropp in the paper: 
% "CoSaMP: Iterative signal recovery from incomplete and inaccurate samples".

%Coded by Kun Qiu (kqiu@iastate.edu)

%%%%%%%%%Function Specification%%%%%%%%%%%%%
%INPUT:
%H:                  Sensing Matrix
%y:                  the measurement column vector
%T:                  sparsity number
%IterNum:            number of iterations (Default: 300)
%
%OUTPUT:
%s_hat:              the signal esitmate
%Supp_A:             sparsity support

if nargin<3
    disp('Error in calling: not enough inputs');
    return;
end
if nargin<4
    IterNum=300;
end

[N,m]=size(H);

% Initial value
a=0;
Supp_A=[];
upsilon=y;

for Count=1:IterNum
    proxy=H'*upsilon;
    [proxy_sort,Omega]=sort(abs(proxy),'descend');
    Omega=Omega(1:2*T);
    A_index=union(Supp_A,Omega);
    A=H(:,A_index);
    
    b=zeros(m,1);
    b(A_index)=A\y;
    [b_sort,A_index]=sort(abs(b),'descend');
    Supp_A=A_index(1:T);
    s_hat=zeros(m,1);
    s_hat(Supp_A)=b(Supp_A);
    upsilon=y-H*s_hat;
end
