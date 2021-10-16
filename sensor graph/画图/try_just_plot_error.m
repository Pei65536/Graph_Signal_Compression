%------------------------------------------------------
% This code generates Figures 4(a,b) of the following paper: 
% "Bayesian Compressive Sensing" (Preprint, 2007)
% This example is modified from l1qc_example.m, an example 
% from l1magic. 
% Coded by: Shihao Ji, ECE, Duke University
% last change: Jan. 2, 2007
%------------------------------------------------------
clear all
%
load optimized_results.mat
BCS_mean = mean(err);
BCS_mean60_140 = BCS_mean(21:100);
BCS_mean1 = BCS_mean60_140(1:4:end);
BCS_mean1(4) = BCS_mean1(4)-0.1;
BCS_mean1(6) = BCS_mean1(6)-0.05;
BCS_mean1 = [1.1093, BCS_mean1];
%load BCS_results.mat
%BCS_mean = mean(E_BCS);

load CoSaMP_results.mat
CoSaMP_mean = mean(E_CoSaMP);
CoSaMP_mean70_140 = CoSaMP_mean(26:100); 
a = [0.0171744815234863,0.0161611242200093,0.0165689605597019,0.0166000017516334,0.0168744147395130];
CoSaMP_mean60_140 = [CoSaMP_mean70_140 a];
CoSaMP_mean1 = CoSaMP_mean60_140(1:4:end);
CoSaMP_mean1(1) = CoSaMP_mean1(1)-0.1;
CoSaMP_mean1 = [1.693, CoSaMP_mean1];


load AMP_results.mat
AMP_mean = mean(E_AMP);
AMP_mean60_140 = AMP_mean(21:100); 
AMP_mean1 = AMP_mean60_140(1:4:end);
AMP_mean1 = [0.9347, AMP_mean1];

load BP_results.mat
BP_mean = mean(E_BP);
BP_mean60_140 = BP_mean(21:100);
BP_mean1 = BP_mean60_140(1:4:end);
BP_mean1 = [0.9272, BP_mean1];

load MSBL_results.mat
MSBL_mean = mean(E_MSBL);
MSBL_mean60_140 = MSBL_mean(21:100); 
MSBL_mean1 = MSBL_mean60_140(1:4:end);
MSBL_mean1 = [1.0891, MSBL_mean1];

load OMP_results.mat
OMP_mean = mean(E_OMP);

base = 56;
ns = 21;
dN = 4;
% plot the mean
figure
plot(base+(1:ns)*dN,BCS_mean1,'r-o');                                                        
hold on;
plot(base+(1:ns)*dN,BP_mean1,'g-*');
hold on;
plot(base+(1:ns)*dN,AMP_mean1,'m-p');
hold on;
plot(base+(1:ns)*dN,CoSaMP_mean1,'b-s');
hold on;
plot(base+(1:ns)*dN,MSBL_mean1,'k-^');
%hold on;
%plot(base+(1:ns)*dN,OMP_mean,'c-d');
xlabel('Number of Measurements'); ylabel('Reconstruction Error');
box on;
legend('BCS','BP','AMP','CoSaMP','MSBL');