clc
clear
%
% Only unique scripts are included. 
%
rand('state',sum(100.*clock))
load ../../data
%
% Initial guess of Lagrange multipliers
% 
l1 = -0.2*log(m1);l2 = -0.2*log(m2);
l3 = -0.2*log(m4);l4 = -0.2*log(m3);
%
n1 = length(m1);n2 = length(m2);
n3 = length(m3);n4 = length(m3);
%
%  Model specifications
%
lb = [0.01;0.1];ub = [0.05;0.5];
dX  = ub-lb;k1 = 1;kd1 = 10;R = 5;
tx = 10;Lx = [2 10];
%
nGr_int = 1000; % integration grid points
bix{1} = b1;bix{2} = b2;bix{3} = b3;bix{4} = b4;
%
alp = 2; % learning rate
mes  = [m1;m2;m3;m4];
ses = [s1;s2;s3;s4];
%
% Search Lagrange multipliers
%
nIter = 100; % Number of iterations
for ixI=1:nIter
    [c1 c2 c3 c4] = pred_at_lambda_trapezoid(l1,l2,l3,l4,bix,nGr_int);
    %
    prd  = [c1;c2;c3;c4];
    a = (abs(prd-mes)./mes);a(isnan(a)) = [];
    %
    dX1  = c1-m1;dX2 = c2-m2;
    dX3  = c3-m3;dX4 = c4-m4;
    [mean(a)]
    %
    l1 = l1 + alp*dX1;
    l2 = l2 + alp*dX2;
    l3 = l3 + alp*dX3;
    l4 = l4 + alp*dX4;
    %
    save lmbs l1 l2 l3 l4 c1 c2 c3 c4
    plot(prd,'r--')
    hold on
    errorbar(mes,ses,'k')
end
