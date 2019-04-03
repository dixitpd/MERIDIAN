clc
clear
rand('state',sum(100.*clock))
%
load ../../data
load('../Find optimal Lambda/lmbs.mat')
load data_bin_occupancy
%
% data_bin_occupancy is a sparse binary matrix. The entries in the matrix
% indicate the 4 bin locations occupied by individual sample points
% sample points are generated using the optimal ensemble
%
% cMat is the covariance matrix among the 44 bin fraction indicators
%
n1 = length(l1);n2 = length(l2);
n3 = length(l3);n4 = length(l4);
%
%
mvec = [m1;m2;m3;m4];svec = [s1;s2;s3;s4]/sqrt(5);smat = diag(svec.*svec);
%
% pseudonverse of the covariance matrix
%
cMatI = pinv(full(cMat),1e-4);
l0  = [l1;l2;l3;l4];
mu0 = full(mean(data_bin_occupancy))';
CC  = cMatI*mu0 + l0;
BB  = -cMatI;
%
% Mean locations and standard deviations of the Lagrange multiplier
% Gaussian distribution
%
lmb_mu = CC + BB*mvec;lmb_ss = BB*smat*BB';
% 
% Model specifications
%
bix{1} = b1;bix{2} = b2;bix{3} = b3;bix{4} = b4;acp = 0;
lb = [0.01;0.1];ub = [0.05;0.5];
dX  = ub-lb;k1 = 1;kd1 = 10;R = 5;
tx = 10;Lx = [2 10];
% 
tBig = [];
%
% Sample Lagrange multipliers according to the Gaussian distribution and
% then sample parameter points using MCMC according to the Lagrange
% multipliers
%
for iter=1:5000
    iter
    ll = mvnrnd(lmb_mu,lmb_ss,1);
    lX(:,iter) = ll';
    lx1=ll(1:n1);lx2=ll(n1+1:n1+n2);lx3=ll(n1+n2+1:n1+n2+n3);lx4=ll(n1+n2+n3+1:end);
    [thx] = pred_at_lambda_sim(lx1,lx2,lx3,lx4,bix,500000,10000);
    tBig{iter} = thx;
    save tBig tBig
end
save tBig tBig

