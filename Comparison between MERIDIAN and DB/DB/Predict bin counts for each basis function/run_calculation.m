clc
clear
%
load ../../data
%
% Model specifications
%
lb = [0.01;0.1];ub = [0.05;0.5];
dX  = ub-lb;k1 = 1;kd1 = 10;R = 5;
tx = 10;Lx = [2 10];
dX = ub-lb;
% 
mes = [m1;m2;m3;m4]; % all mean values
nS = 10000; % number of samples
beta = 32; % see STAR methods
nGrd = 5; % number of grid points per dimension
nDim = nGrd^2; % total number of basis functions points
% %
bix{1} = b1;bix{2} = b2;bix{3} = b3;bix{4} = b4;
nTot = length(m1) + length(m2) + length(m3) + length(m4);
ctx = 1;
%
% prdMat is the matrix of predictions per basis function. 
% When multiplied by weights of each basis functin, it gives DB-based 
% predictions of the bin fractions
%
prdMat = sparse(nTot,nDim); 
for ix=1:nGrd
    mux = lb(1) + dX(1)*(ix)/(nGrd+1);
    for iy = 1:nGrd
        muy = lb(2) + dX(2)*(iy)/(nGrd+1);
        mu  = [mux;muy];
        sg  = (dX/(nGrd))/cfx;
        data = get_bin_fractions(mu,diag(sg),bix,nS);
        s11 = sparse(length(m1),nS);
        s12 = sparse(length(m2),nS);
        s21 = sparse(length(m3),nS);
        s22 = sparse(length(m4),nS);
        for i=1:nS
            s11(data(i,1),i) = 1;
            s12(data(i,2),i) = 1;
            s21(data(i,3),i) = 1;
            s22(data(i,4),i) = 1;
        end
        dtx = [s11;s12;s21;s22];
        cvx = cov(dtx');cvx(isnan(cvx)) = 0;
        prd = mean(dtx')';
        prdMat(:,ctx) = prd;
        ctx = ctx + 1
    end
end
save(strcat('prdMat_',num2str(nGrd),'_',num2str(beta)),'prdMat')
%
