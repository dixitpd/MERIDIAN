clc
clear
%
% We have only kept uniquely used scripts in this directory. Other scripts
% can be found in other directions
%
rand('state',sum(100.*clock))
%
load ../../data
%
% Model specifications
%
lb = [0.01;0.1];ub = [0.05;0.5];
dX  = ub-lb;k1 = 1;kd1 = 10;R = 5;
tx = 10;Lx = [2 10];
%
mes = [m1;m2;m3;m4]; % mean value of measurements
sgs = [s1;s2;s3;s4]/sqrt(5); % standard errors in the means
% 
bix{1} = b1;bix{2} = b2;bix{3} = b3;bix{4} = b4;
% These will be computed using scripts in the directory
% "Predict bin counts for each basis function"
load prdMat_5_1.mat  
nGrd = 5;beta = 1;nDim = nGrd^2;
%
ct = 1;
for ix=1:nGrd
    mux = lb(1) + dX(1)*(ix)/(nGrd+1);
    for iy = 1:nGrd
        muy = lb(2) + dX(2)*(iy)/(nGrd+1);
        mut(ct,:) = [mux muy];
        sgt(ct,:)= (dX/(nGrd))/beta;
        ct = ct+1;
    end
end
%
acp  = 0;
thx = [];
phiO = rand(nDim,1); % Initialize weights
eO = 1e10; % Initialize energy
%
nIter = 40000000;
%
% Perform MCMC chain in the space of weights to sample the posterior
% distribution
%
for iter = 1:nIter 
    if mod(iter,10000) == 0
        100*[iter/nIter acp/iter]
    end
    phiN = phiO + 0.01*randn(nDim,1).*(rand(nDim,1)>0.9)/nDim;
    if min(phiN) > 0
        phiN = phiN/sum(phiN);
        prN  = prdMat*phiN;
        eN   = (prN-mes)./(sgs);eN = 0.5*sum(eN.*eN);
        if exp(-(eN-eO)) > rand
            eO = eN;
            phiO = phiN;
            acp  = acp + 1;
        end
    end
    exx(iter) = eO;
    px(iter) = phiO(1);
    if mod(iter,10000) == 0 && iter > nIter/2 - 1
        fx  = sample_with_phi(phiO,mut,sgt,50);
        thx = [thx;fx];
    end
end
%
i1 = find(thx(:,1)-lb(1) < 0);i2 = find(ub(1) - thx(:,1) < 0);ii = union(i1,i2);
j1 = find(thx(:,2)-lb(2) < 0);j2 = find(ub(2) - thx(:,2) < 0);jj = union(j1,j2);
badI = union(ii,jj);
thx(badI,:) = [];
%
exx = exx((nIter/2)+1:end);
exx = exx(1:10000:end);
save simdata_5_1 exx eO phiO thx

