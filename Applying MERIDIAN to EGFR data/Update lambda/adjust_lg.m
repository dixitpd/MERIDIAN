clc
clear
%
% This script will update the Lagrange multipliers
%
dsx = [0 0.01 0.03 0.1 0.32 1 3.16 10 31.62 100];
tme = [0 5 10 15 30 45 90 180];
%
load akt_bindata
load segfr_bindata
indI = 24;
stx = num2str(indI);
load(strcat('Iter',stx,'/lmbs_notail.mat'));
lst = strcat('Iter',stx,'/data*.mat');myFiles = dir(lst);
cmx = 0;paP = zeros(21,11);seP = zeros(3,11);
cvx = [];kns = [];
for k = 1:length(myFiles)
    baseFileName = myFiles(k).name;
    file_stats = dir(strcat('Iter',stx,'/',baseFileName));
    if file_stats.bytes > 256
        load(strcat('Iter',stx,'/',baseFileName))
        paP = paP + paktbns;
        seP = seP + segfrbns;
        cvx = [cvx crx];
        kns = [kns knx];
        cmx = cmx + ctbin;
    end
end
%
cvx = cvx';
cvx = cov(cvx);
paP = paP/cmx;seP = seP/cmx;
prd = [reshape(paP,21*11,1);reshape(seP,3*11,1)];
mes = [reshape(mbins_akt,21*11,1);reshape(mbins_egfr,3*11,1)];
ses = [reshape(sbins_akt,21*11,1);reshape(sbins_egfr,3*11,1)];
%
ev = eig(cvx);mnx = min(ev(ev>1e-8));
cplus = pinv(full(cvx),mnx);
%
sx = cplus*diag(ses.*ses)*cplus';
for i=1:264
    g = cvx(i,:);
    trm = g*sx*g';
    sprdx(i) = sqrt(trm);
end
%
for iter=1:500
    dxx = (prd-mes);
    dL = 0.01*iter*dxx;
    prx = prd-cvx*dL;
    dlx(:,iter) = dL;
    rprx = mean(abs(prx-mes));
    ex(iter) = rprx;
end
%
ii = find(ex==min(ex));
dL1 = dlx(:,ii);
prx = prd - cvx*dL1;
[mean(abs(mes-prd)./mes) median(abs((prd-mes)./mes))]
[mean(abs(prx-mes)./mes) median(abs(prx-mes)./mes)]
dLA = dL1(1:21*11);dLS = dL1(21*11+1:end);
dLA = reshape(dLA,21,11);dLS = reshape(dLS,3,11);
%
lambda_akt   = lambda_akt    + dLA;
lambda_egfr  = lambda_egfr   + dLS;

lst = strcat('Iter',num2str(indI+1),'/lmbs_notail.mat');
save(lst,'lambda_akt','lambda_egfr','kns')
%