clc
clear
rand('state',sum(100.*clock))
%
doses = [0.1 0.316 1 3.16 10 100];
tmes  = 60*[0 5 15 30 45 180];
%
load param_bounds.dat
lb = param_bounds(:,1);
nPr = length(lb)
ub = param_bounds(:,2);
load cutoffs.dat
%
load akt_bindata
load egfr_bindata
load lmbs_notail
[a nS] = size(kns);
%
% Parameters from the previous lambda iteration are chosen randomly chosen
% as starting point
%
flg = 0;
while flg < 1
   fl1 = 0;fl2 = 0;
   KOld = kns(:,randsample(nS,1));
   d1  = KOld - lb;
   if min(d1) < -1e-5
      fl1 = 1;
   end
   d2  = ub - KOld;
   if min(d2) < -1e-5
      fl2 = 1;
   end
   if fl1==0 && fl2==0
      flg = 1;
   end
end
'Found initial parameters'
%
eOld = 1000;
acc  = 0;
bimbo = 1;
knx  = [];crx = [];
nIter = 50000;
kn = KOld;
paktbns = sparse(21,11);segfrbns = sparse(3,11);ctbin = 0;
%
% MCMC in the parameter space
%
for iter=1:nIter
    %
    [kn eOld flg] = monte_carlo_step(kn,lb,ub,ub-lb,20,eOld,doses,tmes,lambda_akt,lambda_egfr,bins_akt,bins_egfr,cutoffs);
    if mod(iter,50) == 0
        [100*iter/nIter 100*acc/iter]
        if iter > 4999
            prdx = modelPreds(kn,doses,tmes);
            [d_a d_e] = extract_preds(prdx,kn);
            pa = decide_which_bin_akt(bins_akt,d_a);
            sa = decide_which_bin_egfr(bins_egfr,d_e);
            tt = [reshape(pa,21*11,1);reshape(sa,3*11,1)];
            crx = [crx tt];
            knx = [knx kn];
            ctbin     = ctbin + 1
            paktbns   = paktbns   + pa;
            segfrbns  = segfrbns  + sa;
            save datax knx paktbns segfrbns ctbin crx
        end
    end
    acc = acc + flg;
end
% 
save datax knx paktbns segfrbns ctbin crx
