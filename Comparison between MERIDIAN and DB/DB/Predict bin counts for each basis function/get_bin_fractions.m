function [dat] = get_bin_fractions(mu,sg,bix,nS)
%

%
% Given the location and the width of the multivariate Gaussian basis function
% this script samples parameter points and computes which bins the
% individual data point belong to
%
ct = 1;
%
% Model specifications (hard coded)
%
lb = [0.01;0.1];ub = [0.05;0.5];
k1 = 1;kd1 = 10;R = 5;
tx = 10;Lx = [2 10];
%
while ct < nS+1
    prx = mvnrnd(mu,sg.*sg,1)';
    d1  = min(prx-lb);d2 = min(ub-prx);
    if d1 > 0 && d2 > 0
        paramx = [k1;kd1;prx(1);prx(2);R];
        [sl1 ss1] = solve_diffEq(Lx(1),tx,paramx);
        [sl2 ss2] = solve_diffEq(Lx(2),tx,paramx);
        dt = [sl1(2);ss1(2);sl2(2);ss2(2)];
        cx  = decide_which_bin(bix,dt);
        %
        dat(ct,:) = decide_which_bin(bix,dt);
        ct = ct + 1;
    end
end
%


end