function [ct11,ct12,ct21,ct22] = pred_at_lambda_trapezoid(l11,l12,l21,l22,bix,nGr_int)
%
% This script will estimate predicted bin fractions using 
% numerical integration (trapezoidal rule)
%
n11 = length(l11);n12 = length(l12);
n21 = length(l21);n22 = length(l22);
%
lb = [0.01;0.1];ub = [0.05;0.5];
dX  = ub-lb;k1 = 1;kd1 = 10;R = 5;
tx = 10;Lx = [2 10];
%
ct11 = zeros(n11,1);ct12 = zeros(n12,1);
ct21 = zeros(n21,1);ct22 = zeros(n22,1);
%
% Coefficient for each grid pointthe trapezoidal integration rule in 2d
fcx = 2*ones(nGr_int+1,1);fcx(1) = 1;fcx(end) = 1;
fcx = fcx'.*fcx;
%
prt = 0;
for iX=1:nGr_int+1
    kdeg = lb(1) + dX(1)*(iX-1)/nGr_int;
    for jX=1:nGr_int+1
        kdegs = lb(2) + dX(2)*(jX-1)/nGr_int;
        %   Solve
        paramx = [k1;kd1;kdeg;kdegs;R];
        [sl1 ss1] = solve_diffEq(Lx(1),tx,paramx);
        [sl2 ss2] = solve_diffEq(Lx(2),tx,paramx);
        dat = [sl1(2);ss1(2);sl2(2);ss2(2)];
        %
        cx  = decide_which_bin(bix,dat);
        %
        pr  = l11(cx(1)) + l12(cx(2)) + l21(cx(3)) + l22(cx(4));
        pr  = fcx(iX,jX)*exp(-pr);prt = prt  + pr;
        ct11(cx(1)) = ct11(cx(1)) + pr;
        ct12(cx(2)) = ct12(cx(2)) + pr;
        ct21(cx(3)) = ct21(cx(3)) + pr;
        ct22(cx(4)) = ct22(cx(4)) + pr;
    end
end
ct11 = ct11/prt;ct12 = ct12/prt;
ct21 = ct21/prt;ct22 = ct22/prt;
%
end




