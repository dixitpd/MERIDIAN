function [thx] = pred_at_lambda_sim(l1,l2,l3,l4,bix,nS,nBurn)
%
rand('state',sum(100.*clock))
%
lb = [0.01;0.1];ub = [0.05;0.5];
dX  = ub-lb;k1 = 1;kd1 = 10;R = 5;
tx = 10;Lx = [2 10];
%
prt = 0;eO = 1e5;kiO = 0.01 + rand*0.04;kisO = 0.1 + rand*0.4;thx = [];

for iter=1:nS+nBurn
    if mod(iter,50000) == 0
        100*iter/(nS+nBurn)
    end
    %
    kiN = kiO + 0.005*randn;kisN = kisO + 0.05*randn;
    if kiN < ub(1) && kiN > lb(1) && kisN < ub(2) && kisN > lb(2)
        paramx = [k1;kd1;kiN;kisN;R];
        [sl1 ss1] = solve_diffEq(Lx(1),tx,paramx);
        [sl2 ss2] = solve_diffEq(Lx(2),tx,paramx);
        dat = [sl1(2);ss1(2);sl2(2);ss2(2)];
        cx  = decide_which_bin(bix,dat);
        eN  = l1(cx(1)) + l2(cx(2)) + l3(cx(3)) + l4(cx(4));
        if exp(-(eN-eO)) > rand
            kiO = kiN;kisO = kisN;eO = eN;
        end
    end
    if mod(iter,100) == 0 && iter > nBurn-1
        thx = [thx;[kiO kisO]];
    end
end


end