function [ bI ] = decide_which_bin_akt( bns,dat )


% Assign bins to Akt data

%
% Data arranged as [t=0 L1t1 L1t2 ... L2t1 L2t2 ...]
%
[a b] = size(bns);
bI = sparse(a,b+1);
%
% Akt background bin (t=0)
%
akt0 = dat(1,1);
bn0  = akt0 - bns(1,:);
ii   = find(bn0 > 0);
if isempty(ii)
    bI(1,1) = 1;
else
    bI(1,max(ii)+1) = 1;
end
%
dat = dat(2:end,:);dat = reshape(dat,20,1);
for i=2:21
    bnHere = bns(i,:);
    dt     = dat(i-1) - bnHere;
    ii = find(dt>0);
    if isempty(ii)
        bI(i,1) = 1;
    else
        bI(i,max(ii)+1) = 1;
    end
end



end

