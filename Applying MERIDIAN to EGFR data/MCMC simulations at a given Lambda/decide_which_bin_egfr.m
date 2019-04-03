function [ bI ] = decide_which_bin_egfr( bns,dat )


% This script assigns EGFR bins

%
% Data arranged as [0 1 100] ng/ml
%
bI = sparse(3,11);
%
d0 = dat(1) - bns(1,:);
ii = find(d0 > 0);
if isempty(ii)
    bI(1,1) = 1;
else
    bI(1,max(ii) + 1) = 1;
end
%
d0 = dat(2) - bns(2,:);
ii = find(d0 > 0);
if isempty(ii)
    bI(2,1) = 1;
else
    bI(2,max(ii) + 1) = 1;
end
%
d0 = dat(3) - bns(3,:);
ii = find(d0 > 0);
if isempty(ii)
    bI(3,1) = 1;
else
    bI(3,max(ii) + 1) = 1;
end


end

