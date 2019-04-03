function [ cx ] = decide_which_bin( bi,dat )

% Given 4 data points and bin locations for the 4 experimental conditions
% this script decides which bin the data belongs to

% If P is the number of  active receptors, 
% Data arranged as: [P(L1,t1) P(L1,t2) P(L2,t1) P(L2,t2)]


%
for i = 1:4
    d1 = dat(i);
    bn0 = d1-bi{i};
    ii = find(bn0>0);
    if isempty(ii)
        cx(i) = 1;
    else
        cx(i) = max(ii) + 1;
    end
end
%



end

