function [ flg ] = is_within_cutoffs( akt,egfr,cutoffs )
%
% This script figures out whether akt and egfr values are within their ranges
%


% by default assume that everything is bad
flg = 0;
dat = akt;d0 = dat(1,1);
%
dat = dat(2:end,:);dat = reshape(dat,20,1);
dat = [d0;dat;egfr'];
%
% check boundaries
d1  = dat - cutoffs(:,1);
fl1 = 0;
if min(d1) > -1e-5
    fl1 = 1;
end
d2  = cutoffs(:,2) - dat;
fl2 = 0;
if min(d2) > -1e-5
    fl2 = 1;
end
%
flg = fl1*fl2;


end


