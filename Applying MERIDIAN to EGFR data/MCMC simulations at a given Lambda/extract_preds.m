function [ aktD segfrD ] = extract_preds( yy,Kt )

%
% This script extracts Akt and EGFR data from the solved differnetial equations
% Sadly, this file is hard coded. 
%

% Doses are (0.1 0.316 1 3.16 10 100) ng/ml

%
% The doses used in the stimulation are (0.1 0.316 1 3.16 10 100 )
% for Akt
% EGFR will come as an array 0 min (whatever dose), 1 (180 min)
% and 100 (180 min)
%
aktI   = [1 2 4 5 6];
aktD   = [];
%
for i=1:5
    sl = yy{aktI(i)};
    pakt = sl(1:end-1,15);
    aktD = [aktD pakt];
end
%
aktD = aktD + Kt(19);
%
% EGFR for 0 min
%
sl = yy{1};
sl = sl(1,:);
%
r0 = sl(1) + sl(3) + 2*(sl(5) + sl(7) + sl(9) + sl(11) + sl(13) + sl(14)); 
%
% EGFR for 1 ng/ml
%
sl = yy{3};
sl = sl(end,:);
r1 = sl(1) + sl(3) + 2*(sl(5) + sl(7) + sl(9) + sl(11) + sl(13) + sl(14));
%
% EGFR for 100 ng/ml
%
sl   = yy{6};
sl   = sl(end,:);
r100 = sl(1) + sl(3) + 2*(sl(5) + sl(7) + sl(9) + sl(11) + sl(13) + sl(14)); 
segfrD = [r0 r1 r100];
segfrD = segfrD + Kt(20);


end

