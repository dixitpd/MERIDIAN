function [ frc ] = binProp( x0,bins )

%
% This script assigns bin fractions using data and bin locations
%


frc = [];
nb = length(bins);
frc = [frc;sum(x0<bins(1))];

for i=1:nb-1
    s1 = find(x0 > bins(i));s2 = find(x0 < bins(i+1));
    frc = [frc;length(intersect(s1,s2))];
end

frc = [frc;sum(x0>bins(nb))];
frc = frc/sum(frc);

end