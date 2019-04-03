function [ yy ] = modelPreds( K,doses,tms )

% This script solves the differential equation using ode15

Ktrial = K;
nDose  = length(doses);
yy = {};
%
options = odeset('RelTol',1e-5);
for doseI = 1:nDose
    L          = doses(doseI);
    y0         = initial_conditions(Ktrial);
    [t y]      = ode15s(@(t,y)all_derivatives(t,y,Ktrial,L),tms,y0,options);
    yy{doseI} = y;
end