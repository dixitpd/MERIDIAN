function [solx sols] = solve_diffEq(L,t,paramx)


% Given model parameters, this script solves the differential equations
% and computes the number activated at time "t" and at steady state.


k1 = paramx(1);kd1 = paramx(2);ki = paramx(3);kis = paramx(4);
ksyn = ki*paramx(5);

AMat = [-k1*L-ki kd1;k1*L -kis-kd1];
bVec = [ksyn;0];
x0   = [paramx(5);0];
xs   = -inv(AMat)*bVec;
solx = xs + expm(AMat*t)*(x0-xs);
sols = xs;


end

