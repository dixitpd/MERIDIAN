function [ Ktx etx flg_ac ] = monte_carlo_step( KOld,k_l,k_u,dk,npr,eOld,doses,tmes,l_a,l_e,b_a,b_e,cutoff )


%
% This script performs one monte carlo step
%

% accept flag
flg_ac = 0;

% change parameters
ids = floor(10*rand)+1;
ids = randsample(npr,ids);
flg = sparse(npr,1);flg(ids) = 1;
alp = 0.2; % change up to 5%
KNew = KOld + alp*(2*rand(npr,1)-1).*flg.*dk;


% check boundaries
d1  = KNew - k_l;
fl1 = 0;
if min(d1) < -1e-5
    fl1 = 1;
end
d2  = k_u - KNew;
fl2 = 0;
if min(d2) < -1e-5
    fl2 = 1;
end
% flag for Kd
fl3 = 0;
Kd1 = 10^(KNew(2)-KNew(1));
Kd2 = 10^(KNew(4)-KNew(3));
if  Kd1 > 5 && Kd1 < 80 && Kd2 > 50 && Kd2 < 600
   fl3 = 1;
end

% only if fl1 and fl2 are good otherwise no point simulating!
if fl1 == 0 && fl2 == 0 && fl3 == 1
   % evaluate new solution
   kk = modelPreds(KNew,doses,tmes);
   % extract predictions
   [d_a d_e] = extract_preds(kk,KNew);
   % is there a decrease in EGFR levels
   flEGFR = ( d_e(3)-KNew(20) )/( d_e(1) - KNew(20) );
   % has Akt decreased over time at the highest dose?
   flAkt = (d_a(5,5)-KNew(19))/(d_a(2,5)-KNew(19));
   % evaluate energy
   eNew = evaluate_energy( l_a,l_e,b_a,b_e,d_a,d_e );
   if exp(-(eNew-eOld)) > rand && is_within_cutoffs(d_a,d_e,cutoff) == 1 && flEGFR < 1 && flAkt <1
      eOld = eNew;
      KOld = KNew;
      flg_ac  = 1;
   end
end

% return variable!
Ktx = KOld;
etx = eOld;


end

