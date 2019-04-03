function [ y0s ] = initial_conditions( K )

% 
% Initial conditions for individual variables
%

K1 = K(1:17);
K2 = K(18:end);
K1 = 10.^(K1);
K  = [K1;K2];

%
% define rates
%

% EGF binding to EGFR monomer
k1      = K(1);
% EGF unbinding from EGFR
kn1     = K(2);
% EGFR EGF-EGFR dimerization
k2      = K(3);
% EGFR-EGF-EGFR undimerization
kn2     = K(4);
% receptor phosphorylation
kap     = K(5);
% receptor dephosphorylation
kdp     = K(6);
% degradation of inactive
kdeg    = K(7);
% degradation of active
kdegs   = K(8);
% internalization of inactive
ki      = K(9);
% internalization of active
kis     = K(10);
% recycling of inactive
krec    = K(11);
% recycling of active
krecs   = K(12);
% rate of pEGFR binding to Akt
kbind   = K(13);
% rate of pEGFR-Akt unbinding
kunbind = K(14);
% Rate of Akt phosphorylation
kpakt   = K(15);
% rate of pAkt dephosphorylation
kdpakt  = K(16);
% EGFR delivery rate
ksyn    = K(17);

% Total Akt abundance
Akt0   = K(18);
% surface receptors
R0  = (kdeg + krec)*ksyn/(kdeg*ki);
% endosomal receptors
R0i =  ksyn/kdeg;

%
% define concentrations
%

y0s = zeros(16,1);

% ligand free receptors, plasma membrane
y0s(1)    = R0;
% ligand free receptors, endosomes
y0s(2)    = R0i;
% free akt
y0s(16)   = Akt0;



end

