function [ dys ] = all_derivatives( t,y,K,L )

%
% This script computes all derivatives according to model equations
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

%
% define concentrations
%

% ligand free receptors, plasma membrane
r    = y(1);
% ligand free receptors, endosomes
ri   = y(2);
% ligand bound receptors, plasma membrane
b    = y(3);
% ligand bound receptors, endosomes
bi   = y(4);
% 1 ligand bound dimers, plasma membrane
d1   = y(5);
% 1 ligand bound dimers, endosomes
d1i  = y(6);
% 2 ligand bound dimers, plasma membrane
d2   = y(7);
% 2 ligand bound dimers, endosomes
d2i  = y(8);
% 1 ligand bound phosphorylated dimers, plasma membrane
p1   = y(9);
% 1 ligand bound phosphorylated dimers, endosomes
p1i  = y(10);
% 2 ligand bound phosphorylated dimers, plasma membrane
p2   = y(11);
% 2 ligand bound phosphorylated dimers, endosomes
p2i  = y(12);
% 1L dimer bound to Akt
p1a  = y(13);
% 2L dimer bound to Akt
p2a  = y(14);
% pakt
pakt = y(15);
% free akt
akt  = y(16);

%
% Need to set one of the rate constants for thermodynamic consistency
%
% free receptors, plasma membrane
dys(1)  = ksyn - k1*L*r + kn1*b - ki*r + krec*ri - k2*r*b + kn2*d1;
% free receptors, endosomes
dys(2)  = ki*r - krec*ri - kdeg*ri;
% bound receptors, plasma membrane
dys(3)  = k1*L*r - kn1*b - k2*r*b + kn2*d1 - 2*k2*b*b + 2*kn2*d2 - ki*b + krec*bi;
% bound receptors, endosomes
dys(4)  = ki*b - krec*bi - kdeg*bi;
% 1 ligand bound dimer, plasma membrane
dys(5)  = k2*r*b - kn2*d1 - kap*d1 + kdp*p1 - k1*L*d1 + kn1*d2 - ki*d1 + krec*d1i;
% 1 ligand bound dimer, endosomes
dys(6)  = ki*d1 - krec*d1i - kdeg*d1i + kdp*p1i - kap*d1i;
% 2 ligand bound dimer, plasma membrane
dys(7)  = k2*b*b - kn2*d2 - ki*d2 + krec*d2i - kap*d2 + kdp*p2 + k1*L*d1 - kn1*d2;
% 2 ligand bound dimer, endosomes
dys(8)  = ki*d2 - krec*d2i - kdeg*d2i + kdp*p2i - kap*d2i;
% 1 ligand bound phosphorylated dimer, plasma membrane
dys(9)  = kap*d1 - kdp*p1  - kis*p1 + krecs*p1i - k1*L*p1 + kn1*p2 - kbind*akt*p1 + kunbind*p1a + kpakt*p1a;
% 1 ligand bound phosphorylated dimer, endosomes
dys(10) = kis*p1 - krecs*p1i - kdegs*p1i + kap*d1i - kdp*p1i;
% 2 ligand bound phosphorylated dimer, plasma membrane
dys(11) = kap*d2 - kdp*p2 - kis*p2 + krecs*p2i + k1*L*p1 - kn1*p2 - kbind*akt*p2 + kunbind*p2a + kpakt*p2a;
% 2 ligand bound phosphorylated dimer, endosomes
dys(12) = kis*p2 - krecs*p2i - kdegs*p2i - kdp*p2i + kap*d2i;
% p1 bound to Akt
dys(13) = kbind*p1*akt - kpakt*p1a - kunbind*p1a;
% p2 bound to Akt
dys(14) = kbind*p2*akt - kpakt*p2a - kunbind*p2a;
% pAkt
dys(15) = kpakt*(p1a+p2a) - kdpakt*pakt;
% free Akt
dys(16) = -kbind*akt*(p1+p2) + kdpakt*pakt + kunbind*(p1a+p2a);
%

%
dys = dys';

end

