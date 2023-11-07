global dt
%Microtubule Parameters
gv_astral = 0.183; %growth velocity of unbound astral MTs (um/s) (M. Peihl, L. Cassimeris 2003 MBoC, J. Wu et al 2011 MBoC)
sv_astral = 0.3; %shrink velocity of unbound astral MTs (um/s) (J. Wu et al 2011 MBoc)
bsv_astral = 0.057; %shrink velocity of bound astral MTs (um/s) (Laan et al, Cell 2012)

k1 = 0.167; %rescue rate of MTs (Komarova et al. JCB, J. Wu et al. 2011 MBoc)

vf = 0.0083; %microns per second, poleward flux (Mitchison JCB 1989, Reviewed in Rogers et al. JCS 2005)
mt_nuc = dt*2*nc; %nucleation rate of MTs (Piehl et al 2004, PNAS)
kappa=3.3; %bending rigidity of MTs (pN*um^2) (Laan et al. Cell 2012...listed as 3.3e-23 N*m^2)
%kappa = 33.12; %bending rigidity of MTs (pN*um^2) (J. Li, H. Jiang Biophysical Journal
%2017, Pavin et al New Journal of Physics 2012)
f_stall = 5; % (pN) stall force of MTs (R. Ma et al, New Journal of Physics 2014)

%Updated Motor Protein Parameters
%f0_dynein=2;
f0_dynein=3.6; %Elshenawy et al. Nature (2019) ...OLD ->%Laan et al 2012... "dynein can maintain the connection to a shrinking MT for several seconds, against pulling forces up to ~5pN (the connection was lost at a mean force of 2.0 +/- 1.4pN)
f0_Eg5 = 1.5; %stall force of tetrameric Eg5 (Y. Shimamoto et al. Dev Cell 2015)
f0_HSET = 1.1; %stall force of HSET (Reinemann et al. Current Biology 2018)
%v0_dynein = 0.8; %(um/s) unloaded velocity of dynein (J. Wu et al MBoC
%2011)
%v0_dynein=0.54;
v0_dynein = 0.86; %L. Urnavicius, et al. Nature (2018) [Carter Lab], Elshenawy et al. Nature (2019) ...OLD ->"dynein gliding velocity" (Laan et al Cell 2012)
v0_Eg5 = 0.2;%0.025; %unloaded velocity of tetrameric Eg5 (bound to one MT) (Shamimoto et al. Dev Cell 2015)
v0_HSET = 0.2;%0.03; %velocity of HSET (um/s) (Roostalu et al. Cell 2018)
%kid motor forces
f0_chk = 1;  % <1 pN in Yajima 1 N - 12, 1-3 in Brouhard
v0_chk = .063;   %160 nm by Yajima, 1 n - 9 ; 0.063 Brouhard (2005)

%Other
LengthFac=0.25;
CrossLink=0.1;
incDist=gv_astral*dt;
MinDBind=4*gv_astral*dt; %1000kDa is ~0.1um.. Dynein is about 3mDa=3000kDa 0.3um
MinDBind_Int=0.06;%Size of Eg5
dtheta=(10*pi)/180;
probD=0.5;
probChr=0.75;
probEg5=0.5;
probHSET=0.5;
probKT=0.5;
repd=2*c_rad; %repulsive distance
eta = 0.07; %cytoplasmic pulling force per unit MT length (equation from J. Li, H. Jiang Phys Rev E 2018, calculation = 6*pi*viscocity*cargo radius*motor velocity... cargo radius=0.005um, motor velocity=1um/s)
drag=20.6;%viscocity of water/sqrt(permeability) 
%C=5;
drag_chr=30;
drag_alpha=100;
C=10; %((pN*um)from J. Li, H. Jiang Biphysical J. 2017)... %C=2e14; %value from J. Li, H. Jiang Phys Rev E 2018; 200 N*um = 2e14 pN*um)
