% *************************************************************
% defining the helical field 
% *************************************************************

BHpol_map=zeros(NR,NZ);
BHtor_map=zeros(NR,NZ);

BHpol_X_PR_map=zeros(NP,Nradial);
BHpol_Z_PR_map=zeros(NP,Nradial);

BHpol_PR_map=q_PR_map.*Bpol_PR_map;


BHpol_X_PR_map=q_PR_map.*BX_PR_map;
BHpol_Z_PR_map=q_PR_map.*BZ_PR_map;



% *************************************************************
% Helical field curvature
% *************************************************************


BHpol_X_map=q_XZ_map.*BR_XZ_map.*mask_XZ_map;
BHpol_Z_map=q_XZ_map.*BZ_XZ_map.*mask_XZ_map;



BHpol_tot=BHpol_X_map.^2+BHpol_Z_map.^2;
