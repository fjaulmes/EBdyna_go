filename=strcat(DATA_FOLDER,'flux_geometry.mat');
load(filename);
filename=strcat(DATA_FOLDER,'tokamak_PR_map.mat');
load(filename);
filename=strcat(DATA_FOLDER,'B_fields.mat');
load(filename);
filename=strcat(DATA_FOLDER,'q_profile.mat');
load(filename);
filename=strcat(DATA_FOLDER,'pressure_profile.mat');
load(filename);

Btot_PR_map=sqrt(Bpol_PR_map.^2+Btor_PR_map.^2);

NB_THETA=NP;

finesse_data_X=reshape((Rpos_PR_map(:,:)-R0),NB_THETA*Nradial,1);
finesse_data_Z=reshape(Z_PR_map(:,:),NB_THETA*Nradial,1);
[XX ZZ]=meshgrid(X_scale,Z_scale);
% original mapping parameters

%
simulation_size_r=Nradial-2;
mid_Xaxis_large=mid_X_large+round(X_axis/DX);


minX=max(mid_X_large-round(6.2*simulation_size_r),3);
maxX=min(mid_X_large+round(6.2*simulation_size_r),NZ-2);
minZ=max(mid_Z-round(6*elongation*simulation_size_r),3);
maxZ=min(mid_Z+round(6*elongation*simulation_size_r),NZ-2);



size_X=2*ceil(0.5*(maxX-minX));
size_Z=2*ceil(0.5*(maxZ-minZ));

mid_X=ceil(0.5*(maxX-minX))+1;
mid_Z=ceil(0.5*(maxZ-minZ))+1;


% centering the small map on the magnetic axis
X_pos_axis=round(X_axis/DX);
minX=minX+X_pos_axis;
maxX=maxX+X_pos_axis;
size_X=2*ceil(0.5*(maxX-minX));
mid_X=ceil(0.5*(maxX-minX))+1;
mid_Xzero=mid_X;

scale_X=DX*((1:size_X)-mid_X);
scale_Z=DX*((1:size_Z)-mid_Z);





rescaling_to_XZsmall_maps;

Fmirror_coef=(eV/mHe);     % expressing mu in eV.T^-1


Btot_XZ_map=sqrt(BpolX_initial_XZsmall_map.^2+BpolZ_initial_XZsmall_map.^2+Bphi_XZsmall_map.^2);
bX_XZ_map=BpolX_initial_XZsmall_map./Btot_XZ_map;
bZ_XZ_map=BpolZ_initial_XZsmall_map./Btot_XZ_map;
bphi_XZ_map=Bphi_XZsmall_map./Btot_XZ_map;

bX_XZ_map(isnan(bX_XZ_map))=0;
bZ_XZ_map(isnan(bZ_XZ_map))=0;
bphi_XZ_map(isnan(bphi_XZ_map))=0;



gradB_X=zeros(size_X,size_Z);
gradB_Z=zeros(size_X,size_Z);

for (x=3:size_X-2)
    for (z=3:size_Z-2)
        if radial_XZsmall_map(x,z)<Nradial-1
            gradB_X(x,z)=(1/12)*(-Btot_XZ_map(x+2,z)+Btot_XZ_map(x-2,z))+(2/3)*(Btot_XZ_map(x+1,z)-Btot_XZ_map(x-1,z));
            gradB_Z(x,z)=(1/12)*(-Btot_XZ_map(x,z+2)+Btot_XZ_map(x,z-2))+(2/3)*(Btot_XZ_map(x,z+1)-Btot_XZ_map(x,z-1));
        end
    end
end
gradB_X(isnan(gradB_X))=0;
gradB_Z(isnan(gradB_Z))=0;

gradB_X=gradB_X/DX;
gradB_Z=gradB_Z/DX;

vD_X_XZ_map=-Bphi_XZsmall_map.*gradB_Z;
vD_Z_XZ_map=Bphi_XZsmall_map.*gradB_X;
vD_phi_XZ_map=BpolX_initial_XZsmall_map.*gradB_Z-BpolZ_initial_XZsmall_map.*gradB_X;

vD_X_XZ_map=vD_X_XZ_map./(Btot_XZ_map.^2);
vD_Z_XZ_map=vD_Z_XZ_map./(Btot_XZ_map.^2);
vD_phi_XZ_map=vD_phi_XZ_map./(Btot_XZ_map.^2);

% gBX_B_XZ_map=gradB_X./(Btot_XZ_map);
% gBZ_B_XZ_map=gradB_Z./(Btot_XZ_map);

vD_X_XZ_map(isnan(vD_X_XZ_map))=0;
vD_Z_XZ_map(isnan(vD_Z_XZ_map))=0;
vD_phi_XZ_map(isnan(vD_phi_XZ_map))=0;
% gBX_B_XZ_map(isnan(gBX_B_XZ_map))=0;
% gBZ_B_XZ_map(isnan(gBZ_B_XZ_map))=0;
%
%
% Fmirror_XZ_map=(gradB_X.*BpolX_initial_XZsmall_map+gradB_Z.*BpolZ_initial_XZsmall_map)./(Btot_XZ_map);
% Fmirror_XZ_map(isnan(Fmirror_XZ_map))=0;
%
% Fmirror_XZ_map=Fmirror_coef*Fmirror_XZ_map;

%Fmirror_coef=eV/mHe

NB_PSI=Nradial;
NB_THETA=NP;
% NB_THETA=257;
DTHETA=2*pi/(NB_THETA-1);

scale_psi=1:size_r;

scale_phi=2*pi*(0:NB_PHI-1)/(NB_PHI-1);
scale_theta=2*pi*(0:NB_THETA-1)/(NB_THETA-1);
[scale_phi_3D scale_theta_3D scale_psi_3D] =meshgrid(scale_theta,scale_phi,1:size_r);

% psi_scale=(psi_scale-1).*psi_global;


% run('calculate_rotB_vDcurv')

MIN_PSI_VALUE=min(psi_scale)
MAX_PSI_VALUE=max(psi_scale)

psi_XZsmall_map=max(psi_XZsmall_map,MIN_PSI_VALUE);

if psi_XZsmall_map(1,1)>MAX_PSI_VALUE
    psi_norm_XZsmall_map=interp1(psi_scale,1:NB_PSI,min(psi_XZsmall_map,MAX_PSI_VALUE));
else
    psi_norm_XZsmall_map=interp1(psi_scale,1:NB_PSI,psi_XZsmall_map);
    psi_norm_XZsmall_map(isnan(psi_norm_XZsmall_map))=1;
end


% corrected theta maps for complete interpolation
QNB_THETA=round(0.25*NB_THETA);
HQNB_THETA=round(0.5*QNB_THETA);
theta_low_XZsmall_map=theta_XZsmall_map;
theta_up_XZsmall_map=theta_XZsmall_map;
% theta_up_XZsmall_map(find(theta_XZsmall_map<QNB_THETA*DTHETA))=theta_XZsmall_map(find(theta_XZsmall_map<QNB_THETA*DTHETA))+2*pi;
% theta_low_XZsmall_map(find(theta_XZsmall_map>(NB_THETA-QNB_THETA-2)*DTHETA))=theta_XZsmall_map(find(theta_XZsmall_map>(NB_THETA-QNB_THETA-2)*DTHETA))-2*pi;


theta_PR_map_up=theta_PR_map;
theta_PR_map_low=theta_PR_map;

theta_PR_map_up(find(theta_PR_map_up<QNB_THETA*DTHETA))=theta_PR_map_up(find(theta_PR_map_up<QNB_THETA*DTHETA))+2*pi;
theta_PR_map_low(find(theta_PR_map_low>(NB_THETA-QNB_THETA-2)*DTHETA))=theta_PR_map_low(find(theta_PR_map_low>(NB_THETA-QNB_THETA-2)*DTHETA))-2*pi;

theta_data=reshape((theta_PR_map_up(:,:)),NB_THETA*Nradial,1);

theta_up_lin_XZ_map=griddata(finesse_data_X,finesse_data_Z,theta_data,XX,ZZ,'linear');
theta_up_lin_XZ_map(isnan(theta_up_lin_XZ_map)) = 0;
theta_up_lin_XZ_map=theta_up_lin_XZ_map';

theta_up_cub_XZ_map=griddata(finesse_data_X,finesse_data_Z,theta_data,XX,ZZ,'cubic');
theta_up_cub_XZ_map(isnan(theta_up_cub_XZ_map)) = 0;
theta_up_cub_XZ_map=theta_up_cub_XZ_map';

theta_up_XZsmall_map=0.5*(theta_up_lin_XZ_map+theta_up_cub_XZ_map);
theta_up_XZsmall_map=theta_up_XZsmall_map(Xinf:Xsup,Zinf:Zsup);


theta_data=reshape((theta_PR_map_low(:,:)),NB_THETA*Nradial,1);

theta_low_lin_XZ_map=griddata(finesse_data_X,finesse_data_Z,theta_data,XX,ZZ,'linear');
theta_low_lin_XZ_map(isnan(theta_low_lin_XZ_map)) = 0;
theta_low_lin_XZ_map=theta_low_lin_XZ_map';

theta_low_cub_XZ_map=griddata(finesse_data_X,finesse_data_Z,theta_data,XX,ZZ,'cubic');
theta_low_cub_XZ_map(isnan(theta_low_cub_XZ_map)) = 0;
theta_low_cub_XZ_map=theta_low_cub_XZ_map';

theta_low_XZsmall_map=0.5*(theta_low_lin_XZ_map+theta_low_cub_XZ_map);
theta_low_XZsmall_map=theta_low_XZsmall_map(Xinf:Xsup,Zinf:Zsup);

Ne_PR_map=zeros(NP,Nradial);
for p=1:NP
    Ne_PR_map(p,:)=Ne_profile;
end
vA_PR_map=Btot_PR_map./sqrt(mu0*mDT*Ne_PR_map);

vA_data=reshape((vA_PR_map(:,1:end)),NP*Nradial,1);
vA_XZ_map=griddata(finesse_data_X(1:end),finesse_data_Z(1:end),vA_data,XX,ZZ,'cubic');
vA_XZ_map(isnan(vA_XZ_map)) = 0; 
vA_XZ_map=vA_XZ_map';

vA_XZsmall_map=zeros(size_X,size_Z);
vA_XZsmall_map(:,:)=vA_XZ_map(Xinf:Xsup,Zinf:Zsup);

gB0_X_XZ_map=gradB_X;
gB0_Z_XZ_map=gradB_Z;

FILENAME=strcat(DATA_FOLDER,'XZsmall_fields_tokamak_TAE.mat')
save (FILENAME,'q_initial_XZsmall_map','bX_XZ_map','bZ_XZ_map','bphi_XZ_map','size_X','size_Z',...
    'vD_X_XZ_map','vD_Z_XZ_map','vD_phi_XZ_map','gB0_X_XZ_map','gB0_Z_XZ_map','vA_XZsmall_map',...
    'Btot_XZ_map','Bphi_XZsmall_map','BpolX_initial_XZsmall_map','BpolZ_initial_XZsmall_map','Rpos_XZsmall_map',...
    'radial_XZsmall_map','theta_XZsmall_map','theta_low_XZsmall_map','theta_up_XZsmall_map','psiH_XZsmall_map',...
    'psi_norm_XZsmall_map','psi_XZsmall_map','psi_global');


Raxis=R0+X_axis;
FILENAME=strcat(DATA_FOLDER,'motions_map_dimensions.mat')
save (FILENAME,'psi_scale','simulation_size_r','scale_phi','scale_theta','scale_psi','mid_Xaxis_large','mid_Xzero','NB_PSI','NB_THETA','NB_PHI','DTHETA','DPHI',...
    'DX','R0','a','XX_small','ZZ_small','Z_PR_map','scale_X','scale_Z','size_r','X_axis','Z_axis','mid_X','mid_Z','Raxis','Nradial');





