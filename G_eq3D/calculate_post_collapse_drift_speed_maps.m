    filename=strcat(DATA_FOLDER,'flux_geometry.mat');
    load(filename);
    filename=strcat(DATA_FOLDER,'tokamak_PR_map.mat');
    load(filename);
    filename=strcat(DATA_FOLDER,'B_fields.mat');
    load(filename);
    filename=strcat(DATA_FOLDER,'psi_star_evol.mat');
    load(filename);
    filename=strcat(FINESSE_FOLDER,'finesse_data.mat');
    load(filename);   
    
initialize_B_fields_final;

simulation_size_r=222;
mid_Xaxis_large=mid_X_large+round(X_axis/DX);


minX=max(mid_X_large-round(2.5*simulation_size_r),3);
maxX=min(mid_X_large+round(2.5*simulation_size_r),NZ-2);
minZ=max(mid_Z-round(2.2*elongation*simulation_size_r),3);
maxZ=min(mid_Z+round(2.2*elongation*simulation_size_r),NZ-2);

simulation_size_r=round(1.02*simulation_size_r)


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
mid_Xzero=mid_X-X_pos_axis;

scale_X=DX*((1:size_X)-mid_X+X_pos_axis);
scale_Z=DX*((1:size_Z)-mid_Z);

rescaling_to_XZsmall_final_maps


Btot_XZ_map=sqrt(BpolX_final_XZsmall_map.^2+BpolZ_final_XZsmall_map.^2+Bphi_XZsmall_map.^2);
bX_XZ_map=BpolX_final_XZsmall_map./Btot_XZ_map;
bZ_XZ_map=BpolZ_final_XZsmall_map./Btot_XZ_map;
bphi_XZ_map=Bphi_XZsmall_map./Btot_XZ_map;

bX_XZ_map(isnan(bX_XZ_map))=0;
bZ_XZ_map(isnan(bZ_XZ_map))=0;
bphi_XZ_map(isnan(bphi_XZ_map))=0;



% gradB_X=zeros(size_X,size_Z);
% gradB_Z=zeros(size_X,size_Z);
% 
% for (x=3:size_X-2)
%     for (z=3:size_Z-2)
%         if radial_XZsmall_map(x,z)<Nradial-1
%             gradB_X(x,z)=(1/12)*(-Btot_XZ_map(x+2,z)+Btot_XZ_map(x-2,z))+(2/3)*(Btot_XZ_map(x+1,z)-Btot_XZ_map(x-1,z));
%             gradB_Z(x,z)=(1/12)*(-Btot_XZ_map(x,z+2)+Btot_XZ_map(x,z-2))+(2/3)*(Btot_XZ_map(x,z+1)-Btot_XZ_map(x,z-1));
%         end
%     end
% end
% gradB_X(isnan(gradB_X))=0;
% gradB_Z(isnan(gradB_Z))=0;
% 
% gradB_X=gradB_X/DX;
% gradB_Z=gradB_Z/DX;
% 
% vD_X_XZ_map=-Bphi_XZsmall_map.*gradB_Z;
% vD_Z_XZ_map=Bphi_XZsmall_map.*gradB_X;
% vD_phi_XZ_map=BpolX_final_XZsmall_map.*gradB_Z-BpolZ_final_XZsmall_map.*gradB_X;
% 
% vD_X_XZ_map=vD_X_XZ_map./(Btot_XZ_map.^3);
% vD_Z_XZ_map=vD_Z_XZ_map./(Btot_XZ_map.^3);
% vD_phi_XZ_map=vD_phi_XZ_map./(Btot_XZ_map.^3);
% 
% gBX_B_XZ_map=gradB_X./(Btot_XZ_map);
% gBZ_B_XZ_map=gradB_Z./(Btot_XZ_map);
% 
% vD_X_XZ_map(isnan(vD_X_XZ_map))=0;
% vD_Z_XZ_map(isnan(vD_Z_XZ_map))=0;
% vD_phi_XZ_map(isnan(vD_phi_XZ_map))=0;
% gBX_B_XZ_map(isnan(gBX_B_XZ_map))=0;
% gBZ_B_XZ_map(isnan(gBZ_B_XZ_map))=0;





NB_PSI=Nradial
NB_THETA=NP;
% NB_THETA=257;
DTHETA=2*pi/(NB_THETA-1);

scale_psi=1:size_r;
% NB_PHI=129;
% DPHI=2*pi/(NB_PHI-1);
% NB_PHI_DATA_HALF=round(0.5*(NB_PHI-1));

scale_phi=2*pi*(0:NB_PHI-1)/(NB_PHI-1);
scale_theta=2*pi*(0:NB_THETA-1)/(NB_THETA-1);
[scale_phi_3D scale_theta_3D scale_psi_3D] =meshgrid(scale_theta,scale_phi,1:size_r);

% psi_scale=(psi_scale-1).*psi_global;

% psi_scale=(psi_scale-1).*psi_global;



FILENAME=strcat(DATA_FOLDER,'XZsmall_fields_tokamak_post_collapse.mat')
save (FILENAME,'q_final_XZsmall_map','bX_XZ_map','bZ_XZ_map','bphi_XZ_map','size_X','size_Z','Btot_XZ_map','Bphi_XZsmall_map','BpolX_final_XZsmall_map','BpolZ_final_XZsmall_map','Rpos_XZsmall_map','radial_XZsmall_map','theta_XZsmall_map','psiH_XZsmall_map','psi_norm_XZsmall_map','psi_XZsmall_map','psi_global');


Raxis=R0+X_axis;
FILENAME=strcat(DATA_FOLDER,'motions_map_dimensions.mat')
save (FILENAME,'psi_scale','simulation_size_r','scale_phi','scale_theta','scale_psi','mid_Xaxis_large','mid_Xzero','NB_PSI','NB_THETA','NB_PHI','DTHETA','DPHI','NB_PHI_DATA_HALF','DX','R0','a','XX_small','ZZ_small','Z_PR_map','scale_X','scale_Z','size_r','X_axis','Z_axis','mid_X','mid_Z','Raxis');




