ikink_kinetic_energy_phi=zeros(NB_PHI,1);
Btot_ini_PR_map=sqrt(BX_PR_map.^2+BZ_PR_map.^2+Btor_PR_map.^2);

% big arrays, be careful with them
vD_X_PR_map_phi=zeros(NB_PHI,NB_THETA,size_r);
vD_Z_PR_map_phi=zeros(NB_PHI,NB_THETA,size_r);
vD_phi_PR_map_phi=zeros(NB_PHI,NB_THETA,size_r);

% bX_XZ_map(isnan(bX_XZ_map))=0;
% bZ_XZ_map(isnan(bZ_XZ_map))=0;
% bphi_XZ_map(isnan(bphi_XZ_map))=0;
psi_star_dot_omega_map(:,:)=psi_star_dot_evol(FRAME_NUMBER,:,:);

Efield_phi_map_phi=Efield_X_map_phi*0;

for phi_rank=1:NB_PHI
    phi_rank
    B_PR_map=squeeze(Btot_map_phi(phi_rank,:,:));
    bX_PR_map=squeeze(bX_map_phi(phi_rank,:,:));
    bZ_PR_map=squeeze(bZ_map_phi(phi_rank,:,:));
	%gps_PR_map=squeeze(grad_psi_star_map_phi(phi_rank,:,:));
	%BtotX_PR_map=BX_PR_map(:,1:size_r)+bX_PR_map;
	%BtotZ_PR_map=BZ_PR_map(:,1:size_r)+bZ_PR_map;
    %EX_PR_map=squeeze(Efield_X_map_phi(phi_rank,:,:));
    %EZ_PR_map=squeeze(Efield_Z_map_phi(phi_rank,:,:));
    %Ephi_PR_map=squeeze(Efield_phi_map_phi(phi_rank,:,:));
	
	%bphi_PR_map=-(BtotX_PR_map.*EX_PR_map+BtotZ_PR_map.*EZ_PR_map)./max(Ephi_PR_map,1);
	%bphi_PR_map(isinf(bphi_PR_map))=0;
	%Btotphi_PR_map=Btor_PR_map(:,1:size_r)+bphi_PR_map;
    %bphi_PR_map=sqrt(1-(bX_PR_map.^2+bZ_PR_map.^2));
    %BpolX_PR_map=bX_PR_map.*B_PR_map-BXini_PR_map(:,1:size_r);
    %BpolZ_PR_map=bZ_PR_map.*B_PR_map-BZini_PR_map(:,1:size_r);
	% delta Bstar is now stored in the maps
    
    % the mapp contain the total perturbed fields
    BpolX_PR_map=bX_PR_map.*B_PR_map;
    BpolZ_PR_map=bZ_PR_map.*B_PR_map;
    dBpolX_PR_map=BpolX_PR_map-BX_PR_map(:,1:size_r);
    dBpolZ_PR_map=BpolZ_PR_map-BZ_PR_map(:,1:size_r);
%     Bphi_PR_map=bphi_PR_map.*B_PR_map-Btor_PR_map(:,1:size_r);
    dB_PR_map=B_PR_map-Btot_ini_PR_map(:,1:size_r);
	% careful to take all terms up to order 1 in deltaB
	% neglecting small toridal component
	%dB_PR_map=sqrt(dBpolX_PR_map.^2+dBpolZ_PR_map.^2);
	%dB_PR_map=0.5*Btot_ini_PR_map(:,1:size_r).*((B_PR_map./Btot_ini_PR_map(:,1:size_r)).^2-1)

%     B_data=reshape(B_PR_map(:,1:NR),NB_THETA*NR,1);
    Btot_ini_XZ_small_map=sqrt(BpolX_initial_XZsmall_map.^2+BpolZ_initial_XZsmall_map.^2+Bphi_XZsmall_map.^2);
    
    B_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,B_PR_map(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
    B_XZ_map(isnan(B_XZ_map))=0;
    B_XZ_map=B_XZ_map';
    dBX_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,dBpolX_PR_map(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
    dBX_XZ_map(isnan(dBX_XZ_map))=0;
    dBX_XZ_map=dBX_XZ_map';
    dBZ_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,dBpolZ_PR_map(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
    dBZ_XZ_map(isnan(dBZ_XZ_map))=0;
    dBZ_XZ_map=dBZ_XZ_map';

%     Bphi_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,Bphi_PR_map(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
%     Bphi_XZ_map(isnan(Bphi_XZ_map))=0;
%     Bphi_XZ_map=Bphi_XZ_map';
[XX,ZZ] = meshgrid(scale_X,scale_Z);


dbdata=reshape(dB_PR_map,NP*size_r,1);
dB_XZ_map=gridfit(finesse_data_X,finesse_data_Z,dbdata,scale_X,scale_Z,'smoothness',0.5);
dB_XZ_map=dB_XZ_map';

% dbdata=reshape(dB_PR_map,NP*size_r,1);
% dB_XZ_map=griddata(finesse_data_X,finesse_data_Z,dbdata,XX,ZZ,'cubic');
% dB_XZ_map(isnan(dB_XZ_map))=0;
% dB_XZ_map=dB_XZ_map';

%     dB_XZ_map=dtinterp(finesse_mesh,finesse_mesh_dtri,dB_PR_map(IDTRI),XX_small,ZZ_small,DT_INTERPOLATION_METHOD);
%     dB_XZ_map(isnan(dB_XZ_map))=0;
%     dB_XZ_map=dB_XZ_map';
    calculate_kink_kinetic_energy;
    ikink_kinetic_energy_phi(phi_rank)=phi_rank_kinetic_energy;

    gradB_X=zeros(sizeX,sizeZ);
    gradB_Z=zeros(sizeX,sizeZ);
    dgradB_X=zeros(sizeX,sizeZ);
    dgradB_Z=zeros(sizeX,sizeZ);
    
    for (x=3:sizeX-2)
%         for (z=)
%             if radial_XZsmall_map(x,z)<Nradial-1
                gradB_X(x,3:sizeZ-2)=(1/12)*(-Btot_ini_XZ_small_map(x+2,3:sizeZ-2)+Btot_ini_XZ_small_map(x-2,3:sizeZ-2))+(2/3)*(Btot_ini_XZ_small_map(x+1,3:sizeZ-2)-Btot_ini_XZ_small_map(x-1,3:sizeZ-2));
                gradB_Z(x,3:sizeZ-2)=(1/12)*(-Btot_ini_XZ_small_map(x,(3:sizeZ-2)+2)+Btot_ini_XZ_small_map(x,(3:sizeZ-2)-2))+(2/3)*(Btot_ini_XZ_small_map(x,(3:sizeZ-2)+1)-Btot_ini_XZ_small_map(x,(3:sizeZ-2)-1));
                dgradB_X(x,3:sizeZ-2)=(1/12)*(-dB_XZ_map(x+2,3:sizeZ-2)+dB_XZ_map(x-2,3:sizeZ-2))+(2/3)*(dB_XZ_map(x+1,3:sizeZ-2)-dB_XZ_map(x-1,3:sizeZ-2));
                dgradB_Z(x,3:sizeZ-2)=(1/12)*(-dB_XZ_map(x,(3:sizeZ-2)+2)+dB_XZ_map(x,(3:sizeZ-2)-2))+(2/3)*(dB_XZ_map(x,(3:sizeZ-2)+1)-dB_XZ_map(x,(3:sizeZ-2)-1));
%             end
%         end
    end

    %remove outer derivative
    gradB_X=gradB_X.*radial_mask;
    gradB_Z=gradB_Z.*radial_mask;
    dgradB_X=dgradB_X.*radial_mask;
    dgradB_Z=dgradB_Z.*radial_mask;

    gradB_X=gradB_X/DX;
    gradB_Z=gradB_Z/DX;
    dgradB_X=dgradB_X/DX;
    dgradB_Z=dgradB_Z/DX;

    gradB_X(isnan(gradB_X))=0;
    gradB_Z(isnan(gradB_Z))=0;
    dgradB_X(isnan(dgradB_X))=0;
    dgradB_Z(isnan(dgradB_Z))=0;
    
    
    vD_X_XZ_map=-Bphi_XZsmall_map.*dgradB_Z;
    vD_Z_XZ_map=Bphi_XZsmall_map.*dgradB_X;
    vD_phi_XZ_map=BpolX_initial_XZsmall_map.*dgradB_Z-BpolZ_initial_XZsmall_map.*dgradB_X+dBX_XZ_map.*gradB_Z-dBZ_XZ_map.*gradB_X;
    
    vD_X_XZ_map=vD_X_XZ_map./(B_XZ_map.^3);
    vD_Z_XZ_map=vD_Z_XZ_map./(B_XZ_map.^3);
    vD_phi_XZ_map=vD_phi_XZ_map./(B_XZ_map.^3);
    
    vD_X_XZ_map(isnan(vD_X_XZ_map))=0;
    vD_Z_XZ_map(isnan(vD_Z_XZ_map))=0;
    vD_phi_XZ_map(isnan(vD_phi_XZ_map))=0;

    vD_X_XZ_map=vD_X_XZ_map.*psi_XZ_map_mask;
    vD_Z_XZ_map=vD_Z_XZ_map.*psi_XZ_map_mask;
    vD_phi_XZ_map=vD_phi_XZ_map.*psi_XZ_map_mask;
    
    
    vD_data=reshape(vD_X_XZ_map(:,:)',sizeX*sizeZ,1);
    vD_X_PR_map=griddata(X_scale_data,Z_scale_data,vD_data,X_PR_map(:,1:size_r),Z_PR_map(:,1:size_r),'cubic');
    vD_data=reshape(vD_Z_XZ_map(:,:)',sizeX*sizeZ,1);
    vD_Z_PR_map=griddata(X_scale_data,Z_scale_data,vD_data,X_PR_map(:,1:size_r),Z_PR_map(:,1:size_r),'cubic');
    vD_data=reshape(vD_phi_XZ_map(:,:)',sizeX*sizeZ,1);
    vD_phi_PR_map=griddata(X_scale_data,Z_scale_data,vD_data,X_PR_map(:,1:size_r),Z_PR_map(:,1:size_r),'cubic');
    vD_X_PR_map(isnan(vD_X_PR_map))=0;
    vD_Z_PR_map(isnan(vD_Z_PR_map))=0;
    vD_phi_PR_map(isnan(vD_phi_PR_map))=0;
    vD_X_PR_map(:,size_r-5:end)=0;
    vD_Z_PR_map(:,size_r-5:end)=0;
    vD_phi_PR_map(:,size_r-5:end)=0;

    vD_X_PR_map_phi(phi_rank,:,:)=vD_X_PR_map;
    vD_Z_PR_map_phi(phi_rank,:,:)=vD_Z_PR_map;
    vD_phi_PR_map_phi(phi_rank,:,:)=vD_phi_PR_map;
    
end