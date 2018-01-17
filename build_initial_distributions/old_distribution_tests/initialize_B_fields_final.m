
Btot_PR_map=sqrt(Bpol_PR_map.^2+Btor_PR_map.^2);

finesse_data_X=reshape((Rpos_PR_map(:,:)-R0),NP*Nradial,1);
finesse_data_Z=reshape(Z_PR_map(:,:),NP*Nradial,1);
[XX,ZZ] = meshgrid(X_scale,Z_scale);

frame_rank=1001;

disp('loading new frame rank#');
disp(frame_rank);
f=frame_rank;

if (f<10)
    frame_name='00';
elseif (f<100)
    frame_name='0';
else
    frame_name='';
end

filename='..\B_maps\B0';
frame_name=strcat(frame_name,num2str(f));
filename=strcat(filename,frame_name,'.mat')
load(filename);

NB_THETA=NP
update_psi_star_frame_rank;


psi_PR_map(:,1:size_r-1)=psiH_PR_map(:,1:size_r-1)+psi_star_PR_map(:,1:size_r-1);

psi_scale=mean(psi_PR_map(1:256,:));

psi_data=reshape(psi_PR_map(:,:),NB_THETA*257,1);
psi_XZ_map=griddata(finesse_data_X,finesse_data_Z,psi_data,XX,ZZ,'cubic');
psi_XZ_map(isnan(psi_XZ_map))=0;
psi_XZ_map=psi_XZ_map';

psi2D=psi_XZ_map;
gpsi_R=zeros(NZ,NZ);
gpsi_Z=zeros(NZ,NZ);

for (x=3:NZ-2)
    for (z=3:NZ-2)
        gpsi_R(x,z)=(1/12)*(-psi2D(x+2,z)+psi2D(x-2,z))+(2/3)*(psi2D(x+1,z)-psi2D(x-1,z));
        gpsi_Z(x,z)=(1/12)*(-psi2D(x,z+2)+psi2D(x,z-2))+(2/3)*(psi2D(x,z+1)-psi2D(x,z-1));
%         gpsi_Z(x,z)=0.5*(psi2D(x,z+1)-psi2D(x,z-1));
    end
end
gpsi_R=gpsi_R/DX;
gpsi_Z=gpsi_Z/DX;

BX_XZ_map=-gpsi_Z./Rpos_map;
BZ_XZ_map=gpsi_R./Rpos_map;
Btot_XZ_map=sqrt(BX_XZ_map.^2+BZ_XZ_map.^2+Bphi_XZ_map.^2);


bX_PR_map(:,:)=bX_map_phi(1,:,:);
bZ_PR_map(:,:)=bZ_map_phi(1,:,:);
B_PR_map_partial(:,:)=Btot_map_phi(1,:,:);

BX_PR_map_partial=bX_PR_map.*B_PR_map_partial;
BZ_PR_map_partial=bZ_PR_map.*B_PR_map_partial;

BX_PR_map(:,1:size_r-1)=BX_PR_map_partial(:,1:size_r-1);
BZ_PR_map(:,1:size_r-1)=BZ_PR_map_partial(:,1:size_r-1);
Btot_PR_map(:,1:size_r-1)=B_PR_map_partial(:,1:size_r-1);

B_data=reshape(BX_PR_map(:,:),NB_THETA*257,1);
BX_XZ_map=griddata(finesse_data_X,finesse_data_Z,B_data,XX,ZZ,'cubic');
BX_XZ_map(isnan(BX_XZ_map))=0;
BX_XZ_map=BX_XZ_map';

B_data=reshape(BZ_PR_map(:,:),NB_THETA*257,1);
BZ_XZ_map=griddata(finesse_data_X,finesse_data_Z,B_data,XX,ZZ,'cubic');
BZ_XZ_map(isnan(BZ_XZ_map))=0;
BZ_XZ_map=BZ_XZ_map';

B_data=reshape(Btot_PR_map(:,:),NB_THETA*257,1);
Btot_XZ_map=griddata(finesse_data_X,finesse_data_Z,B_data,XX,ZZ,'cubic');
Btot_XZ_map(isnan(Btot_XZ_map))=0;
Btot_XZ_map=Btot_XZ_map';

