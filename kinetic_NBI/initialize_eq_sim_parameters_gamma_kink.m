




OMEGA_KINK_REAL=0e4
PERIOD_KINK=2*pi/OMEGA_KINK_REAL;





%simulation options
SAVE_DATA_FILE=1;
USE_LAP_PSI=0;

% only take one frame for the energy calculation
FRAME_NUMBER=11
FRAME_NUMBER_GAMMA=10*(FRAME_NUMBER-1)+1

NB_BINS_OMEGA=32;
DELTA_OMEGA=2*pi/NB_BINS_OMEGA;
bins_omega_bounds=(0:DELTA_OMEGA:2*pi);
bins_omega=bins_omega_bounds(1:end-1)+0.5*DELTA_OMEGA;

NB_BINS_THETA=32;
DELTA_THETA=2*pi/NB_BINS_THETA;
bins_theta_bounds=(0:DELTA_THETA:2*pi);
bins_theta=bins_theta_bounds(1:end-1)+0.5*DELTA_THETA;

DELTA_PSI=8;
bins_psi_bounds=(0:DELTA_PSI:160);
bins_psi=bins_psi_bounds(1:end-1)+0.5*DELTA_PSI;
NB_BINS_PSI=length(bins_psi)



% DISTNAME=strcat('initial_NBI60keV_co_D_distribution',num2str(PROCESS_NUMBER),'.mat');
% %DISTNAME=strcat('initial_alphas_MB_D_distribution',num2str(PROCESS_NUMBER),'.mat')
% %DISTNAME=strcat('initial_flatD_',num2str(PROCESS_VAL),'keV_distribution.mat')
% STATNAME=strcat('initial_NBI60keV_co_precession_stats',num2str(PROCESS_NUMBER),'.mat');
% 
% LOADNAME=strcat('initial_co_NBI60keV_precession',num2str(PROCESS_NUMBER),'.mat');
% %SAVENAME=strcat('initial_MB_D_precession',num2str(PROCESS_NUMBER),'.mat')
% %SAVENAME=strcat('initial_flatD_',num2str(PROCESS_VAL),'keV_precession.mat')


disp('no graphics will be displayed during this simulation')

MAX_GAMMA=2e6
MIN_GAMMA=-2e6


DT_INTERPOLATION_METHOD='quadratic'  
finesse_data_X=reshape((Rpos_PR_map(:,1:size_r)-R0),NP*size_r,1);
finesse_data_Z=reshape(Z_PR_map(:,1:size_r),NP*size_r,1);

minX=min(finesse_data_X);
maxX=max(finesse_data_X);
minZ=min(finesse_data_Z);
maxZ=max(finesse_data_Z);
    
sizeX=round(1.2*(maxX-minX)/DX);
sizeZ=round(1.2*(maxZ-minZ)/DX);
    
mid_X=ceil(0.5*(maxX-minX)/DX)+1;
mid_Z=ceil(0.5*(maxZ-minZ)/DX)+1;
scale_X=DX*((1:sizeX)-mid_X);
scale_Z=DX*((1:sizeZ)-(mid_Z));

finesse_mesh=[finesse_data_X  finesse_data_Z];
[finesse_mesh, IDTRI, J] = unique(finesse_mesh,'last','rows');
finesse_mesh_dtri=DelaunayTri(finesse_mesh);

[XX_small ZZ_small]=meshgrid(scale_X,scale_Z);
X_scale_data=reshape(XX_small,sizeX*sizeZ,1);
Z_scale_data=reshape(ZZ_small,sizeX*sizeZ,1);
XZ_mesh=[X_scale_data Z_scale_data];
XZ_mesh_dtri=DelaunayTri(XZ_mesh);




[residue Xinf]=(min(abs(X_scale-scale_X(1))));
[residue Xsup]=(min(abs(X_scale-scale_X(end))));
[residue Zinf]=(min(abs(Z_scale-scale_Z(1))));
[residue Zsup]=(min(abs(Z_scale-scale_Z(end))));

size_X=length(scale_X);
size_Z=length(scale_Z);

Btot_ini=sqrt(Bphi_XZ_map.^2+BpolX_initial_map.^2+BpolZ_initial_map.^2);

Bphi_XZsmall_map=zeros(size_X,size_Z);
Bphi_XZsmall_map(:,:)=Bphi_XZ_map(Xinf:Xsup,Zinf:Zsup);

BpolX_initial_XZsmall_map=zeros(size_X,size_Z);
BpolX_initial_XZsmall_map(:,:)=BpolX_initial_map(Xinf:Xsup,Zinf:Zsup);

BpolZ_initial_XZsmall_map=zeros(size_X,size_Z);
BpolZ_initial_XZsmall_map(:,:)=BpolZ_initial_map(Xinf:Xsup,Zinf:Zsup);

Rpos_XZsmall_map=zeros(size_X,size_Z);
Rpos_XZsmall_map(:,:)=Rpos_map(Xinf:Xsup,Zinf:Zsup);
	
Btot_XZ_map=zeros(size_X,size_Z);
Btot_XZ_map(:,:)=Btot_ini(Xinf:Xsup,Zinf:Zsup);

psi_XZsmall_map=zeros(size_X,size_Z);
psi_XZsmall_map(:,:)=psi_XZ_map(Xinf:Xsup,Zinf:Zsup);

ion_density_XZ_map=interp1(psi_scale,Ne_profile,psi_XZsmall_map);
ion_density_XZ_map(isnan(ion_density_XZ_map))=0;

radial_XZsmall_map=zeros(size_X,size_Z);
radial_XZsmall_map(:,:)=radial_XZ_map(Xinf:Xsup,Zinf:Zsup);
radial_mask=(radial_XZsmall_map<size_r-5);


psi_q1=interp1(q_initial_profile,1:Nradial,1)
psi_max=round(psi_q1)+3

NR=size_r+10
finesse_data_X_extended=reshape((Rpos_PR_map(:,1:NR)-R0),NP*NR,1);
finesse_data_Z_extended=reshape(Z_PR_map(:,1:NR),NP*NR,1);
finesse_mesh_extended=[finesse_data_X_extended  finesse_data_Z_extended];
[finesse_mesh_extended, IDTRI_EXT, J] = unique(finesse_mesh_extended,'last','rows');
finesse_mesh_extended_dtri=DelaunayTri(finesse_mesh_extended);

psi_max_value=interp1(1:Nradial,psi_scale,psi_max+2)
psi_data=reshape(psi_PR_map(:,1:NR),NB_THETA*NR,1);
psi_XZ_map_ini=dtinterp(finesse_mesh_extended,finesse_mesh_extended_dtri,psi_PR_map(IDTRI_EXT),XX_small,ZZ_small,'quadratic');
psi_XZ_map_ini(isnan(psi_XZ_map_ini))=0;
psi_XZ_map_ini=psi_XZ_map_ini';
psi_XZ_map_mask=psi_XZ_map_ini*0;

for (x=3:sizeX-2)
    for (z=3:sizeZ-2)
        if psi_max_value<0
            if psi_XZ_map_ini(x,z)<psi_max_value
                psi_XZ_map_mask(x,z)=1;
            end
        else
            if psi_XZ_map_ini(x,z)>psi_max_value
                psi_XZ_map_mask(x,z)=1;
            end
        end
    end
end