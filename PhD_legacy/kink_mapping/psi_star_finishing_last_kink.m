
%clear all;
%close all;

psi_star_2D_transition=psi_star_2D;


%alpha=atan(2);

disp('----------------');
clear S3 S3_filt;



size_r=round(pos_psi_rx/PRECISE_MESH+4)



pos_rmix=xPsih_zero;
pos_psi_rx=xPsih_zero;
rmix=rPsih_zero;

rx=rPsih_zero;
ksi0=-rx;
sign_ksi0=sign(ksi0);




% F=getframe;
% [im,map] = frame2im(F);    %Return associated image data 
% 
% imwrite(im,filename,'bmp');

disp('...done...');
%pause



time_evol(f_rank)=2*time_evol(f_rank-1)-time_evol(f_rank-2);
psi_star_2D_evol(f_rank,:,:)=psi_star_2D;
rx_evol(f_rank)=rx_evol(f_rank-1);
ksi0_evol(f_rank)=ksi0_evol(f_rank-1);

% f=99;
% f_rank=100;
% time_value=f*0.01
% time_evol(f_rank)=time_value;
% psi_star_2D_evol(f_rank,:,:)=psi_star_2D;
% rx_evol(f_rank)=rmix;
% ksi0_evol(f_rank)=rmix;

f_rank=f_rank+1;
time_value=1
time_evol(f_rank)=2*time_evol(f_rank-1)-time_evol(f_rank-2);
psi_star_2D_evol(f_rank,:,:)=psi_star_2D;
rx_evol(f_rank)=rx_evol(f_rank-1);
ksi0_evol(f_rank)=ksi0_evol(f_rank-1);

time_evol=time_evol/max(time_evol);

psi_star_2D_evol=SIGN_PSIH*psi_star_2D_evol;

Ntime=length(time_evol)

FILENAME=strcat(DATA_FOLDER,'psi_star_evol.mat')
save (FILENAME,'psi_star_2D_evol','ksi0_evol','rx_evol','time_evol','NR','Nomega','Ntime')

FILENAME=strcat(DATA_FOLDER,'tokamak_map_dimensions.mat')
save (FILENAME,'NX','NZ','sup_X','sup_Z','inf_X','inf_Z','mid_X','mid_Z','DX','Rpos_map','Rpos','X_scale','Z_scale','Domega');

psi_star_final_map(:,:)=psi_star_2D_evol(end,:,:);
% psi_star_final_map(:,:)=psi_region3;
psi_star_final_profile=mean(psi_star_final_map(:,1:end-1),2);
% psi_star_final_profile(end+1)=psi_star_final_profile(end);

psi_star_final_profile_result=psi_star_initial;
psi_star_final_profile_result(1:round(xPsih_zero/PRECISE_MESH))=interp1(scale_r(1:length(psi_star_final_profile)),psi_star_final_profile,radial_r_value(1:round(xPsih_zero/PRECISE_MESH)),'cubic');

psi_star_final_profile=psi_star_final_profile_result;
clear psi_final_profile;
psi_final_profile=mean(psiH_PR_map(1:end-1,1:size_r),1)+psi_star_final_profile(1:size_r);
psi_final_profile=[ psi_final_profile mean(psi_PR_map(1:end-1,size_r+1:end),1) ];


FILENAME=strcat(DATA_FOLDER,'psi_profiles_kadomtsev.mat')
save(FILENAME,'-append','psi_star_final_profile','psi_final_profile');

%%
psi_final_PR_map=zeros(Nradial,NP);
for om=1:Nomega
    psi_final_PR_map(size_r+1:end,om)=psi_final_profile(size_r+1:end);
end
t=size(psi_star_2D_evol,1)
for om=1:Nomega
    psi_final_PR_map(1:size_r,om)=interp1(scale_r,psi_star_final_map(:,om),radial_r_value(1:size_r),'cubic');
end
psi_final_PR_map=psi_final_PR_map';

%%
FILENAME=strcat(FINESSE_FOLDER,'finesse_data.mat');
load(FILENAME);

psi_final_data=reshape((psi_final_PR_map(:,:)),NP*Nradial,1);
finesse_data_X=reshape((Rpos_PR_map(:,:)-R0),NP*Nradial,1);
finesse_data_Z=reshape(Z_PR_map(:,:),NP*Nradial,1);

psi_XZ_map_final=griddata(finesse_data_X,finesse_data_Z,psi_final_data,X_scale,Z_scale);
% psi_XZ_map_final=gridfit(finesse_data_X,finesse_data_Z,psi_final_data,X_scale,Z_scale);
psi_XZ_map_final=psi_XZ_map_final';
FILENAME=strcat(DATA_FOLDER,'psi_profiles_kadomtsev.mat')
save(FILENAME,'-append','psi_XZ_map_final');

%%
figure(9)

set(gca,'FontSize',22);
plot(radial_r_value,psi_star_initial,'b--','LineWidth',1.5);
grid on;
hold on;
plot(radial_r_value(1:size_r),psi_star_final_profile(1:size_r),'r','LineWidth',1.5);
xlim([0 r_value_q1_mean*1.6]);
xlabel('r (m)');
ylabel('\Psi_* (T.m^{-2})');
legend('\psi_{*-}','\psi_{*+}');
ylim([-0.006 0.006])
