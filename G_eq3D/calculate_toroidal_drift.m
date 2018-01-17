
alphas_ejected_collapse=alphas_ejected;

load('initialG_alphas_precession.mat', 'alphas_ejected')
load('initialG_alphas_precession.mat', 'posX_output')
load('initialG_alphas_precession.mat', 'posZ_output')

alphas_ejected=alphas_ejected_collapse+alphas_ejected;
alphas_ejected(alphas_ejected>1)=1;

calculate_vD_XZ_map;
vDphi_output=interp2(scale_X,scale_Z,vD_phi_XZ_map',posX_output,posZ_output,'*linear');
r_output=sqrt((posX_output-X_axis).^2+posZ_output.^2);
PART_POP=(~alphas_ejected);
vDphi_output_corr=vDphi_output;

end_ts=size(posX_output,1);
for ts=1:end_ts
    vDphi_output_corr(ts,:)=(vDphi_output(ts,:)').*(2*alphas_Ekin-Bavg*alphas_mm)/(ZHe*Bavg);
end

vDphi_avg=mean(vDphi_output_corr(:,:),1);
omega_precess_avg=-vDphi_avg/R0;

r_avg=mean(r_output(:,:),1)';
alphas_kappa=sqrt((alphas_Ekin*R0+Bavg*alphas_mm.*(r_avg-R0))./(2*alphas_mm.*r_avg*Bavg));
