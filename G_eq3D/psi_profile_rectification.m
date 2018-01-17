close all

NR_CORR=size_r+5;
psi_star_final_profile=mean(psi_PR_map(1:end-1,1:NR_CORR+1),1);
psi_star_final_profile=psi_star_final_profile-mean(psiH_PR_map(1:end-1,1:NR_CORR+1),1);
psi_final_profile_recalc=psi_star_final_profile;
correction_psi_dr=(1:NR_CORR)/NR_CORR;
correction_psi_dr=correction_psi_dr.^8;
psi_final_profile_der=zeros(size_r,1);
psi_der_recalc=zeros(size_r,1);
dr=1/Nradial;

for (r_pos=2:NR_CORR)
    psi_final_profile_der(r_pos)=(psi_star_final_profile(r_pos+1)-psi_star_final_profile(r_pos-1))/(2*dr);
end

correction_psi_dr=correction_psi_dr*psi_final_profile_der(NR_CORR);

psi_der_recalc=psi_final_profile_der;

correction_psi_dr=correction_psi_dr';
correction_psi=(1:NR_CORR)/(NR_CORR);
correction_psi=(correction_psi.^4)';
psi_der_recalc=psi_final_profile_der.*(1-correction_psi)+correction_psi.*correction_psi_dr;

% for (r_pos=size_r+1:-1:1)
%     psi_der_recalc(r_pos)=0.2*psi_final_profile_der(r_pos)+0.8*correction_psi_dr(r_pos);
% end


figure(1)
plot(psi_final_profile_der)
hold on
plot(psi_der_recalc,'r')
plot(correction_psi_dr,'k')


psi_star_profile=psi_star_final_profile;
clear psi_star_final_profile
psi_star_final_profile=psi_star_profile(1:NR_CORR);

psi_final_profile_recalc=psi_star_final_profile;
for (r_pos=NR_CORR-1:-1:2)
    psi_final_profile_recalc(r_pos)=psi_final_profile_recalc(r_pos+1)-psi_der_recalc(r_pos)*dr;
end
psi_final_profile_recalc(end)=psi_star_final_profile(end);

DPSI=psi_final_profile_recalc(end-1)-psi_star_final_profile(end-1);
DPSI_GLOBAL_INI=psi_star_final_profile(end-1)-psi_star_final_profile(2);

psi_final_profile_recalc=psi_final_profile_recalc-DPSI;
DPSI_GLOBAL_RECALC=psi_final_profile_recalc(end-1)-psi_final_profile_recalc(2);

psi_final_profile_recalc=(DPSI_GLOBAL_INI/DPSI_GLOBAL_RECALC)*psi_final_profile_recalc;
psi_final_profile_recalc(end)=psi_star_final_profile(end);


correction_psi=(1:NR_CORR)/NR_CORR;
correction_psi=correction_psi.^8;
correction_psi=1-correction_psi;
psi_final_profile_recalc=correction_psi.*psi_star_final_profile+(1-correction_psi).*psi_final_profile_recalc;

figure(2)
plot(psi_star_final_profile);
hold on
plot(psi_final_profile_recalc,'r');

psi_final_profile_recalc=psi_final_profile_recalc+mean(psiH_PR_map(1:end-1,1:NR_CORR),1);

psi_PR_map_recalc=psi_PR_map;
for (r_pos=2:NR_CORR)
    psi_PR_map_recalc(:,r_pos)=ones(257,1)*psi_final_profile_recalc(r_pos);
end
psi_PR_map=psi_PR_map_recalc;