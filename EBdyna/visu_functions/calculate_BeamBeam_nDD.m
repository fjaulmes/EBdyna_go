% Bosch-Hale coefficients for sigma of DD from beam target
A1_NDD =   53701;
A2_NDD =   330.2700;
A3_NDD =  -0.1271;
A4_NDD =   2.9327e-05;
A5_NDD =  -2.5151e-09;

summary.PTOT_NBI=4;
summary.PSIM_NBI=1;
mu=1.66053906660e-27;
mHe3=3.0160293.*mu;

%%
BINVALS=nan(length(ejected),NB_TS_AVG+1);

output.rhot=interp1(psi_scale_ext,rho_tor_scale_ext,output.psi(:,:));
PROFILES_POP=~ejected;
[NMARKERS,BINVALS(PROFILES_POP,end)] = histc(output.rhot(PROFILES_POP,end),rho_tor_scale_adjust);





% info for previous time stamp (1 before last)
TS=1;
BORN_MARKERS=logical(birth_matrix(end-TS,:)');
n_ej=logical(logical(output.time_stamp_loss>=par.NB_TIME_STAMPS-TS)+(isnan(output.time_stamp_loss)));
POP_TS=logical(n_ej.*BORN_MARKERS);
[NMARKERS,BINVALS(POP_TS,end-TS)] = histc(output.rhot(POP_TS,end-TS),rho_tor_scale_adjust);

% test population used for Beam-Beam statistics
TEST_POP=find(PROFILES_POP);
nDD_BeamBeam_approx=ejected*0;
Ecom_BeamBeam_approx=ejected*0;
Erel_BeamBeam_approx=ejected*0;
sv_BeamBeam_approx=ejected*0;
% increasing stats by 4
nDD_BeamBeam=zeros(length(TEST_POP),4);
Ecom_BeamBeam=zeros(length(TEST_POP),4);
Erel_BeamBeam=zeros(length(TEST_POP),4);
sv_BeamBeam=zeros(length(TEST_POP),4);

v1=squeeze(output.v(:,end,:));
v2=squeeze(output.v(:,end-TS,:));

for ii=1:length(TEST_POP)
    % calculate Beam-Beam interaction for 1 ion!
    if mod(ii,1000)==0
        disp(['processing markers: ' num2str(round(ii)*100/length(TEST_POP)) ' %'])
    end
    BIN_VAL=BINVALS(TEST_POP(ii),end);
    nFast_beam=summary.density_markers(BIN_VAL);
    vR=v1(TEST_POP(ii),1)-v2(BINVALS(POP_TS,end-TS)==BIN_VAL,1);
    vZ=v1(TEST_POP(ii),2)-v2(BINVALS(POP_TS,end-TS)==BIN_VAL,2);
    vphi=v1(TEST_POP(ii),3)-v2(BINVALS(POP_TS,end-TS)==BIN_VAL,3);
    vrel2=(vR.^2+vZ.^2+vphi.^2);
    
    vR=0.5.*(v1(TEST_POP(ii),1)+v2(BINVALS(POP_TS,end-TS)==BIN_VAL,1));
    vZ=0.5.*(v1(TEST_POP(ii),2)+v2(BINVALS(POP_TS,end-TS)==BIN_VAL,2));
    vphi=0.5.*(v1(TEST_POP(ii),3)+v2(BINVALS(POP_TS,end-TS)==BIN_VAL,3));
    vCOM2=(vR.^2+vZ.^2+vphi.^2);

    Erel=0.5.*(mD/eV).*vrel2; 
    vrel=sqrt(vrel2);
    Ecom=0.5.*(mD/eV).*vCOM2; 
    vCOM=sqrt(vCOM2);
    E_sigv=0.5*Erel*1e-3; % D has 2x mass of H
    % cross section for DD
    s_E=(A1_NDD+E_sigv.*(A2_NDD+E_sigv.*(A3_NDD+E_sigv.*(A4_NDD+E_sigv.*A5_NDD))));
    sv_val=s_E./(E_sigv.*exp(31.397./sqrt(E_sigv)));  % in millibarns (10^-31 m^2)
    sv_val=sv_val.*vrel*1e-31;
    Ecom_BeamBeam_approx(TEST_POP(ii))=mean(Ecom);
    Erel_BeamBeam_approx(TEST_POP(ii))=mean(Erel);
    nDD_BeamBeam_approx(TEST_POP(ii))=mean(par.MARKER_WEIGHT.*sv_val.*(summary.PTOT_NBI/summary.PSIM_NBI).*nFast_beam);
    sv_BeamBeam_approx(TEST_POP(ii))=mean(sv_val);
    Ecom_BeamBeam(ii,:)=[mean(Ecom(1:4:end)) mean(Ecom(2:4:end)) mean(Ecom(3:4:end)) mean(Ecom(4:4:end))];
    Erel_BeamBeam(ii,:)=[mean(Erel(1:4:end)) mean(Erel(2:4:end)) mean(Erel(3:4:end)) mean(Erel(4:4:end))];
    sv_BeamBeam(ii,:)=[mean(sv_val(1:4:end)) mean(sv_val(2:4:end)) mean(sv_val(3:4:end)) mean(sv_val(4:4:end))];
    nDD_BeamBeam(ii,:)=(par.MARKER_WEIGHT.*sv_BeamBeam(ii,:).*(summary.PTOT_NBI/summary.PSIM_NBI).*nFast_beam)/4;
end

% vectorizing and cleaning
nDD_BeamBeam=nDD_BeamBeam(:);
sv_BeamBeam=sv_BeamBeam(:);
Erel_BeamBeam=Erel_BeamBeam(:);
Ecom_BeamBeam=Ecom_BeamBeam(:);
nDD_BeamBeam=nDD_BeamBeam(nDD_BeamBeam>0);
sv_BeamBeam=sv_BeamBeam(sv_BeamBeam>0);
Erel_BeamBeam=Erel_BeamBeam(Erel_BeamBeam>0);
Ecom_BeamBeam=Ecom_BeamBeam(Ecom_BeamBeam>0);

output.nDD_BeamBeam=nDD_BeamBeam_approx;
output.Ecom_BeamBeam=Ecom_BeamBeam_approx;
output.Erel_BeamBeam=Erel_BeamBeam_approx;
output.sv_BeamBeam=sv_BeamBeam_approx;

%%
values_markers_ndd=nDD_BeamBeam_approx(PROFILES_POP);

ndd_BB_profile=rho_tor_bins*0;

for ii=1:length(rho_tor_bins)
    ndd_BB_profile(ii)=sum(values_markers_ndd(BINVALS(TEST_POP,end)==ii));
end

ndd_BB_profile=ndd_BB_profile./dv_flux;

summary.ndd_BB_profile=ndd_BB_profile;


%%

% now the spectrum of the emitted neutrons

ndd_BB_Ekin_rel_dist=Ekin_values_ndd*0;
% Erel_BeamBeam(ejected)=nan;
[EKIN_DIST_vals EKIN_DIST_POP]=histc(Erel_BeamBeam,Ekin_bins_ndd);
for ek=1:length(EKIN_DIST_vals)-1
    ndd_BB_Ekin_rel_dist(ek)=sum(nDD_BeamBeam(EKIN_DIST_POP==ek));
end

summary.ndd_BB_Ekin_rel_dist=ndd_BB_Ekin_rel_dist;




ndd_BB_Ekin_COM_dist=Ekin_values_ndd*0;
% Erel_BeamBeam(ejected)=nan;
[EKIN_DIST_vals EKIN_DIST_POP]=histc(Ecom_BeamBeam,Ekin_bins_ndd);
for ek=1:length(EKIN_DIST_vals)-1
    ndd_BB_Ekin_COM_dist(ek)=sum(nDD_BeamBeam(EKIN_DIST_POP==ek));
end

summary.ndd_BB_Ekin_COM_dist=ndd_BB_Ekin_COM_dist;

%%

try
    save(FNAME,'-append','summary')
    save(FNAME,'-append','output')
catch
    warning('could not save summary to file!!!')
end