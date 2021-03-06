
close all;
pause(0.2);

REINIT_ALL_MAPS=0;
REINIT_SIMULATION_DATA=1;

if REINIT_ALL_MAPS==1
    clear all
    load flux_geometry.mat
    load B_fields.mat
    load tokamak_map_dimensions.mat
    load tokamak_PR_map.mat
    mid_X_large=mid_X;
    mid_Z_large=mid_Z;
    calculate_lap_psi;
    calculate_lap_psiH;
    calculate_lap_B;
    save large_equilibrium_XZ_maps.mat lap_B DX X_axis Z_axis mid_X mid_X_large mid_Z mid_Z_large elongation Z_scale X_scale radial_XZ_map theta_XZ_map Rpos_map NZ radial_r_value_flux
    clear all
    load physics_constants.mat
    load tokamak_PR_map.mat
    calculate_pre_collapse_drift_speed_maps;
    REINIT_SIMULATION_DATA=1
end

if REINIT_SIMULATION_DATA==1
    clear all;
    clc
    format compact
    
    load physics_constants.mat
    load XZsmall_fields_tokamak_pre_collapse.mat
    load motions_map_dimensions.mat
%     load psi_star_evol.mat
    load pre_collapse_XZsmall_maps.mat
    load('flux_geometry.mat', 'dr_X_PR_map')
    load('flux_geometry.mat', 'dr_Z_PR_map')

    REINIT_SIMULATION_DATA=1;
   
end


% initialize completely the maps for pre-collapse calculations
FAST_SAWTOOTH=1;
tau_sim=0.1e-4;

TPRECISE=1;
TIME_STAMP_PRECISION=1;
DELTA_TIME=(1e-9)/TPRECISE;
h=DELTA_TIME;
NB_TIME_STEPS=round(0.1*(tau_sim/FAST_SAWTOOTH)/DELTA_TIME)
% one time stamp in one loop
NB_TIME_STAMPS=round(NB_TIME_STEPS/TIME_STAMP_PRECISION);
% ten time steps in one loop
time_scale=(1:NB_TIME_STAMPS)*DELTA_TIME*TIME_STAMP_PRECISION*10;

%simulation options
SAVE_DATA_FILE=1;
USE_LAP_PSI=0;
SAVENAME=strcat('many_particles_eq_Gloca200213_h',num2str(TPRECISE));



if SAVE_DATA_FILE==0
    disp('no data file for this simulation')
end

disp('no graphics will be displayed during this simulation')



%run('calculate_collapse_frozen_drift_speed_maps')


% corrected theta maps for complete interpolation
QNB_THETA=round(0.25*NB_THETA);
HQNB_THETA=round(0.5*QNB_THETA);
theta_low_XZsmall_map=theta_XZsmall_map;
theta_up_XZsmall_map=theta_XZsmall_map;
theta_up_XZsmall_map(find(theta_XZsmall_map<QNB_THETA*DTHETA))=theta_XZsmall_map(find(theta_XZsmall_map<QNB_THETA*DTHETA))+2*pi;
theta_low_XZsmall_map(find(theta_XZsmall_map>(NB_THETA-QNB_THETA-2)*DTHETA))=theta_XZsmall_map(find(theta_XZsmall_map>(NB_THETA-QNB_THETA-2)*DTHETA))-2*pi;

% for correction of Eperp
BMAX=max(max(Btot_XZ_map))
BMIN=min(min(Btot_XZ_map(Btot_XZ_map>1)))


load initial_alphas_distributions.mat

alphas_Bfield=interp2(scale_X,scale_Z,Btot_XZ_map',alphas_pos_x,alphas_pos_z,'*linear');

alphas_mm=(alphas_Ekin-0.5*((mHe/eV)*alphas_vpll.^2))./alphas_Bfield;
alphas_ejected=zeros(Nalphas_simulated,1);
alphas_Eperp=alphas_Bfield.*alphas_mm;



Bavg=mean(Btot_XZ_map(round(0.8*mid_X),:));



alphas_Epll=0.5*(mHe/eV)*alphas_vpll.^2;

%initial position values
alphas_Ekin0=alphas_Ekin;
alphas_mm0=alphas_mm;
X0=alphas_pos_x;
Z0=alphas_pos_z;
phi0=alphas_pos_phi;


alphas_Eperp=alphas_Ekin-alphas_Epll;
alphas_vperp=sqrt(2*alphas_Eperp*eV/mHe);
alphas_rhoL=(mHe/eV)*(1/ZHe)*alphas_vperp./alphas_Bfield;


disp('**************************************************************');

%initializing simulation arrays

alphas_psi=zeros(Nalphas_simulated,1);
alphas_theta=zeros(Nalphas_simulated,1);
bX=zeros(Nalphas_simulated,1);
bZ=zeros(Nalphas_simulated,1);
bphi=zeros(Nalphas_simulated,1);
alphas_psi_value=zeros(Nalphas_simulated,1);
alphas_psi_value_corr=zeros(Nalphas_simulated,1);




% initial field and position
alphas_psi_value=interp2(scale_X,scale_Z,psi_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear');
alphas_psi_value=max(alphas_psi_value,-psi_global);
alphas_psi=interp1(psi_scale,1:257,alphas_psi_value);

alphas_theta=interp2(scale_X,scale_Z,theta_XZsmall_map',alphas_pos_x,alphas_pos_z,'*nearest');



time=0;





% estimating psi and theta values
interpolate_theta_psi_fromXZ;
alphas_pos_phi_wrapped=wrap2pi(alphas_pos_phi);


% direction of the B field at local positions
bX=interp2_XZ(interp_x,interp_z,bX_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
bZ=interp2_XZ(interp_x,interp_z,bZ_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
bphi=sqrt(1-(bX.^2+bZ.^2));

norm_angle=rand(Nalphas_simulated,1)*2*pi;
uX=sqrt(1./(1+(bX./bphi).^2));
uZ=0*uX;
uphi=-(bX./bphi).*uX;
unorm=sqrt(uX.^2+uZ.^2+uphi.^2);

wX=1./sqrt((uX./uphi).^2+1+((1./bZ).*(uX./uphi).*(bphi-bX)).^2);
wZ=(wX./bZ).*((uX./uphi).*bphi-bX);
wphi=-(uX./uphi).*wX;
wnorm=sqrt(wX.^2+wZ.^2+wphi.^2);
wX=wX./wnorm;
wZ=wZ./wnorm;
wphi=wphi./wnorm;

% normal vector N = (cos) u + (sin) w
% that N.b=0 precisely
NX=(cos(norm_angle).*uX+sin(norm_angle).*wX);
NZ=((cos(norm_angle).*uZ+sin(norm_angle).*wZ));
Nphi=(cos(norm_angle).*uphi+sin(norm_angle).*wphi);


v0_X=alphas_vpll.*bX+alphas_vperp.*NX;
v0_Z=alphas_vpll.*bZ+alphas_vperp.*NZ;
v0_phi=alphas_vpll.*bphi+alphas_vperp.*Nphi;

v_X=v0_X;
v_Z=v0_Z;
v_phi=v0_phi;

alphas_Omega=(ZHe*eV/mHe)*alphas_Bfield;
epsilon=0.5*h*alphas_Omega;
alphas_determinant=1+(epsilon.^2);

% previous mid point speed values
v_X_prev=v0_X+epsilon.*(-bphi.*v0_Z+bZ.*v0_phi);
v_Z_prev=v0_Z+epsilon.*(bphi.*v0_X-bX.*v0_phi);
v_phi_prev=v0_phi+epsilon.*(-bZ.*v0_X+bX.*v0_Z);

v_phi_prev_prev=v_phi_prev;
% v_phi_prev=v_phi;

%initialize Bfield properly
alphas_mm_part=alphas_mm;
% update_GT_3D_collapse;

interpolate_theta_psi_fromXZ;

% amplitude of the B field and potential at half time step local positions
alphas_Bfield=interp2_XZ(interp_x,interp_z,Btot_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);

% Canonical angular momentum evolution
alphas_omega=wrap2pi(alphas_theta-alphas_pos_phi_wrapped);

alphas_psi_value_corr=alphas_psi_value;

alphas_pphi0=(mHe/eV)*alphas_Rpos.*v_phi-(ZHe)*alphas_psi_value_corr;



%trapping parameter
radial_pos=(a/257)*interp2(scale_X,scale_Z,radial_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear');
Bavg=mean(alphas_Bfield);
alphas_kappa=sqrt((alphas_Ekin*R0+Bavg*alphas_mm.*(radial_pos-R0))./(2*alphas_mm.*radial_pos*Bavg));
alphas_lambda=Bavg*alphas_mm./alphas_Ekin;

alphas_Eperp=alphas_Ekin-alphas_Epll;
alphas_vperp=sqrt(2*alphas_Eperp*eV/mHe);
alphas_rhoL=(mHe/eV)*(1/ZHe)*alphas_vperp./alphas_Bfield;




alphas_Ekin_half=alphas_Ekin;
alphas_Ekin_prev=alphas_Ekin;
alphas_pphi0_prev=alphas_pphi0;



interpolate_theta_psi_fromXZ;


    
% benchmarking time between two data recordings
disp('**************************************************************');
tic

time=0;


for time_step=1:NB_TIME_STEPS
   
    time_step_integration_GT_eq;
    time_step_integration_GT_eq;

    time_step_integration_GT_eq;
    time_step_integration_GT_eq;
    
    time_step_integration_GT_eq;
    time_step_integration_GT_eq;
    
    time_step_integration_GT_eq;
    time_step_integration_GT_eq;

    time_step_integration_GT_eq;
%     time_step_integration_GT_collapse;
    time_step_integration_GT_eq;

    
%    estimate_half_time_step_positions_eq_G;

%    adapt_speed_pphi_G;
%    v_phi=v_phi_tilde;


    if (mod(time_step,TIME_STAMP_PRECISION)==0)
        % making all values up to time step
%         estimate_half_time_step_positions_eq_G;
        adapt_speed_Ekin_G;
%         v_X_prev=v_X;
%         v_Z_prev=v_Z;
% %         v_phi=v_phi_tilde;
%         v_phi_prev_prev=v_phi_prev;
%         v_phi_prev=v_phi;
        
        time_stamp=ceil((time_step)/TIME_STAMP_PRECISION);
        alphas_omega=wrap2pi(alphas_theta-alphas_pos_phi_wrapped);
        alphas_psi_value_corr=alphas_psi_value;
        
        % rescaling the value of pphi0 to our knowledge of the toroidal
        % speed ans psi value
        % alphas_pphi0=(mHe/eV)*alphas_Rpos.*v_phi-(ZHe)*alphas_psi_value_corr;

        % amplitude of the B field and potential at time step local positions

        alphas_Bfield=interp2_XZ(interp_x,interp_z,Btot_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);

		% direction of the B field at local positions

		bX=interp2_XZ(interp_x,interp_z,bX_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
		bZ=interp2_XZ(interp_x,interp_z,bZ_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);

		bphi=sqrt(1-(bX.^2+bZ.^2));

        alphas_vpll=v_X_step.*bX+v_Z_step.*bZ+v_phi_step.*bphi;
        alphas_Epll=0.5*(mHe/eV)*alphas_vpll.^2;
        alphas_Eperp=max(alphas_Ekin-alphas_Epll,0);
        outcast=find(alphas_psi>simulation_size_r+22);
        alphas_vperp=sqrt(2*alphas_Eperp*eV/mHe);
        alphas_rhoL=(mHe/eV)*(1/ZHe)*alphas_vperp./alphas_Bfield;
        if (~isempty(outcast))
            recast=randi(Nalphas_simulated,size(outcast,1),1);
            reposition_lost_particles_3DG;
            alphas_ejected(outcast)=1;
            disp(strcat('number of ejected particles =  ',num2str(size(outcast,1))));
        end
        Xpos_output(time_stamp,:)=alphas_pos_x;
        Zpos_output(time_stamp,:)=alphas_pos_z;
        phipos_output(time_stamp,:)=alphas_pos_phi;
        psipos_output(time_stamp,:)=alphas_psi;
        vparallel_output(time_stamp,:)=alphas_vpll;
        Eperp_output(time_stamp,:)=alphas_Eperp;
        theta_output(time_stamp,:)=alphas_theta;
        Ekin_output(time_stamp,:)=alphas_Ekin;
        vphi_output(time_stamp,:)=v_phi_step;
        pphi_output(time_stamp,:)=alphas_pphi0;
        psi_value_output(time_stamp,:)=alphas_psi_value_corr;
        rhoL_output(time_stamp,:)=alphas_rhoL;
        Bfield_output(time_stamp,:)=alphas_Bfield;

    end

    if (SAVE_DATA_FILE==1) && (mod(time_step,250)==0)
        toc
        tic
        clear outcast recast
        % correct the particles that are out of simulation domain
        % by giving them a position randomly in the initial distribution
        outcast=find(alphas_psi>simulation_size_r+22);
        recast=randi(Nalphas_simulated,size(outcast,1),1);
        if (~isempty(outcast))
            reposition_lost_particles_3DG;
            alphas_ejected(outcast)=1;
        end

        disp('------------------------------------------------')
        disp(strcat('time_step # ',num2str(time_step),' ; time = ',num2str(time)));
        disp(strcat('time_stamp # ',num2str(time_stamp),' ; time = ',num2str(time)));
        disp(strcat('number of repositionned particles =  ',num2str(size(outcast,1))));
        disp(strcat('total number of ejected particles =  ',num2str(size(find(alphas_ejected),1))));

        filename=strcat('alphas_pre_Grecord_',num2str(time_step),'.mat');
        save(filename,'time','alphas_pos_x','alphas_pos_z','alphas_pos_phi','alphas_vpll','alphas_mm_part','alphas_psi','alphas_ejected');
        disp('------------------------------------------------')
    end



end

toc




if SAVE_DATA_FILE==1
    save ('initial_Galphas_pre_collapse.mat','alphas_pos_x','alphas_pos_z','alphas_pos_phi','alphas_vpll','alphas_mm_part','alphas_psi','alphas_ejected');
end
disp('**************************************************************');
disp('done');
Nalphas_simulated
