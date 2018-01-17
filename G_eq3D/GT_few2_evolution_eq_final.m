
close all;
pause(0.2);

REINIT_ALL_MAPS=1;
REINIT_SIMULATION_DATA=1;

if REINIT_ALL_MAPS==1
    clear all
    initialize_folder_names;
    filename=strcat(DATA_FOLDER,'physics_constants.mat');
    load(filename);
    filename=strcat(DATA_FOLDER,'tokamak_map_dimensions.mat');
    load(filename);

    mid_X_large=mid_X;
    mid_Z_large=mid_Z;    
    NB_PHI=129;
    DPHI=2*pi/(NB_PHI-1);
    NB_PHI_DATA_HALF=round(0.5*(NB_PHI-1));
    
    calculate_post_collapse_drift_speed_maps;
    REINIT_SIMULATION_DATA=1
end

if REINIT_SIMULATION_DATA==1
    clear all;
    clc
    format compact
    initialize_folder_names;
    
    filename=strcat(DATA_FOLDER,'physics_constants.mat');
    load(filename);
    filename=strcat(DATA_FOLDER,'XZsmall_fields_tokamak_post_collapse.mat');
    load(filename);
    filename=strcat(DATA_FOLDER,'motions_map_dimensions.mat');
    load(filename);
%     filename=strcat(DATA_FOLDER,'psi_star_evol.mat');
%     load(filename, 'psi_star_2D_evol_interp');
%     psi_star_map_phi_evol=interp1(1:1001,psi_star_2D_evol_interp,1:101,'cubic');

    filename=strcat(DATA_FOLDER,'flux_geometry.mat');
    load(filename, 'dr_X_PR_map');
    load(filename, 'dr_Z_PR_map');

    REINIT_SIMULATION_DATA=1;
   
end

TOROIDAL_FIELD_DIRECTION=sign(mean(mean(Bphi_XZsmall_map)))


% initialize completely the maps for pre-collapse calculations
FAST_SAWTOOTH=1;
tau_cr=4e-4;

TPRECISE=2;
TIME_STAMP_PRECISION=1;
DELTA_TIME=(1e-9)/TPRECISE;
h=DELTA_TIME;
NB_TIME_STEPS=round(0.01*(tau_cr/FAST_SAWTOOTH)/DELTA_TIME)
% one time stamp in one loop
NB_TIME_STAMPS=round(NB_TIME_STEPS/TIME_STAMP_PRECISION);
% ten time steps in one loop
time_scale=(1:NB_TIME_STAMPS)*DELTA_TIME*TIME_STAMP_PRECISION*10;

%simulation options
SAVE_DATA_FILE=1;
DISPLAY_OUTPUTS=0;
EVOLVE_VEXB=1;
USE_LAP_PSI=0;
SAVENAME=strcat('few_particles_nocorr_eq_final_G250413_h',num2str(TPRECISE))



if SAVE_DATA_FILE==0
    disp('no data file for this simulation')
end

if DISPLAY_OUTPUTS==0
    disp('no graphics will be displayed during this simulation')
else
    filename=strcat(DATA_FOLDER,'tokamak_PR_map.mat');
    load(filename);
    [XXsmall ZZsmall]=meshgrid(scale_X,scale_Z);
    finesse_data_X=reshape((Rpos_PR_map(:,1:size_r)-R0),NP*size_r,1);
    finesse_data_Z=reshape(Z_PR_map(:,1:size_r),NP*size_r,1);
end


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


%Nalphas_simulated=128;
Nalphas_simulated=64;

alphas_pos_z=[-0.38 -0.3 -0.05 -0.02 0.02 0.05 0.3 0.38 -0.38 -0.3 -0.12 -0.02 0.02 0.12 0.3 0.38 0.22 0.22 0.22 0.22 0.22 0.22 0.22 0.22 -0.02 0.02  0 0 -0.2 -0.2 0.2 0.2];
alphas_pos_z=[alphas_pos_z alphas_pos_z];
alphas_pos_x=alphas_pos_z*0;

%alphas_pos_x=[alphas_pos_x alphas_pos_z];
%alphas_pos_z=[alphas_pos_z alphas_pos_z*0];

alphas_pos_z=alphas_pos_z' + Z_axis;
alphas_pos_x=alphas_pos_x' + Raxis-R0 - 0.01;

% alphas_pos_x=[-0.26 -0.2 -0.12 -0.01 0.01 0.12 0.2 0.26 -0.26 -0.2 -0.12 -0.01 0.01 0.12 0.2 0.26 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 -0.01 0.01  0 0 -0.2 -0.2 0.2 0.2];
% alphas_pos_x=[alphas_pos_x alphas_pos_x]';
% alphas_pos_x=alphas_pos_x+Raxis-R0;
% alphas_pos_z=0*alphas_pos_x + Z_axis;

alphas_pos_phi=0*alphas_pos_x+1.25;
alphas_Bfield=interp2(scale_X,scale_Z,Btot_XZ_map',alphas_pos_x,alphas_pos_z,'*linear');

alphas_Ekin=2*[ 800 800 800 800 800 800 800 800 600 600 600 600 600 600 600 600 200 300 400 500 600 800 8 22 8 22 8 22 8 22 8 22]*1e3;  % keV
%alphas_Ekin=[alphas_Ekin 2*alphas_Ekin alphas_Ekin 2*alphas_Ekin]
alphas_Ekin=[alphas_Ekin 2*alphas_Ekin ]
alphas_Ekin=alphas_Ekin';


alphas_vpll=(ones(Nalphas_simulated,1)*0.1);

alphas_vpll=alphas_vpll.*sqrt(alphas_Ekin*2*eV/mHe);
alphas_vpll(33:64)=2*alphas_vpll(33:64);
alphas_mm=(alphas_Ekin-0.5*((mHe/eV)*alphas_vpll.^2))./alphas_Bfield;
alphas_ejected=zeros(Nalphas_simulated,1);
alphas_Eperp=alphas_Bfield.*alphas_mm;
Ekin_adapt=zeros(Nalphas_simulated,1);



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
Xpos_output=zeros(NB_TIME_STAMPS,Nalphas_simulated);
Zpos_output=zeros(NB_TIME_STAMPS,Nalphas_simulated);
phipos_output=zeros(NB_TIME_STAMPS,Nalphas_simulated);
psipos_output=zeros(NB_TIME_STAMPS,Nalphas_simulated);
vD_output=zeros(NB_TIME_STAMPS,Nalphas_simulated);
vparallel_output=zeros(NB_TIME_STAMPS,Nalphas_simulated);
Eperp_output=zeros(NB_TIME_STAMPS,Nalphas_simulated);
theta_output=zeros(NB_TIME_STAMPS,Nalphas_simulated);
Ekin_output=zeros(NB_TIME_STAMPS,Nalphas_simulated);
vphi_output=zeros(NB_TIME_STAMPS,Nalphas_simulated);
pphi_output=zeros(NB_TIME_STAMPS,Nalphas_simulated);
psi_value_output=zeros(NB_TIME_STAMPS,Nalphas_simulated);
rhoL_output=zeros(NB_TIME_STAMPS,Nalphas_simulated);
Bfield_output=zeros(NB_TIME_STAMPS,Nalphas_simulated);




alphas_psi=zeros(Nalphas_simulated,1);
alphas_theta=zeros(Nalphas_simulated,1);
bX=zeros(Nalphas_simulated,1);
bZ=zeros(Nalphas_simulated,1);
bphi=zeros(Nalphas_simulated,1);
alphas_psi_value=zeros(Nalphas_simulated,1);
alphas_psi_value_corr=zeros(Nalphas_simulated,1);



% intermediate calculations arrays
delta_E=zeros(Nalphas_simulated,1);
alphas_vpll_int=alphas_vpll;
vparallel_tilde=alphas_vpll;


% initial field and position
alphas_psi_value=interp2(scale_X,scale_Z,psi_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear');
alphas_psi_value=max(alphas_psi_value,-psi_global);
alphas_psi=interp1(psi_scale,1:NB_PSI,alphas_psi_value);

alphas_theta=interp2(scale_X,scale_Z,theta_XZsmall_map',alphas_pos_x,alphas_pos_z,'*nearest');



time=0;





% estimating psi and theta values
interpolate_theta_psi_fromXZ;
alphas_pos_phi_wrapped=wrap2pi(alphas_pos_phi);


% direction of the B field at local positions
bX=interp2_XZ(interp_x,interp_z,bX_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
bZ=interp2_XZ(interp_x,interp_z,bZ_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
bphi=TOROIDAL_FIELD_DIRECTION*sqrt(1-(bX.^2+bZ.^2));

% norm_angle=rand(Nalphas_simulated,1)*2*pi;
norm_angle=ones(Nalphas_simulated,1)*0.5*pi;
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
alphas_pphi0_half=alphas_pphi0;



%trapping parameter
radial_pos=(a/NB_PSI)*interp2(scale_X,scale_Z,radial_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear');
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

    
%     estimate_half_time_step_positions_eq_G;
% 
%     adapt_speed_pphi_G;
%    v_phi=v_phi_tilde;


    if (mod(time_step,TIME_STAMP_PRECISION)==0)
        % making all values up to time step
%         estimate_half_time_step_positions_eq_G;


        adapt_speed_Ekin_G;
%         v_phi=0.5*(v_phi+v_phi_tilde);
%         adapt_speed_Ekin_G;

        v_X_step=(0.5*v_X+0.5*v_X_prev);
        v_Z_step=(0.5*v_Z+0.5*v_Z_prev);
        v_phi_step=(0.5*v_phi+0.5*v_phi_prev);
%         v_X_step=(1.5*v_X-0.5*v_X_prev);
%         v_Z_step=(1.5*v_Z-0.5*v_Z_prev);
%         v_phi_step=(1.5*v_phi-0.5*v_phi_prev);
        
        v_X_prev=v_X;
        v_Z_prev=v_Z;
        v_phi_prev=v_phi;
        
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

		bphi=TOROIDAL_FIELD_DIRECTION*sqrt(1-(bX.^2+bZ.^2));

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
        Ekin_output(time_stamp,:)=0.5*(mHe/eV)*(v_X_step.^2+v_Z_step.^2+v_phi_step.^2);
        vphi_output(time_stamp,:)=v_phi_step;
        pphi_output(time_stamp,:)=alphas_pphi0;
        psi_value_output(time_stamp,:)=alphas_psi_value_corr;
        rhoL_output(time_stamp,:)=alphas_rhoL;
        Bfield_output(time_stamp,:)=alphas_Bfield;

    end

    if (SAVE_DATA_FILE==1) && (mod(time_step,1000)==0)
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

        save (SAVENAME ,'Bfield_output','Xpos_output','Zpos_output','phipos_output','theta_output',...
            'psi_value_output','vparallel_output','Ekin_output','Eperp_output','rhoL_output','psipos_output',...
            'vphi_output','pphi_output','time_scale','alphas_ejected');
        disp('------------------------------------------------')
    end



end

toc

Xpos_output=Xpos_output';
Zpos_output=Zpos_output';
phipos_output=phipos_output';
vparallel_output=vparallel_output';
Eperp_output=Eperp_output';
theta_output=theta_output';
psipos_output=psipos_output';
Ekin_output=Ekin_output';
vphi_output=vphi_output';
pphi_output=pphi_output';
psi_value_output=psi_value_output';
rhoL_output=rhoL_output';
Bfield_output=Bfield_output';


Xpos_outputG=Xpos_output;
Zpos_outputG=Zpos_output;
phipos_outputG=phipos_output;
vparallel_outputG=vparallel_output;
Eperp_outputG=Eperp_output;
theta_outputG=theta_output;
psipos_outputG=psipos_output;
Ekin_outputG=Ekin_output;
vphi_outputG=vphi_output;
pphi_outputG=pphi_output;
psi_value_outputG=psi_value_output;
rhoL_outputG=rhoL_output;
Bfield_outputG=Bfield_output;


alphas_ejected_G=alphas_ejected;
time_scale_G=time_scale;

if SAVE_DATA_FILE==1
     save (SAVENAME ,'Bfield_outputG','Xpos_outputG','Zpos_outputG','phipos_outputG','theta_outputG',...
         'psi_value_outputG','vparallel_outputG','Ekin_outputG','Eperp_outputG','rhoL_outputG','psipos_outputG',...
         'vphi_outputG','pphi_outputG','time_scale_G','alphas_ejected_G');
end
disp('**************************************************************');
disp('done');
Nalphas_simulated


