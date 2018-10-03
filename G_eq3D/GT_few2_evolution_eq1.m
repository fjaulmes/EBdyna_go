%GT_few_evolution_eq || Simulate ~64 (fast) ions in sawtooth pre-crash equilibrium
%   Simulate ~64 (fast) ions in sawtooth pre-crash equilibrium
% First we initialize maps, define parameters (e.g. simulation time) and
% interpolate the magnetic fields in equilibrium. Next the sawtooth crash
% is simulated using a leap-frog integration scheme.


%% Cleaning
close all;
drawnow;

%% Initialization maps
REINIT_ALL_MAPS=0;
REINIT_SIMULATION_DATA=1;

%% Reinitialization
% Maps generation
if REINIT_ALL_MAPS==1
    clear all
    initialize_folder_names;    %Find the folders and load:
    filename=strcat(DATA_FOLDER,'physics_constants.mat');
    load(filename);
    filename=strcat(DATA_FOLDER,'tokamak_map_dimensions.mat');
    load(filename);

    mid_X_large=mid_X;
    mid_Z_large=mid_Z;    
    NB_PHI=129; % Number of toroidal positions incl. 0 and 2pi
    DPHI=2*pi/(NB_PHI-1);   %Angle between toroidal positions
    NB_PHI_DATA_HALF=round(0.5*(NB_PHI-1)); %Number of half the toroidal positions
    calculate_pre_collapse_drift_speed_maps; %Calculations and storage of maps
    clear all;
    REINIT_SIMULATION_DATA=1;
end

% Simulation data
if REINIT_SIMULATION_DATA==1
    clear all;
    clc
    format compact
    initialize_folder_names;
    
    filename=strcat(DATA_FOLDER,'physics_constants.mat');
    load(filename);
    filename=strcat(DATA_FOLDER,'XZsmall_fields_tokamak_pre_collapse.mat');
    load(filename);
    filename=strcat(DATA_FOLDER,'motions_map_dimensions.mat');
    load(filename);
    
    APPLY_RMP_FIELD=1
    POINCARE_SIMULATION=1
    CALCULATE_TRUE_PPHI=1

    if APPLY_RMP_FIELD==1
        
        load('BS_AUG_theta_psi_phi_2016-03-08.mat');
        BR_RMP=BR;
        BZ_RMP=BZ;
        Bphi_RMP=Bphi;
        if CALCULATE_TRUE_PPHI==1
            dAphi_dphi=permute(dAphi_dphi,[3 1 2]);
        end
        
        BR_RMP=permute(BR_RMP,[3 1 2]);
        BZ_RMP=permute(BZ_RMP,[3 1 2]);
        Bphi_RMP=permute(Bphi_RMP,[3 1 2]);
        
        %BR_RMP=BR_RMP*0;
        %BZ_RMP=BR_RMP*0;
        %Bphi_RMP=Bphi_RMP*0;
        
        NB_PHI=size(BR_RMP,1)       % Number of toroidal positions incl. 0 and 2pi
        NB_THETA_RMP=size(BR_RMP,2) % Number of poloidal positions incl. 0 and 2pi
        NB_PSI=size(BR_RMP,3)       % Number of radial positions
        DPSIVAL=round(Nradial/NB_PSI)  % estimating radial interval for 3D fields
		
        DTHETA=2*pi/(NB_THETA_RMP-1)
        DPHI=2*pi/(NB_PHI-1);   %Angle between toroidal positions

        scale_phi=2*pi*((0:(NB_PHI-1))/(NB_PHI-1));
        scale_theta=2*pi*((0:(NB_THETA_RMP-1))/(NB_THETA_RMP-1));
        scale_psi=1:DPSIVAL:Nradial;
        NB_PHI_DATA_HALF=round(0.5*(NB_PHI-1)); %Number of half the toroidal positions (not used here)

        clear BR BZ Bphi;
       
    end
    simulation_size_r=Nradial-1

    REINIT_SIMULATION_DATA=1;
end

%% Simulation parameters
TOROIDAL_FIELD_DIRECTION=sign(mean(mean(Bphi_XZsmall_map)))


% initialize completely the maps for pre-collapse calculations
FAST_SAWTOOTH=1;        % 0.5 -> 8 micro s. | 1 -> 4 micro s.
ZHe=1                   % Deuterium simulation
mHe=mD


if POINCARE_SIMULATION==0
    tau_cr=10e-4;            % Simulation time (typical crash time)
    TPRECISE=8;             % Precision of timesteps in a.u. larger means smaller time steps
    TIME_STAMP_PRECISION=2*TPRECISE; % stamping of time stamps
    DELTA_TIME=(1e-9)/TPRECISE; % size of time stamps (1e-9 timestamp)
    RECORD_TIME_STAMP=10000
else
    tau_cr=0.005;                % Long simulation necessary for field line tracing
    TPRECISE=4;                  % Precision of timesteps in a.u. larger means smaller time steps
    TIME_STAMP_PRECISION=10*TPRECISE;     % stamping of time steps
    DELTA_TIME=(1e-9)/TPRECISE;       % size of time stamps (1e-6 timestamp)
    RECORD_TIME_STAMP=10000
end

% coefficients used to correct the integration 
% when potentials are known
ECOEF=0.0/TPRECISE
PPHI_CORR_FACTOR=0.0000/TPRECISE
PPHI_CORR_INTEG_FACTOR=0.0000


h=DELTA_TIME;
NB_TIME_STEPS=round(0.1*(tau_cr/FAST_SAWTOOTH)/DELTA_TIME)  % number of time steps within a loop
% one time stamp in one loop
NB_TIME_STAMPS=round(NB_TIME_STEPS/TIME_STAMP_PRECISION);   % number of time stamps during simulation
% ten time steps in one loop
time_scale=(1:NB_TIME_STAMPS)*DELTA_TIME*TIME_STAMP_PRECISION*10;

%% Simulation options
SAVE_DATA_FILE=1;
DISPLAY_OUTPUTS=0;
if POINCARE_SIMULATION==0
    SAVENAME=strcat('few_particles_eq_8keV_Ecoef0p2_pphi0p0002_G230216_h',num2str(TPRECISE))
else
    SAVENAME=strcat('few_particles_eq_poincarre_noRMP_25keV_p0p1_G090315_h',num2str(TPRECISE))
end


if SAVE_DATA_FILE==0
    disp('no data file for this simulation')
end

if DISPLAY_OUTPUTS==0
    disp('**************************************************************');
    disp('no performance outputs will be displayed during this simulation')
    disp('**************************************************************');
else
    disp('**************************************************************');
end

%% Correct theta maps for interpolation
% corrected theta maps for complete interpolation (correct for 2 pi crossing)
QNB_THETA=round(0.25*NB_THETA); %quarter of NB_THETA
HQNB_THETA=round(0.5*QNB_THETA); %half of quarter of NB_THETA

%Low map
theta_low_XZsmall_map=theta_XZsmall_map;
theta_low_XZsmall_map(theta_XZsmall_map>(NB_THETA-QNB_THETA-2)*DTHETA)=theta_XZsmall_map(theta_XZsmall_map>(NB_THETA-QNB_THETA-2)*DTHETA)-2*pi;

% High map
theta_up_XZsmall_map=theta_XZsmall_map;
theta_up_XZsmall_map(theta_XZsmall_map<QNB_THETA*DTHETA)=theta_XZsmall_map(theta_XZsmall_map<QNB_THETA*DTHETA)+2*pi;

%% Initialization of ions

% POSITIONS AND INITIALIZATION OF VELOCITY SPACE

if POINCARE_SIMULATION==0
    Nalphas_simulated=64; %number of (fast) ions
    
    alphas_pos_z=1.2*((0:31)/31-0.6)
    alphas_pos_z=repmat(alphas_pos_z',2,1); %Redone this ans removed transpose in next lines
    alphas_pos_x=zeros(size(alphas_pos_z));
    alphas_pos_z=alphas_pos_z + Z_axis;
    alphas_pos_x=alphas_pos_x + Raxis-R0; % magnetic axis and major axis
    alphas_pos_phi=zeros(size(alphas_pos_x))+1.25; % Arbitraty toroidal angle
    
    % MAGNETIC FIELD
    alphas_Bfield=interp2(scale_X,scale_Z,Btot_XZ_map',alphas_pos_x,alphas_pos_z,'*linear');
    
    % ENERGY AND VELOCITY
    alphas_Ekin=1*repmat(8*1e3,1,16);  % keV
    alphas_Ekin = repmat(alphas_Ekin',2,1); %(eV)
    alphas_vtot=sqrt(alphas_Ekin*2*eV/mHe); %Absolute value for speed |v|
    alphas_vpll=alphas_vtot*0.15;
    % Alter velocities for distribution
    alphas_vpll(33:64)=4*alphas_vpll(33:64); %Increase last half by 4
    alphas_vpll(33:2:64)=-alphas_vpll(33:2:64);
    alphas_Epll=0.5*(mHe/eV)*alphas_vpll.^2; % Parallel kinetic energy (eV)
    alphas_Eperp=alphas_Ekin-alphas_Epll;    % Perpendicular kinetic energy (eV)
    alphas_mm=alphas_Eperp./alphas_Bfield;   % magnetic moment (eV/T)
    
else
    Nalphas_simulated=128;
    Nalphas_half=round(0.5*Nalphas_simulated)-1;
    Nalphas_fourth=round(0.25*Nalphas_simulated)-1;
	alphas_pos_z=zeros(1,Nalphas_simulated);
    
    alphas_pos_z(1:Nalphas_fourth+1)=0.38*((0:Nalphas_fourth)/Nalphas_fourth);
	alphas_pos_z(Nalphas_fourth+2:Nalphas_half+1)=-alphas_pos_z(1:Nalphas_fourth+1)-0.3;
	alphas_pos_z(1:Nalphas_fourth+1)=alphas_pos_z(1:Nalphas_fourth+1)+0.3;
	alphas_pos_z=alphas_pos_z+0.001;

    alphas_pos_x=alphas_pos_z*0;
    alphas_pos_x(round(0.5*Nalphas_simulated)+1:Nalphas_simulated)=0.9*((0:Nalphas_half)/Nalphas_half-0.5)
	alphas_pos_x=alphas_pos_x-0.001;
    alphas_pos_z=alphas_pos_z' + Z_axis;
    alphas_pos_x=alphas_pos_x';
	%alphas_pos_x=alphas_pos_x' + Raxis-R0;
    
    alphas_pos_phi=0*alphas_pos_x;
    alphas_Bfield=interp2(scale_X,scale_Z,Btot_XZ_map',alphas_pos_x,alphas_pos_z,'*linear');
    
    alphas_Ekin=(1e3)*repmat(25,Nalphas_simulated,1);  % keV
	alphas_vtot=sqrt(alphas_Ekin*2*eV/mHe); %Absolute value for speed |v|
    alphas_vpll=(alphas_vtot*0.1);
% 	alphas_vpll(1:2:Nalphas_simulated)=alphas_vpll(1:2:Nalphas_simulated)*0.6;
	alphas_vpll(1:2:Nalphas_simulated)=-alphas_vpll(1:2:Nalphas_simulated);
%	alphas_vpll(1:4:Nalphas_simulated)=alphas_vpll(1:4:Nalphas_simulated)*0.5;
% 	alphas_vpll(2:4:Nalphas_simulated)=alphas_vpll(2:4:Nalphas_simulated)*0.15;
%   alphas_vpll=alphas_vpll.*sqrt(alphas_Ekin*2*eV/mHe);
    alphas_mm=(alphas_Ekin-0.5*((mHe/eV)*alphas_vpll.^2))./alphas_Bfield;
end
alphas_Epll=0.5*(mHe/eV)*alphas_vpll.^2;
alphas_Eperp=max(alphas_Ekin-alphas_Epll,0);

alphas_ejected=zeros(Nalphas_simulated,1); %?

%% Storing initial values and pre-allocating
alphas_Ekin0=alphas_Ekin;
alphas_mm0=alphas_mm;
X0=alphas_pos_x;
Z0=alphas_pos_z;
phi0=alphas_pos_phi; 

alphas_vperp=sqrt(2*alphas_Eperp*eV/mHe); %Perpendicular velocity
alphas_rhoL=(mHe/eV)*(1/ZHe)*alphas_vperp./alphas_Bfield; %Gyroradius ions

disp('**************************************************************');

%initializing simulation arrays
Xpos_outputG=zeros(Nalphas_simulated,NB_TIME_STAMPS);
Zpos_outputG=zeros(Nalphas_simulated,NB_TIME_STAMPS);
phipos_outputG=zeros(Nalphas_simulated,NB_TIME_STAMPS);
psipos_outputG=zeros(Nalphas_simulated,NB_TIME_STAMPS);
vD_outputG=zeros(Nalphas_simulated,NB_TIME_STAMPS);
vparallel_outputG=zeros(Nalphas_simulated,NB_TIME_STAMPS);
Eperp_outputG=zeros(Nalphas_simulated,NB_TIME_STAMPS);
theta_outputG=zeros(Nalphas_simulated,NB_TIME_STAMPS);
Ekin_outputG=zeros(Nalphas_simulated,NB_TIME_STAMPS);
vphi_outputG=zeros(Nalphas_simulated,NB_TIME_STAMPS);
pphi_outputG=zeros(Nalphas_simulated,NB_TIME_STAMPS);
pphi_recalc_outputG=pphi_outputG;
psi_value_outputG=zeros(Nalphas_simulated,NB_TIME_STAMPS);
rhoL_outputG=zeros(Nalphas_simulated,NB_TIME_STAMPS);
Bfield_outputG=zeros(Nalphas_simulated,NB_TIME_STAMPS);
Ekin_half_outputG=zeros(Nalphas_simulated,NB_TIME_STAMPS);
Xpos_gc_outputG=Xpos_outputG*0;
Zpos_gc_outputG=Zpos_outputG*0;


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


% initial field and position (psi and theta)
alphas_psi_value=interp2(scale_X,scale_Z,psi_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear'); %Flux coordinate
alphas_psi_value=max(alphas_psi_value,-psi_global);
alphas_psi=interp1(psi_scale,1:NB_PSI,alphas_psi_value);

alphas_theta=interp2(scale_X,scale_Z,theta_XZsmall_map',alphas_pos_x,alphas_pos_z,'*nearest'); % Poloidal angle for initialization (nearest to prevent errors caused by 0 and 2pi border.)

% estimating psi and theta values (redone using interpolation)
interpolate_theta_psi_fromXZ;
alphas_pos_phi_wrapped=wrap2pi(alphas_pos_phi);


%% Direction of B
% direction of the B field at local positions
bX=interp2_XZ(interp_x,interp_z,bX_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
bZ=interp2_XZ(interp_x,interp_z,bZ_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
bphi=sqrt(1-(bX.^2+bZ.^2));

% norm_angle=rand(Nalphas_simulated,1)*2*pi;
norm_angle=ones(Nalphas_simulated,1)*0.5*pi; %Gyro-angle

% Vectors perpendicular to B-field (arbitrary sign / relation)
uX=sqrt(1./(1+(bX./bphi).^2));  % Makes the u-vector normalized
uZ=0*uX;                        % Set to 0 for convenience
uphi=-(bX./bphi).*uX;           % Makes the u-vector perpendicular to B-field

% Second orthogonal vector by cross product
wX=uZ.*bphi-uphi.*bZ; 
wZ=uphi.*bX-uX.*bphi;
wphi=uX.*bZ-uZ.*bX;

% normal vector N = (cos) u + (sin) w. N is the direction of perpendicular
% velocity (b of parallel velocity). 
% ! Note the gyro-angle should be random since u and w and in a random direction!
NX=(cos(norm_angle).*uX+sin(norm_angle).*wX);
NZ=((cos(norm_angle).*uZ+sin(norm_angle).*wZ));
Nphi=(cos(norm_angle).*uphi+sin(norm_angle).*wphi);

v0_X=alphas_vpll.*bX+alphas_vperp.*NX;
v0_Z=alphas_vpll.*bZ+alphas_vperp.*NZ;
v0_phi=alphas_vpll.*bphi+alphas_vperp.*Nphi;

v_X=v0_X;
v_Z=v0_Z;
v_phi=v0_phi;

alphas_Omega=(ZHe*eV/mHe)*alphas_Bfield;    % Larmor frequency
epsilon=0.5*h*alphas_Omega;                 % Parameter for matrix evolution (see 4.9 Thesis Fabien)
alphas_determinant=1+(epsilon.^2);          

% previous mid point speed values
v_X_prev=v0_X+epsilon.*(-bphi.*v0_Z+bZ.*v0_phi);
v_Z_prev=v0_Z+epsilon.*(bphi.*v0_X-bX.*v0_phi);
v_phi_prev=v0_phi+epsilon.*(-bZ.*v0_X+bX.*v0_Z);

v_phi_prev_prev=v_phi_prev;
v_phi_corr_integ=v_phi_prev*0;
% v_phi_prev=v_phi;

%initialize Bfield properly
alphas_mm_part=alphas_mm;
% update_GT_3D_collapse;

interpolate_theta_psi_fromXZ;

% amplitude of the B field and potential at half time step local positions
alphas_Bfield=interp2_XZ(interp_x,interp_z,Btot_XZ_map,INDEX_LIST_1,INDEX_LIST_2,INDEX_LIST_3,INDEX_LIST_4);
if APPLY_RMP_FIELD==1
    BR_tilde=alphas_Bfield*0;
    BZ_tilde=alphas_Bfield*0;
    Bphi_tilde=alphas_Bfield*0;
    [IL3D_1 IL3D_2 IL3D_3 IL3D_4 IL3D_5 IL3D_6 IL3D_7 IL3D_8 slopex slopey slopez] = ...
        build_3Dinterp_indexarrays(scale_phi, scale_theta, scale_psi, DPHI,DTHETA,DPSIVAL,alphas_pos_phi_wrapped, alphas_theta, alphas_psi);
    BR_tilde=lininterp3( BR_RMP,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
    BZ_tilde=lininterp3( BZ_RMP,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
    Bphi_tilde=lininterp3( Bphi_RMP,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
    if CALCULATE_TRUE_PPHI==1
        alphas_dAphi_dphi=alphas_mm*0;
        alphas_Delta_pphi=alphas_mm*0;
        alphas_dAphi_dphi=lininterp3( dAphi_dphi,IL3D_1,IL3D_2,IL3D_3,IL3D_4,IL3D_5,IL3D_6,IL3D_7,IL3D_8, slopex,slopey,slopez);
    end
end

% Canonical angular momentum evolution
alphas_omega=wrap2pi(alphas_theta-alphas_pos_phi_wrapped);

alphas_psi_value_corr=alphas_psi_value;

alphas_pphi0=(mHe/eV)*alphas_Rpos.*v_phi-(ZHe)*alphas_psi_value_corr;
alphas_pphi0_half=alphas_pphi0;

alphas_Ekin_half=alphas_Ekin;
alphas_Ekin_prev=alphas_Ekin;
alphas_pphi0_prev=alphas_pphi0;

    
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
    time_step_integration_GT_eq;
    
    
    %calculating integer time step velocity values
    if CALCULATE_TRUE_PPHI==0
        update_GT_3D_eq;
        v_X_step=(0.5*v_X+0.5*v_X_prev);
        v_Z_step=(0.5*v_Z+0.5*v_Z_prev);
        v_phi_step=(0.5*v_phi+0.5*v_phi_prev);
        v_X=v_X_prev;
        v_Z=v_Z_prev;
        v_phi=v_phi_prev;
    else
        alphas_pphi0=alphas_pphi0_prev+alphas_Delta_pphi;
        alphas_pphi0_prev=alphas_pphi0;
        alphas_Delta_pphi=zeros(Nalphas_simulated,1);
    end
	
    %applying vphi correction
    if PPHI_CORR_FACTOR~=0
        v_phi_step_recalc=(alphas_pphi0+(ZHe)*alphas_psi_value)./(alphas_Rpos);
        v_phi_step_recalc=(eV/mHe)*v_phi_step_recalc;
        adapt_speed_pphi_G;
    end
    %applying Ekin correction
    if ECOEF~=0
        adapt_speed_Ekin_G;
    end
    %applying vphi correction
    if PPHI_CORR_FACTOR~=0
        update_GT_3D_eq;
        v_X_step=(0.5*v_X+0.5*v_X_prev);
        v_Z_step=(0.5*v_Z+0.5*v_Z_prev);
        v_phi_step=(0.5*v_phi+0.5*v_phi_prev);
        v_X=v_X_prev;
        v_Z=v_Z_prev;
        v_phi=v_phi_prev;
        adapt_speed_pphi_G;
    end
    

	% checking for ejected particles
    alphas_vpll=v_X_step.*bX+v_Z_step.*bZ+v_phi_step.*bphi;
    alphas_Epll=0.5*(mHe/eV)*alphas_vpll.^2;
    alphas_Eperp=max(alphas_Ekin-alphas_Epll,0);
    outcast=find(alphas_psi>simulation_size_r);
    alphas_vperp=sqrt(2*alphas_Eperp*eV/mHe);
    if (~isempty(outcast))
        recast=randi(Nalphas_simulated,size(outcast,1),1);
        reposition_lost_particles_3DG;
        alphas_ejected(outcast)=1;
        disp(strcat('number of ejected particles =  ',num2str(size(outcast,1))));
    end

    if (mod(time_step,TIME_STAMP_PRECISION)==0)
        % recalculating velocity time step values
        update_GT_3D_eq;
        v_X_step=(0.5*v_X+0.5*v_X_prev);
        v_Z_step=(0.5*v_Z+0.5*v_Z_prev);
        v_phi_step=(0.5*v_phi+0.5*v_phi_prev);
        v_X=v_X_prev;
        v_Z=v_Z_prev;
        v_phi=v_phi_prev;
        
        time_stamp=ceil((time_step)/TIME_STAMP_PRECISION);
        alphas_omega=wrap2pi(alphas_theta-alphas_pos_phi_wrapped);
        
        % precise value of psi for pphi comparison
        alphas_psi_value=interp2(scale_X,scale_Z,psi_XZsmall_map',alphas_pos_x,alphas_pos_z,'*cubic');
        alphas_psi_value_corr=alphas_psi_value;
        
        % corrections to recover guiding center position
        pos_X_gc_corr=(mHe/eV)*(1/ZHe)*(v_Z_step.*bphi-v_phi_step.*bZ)./alphas_Bfield;
        pos_Z_gc_corr=(mHe/eV)*(1/ZHe)*(v_phi_step.*bX-v_X_step.*bphi)./alphas_Bfield;
        

        alphas_vpll=v_X_step.*bX+v_Z_step.*bZ+v_phi_step.*bphi;
        alphas_Epll=0.5*(mHe/eV)*alphas_vpll.^2;
        alphas_Eperp=max(alphas_Ekin-alphas_Epll,0);
        alphas_vperp=sqrt(2*alphas_Eperp*eV/mHe);
        alphas_rhoL=(mHe/eV)*(1/ZHe)*alphas_vperp./alphas_Bfield;
		
        Xpos_outputG(:,time_stamp)=alphas_pos_x;
        Zpos_outputG(:,time_stamp,:)=alphas_pos_z;
        phipos_outputG(:,time_stamp,:)=alphas_pos_phi;
        psipos_outputG(:,time_stamp,:)=alphas_psi;
        vparallel_outputG(:,time_stamp,:)=alphas_vpll;
        Eperp_outputG(:,time_stamp,:)=alphas_Eperp;
        theta_outputG(:,time_stamp,:)=alphas_theta;
        Ekin_outputG(:,time_stamp,:)=0.5*(mHe/eV)*(v_X_step.^2+v_Z_step.^2+v_phi_step.^2);
		Ekin_half_outputG(:,time_stamp,:)=0.5*(mHe/eV)*(v_X.^2+v_Z.^2+v_phi.^2);
        vphi_outputG(:,time_stamp,:)=v_phi_step;
        pphi_outputG(:,time_stamp,:)=alphas_pphi0;
        pphi_recalc_outputG(:,time_stamp,:)=(mHe/eV)*(alphas_pos_x+R0).*v_phi_step-(ZHe)*alphas_psi_value_corr;
        psi_value_outputG(:,time_stamp,:)=alphas_psi_value_corr;
        rhoL_outputG(:,time_stamp,:)=alphas_rhoL;
        Bfield_outputG(:,time_stamp,:)=alphas_Bfield;
        Xpos_gc_outputG(:,time_stamp)=alphas_pos_x+pos_X_gc_corr;
        Zpos_gc_outputG(:,time_stamp,:)=alphas_pos_z+pos_Z_gc_corr;
        
        % displaying the difference with the initial values of the COMs
        if DISPLAY_OUTPUTS==1
            if time_step>RECORD_TIME_STAMP
                pphi_recalc=(mHe/eV)*(alphas_pos_x+R0).*v_phi_step-(ZHe)*alphas_psi_value_corr;
                disp('------------------------------')
                disp('pphi error')
                mean(abs(pphi_recalc-alphas_pphi0))
                disp('------------------------------')
                Ekin_recalc=0.5*(mHe/eV)*(v_X_step.^2+v_Z_step.^2+v_phi_step.^2);
                Ekin_recalc2=0.5*(mHe/eV)*(v_X.^2+v_Z.^2+v_phi.^2);
                disp('------------------------------')
                disp('Ekin error')
                mean(abs(alphas_Ekin-Ekin_recalc))
                mean(abs(alphas_Ekin-Ekin_recalc2))
                disp('------------------------------')
            end
        end
    end

    if (SAVE_DATA_FILE==1) && (mod(time_step,RECORD_TIME_STAMP)==0)
        toc
        tic
        clear outcast recast
        % correct the particles that are out of simulation domain
        % by giving them a position randomly in the initial distribution
        outcast=find(alphas_psi>simulation_size_r);
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

        save (SAVENAME ,'alphas_Ekin','alphas_pphi0','Bfield_outputG','Xpos_outputG','Zpos_outputG','Xpos_gc_outputG','Zpos_gc_outputG','phipos_outputG','theta_outputG',...
            'psi_value_outputG','vparallel_outputG','Ekin_outputG','Ekin_half_outputG','Eperp_outputG','rhoL_outputG','psipos_outputG',...
            'vphi_outputG','pphi_outputG','pphi_recalc_outputG','time_scale','alphas_ejected');
        disp('------------------------------------------------')
    end



end