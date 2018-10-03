

close all;
pause(0.2);

REINIT_SIMULATION_DATA=1;


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
    filename=strcat(DATA_FOLDER,'tokamak_PR_map.mat');
    load(filename);
    filename=strcat(DATA_FOLDER,'B_fields.mat');
    load(filename);
    filename=strcat(DATA_FOLDER,'tokamak_map_dimensions.mat');
    load(filename);
    filename=strcat(DATA_FOLDER,'flux_geometry.mat');
    load(filename);
    filename=strcat(DATA_FOLDER,'pressure_profile.mat');
    load(filename);
    filename=strcat(DATA_FOLDER,'q_profile.mat');
    load(filename);
    filename=strcat(DATA_FOLDER,'Epot_psi_star_dot_evol.mat');
    load(filename);
    
    filename=strcat(DATA_FOLDER,'psi_star_evol.mat');
    load(filename,'ksi0_evol_lin');
    %     load(filename, 'psi_star_2D_evol_interp');
    %     psi_star_map_phi_evol=interp1(1:1001,psi_star_2D_evol_interp,1:101,'cubic');
    
    
    REINIT_SIMULATION_DATA=1;
end
%theta_scale=(0:NP-1)*2*pi/(NP-1);
load('Wkink_evol_global.mat', 'Wkink_deltaK_evol');


Bini_PR_map=sqrt(Btor_PR_map.^2+Bpol_PR_map.^2);
BXini_PR_map=BX_PR_map;
BZini_PR_map=BZ_PR_map;

TOROIDAL_FIELD_DIRECTION=sign(mean(mean(Bphi_XZsmall_map)))
initialize_eq_sim_parameters_gamma_kink;


%evaluate MHD energy
derq=gradient(q_initial_profile,radial_r_value_flux);
Btot_PR_map=sqrt(BX_PR_map.^2+BZ_PR_map.^2+Btor_PR_map.^2);
B1=mean(Btot_PR_map(1:end-1,psi_rank_q1))
P1=P_initial_profile(psi_rank_q1)
n1=Ne_profile(psi_rank_q1)
vA1=B1/sqrt(mu0*mD*n1);
derq1=derq(psi_rank_q1);
omegaA1=vA1/(sqrt(3)*r_value_q1_mean*R0*derq1)
tauA1=1/omegaA1

Bpol1sq=mean(Bpol_PR_map(1:end-1,psi_rank_q1).^2);
dr_avg=mean(dr_PR_map(1:end-1,:));
integ1=0;
for r=2:psi_rank_q1
    integ1=integ1+dr_avg(r)*(P_initial_profile(r)-P1)*2*radial_r_value_flux(r);
end
beta_pol1=(2*mu0/Bpol1sq)*integ1/r_value_q1_mean^2

%  gamma Bussac 
%  ksi=1
MHD_raw=(6*pi^2*Bphi0^2*r_value_q1_mean^4)/(mu0*R0^3)*(1-q_initial_profile(1))*(beta_pol1^2-13/144);
gammaI=(3*pi*omegaA1*r_value_q1_mean^2)/(R0^2)*(1-q_initial_profile(1))*(beta_pol1^2-13/144)


for f=14:2:20
    
    ksi0=ksi0_evol_lin(f)
    %WMHD_KINK=(MHD_raw*ksi0^2)
    WMHD_KINK=Wkink_deltaK_evol(f)
    
    
	FRAME_NUMBER=f
	FRAME_NUMBER_GAMMA=10*(FRAME_NUMBER-1)+1

	frame_rank=FRAME_NUMBER_GAMMA;
	load_Bmaps_frame_rank;
	load_Emaps_frame_rank;

	%build_dvD_map_frame_rank;
	%save (strcat('dvD_map_frame_rank',num2str(FRAME_NUMBER),'.mat'),'vD_X_PR_map_phi','vD_Z_PR_map_phi','vD_phi_PR_map_phi','ikink_kinetic_energy_phi','Efield_phi_map_phi');
	build_Ephi_map_frame_rank;

	for PROCESS_NUMBER=1:2
		close all
		
		% initialize completely the maps for pre-collapse calculations
		if FRAME_NUMBER==14
		DISTNAME=strcat('../../AUG_30382_2p5/particles_equilibrium/initial_NBI60keV_transpMS_D_distribution',num2str(PROCESS_NUMBER),'.mat')
		STATNAME=strcat('../../AUG_30382_2p5/particles_equilibrium/initial_NBI60keV_precession_stats',num2str(PROCESS_NUMBER),'.mat')
		LOADNAME=strcat('../../AUG_30382_2p5/particles_equilibrium/initial_NBI60keV_precession',num2str(PROCESS_NUMBER),'.mat')
		
		load(DISTNAME);
		load(LOADNAME);
		load(STATNAME);
        end
		
        NB_SPLITS=4;

        % 1 million particles in this one!
        particles_weight=2.82871e+13 
        particles_weight=particles_weight/NB_SPLITS

		initialize_eq_sim_maps_kink;
		
		NB_TIME_STAMP_PRECESS_DATA=size(Xpos_gc_output,1)
		PART_POP=find((~alphas_ejected).*(~isnan(alphas_Ekin)));
		
		% co current NBI means counter field passing ions
% 		PART_POP=PART_POP(COUNTER_PASSING);
		disp('number of ions considered=')
		disp(length(PART_POP));
		
		if NB_TIME_STAMP_PRECESS_DATA~=NB_TIME_STAMPS
			disp('Warning: the number of time points isnt consistent with the initialize_eq_sim_parameters_gamma_kink paramteters!')
		end
		
	%     clear Bphi_XZ_map BpolX_initial_map BpolZ_initial_map Rpos_map Rpos_map theta_XZ_map psi_XZ_map
	%     clear Btot_ini
		
		% Initial field value
		% frame_rank=1
		% load_Bmaps_frame_rank;
		% Btot_map_phi_ini=Btot_map_phi;
		% bX_map_phi_ini=bX_map_phi;
		% bZ_map_phi_ini=bZ_map_phi;
		% grad_psi_star_map_phi_ini=grad_psi_star_map_phi;
		
		
		% Initialize perturbed field values at considered time point
		frame_rank=FRAME_NUMBER_GAMMA;
		FRAME_NUMBER
% 		load('energy_mag_evol.mat')
% 		WMHD_KINK=deltaB2_tot_evol_global(FRAME_NUMBER);
		
		frame_rank=FRAME_NUMBER_GAMMA-20;
		load_Bmaps_frame_rank;
		Btot_map_phi_prev_prev=Btot_map_phi;
		frame_rank=FRAME_NUMBER_GAMMA-10;
		load_Bmaps_frame_rank;
		Btot_map_phi_prev=Btot_map_phi;
		frame_rank=FRAME_NUMBER_GAMMA+10;
		load_Bmaps_frame_rank;
		Btot_map_phi_next=Btot_map_phi;
		frame_rank=FRAME_NUMBER_GAMMA+20;
		load_Bmaps_frame_rank;
		Btot_map_phi_next_next=Btot_map_phi;
		dBtot_map_phi=-(1/12)*Btot_map_phi_next_next+(2/3)*Btot_map_phi_next-(2/3)*Btot_map_phi_prev+(1/12)*Btot_map_phi_prev_prev;
		DTEXP=(4e-6)
		dBtot_map_phi=dBtot_map_phi/DTEXP;
		
		frame_rank=FRAME_NUMBER_GAMMA;
		load_Bmaps_frame_rank;
		
		% save (strcat('gB_map_frame_rank',num2str(FRAME_NUMBER),'.mat'),'gB_X_PR_map_phi','gB_Z_PR_map_phi');
		
		%%
		
		load(strcat('dvD_map_frame_rank',num2str(FRAME_NUMBER),'.mat'));
% 		EK_KINK=2*sum(ikink_kinetic_energy_phi(1:end-1))
% 		WKINK=EK_KINK+WMHD_KINK
        WKINK=WMHD_KINK
		
		init_gamma_kink_omega_maps;
		
		% build_E_map_frame_rank;
		
		alphas_Ekin_pop=alphas_Ekin(PART_POP);
		alphas_mm_pop=alphas_mm(PART_POP);
		
		
		%delta_E_ions_evol=zeros(NB_TIME_STAMPS,1);
		%gamma_ions_evol=zeros(NB_TIME_STAMPS,1);
		ions_gamma=zeros(length(PART_POP),1);
		ions_gamma_vD=zeros(length(PART_POP),1);
		ions_gamma_mu=zeros(length(PART_POP),1);
		ions_deltaE=zeros(length(PART_POP),1);
		
		%%
		TIMEINI=10;
		
		alphas_Xini=Xpos_gc_output(TIMEINI,PART_POP);
		alphas_Zini=Zpos_gc_output(TIMEINI,PART_POP);
		alphas_excursion_output=Xpos_gc_output(:,PART_POP)*0;
		for n=1:length(PART_POP)
			alphas_excursion_output(:,n)=(Xpos_gc_output(:,PART_POP(n))-alphas_Xini(n)).^2+(Zpos_gc_output(:,PART_POP(n))-alphas_Zini(n)).^2;
		end
		alphas_excursion_output=sqrt(alphas_excursion_output);
		alphas_time_orbit=alphas_Xini*0;
		alphas_delta_time=alphas_Xini*0;
		
		
		
		for n=1:length(PART_POP)
			time_orbit=0;
			DIST_CLOSED_LOOP=5e-3;
			%     [min_dist time_orbit]=min(alphas_excursion_output(TIMEINI+12:end,n));
			CLOSED_ORBIT=find(alphas_excursion_output(TIMEINI+12:end,n)<DIST_CLOSED_LOOP);
			
			
			while ((isempty(CLOSED_ORBIT))&&(DIST_CLOSED_LOOP<0.1))
				DIST_CLOSED_LOOP=DIST_CLOSED_LOOP+5e-3;
				CLOSED_ORBIT=find(alphas_excursion_output(TIMEINI+12:end,n)<DIST_CLOSED_LOOP);
				if DIST_CLOSED_LOOP>0.095
					disp('found orbit difficult to close at ')
					n
					DIST_CLOSED_LOOP
				end
			end
			if (~isempty(CLOSED_ORBIT))
				time_orbit=CLOSED_ORBIT(1);
				time_orbit=time_orbit+11+TIMEINI;
				if (time_orbit<length(time_scale)-1)
					while (alphas_excursion_output(time_orbit+1,n)<=alphas_excursion_output(time_orbit,n))&&(time_orbit<length(time_scale)-1)
						time_orbit=time_orbit+1;
					end
				end
			end
			
			%     if (~isempty(CLOSED_ORBIT)) && (time_orbit==0)
			%         time_orbit=CLOSED_ORBIT(1);
			%         time_orbit=time_orbit+11+TIMEINI;
			%     end
			if (isempty(CLOSED_ORBIT))
				time_orbit=0
			end
			if (time_orbit>12+TIMEINI)&&(time_orbit<length(time_scale))
				% as many orbits as possible
				%while time_orbit<length(time_scale)-time_orbit
				%			time_orbit=time_orbit+time_orbit;
				%end
				alphas_time_orbit(n)=time_orbit;
				alphas_delta_time(n)=time_scale(time_orbit)-time_scale(TIMEINI);
			else
				time_orbit=0;
				alphas_time_orbit(n)=0;
				alphas_delta_time(n)=0;
			end
			
		end
		%%
		EJECTED_PART=(alphas_ejected);
		%DELTA_TIME=time_scale(end)-time_scale(1)
		figure(2)
		hold on
		imagesc((0:NP-1)*2*pi/(NP-1),1:size(vD_dot_E_map_theta,2),-vD_dot_E_map_theta');
		xlim([0 2*pi])
		ylim([1 size(vD_dot_E_map_theta,2)])
		
		for n=1:length(PART_POP)
			calculate_ion_gamma_avg;
		end
		
		% gamma_mu_omega_psi_map1=gamma_mu_omega_psi_map;
		% gamma_mu_theta_psi_map1=gamma_mu_theta_psi_map;
		% gamma_vD_omega_psi_map1=gamma_vD_omega_psi_map;
		% gamma_vD_theta_psi_map1=gamma_vD_theta_psi_map;
		% gamma_vD_phi_psi_map1=gamma_vD_phi_psi_map;
		
		save (strcat('dgamma_ions_coNBI',num2str(PROCESS_NUMBER),'_f',num2str(FRAME_NUMBER),'.mat'),'ions_gamma','ions_deltaE',...
			'ions_gamma_vD','ions_gamma_mu',...
			'gamma_mu_omega_psi_map','gamma_mu_theta_psi_map','gamma_vD_omega_psi_map','gamma_vD_theta_psi_map','gamma_vD_phi_psi_map');
		
		disp('**************************************************************');
		disp('done');
		Nalphas_simulated
		% gamma_fast_ions=sum(gamma_ions_evol)
		ions_gamma=min(ions_gamma,MAX_GAMMA);
		ions_gamma=max(ions_gamma,MIN_GAMMA);
		gamma_fast_ions=sum(ions_gamma)
		disp('**************************************************************');
		
	 end   
end
% exit
