

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
    
    %     filename=strcat(DATA_FOLDER,'psi_star_evol.mat');
    %     load(filename, 'psi_star_2D_evol_interp');
    %     psi_star_map_phi_evol=interp1(1:1001,psi_star_2D_evol_interp,1:101,'cubic');
    
    
    REINIT_SIMULATION_DATA=1;
end
Bini_PR_map=sqrt(Btor_PR_map.^2+Bpol_PR_map.^2);
BXini_PR_map=BX_PR_map;
BZini_PR_map=BZ_PR_map;

TOROIDAL_FIELD_DIRECTION=sign(mean(mean(Bphi_XZsmall_map)))
initialize_eq_sim_parameters_gamma_kink;

for f=23:2:25

	FRAME_NUMBER=f
	FRAME_NUMBER_GAMMA=10*(FRAME_NUMBER-1)+1

	frame_rank=FRAME_NUMBER_GAMMA;
	load_Bmaps_frame_rank;
	load_Emaps_frame_rank;

	%build_dvD_map_frame_rank;
	%save (strcat('dvD_map_frame_rank',num2str(FRAME_NUMBER),'.mat'),'vD_X_PR_map_phi','vD_Z_PR_map_phi','vD_phi_PR_map_phi','ikink_kinetic_energy_phi');

	for PROCESS_NUMBER=1:8
		close all
		
		% initialize completely the maps for pre-collapse calculations
		
		DISTNAME=strcat('initial_sNBI60keV_counter_D_distribution',num2str(PROCESS_NUMBER),'.mat')
		STATNAME=strcat('initial_NBI60keV_counter_precession_stats',num2str(PROCESS_NUMBER),'.mat')
		LOADNAME=strcat('initial_counter_NBI60keV_precession',num2str(PROCESS_NUMBER),'.mat')
		
		load(DISTNAME);
		load(LOADNAME);
		load(STATNAME);
		
		transp_weight=6.28259e+13
		NBSPLITS=2
		initialize_eq_sim_maps_kink;
		
		NB_TIME_STAMP_PRECESS_DATA=size(Xpos_gc_output,1)
		PART_POP=find((~alphas_ejected).*(~isnan(alphas_Ekin)));
		
		% counter current NBI means co field stabilizing passing ions
		%PART_POP=PART_POP(CO_PASSING);
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
		load('energy_mag_evol.mat')
		WMHD_KINK=deltaB2_tot_evol_global(FRAME_NUMBER);
		
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
		EK_KINK=2*sum(ikink_kinetic_energy_phi(1:end-1))
		WKINK=2*EK_KINK
		
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
			
			
			while ((isempty(CLOSED_ORBIT))&&(DIST_CLOSED_LOOP<0.05))
				DIST_CLOSED_LOOP=DIST_CLOSED_LOOP+5e-3;
				CLOSED_ORBIT=find(alphas_excursion_output(TIMEINI+12:end,n)<DIST_CLOSED_LOOP);
				if DIST_CLOSED_LOOP>3e-2
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
		imagesc(theta_scale,1:size(vD_dot_E_map_theta,2),-vD_dot_E_map_theta');
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
		
		save (strcat('dgamma_ions_counterNBI',num2str(PROCESS_NUMBER),'_f',num2str(FRAME_NUMBER),'.mat'),'ions_gamma','ions_deltaE',...
			'ions_gamma_vD','ions_gamma_mu',...
			'gamma_mu_omega_psi_map','gamma_mu_theta_psi_map','gamma_vD_omega_psi_map','gamma_vD_theta_psi_map','gamma_vD_phi_psi_map');
		
		disp('**************************************************************');
		disp('done');
		Nalphas_simulated
		% gamma_fast_ions=sum(gamma_ions_evol)
		ions_gamma=min(ions_gamma,MAX_GAMMA);
		ions_gamma=max(ions_gamma,MIN_GAMMA);
		gamma_fast_ions=sum(ions_gamma)
		
	 end   
end
% exit
