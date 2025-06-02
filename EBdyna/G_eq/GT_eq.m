function [exitcode]=GT_eq(x,v,input,output,ejected)
global maps dim par time qom const
% % Functions to convert energy [eV] to velocity [m/s] and vica verca
% E_to_v=@(m,E) sqrt(2*E*(const.eV/input.m));
v_to_E=@(m,v) 0.5*(m/const.eV)*sum(v.^2,2);

%GT_eq simulation code for integration of particle motion
%   Processes particles according to parameters, maps, dim etc.
%
%   Works with modes:
% [] - DEFAULT (first 2 then 3)
% [] - DEFAULT (first 2 then 3)
% 1  - TEST
% 2  - PRECESSION
%       Will use a (shorter) simulation in 2D field and categorize particles.
%       Determines bounce / drift frequencies
%       Calls G_eq with same 'box'-nr for full simulation
%       Produces:
%           RAW     - file:     In which par and output are stored
%           STATS   - file:     With precession data.
%               .mat        - succesvol
%               .matprec    - failed. Can be used for new precession or is used if all particles in equilibrium were ejected
% 3/4  - FULL (after precession/without precession)
%       Longer simulation in required field (e.g. RMP-field, with ST)
%       Produces:
%           RAW     - file:     In which par and output are stored
%           FULL    - file:     After evaluation (maybe with precession data)
% 5  - POINCARE
%       Plots the poincare map at phi==0 in a figure if GUI is used. Stores
%       position based on phi=0 crossing
%       Produces:
%           RAW     - file:     In which par and output are stored.
% NOTE: THE RAW FILE COULD BE USED TO CONTINUE FROM THIS POINT. DELETE RAW
% FILES WITH INCORRECT PARAMETERS OR EBdyna_go MIGH CRASH / PRODUCE FAULTY
% OUTPUT
TESTING_CX=0;TESTING_SD=0;
if isfield(par,'TESTING_CX')
    TESTING_CX=par.TESTING_CX;
else
    par.TESTING_CX=0;
end
if isfield(par,'TESTING_SD')
    TESTING_SD=par.TESTING_SD;
else
    par.TESTING_SD=0;
end
if ~isfield(par,'remove_thermalized_markers')
    par.remove_thermalized_markers=0;
end
if ~isfield(par,'full_2D_neutrals')
    par.full_2D_neutrals=0;
end
if isfield(par,'isDTbackground')
    if par.isDTbackground
        disp('DT background simulation, changing bulk plasma isotope mass.')
        const.mBulk=0.5.*(const.mD+const.mT);
    else
        const.mBulk=const.mD;
    end
else
    const.mBulk=const.mD;
end
if isfield(par,'Hisotope_frac')
    if par.Hisotope_frac>0
        disp('Hydrogen minority in the background is changing isotope mass of bulk')
        const.mBulk=(par.Hisotope_frac*const.mH+(1-par.Hisotope_frac)*const.mBulk);
    else
        const.mBulk=const.mD;
    end
else
    const.mBulk=const.mD;
	par.Hisotope_frac=0;
end
% const.amu= 1.66053906660e-27;
% mass numbers for RABBIT Coulomb log calculations
% taken wrt. mH it seems
const.ABulk=const.mBulk./const.mH;
const.Ae=const.me./const.mH;
input.Afast=input.m./const.mH;
input.ABulk=const.ABulk; % for the record

%#ok<*NASGU>  % Removes warning in Matlab editor that some values are not used (they are often saved in the output file(s)).
%% data for ndd / slowing down calculation
if (par.CALCULATE_NDD==1) || (par.COULOMB_COLL) || (par.CALCULATE_CX)
    dim.pn=dim.psi_norm;
    dim.xlength=length(dim.pn);
    dim.pn_max=max(dim.pn);

%     ADJUST_NUI=1.0;  % collision with background ions can be adjusted to match experiments
    psi_xi_pos    = ejected*0;
    dxi           = ejected*0;
    dx            = ejected*0;
    psi_slope     = ejected*0;
	psi_gc_value  = ejected*0;
    psi_value_avg = ejected*0;
    delta_Ekin    = ejected*0;
	Xgc_ind=((x(:,1)-dim.R0)*dim.DX_inv)+dim.mid_Xzero;
    Zgc_ind=( x(:,2)        *dim.DZ_inv)+dim.mid_Z;
    
    sqrt_x_v    = ejected*0;
    x_v         = ejected*0;
    psi_x       = ejected*0;
    psi_prime_x = ejected*0;
    pmpp_x      = ejected*0;
    nu_perp     = ejected*0;
    
end

if par.CALCULATE_CX==1 || par.CALCULATE_NDD==1
    dim.pn=dim.psi_norm;
    dim.xlength=length(dim.pn);
	psi_gc_value  = ejected*0;
    psi_value_avg = ejected*0;
    Ekin_val   = ejected*0;
    E_sigv     = ejected*0;
    vnbiD      = ejected*0;
end

ejected_CX      = ejected*0; % list of ejected due to CX
COLL_GROUP      = ejected*0;
if par.CALCULATE_CX==1 
    COLL_GROUP_N      = ejected*0;
end

if par.CALCULATE_CX==1
    init_GTeq_CX;
end

if par.CALCULATE_NDD==1
    % Bosch-Hale coefficients for sigma of DD from beam target
    A1_NDD =   53701;
    A2_NDD =   330.2700;
    A3_NDD =  -0.1271;
    A4_NDD =   2.9327e-05;
    A5_NDD =  -2.5151e-09;
    s_E        = ejected*0;
    sv_val     = ejected*0;
    ndd_cum    = ejected*0; % cumulated ndd
    ndd_val    = ejected*0; % instantaneous ndd Beam-Plasma
	
	const.m_neutron  = 1.674927498e-27; % kg
	const.m_r_nDD    = (input.m*const.mBulk)/(input.m+const.mBulk);
	const.E0_nDD     = 2.45*1e6;    % 2.45 MeV neutrons from DD
	const.v0_nDD     = sqrt(2*const.E0_nDD*(const.eV/const.m_neutron));    % velocity of neutrons from DD
			
	if isfield(par,'RECORD_VSPACE_NDD')
		if par.RECORD_VSPACE_NDD
		    sigma_MB         = ejected*0;
			
			theta_iso        = ejected*0;  % angle for isotropic distribution of velocities
			alpha_iso        = ejected*0; % angle for isotropic distribution of velocities
			vXYZ_neutrons    = v*0; % isotropic neutron source (in COM referential)
			
			vXYZ_MB_Ti        = v*0; % Maxwellian in Cartesian coordinates, based on Ti
			vXYZ_COM_BP       = v*0; % COM velocity for BP interaction, in cartesians
			vXYZ_COM_BB       = v*0; % COM velocity for BB interaction, in cartesians
			E_COM_BB          = ejected*0; % COM energy for BB interaction, in cartesians
			vXYZ_fast         = v*0; % Beam velocity in cartesian coordinates
			vXYZ_fast_prev    = v*0; % Beam velocity in cartesian coordinates (previous time stamp)
			% psi_gc_value_prev = ejected*0;  % radial position (fast ion or gc of fast ion) at previous time stamp
			BINVALS           = ejected*0;  % radial bin values on nDD bins
			BINVALS_PREV      = ejected*0;  % radial bin values on nDD bins at previous time stamp
			
            nfast         = ejected*0; % density of fast ions from previous run
			ndd_val_BP    = ejected*0; % instantaneous ndd Beam-Plasma
			ndd_val_BB    = ejected*0; % instantaneous ndd Beam-Beam
			v_n_BP     = v*0; % neutron velocity [cartesians!] from BP (1 interaction with random velocity pick in Ti Maxwellian) 
			v_n_BB     = v*0; % neutron  velocity [cartesians!] from BB (1 interaction with random velocity pick in previous step)
			v_rel_BB   = ejected*0;
			v_rel_BP   = ejected*0;
		end
	else
		par.RECORD_VSPACE_NDD=0;
    end   
    if isfield(par,'CALCULATE_NDD_THERMAL')
		if par.CALCULATE_NDD_THERMAL
			x_ni_XYZi          = x*0; % Markers for the background density
			vXYZ_MB_Ti         = v*0; % Maxwellian velocity in cartesian coordinates
			vXYZ_MB_Ti_prev    = v*0; % Maxwellian velocity in cartesian coordinates (previous time stamp)
			ndd_val_th         = ejected*0; % instantaneous ndd Beam-Beam
            % fields required from initial maps:            dim.x_ni dim.Ti_x_ni dim.ni_markers_weight
        end
	else
		par.CALCULATE_NDD_THERMAL=0;
    end   
	Ti                  = ejected*0;
    ni                  = ejected*0;
    v_norm              = ejected*0;

end

%% Default start parameters 
output.loss         = zeros(par.NB_TIME_STEPS,1);
PSI_LIMIT_EJECTED=par.PSI_LIMIT_EJECTED;
psi_value    = ejected*0;
RZout_check  = ejected*0;
ejected_wall = ejected*0;
EJECTED_THRESH = 0.5;   % value for ejected flag [compared to interpolation from map between 0 and 1]
dim.Zup_lim     = max(dim.scale_Z);
dim.Zlow_lim    = min(dim.scale_Z);
dim.Rup_lim     = max(dim.scale_X)+dim.R0;
dim.Rlow_lim    = min(dim.scale_X)+dim.R0;
% Rin_lim     = ejected*0;
% Rout_lim    = ejected*0;
% Rout_check  = ejected*0;
% Rin_check   = ejected*0;
% Rup_check   = ejected*0;
% Rlow_check  = ejected*0;
SIM_BORN    = ejected*0;
if par.COULOMB_COLL
    THERMALIZED_MARKERS = ejected*0;
    const.A_e=const.me/const.mH;
    const.A_i=const.mBulk/const.mH;
    input.Afast=input.m/const.mH;
    nu_e                = ejected*0;
    nu_i                = ejected*0;
    log_lambda_Di       = ejected*0;
    log_lambda_De       = ejected*0;
    vth_e               = ejected*0;
    vth_i               = ejected*0;
    nu_tot              = ejected*0;
    Te                  = ejected*0;
    Ti                  = ejected*0;
    ne                  = ejected*0;
    ni                  = ejected*0;
    v_crit              = ejected*0;
    delta_v             = ejected*0;
    delta_v_spread      = ejected*0;
    vpll_sd             = ejected*0;
    v_norm              = ejected*0;
    delta_omega         = ejected*0;
    delta_omega_spread  = ejected*0;
    SLOW_ELEC_GROUP     = ejected*0;
    Ekin_prev           = ejected*0;
    Ekin_new            = ejected*0;
    % this minimum velocity value is taken in deep SOL region: used only
    % for non slowing down simulation (eg. MB dist)
    MIN_COLL_V          = sqrt(2*mean(dim.Te_prof(end-2:end))*(const.eV/input.m)); % E_to_v(input.m,mean(dim.Te_prof(end-1:end)));
    X_SLOW_ELEC_THRESH  = 1.062325;  % intersection for psi - psi' for x<<1 and x>>1
    %if par.CALCULATE_PDEP
    %         dv_e              = ejected*0;
    %         dv_i              = ejected*0;
    Eexch_e           = ejected*0;
    Eexch_i           = ejected*0;
    Etot_e           = ejected*0;
    Etot_i           = ejected*0;
    mom_tor_tot      = ejected*0;
    vprev_phi        = ejected*0;
    %end
    if par.IMPURITY_ZAVG>0
        Eexch_imp           = ejected*0;
        Etot_imp            = ejected*0;
        par.Aimp=par.mimp./const.mH;
    end
    if par.CALCULATE_DEVIATION
        elem_deviation      = ejected*0;
        deviation           = ejected*0;
    end
    output.loss_wall    = zeros(par.NB_TIME_STEPS,1);    
end

ejected_sd          = ejected*0;

if par.calculate_length_trajectory
    Delta_l=ejected*0;
    elem_Delta_l=ejected*0;
end
N_simulated = length(ejected);

if par.PROGRESSIVE_BIRTH==0
    SIM_BORN=ejected*0+logical(1);
else
    SIM_BORN=logical(par.birth_matrix(1,:))';
end
% (Set first time step to 1, but this could be changed when continueing from a RAW-file
start_time_step=1;
time=0;

%% Load PREC-data if requested
if par.GET_PREC_INFO && isfield(input,'prec')
    prec=input.prec; input=rmfield(input,'prec');
	prec_data=true; % For Poincare plotting in particle categories
else
	prec_data=false; % For Poincare plotting in particle categories
end

%% Add data for specific mode
if par.mode==5
    % Alter output for variable saving of parameters
    output.x=NaN(input.N_job,3,1e4);    output.x_gc=output.x;	output.v=output.x;
    ind_time=zeros(input.N_job,1);   % Vector counting number of indexes already stored
    par.time_scale=NaN(input.N_job,size(output.x,3)); % Time of each save
    
    switch par.poincare_type
        case 'sfl'
            % Pre-allocate and determine size of psi's to be plotted
            psi_large  =NaN(size(output.x,1),1,size(output.x,3));
            theta_large=NaN(size(output.x,1),1,size(output.x,3));
        case 'xy'
            theta_large=NaN(size(output.x,1),1,size(output.x,3));
            output.x_tilde  =NaN(size(output.x,1),1,size(output.x,3));
            output.y_tilde=NaN(size(output.x,1),1,size(output.x,3));
    end
    if par.Poincare_plot
        [ha,~,psi_to_psi_overline]=make_poincare;
    end
    par_dt_or=par.dt;
end

%% Add additional output
%Parameters to determine analytic value of pphi
if par.CALCULATE_PPHI_3D
    output.pphi_an_temp=zeros(input.N_job,par.NB_TIME_STAMPS);
    Delta_pphi_an=zeros(input.N_job,1);
end

% call this script (continue_from_raw) if you have to continue from raw files
% it is quite annoying and useless most of the time so now
% I have cpmmented this part out and pushed it to a separate script

% continue_from_raw;


%% Display input / output
disp('**************************************************************');
disp('PARAMETERS:')
disp(par)
disp('**************************************************************');
disp('INPUT PARTICLES')
disp(input)
disp('**************************************************************');
disp('RAW CONSISTS OF:')
disp(output)
% benchmarking time between two data recordings
disp('**************************************************************');
disp('**************************************************************');

%% SET START PARAMETERS (Temporarily stored parameters etc.)
[~,B]=B_interpolation(x);
Bfield_sq=dot(B,B,2);
Bfield=sqrt(Bfield_sq);
vpll_ej=dot(B,v,2)./Bfield;
Eperp=0.5*(input.m/const.eV)*((dot(v,v,2))-vpll_ej.^2);
x_gc=x+bsxfun(@times,1./(qom*Bfield_sq),cross(v,B,2));

v_n=v;
n_ejected=~ejected;
sn_ejected=n_ejected;
x_temp_prev3=x;
x_temp_prev2=x;
x_temp_prev=x;
x_temp=x;
v_temp=v;

sign_Z_nu=sign(x_gc(:,2)-dim.func_Z_cross(x_gc(:,1))); % Add shift for current equilibrium found empirically! Change when using other equilibrium
sign_vpll_nu=sign(dot(v,B,2));

% prelimnar initialization of storage variable
sign_vpll_prev=sign_vpll_nu;
sign_Z_prev=sign_Z_nu;
x_prev=x;
X_ind=((x(:,1)-dim.R0)*dim.DX_inv)+dim.mid_Xzero;
Z_ind=( x(:,2)        *dim.DZ_inv)+dim.mid_Z;

% if par.CALCULATE_PPHI_3D
v_temp_2=v;
% end

part_just_lost=ejected;
part_lost_recently=ejected;
NAN_OUTPUTS=ejected;
IMAG_OUTPUTS=ejected;

process_time=tic;
subtime=tic;

time_stamp=1;

%% MAIN SIMULATION LOOP


% Use matprec-file if code has produced a matprec file previouslly, but
% extract_precession failed.
%if par.mode==2 && exist([par.SAVENAME_RAW,'prec'],'file')
%    disp('prec file already exists in output folder !')
%    TEMP_NAME=[par.SAVENAME_STATS(1:end-4),'_matprec.mat'];
%    copyfile([par.SAVENAME_RAW,'prec'],TEMP_NAME)
%    %load(TEMP_NAME);
%    delete([par.SAVENAME_RAW,'prec']);
%end
    for time_step=start_time_step:par.NB_TIME_STEPS
        %% Update in advance of next iteration loop
        %Store not ejected variable (to compare which particle have been
        %lost in this loop)
        % n_ejected(SIM_BORN)=~ejected(SIM_BORN);
%         if (sum(isnan(v))+sum(isnan(x)))>0
%             disp('error')
%         end
        
        sn_ejected=~ejected;
        sn_ejected=sn_ejected&SIM_BORN;
        if par.CALCULATE_CX
            IONS_GROUP=sn_ejected&~CX_NEUTRALS;
            NEUTRALS_GROUP=sn_ejected&CX_NEUTRALS;
        else
            IONS_GROUP=sn_ejected;
        end
        
%         part_lost_recently=(output.flag_loss_counter>=1)&(output.flag_loss_counter<=par.RECORD_AFTER_LOSS);
        if length(find(part_lost_recently))>0
            sn_ejected(part_lost_recently)=1;
        end
        if length(find(NAN_OUTPUTS))>0
            sn_ejected(NAN_OUTPUTS)=0;
        end
        if length(find(IMAG_OUTPUTS))>0
            sn_ejected(IMAG_OUTPUTS)=0;
        end
        
        % Temporarily store position for recording where particles are lost
        x_temp_prev3(sn_ejected,:)=x_temp_prev2(sn_ejected,:);
        x_temp_prev2(sn_ejected,:)=x_temp_prev(sn_ejected,:);
        x_temp_prev(sn_ejected,:)=x_temp(sn_ejected,:);
        x_temp(sn_ejected,:)=x(sn_ejected,:);
        v_temp(sn_ejected,:)=v(sn_ejected,:);

        
        % Temporarily store signs of current vpll and Z
        sign_vpll_prev(sn_ejected)=sign_vpll_nu(sn_ejected);
        sign_Z_prev(sn_ejected)=sign_Z_nu(sn_ejected);
        		
        %% Time step
		% for covered distance calculation (Poincare plot)
        if par.calculate_length_trajectory
            elem_Delta_l(sn_ejected)=0*Bfield_sq(sn_ejected);
        end

        if ~par.CALCULATE_PPHI_3D


            for i=1:par.NR_FUND_IN_LOOP
                [x LOSS_POP_OUT]=apply_domain_limits(x,sn_ejected,dim);
                if par.calculate_length_trajectory
                    x_prev(IONS_GROUP,:)=x(IONS_GROUP,:);
                end
                v_temp_2(sn_ejected,:)=v(sn_ejected,:);
				if ~isreal(x(sn_ejected,3))||~isreal(v(sn_ejected,3))
					disp('warning : imaginary found after motion operator')
					IMAG_OUTPUTS(sn_ejected)=abs(imag(x(sn_ejected,3)))>0|abs(imag(v(sn_ejected,3)))>0;
					x(IMAG_OUTPUTS,:)=real(x_temp(IMAG_OUTPUTS,:));
					v(IMAG_OUTPUTS,:)=real(v_temp(IMAG_OUTPUTS,:));
					% COLL_GROUP(IMAG_OUTPUTS)=0;
					IONS_GROUP(IMAG_OUTPUTS)=0;
					sn_ejected(IMAG_OUTPUTS)=0;
					ejected(IMAG_OUTPUTS)=1;
					ejected_wall(IMAG_OUTPUTS)=1;
				end
                [x(IONS_GROUP,:),v(IONS_GROUP,:)]=time_step_integration_GT_eq_struct(x(IONS_GROUP,:),v(IONS_GROUP,:),input.Fc_field(IONS_GROUP,:));
                time=time+par.dt;

				if par.calculate_length_trajectory
					Delta_l(IONS_GROUP)=sqrt((x(IONS_GROUP,1).*cos(x(IONS_GROUP,3))-x_prev(IONS_GROUP,1).*cos(x_prev(IONS_GROUP,3))).^2+...
						 (x(IONS_GROUP,1).*sin(x(IONS_GROUP,3))-x_prev(IONS_GROUP,1).*sin(x_prev(IONS_GROUP,3))).^2+...
						 ((x(IONS_GROUP,2))-x_prev(IONS_GROUP,2)).^2);
					elem_Delta_l(IONS_GROUP)=elem_Delta_l(IONS_GROUP)+Delta_l(IONS_GROUP);
                end
                
            end
        else
            for i=1:par.NR_FUND_IN_LOOP
			    [x LOSS_POP_OUT]=apply_domain_limits(x,sn_ejected,dim);
                if par.calculate_length_trajectory
                    x_prev(IONS_GROUP,:)=x(IONS_GROUP,:);
                end
                v_temp_2(IONS_GROUP,:)=v(IONS_GROUP,:);
                field_3D=find_3D_Afield(x(IONS_GROUP,:),{'dAphi_dphi','dAR_dphi','dAZ_dphi'});
                
                % The time step
                [x(IONS_GROUP,:),v(IONS_GROUP,:)]=time_step_integration_GT_eq_struct(x(IONS_GROUP,:),v(IONS_GROUP,:),input.Fc_field(IONS_GROUP,:));
                time=time+par.dt;
                % Estimate velocity at n
                v_n(IONS_GROUP,:)=0.5*(v(IONS_GROUP,:)+v_temp_2(IONS_GROUP,:));
                
                % Evaluate pphi
                Delta_pphi_an(IONS_GROUP)=Delta_pphi_an(IONS_GROUP)...
                    +v_n(IONS_GROUP,1).*field_3D.dAR_dphi...
                    +v_n(IONS_GROUP,2).*field_3D.dAZ_dphi...
                    +v_n(IONS_GROUP,3).*field_3D.dAphi_dphi;
					
				if par.calculate_length_trajectory
					Delta_l(IONS_GROUP)=sqrt((x(IONS_GROUP,1).*cos(x(IONS_GROUP,3))-x_prev(IONS_GROUP,1).*cos(x_prev(IONS_GROUP,3))).^2+...
						 (x(IONS_GROUP,1).*sin(x(IONS_GROUP,3))-x_prev(IONS_GROUP,1).*sin(x_prev(IONS_GROUP,3))).^2+...
						 ((x(IONS_GROUP,2))-x_prev(IONS_GROUP,2)).^2);
					elem_Delta_l(IONS_GROUP)=elem_Delta_l(IONS_GROUP)+Delta_l(IONS_GROUP);
				end
            end
        end
        % moving the neutrals is rather straightforward!
        if par.CALCULATE_CX
            x(NEUTRALS_GROUP,:)=x(NEUTRALS_GROUP,:)+v(NEUTRALS_GROUP,:).*(par.NR_FUND_IN_LOOP.*par.dt);
        end
            
        % security to avoid interpolation errors related to position
		% for ensuring precision on phi values (not for prec calculation)
		if par.mode~=2
		    % checking for NAN_OUTPUTS and imaginary funny guys.....
			% happens when re-ionization are in impossible location


  
            if (sum(isnan(v(sn_ejected,:)))+sum(isnan(x(sn_ejected,:))))>0
 				disp('error : NaN found after motion loop')
                NAN_OUTPUTS(sn_ejected)=isnan(sum(v(sn_ejected,:),2))|isnan(sum(x(sn_ejected,:),2));
                x(NAN_OUTPUTS,:)=x_temp(NAN_OUTPUTS,:);
                v(NAN_OUTPUTS,:)=v_temp(NAN_OUTPUTS,:);
                IONS_GROUP(NAN_OUTPUTS)=0;
                ejected(NAN_OUTPUTS)=1;
                ejected_wall(NAN_OUTPUTS)=1;
            end  
% 		    if sum(sum(isnan(v(sn_ejected,:)))+sum(isnan(x(sn_ejected,:))))>0
% 				disp('error')
%                 NAN_OUTPUTS(sn_ejected)=isnan(sum(v(sn_ejected,:),2))|isnan(sum(x(sn_ejected,:),2));
%                 x(NAN_OUTPUTS,:)=x_temp(NAN_OUTPUTS,:);
%                 v(NAN_OUTPUTS,:)=v_temp(NAN_OUTPUTS,:);
%                 COLL_GROUP(NAN_OUTPUTS)=0;
%                 IONS_GROUP(NAN_OUTPUTS)=0;
%                 sn_ejected(NAN_OUTPUTS)=0;
%                 ejected(NAN_OUTPUTS)=1;
%                 ejected_wall(NAN_OUTPUTS)=1;
%             end
            x(sn_ejected,3)=mod(x(sn_ejected,3),2*pi);

        end
        
        % make sure neutrals are still inside grid
		[x LOSS_POP_OUT]=apply_domain_limits(x,sn_ejected,dim);

		% Estimate velocity at n
        v_n(sn_ejected,:)=0.5*(v(sn_ejected,:)+v_temp_2(sn_ejected,:));
        % Estimate x_gc at n
        [~,B(sn_ejected,:)]=B_interpolation(x(sn_ejected,:));
        Bfield_sq(sn_ejected)=dot(B(sn_ejected,:),B(sn_ejected,:),2);
        x_gc(IONS_GROUP,:)=x(IONS_GROUP,:)+bsxfun(@times,1./(qom*Bfield_sq(IONS_GROUP)),cross(v_n(IONS_GROUP,:),B(IONS_GROUP,:),2));
		% overriding gc position for outsiders
		x_gc(LOSS_POP_OUT,:)=x(LOSS_POP_OUT,:);
		
        if par.CALCULATE_CX
		    % neutrals need to be contained!!!!!!
			x(NEUTRALS_GROUP,2)=max(x(NEUTRALS_GROUP,2),dim.Zlow_lim);
			x(NEUTRALS_GROUP,2)=min(x(NEUTRALS_GROUP,2),dim.Zup_lim);
			x(NEUTRALS_GROUP,1)=max(x(NEUTRALS_GROUP,1),dim.Rlow_lim); 
			x(NEUTRALS_GROUP,1)=min(x(NEUTRALS_GROUP,1),dim.Rup_lim); 
            x_gc(NEUTRALS_GROUP,:)=x(NEUTRALS_GROUP,:);
        end

		
		% for cumulated covered distance calculation (Poincare plot)
		if par.calculate_length_trajectory
			output.Delta_l(sn_ejected)=output.Delta_l(sn_ejected)+elem_Delta_l(sn_ejected);
        end
        
		X_ind(sn_ejected)=((x(sn_ejected,1)-dim.R0)*dim.DX_inv)+dim.mid_Xzero;
        Z_ind(sn_ejected)=( x(sn_ejected,2)        *dim.DZ_inv)+dim.mid_Z;
        
        % calculate psi values for profile estimates
        if par.COULOMB_COLL || par.CALCULATE_NDD || par.CALCULATE_CX
            if par.USE_1DGC_POS
                Xgc_ind(sn_ejected)=((x_gc(sn_ejected,1)-dim.R0)*dim.DX_inv)+dim.mid_Xzero;
                Zgc_ind(sn_ejected)=( x_gc(sn_ejected,2)        *dim.DZ_inv)+dim.mid_Z;
            end
            v_norm(sn_ejected)=sqrt(dot(v(sn_ejected,:),v(sn_ejected,:),2));
			if ~isreal(v_norm(sn_ejected,:))
 				disp('WARNING : imaginary found in v_norm ... fixing for velocity and position')
                IMAG_OUTPUTS(sn_ejected)=abs(imag(v_norm(sn_ejected)))>0;
                x(IMAG_OUTPUTS,:)=real(x_temp(IMAG_OUTPUTS,:));			
				v(IMAG_OUTPUTS,:)=real(v_temp(IMAG_OUTPUTS,:));
                IONS_GROUP(IMAG_OUTPUTS)=0;
                sn_ejected(IMAG_OUTPUTS)=0;
                ejected(IMAG_OUTPUTS)=1;
                ejected_wall(IMAG_OUTPUTS)=1;
				v_norm(IMAG_OUTPUTS)=real(v_norm(IMAG_OUTPUTS));
            end
			if par.COULOMB_COLL 
                if ~par.CALCULATE_CX
				COLL_GROUP=(v_norm>MIN_COLL_V)&sn_ejected;
                COLL_GROUP_N=COLL_GROUP;
                else
 				COLL_GROUP=(v_norm>MIN_COLL_V)&IONS_GROUP;
                % added to update ion density for neutrals reionizations
 				COLL_GROUP_N=COLL_GROUP|NEUTRALS_GROUP;
               end
			else
				COLL_GROUP=sn_ejected;
                COLL_GROUP_N=COLL_GROUP;
            end
            if par.USE_1DGC_POS
                psi_gc_value(COLL_GROUP_N)=ba_interp2(maps(1).psi_norm1_XZ,Zgc_ind(COLL_GROUP_N),Xgc_ind(COLL_GROUP_N),'linear'); % find psi-gc-value
            else
                psi_gc_value(COLL_GROUP_N)=ba_interp2(maps(1).psi_norm1_XZ,Z_ind(COLL_GROUP_N),X_ind(COLL_GROUP_N),'linear'); % find psi-value at particle [better!]
                psi_gc_value(COLL_GROUP_N)=min(psi_gc_value(COLL_GROUP_N),dim.pn_max);  % security because of poor rounding of interpolation function!
            end
            psi_value_avg(COLL_GROUP_N)=psi_value_avg(COLL_GROUP_N)+psi_gc_value(COLL_GROUP_N);
            % index vector
            [~,psi_xi_pos(COLL_GROUP_N)] = histc(psi_gc_value(COLL_GROUP_N),dim.pn);
            psi_xi_pos(COLL_GROUP_N) = max(psi_xi_pos(COLL_GROUP_N),1);     % To avoid index=0 when xi < x(1)
            psi_xi_pos(COLL_GROUP_N) = min(psi_xi_pos(COLL_GROUP_N),dim.xlength-1);   % To avoid index=m+1 when xi > x(end).
            % 't' slope vector [p x 1]
            dxi(COLL_GROUP_N) = psi_gc_value(COLL_GROUP_N)-dim.pn(psi_xi_pos(COLL_GROUP_N));
            dx(COLL_GROUP_N) = dim.pn(psi_xi_pos(COLL_GROUP_N)+1)-dim.pn(psi_xi_pos(COLL_GROUP_N));
            psi_slope(COLL_GROUP_N) = dxi(COLL_GROUP_N)./dx(COLL_GROUP_N);
			
			ni(COLL_GROUP_N)=interp1qr(dim.pn,dim.ni_prof,psi_gc_value(COLL_GROUP_N),psi_xi_pos(COLL_GROUP_N),psi_slope(COLL_GROUP_N));
            Ti(COLL_GROUP_N)=interp1qr(dim.pn,dim.Ti_prof,psi_gc_value(COLL_GROUP_N),psi_xi_pos(COLL_GROUP_N),psi_slope(COLL_GROUP_N));               
			vth_i(COLL_GROUP)=sqrt(2*Ti(COLL_GROUP)*(const.eV/const.mBulk));
            
            if par.CALCULATE_CX
                if ~par.full_2D_neutrals
                    neutral_density(COLL_GROUP)=interp1qr(dim.pn,dim.n0_prof,psi_gc_value(COLL_GROUP),psi_xi_pos(COLL_GROUP),psi_slope(COLL_GROUP));
                else
                    if par.N0_FAC_D2>0
                        [neutral_density(COLL_GROUP),neutral_density_D2(COLL_GROUP)]=neutrals_2D_interp(x(COLL_GROUP,:));
                    else
                        [neutral_density(COLL_GROUP),~]=neutrals_2D_interp(x(COLL_GROUP,:));
                    end
                end
                if par.USE_T0_TABLE
                    neutral_temp(COLL_GROUP)=interp1qr(dim.pn,dim.T0_prof,psi_gc_value(COLL_GROUP),psi_xi_pos(COLL_GROUP),psi_slope(COLL_GROUP));
                end
            end
        end
        if par.COULOMB_COLL
            
            Ekin_prev(COLL_GROUP)=0.5*(input.m/const.eV).*sum(v(COLL_GROUP,:).^2,2);
%             v(COLL_GROUP,:)=bsxfun(@times,(v_norm(COLL_GROUP)-delta_v(COLL_GROUP))./v_norm(COLL_GROUP),v(COLL_GROUP,:));

            Bfield(COLL_GROUP)=sqrt(Bfield_sq(COLL_GROUP));
			vprev_phi(COLL_GROUP)=v(COLL_GROUP,3);
            if par.IMPURITY_ZAVG>0
                [ v(COLL_GROUP,:),Eexch_e(COLL_GROUP),Eexch_i(COLL_GROUP),Eexch_imp(COLL_GROUP),Etot_e(COLL_GROUP),Etot_i(COLL_GROUP),Etot_imp(COLL_GROUP),delta_omega(COLL_GROUP)] = ...
                    calculate_collisions_imp(const,par,dim,X_SLOW_ELEC_THRESH,input,Ti(COLL_GROUP),vth_i(COLL_GROUP),ni(COLL_GROUP),psi_gc_value(COLL_GROUP),psi_xi_pos(COLL_GROUP),psi_slope(COLL_GROUP),...
                    v_norm(COLL_GROUP),Ekin_prev(COLL_GROUP),v(COLL_GROUP,:),B(COLL_GROUP,:),Bfield(COLL_GROUP),Eexch_e(COLL_GROUP),Eexch_i(COLL_GROUP),Eexch_imp(COLL_GROUP),...
                    Etot_e(COLL_GROUP),Etot_i(COLL_GROUP),Etot_imp(COLL_GROUP),TESTING_SD);
            else
                [ v(COLL_GROUP,:),Eexch_e(COLL_GROUP),Eexch_i(COLL_GROUP),Etot_e(COLL_GROUP),Etot_i(COLL_GROUP),delta_omega(COLL_GROUP)] = ...
                    calculate_collisions(const,par,dim,X_SLOW_ELEC_THRESH,input,Ti(COLL_GROUP),vth_i(COLL_GROUP),ni(COLL_GROUP),psi_gc_value(COLL_GROUP),psi_xi_pos(COLL_GROUP),psi_slope(COLL_GROUP),...
                    v_norm(COLL_GROUP),Ekin_prev(COLL_GROUP),v(COLL_GROUP,:),B(COLL_GROUP,:),Bfield(COLL_GROUP),Eexch_e(COLL_GROUP),Eexch_i(COLL_GROUP),...
                    Etot_e(COLL_GROUP),Etot_i(COLL_GROUP),TESTING_SD);
            end

%             if (sum(isnan(v(sn_ejected,:)))+sum(isnan(x(sn_ejected,:))))>0
%  				disp('error : NaN found after collision operator')
%                 NAN_OUTPUTS(sn_ejected)=isnan(sum(v(sn_ejected,:),2))|isnan(sum(x(sn_ejected,:),2));
%                 x(NAN_OUTPUTS,:)=x_temp(NAN_OUTPUTS,:);
%                 v(NAN_OUTPUTS,:)=v_temp(NAN_OUTPUTS,:);
%                 COLL_GROUP(NAN_OUTPUTS)=0;
%                 IONS_GROUP(NAN_OUTPUTS)=0;
% %                 sn_ejected(NAN_OUTPUTS)=0;
%                 ejected(NAN_OUTPUTS)=1;
%                 ejected_wall(NAN_OUTPUTS)=1;
%             end
%            if sum(isnan(delta_omega(sn_ejected)))>0
%  				disp('error : NaN found in delta_omega');
% 				NAN_OUTPUTS(sn_ejected)=isnan(delta_omega(sn_ejected));
% 				delta_omega(NAN_OUTPUTS)=0;
% 		   end
 			
            Ekin_new(COLL_GROUP)=0.5*(input.m/const.eV).*sum(v(COLL_GROUP,:).^2,2);
			
			%added security to avoid pesky residual imaginaries....
            delta_Ekin(COLL_GROUP)=delta_Ekin(COLL_GROUP)+real(Ekin_new(COLL_GROUP)-Ekin_prev(COLL_GROUP));
            

            [v(COLL_GROUP,:) ]=deviate_velocity_vector(v(COLL_GROUP,:),v_norm(COLL_GROUP),delta_omega(COLL_GROUP));
            % total toroidal angular momentum given to plasma through collision
            mom_tor_tot(COLL_GROUP)=mom_tor_tot(COLL_GROUP)+par.MARKER_WEIGHT.*const.mBulk.*(vprev_phi(COLL_GROUP)-v(COLL_GROUP,3)).*squeeze(x(COLL_GROUP,1));

            % security check for errors : usually caused by incorrect maps / profiles

%             if ~isreal(x(sn_ejected,3))||~isreal(v(sn_ejected,3))
%  				disp('warning : imaginary found after collision operator')
%                 IMAG_OUTPUTS(sn_ejected)=abs(imag(x(sn_ejected,3)))>0|abs(imag(v(sn_ejected,3)))>0;
%                 x(IMAG_OUTPUTS,:)=real(x_temp(IMAG_OUTPUTS,:));
%                 v(IMAG_OUTPUTS,:)=real(v_temp(IMAG_OUTPUTS,:));
%                 COLL_GROUP(IMAG_OUTPUTS)=0;
%                 IONS_GROUP(IMAG_OUTPUTS)=0;
% %                 sn_ejected(IMAG_OUTPUTS)=0;
%                 ejected(IMAG_OUTPUTS)=1;
%                 ejected_wall(IMAG_OUTPUTS)=1;
%             end
            
            if par.CALCULATE_DEVIATION
				deviation(COLL_GROUP)=deviation(COLL_GROUP)+delta_omega(COLL_GROUP);
			end
        end
        if par.CALCULATE_NDD || par.CALCULATE_CX
            Ekin_val(COLL_GROUP)=0.5.*(input.m/const.eV).*(v_norm(COLL_GROUP)).^2;
            % cross section correted to D mass
            vnbiD(COLL_GROUP)=v_norm(COLL_GROUP); % can be adjusted later on according to toroidal rotation
            if par.CALCULATE_CX
                % isotope (D2 instead of H2) effect in the cross section formula
                E_sigv(COLL_GROUP)=0.5*Ekin_val(COLL_GROUP)*1e-3;
                % cross section wrt. D2 background : no D2 inside separatrix
                COLL_GROUP_D2=(psi_gc_value>1)&COLL_GROUP;
                sigma_cx_D2(COLL_GROUP_D2)=(1e-20).*b1_CX.*exp(-b2_CX./E_sigv(COLL_GROUP_D2).^b3_CX)./(b4_CX.*E_sigv(COLL_GROUP_D2).^b5_CX+b6_CX.*E_sigv(COLL_GROUP_D2).^b7_CX+b8_CX.*E_sigv(COLL_GROUP_D2).^b9_CX+b10_CX.*E_sigv(COLL_GROUP_D2).^b11_CX);
                % the D CX cross section was thought to be in center of mass referential but actually is NOT
                % E_sigv(COLL_GROUP)=0.5*E_sigv(COLL_GROUP);
                if ~par.USE_T0_TABLE
                    sigma_cx(COLL_GROUP)=10^-20 .* (A1_CX .*log (A2_CX./E_sigv(COLL_GROUP) + A6_CX))./(1+A3_CX.*E_sigv(COLL_GROUP)+A4_CX.*E_sigv(COLL_GROUP).^3.5 + A5_CX.*E_sigv(COLL_GROUP).^(5.4));
                    if TESTING_CX
                        sigma_cx(COLL_GROUP)=100*sigma_cx(COLL_GROUP); % (for testing)
                        sigma_cx_D2(COLL_GROUP)=100*sigma_cx_D2(COLL_GROUP);
                    end
                    CX_rate(COLL_GROUP)=sigma_cx(COLL_GROUP).*neutral_density(COLL_GROUP).*vnbiD(COLL_GROUP);
                else
                    % Indexes 2D using 'sigma_cxD_T0_loglog' (careful here, tabulated using deuterium "real" kin. energy)
                    % 'Elog_scaling', 'Emin_eV', 'T0_values_log_lin', 'T0log_scaling', 'T0min_eV', 'sigma_cxD_T0_loglog'
                    E_sigv_log(COLL_GROUP)  = dim.Elog_scaling.*log(Ekin_val(COLL_GROUP)./dim.Emin_eV);
                    Elog_ind(COLL_GROUP)    = E_sigv_log(COLL_GROUP)*dim.DElog_inv;
                    neutral_temp(COLL_GROUP)= max(neutral_temp(COLL_GROUP),dim.T0min_eV);
                    T0_ind(COLL_GROUP)      = dim.T0log_scaling.*log10(neutral_temp(COLL_GROUP)./dim.T0min_eV)*dim.DT0_inv;
                    CX_rate(COLL_GROUP)     = ba_interp2(dim.sigma_cxD_T0_loglog ,T0_ind(COLL_GROUP),Elog_ind(COLL_GROUP),'linear');
                    CX_rate(COLL_GROUP)     = CX_rate(COLL_GROUP).*neutral_density(COLL_GROUP);
                    if TESTING_CX
                        CX_rate(COLL_GROUP)=100*CX_rate(COLL_GROUP); % (for testing)
                        sigma_cx_D2(COLL_GROUP)=100*sigma_cx_D2(COLL_GROUP);
                    end
                end
                if ~par.full_2D_neutrals
                    neutral_density_D2(COLL_GROUP_D2)=interp1qr(dim.pn,dim.n0_prof_D2,psi_gc_value(COLL_GROUP_D2),psi_xi_pos(COLL_GROUP_D2),psi_slope(COLL_GROUP_D2));
                end
                CX_rate(COLL_GROUP_D2)=CX_rate(COLL_GROUP_D2)+sigma_cx_D2(COLL_GROUP_D2).*neutral_density_D2(COLL_GROUP_D2).*vnbiD(COLL_GROUP_D2);
                % decreased cumulative rate to reflect extended marker life
				CX_val(COLL_GROUP)=CX_val(COLL_GROUP).*exp(-(par.NR_FUND_IN_LOOP.*par.dt).*CX_rate(COLL_GROUP));
				CX_cum(COLL_GROUP)=1 - CX_val(COLL_GROUP);
                % CX_cum(COLL_GROUP)=CX_cum(COLL_GROUP)+CX_val(COLL_GROUP) .* CX_rate(COLL_GROUP).*(par.NR_FUND_IN_LOOP.*par.dt);
                % CX_val(COLL_GROUP)=1-CX_cum(COLL_GROUP);
                
                if sum(NEUTRALS_GROUP)>0
                    % impact ionization rates then CX processses
                    if ~par.USE_T0_TABLE
                        % first consider impact ionization rate
                        sigma_cx(NEUTRALS_GROUP)=interp1qr(Eb_ionization, sigma_ionization_val ,Ekin_val(NEUTRALS_GROUP));
                        NCX_rate(NEUTRALS_GROUP)=sigma_cx(NEUTRALS_GROUP).*ni(NEUTRALS_GROUP).*v_norm(NEUTRALS_GROUP);
                        sigma_cx(NEUTRALS_GROUP)=10^-20 .* (A1_CX .*log (A2_CX./E_sigv(NEUTRALS_GROUP) + A6_CX))./(1+A3_CX.*E_sigv(NEUTRALS_GROUP)+A4_CX.*E_sigv(NEUTRALS_GROUP).^3.5 + A5_CX.*E_sigv(NEUTRALS_GROUP).^(5.4));
                        if TESTING_CX 
                            sigma_cx(NEUTRALS_GROUP)=40*sigma_cx(NEUTRALS_GROUP); % (for testing)
                        end
                        NCX_rate(NEUTRALS_GROUP)=NCX_rate(NEUTRALS_GROUP)+sigma_cx(NEUTRALS_GROUP).*ni(NEUTRALS_GROUP).*v_norm(NEUTRALS_GROUP);
                    else
                        % Indexes 2D using 'sigma_cxD_T0_log' : note ion
                        % temp used for both neutrals and ions (bad approx. in SOL but ion density small there)
                        E_sigv_log(NEUTRALS_GROUP)  = dim.Elog_scaling.*log(Ekin_val(NEUTRALS_GROUP)./dim.Emin_eV);
                        Elog_ind(NEUTRALS_GROUP)    = E_sigv_log(NEUTRALS_GROUP)*dim.DElog_inv;
    %                     Ti_ind(NEUTRALS_GROUP)     = Ti(NEUTRALS_GROUP)*dim.DT0_inv;
                        Ti_ind(NEUTRALS_GROUP)      = dim.T0log_scaling.*log10(max(Ti(NEUTRALS_GROUP),dim.T0min_eV)./dim.T0min_eV)*dim.DT0_inv;
                         % first consider impact ionization rate
                        NCX_rate(NEUTRALS_GROUP)=ba_interp2(dim.sigma_ionD_Ti_loglog ,Ti_ind(NEUTRALS_GROUP),Elog_ind(NEUTRALS_GROUP),'linear').*ni(NEUTRALS_GROUP);
                        NCX_rate(NEUTRALS_GROUP)    = NCX_rate(NEUTRALS_GROUP)+ba_interp2(dim.sigma_cxD_T0_loglog ,Ti_ind(NEUTRALS_GROUP),Elog_ind(NEUTRALS_GROUP),'linear').*ni(NEUTRALS_GROUP);
                        % NCX_rate(NEUTRALS_GROUP)    = NCX_rate(NEUTRALS_GROUP).*ni(NEUTRALS_GROUP);
                    end
                     NCX_val(NEUTRALS_GROUP)=NCX_val(NEUTRALS_GROUP) .* exp(-(par.NR_FUND_IN_LOOP.*par.dt).*NCX_rate(NEUTRALS_GROUP));
                    NCX_cum(NEUTRALS_GROUP)=1 - NCX_val(NEUTRALS_GROUP);
                end
				
                % NCX_cum(NEUTRALS_GROUP)=NCX_cum(NEUTRALS_GROUP)+NCX_val(NEUTRALS_GROUP) .* NCX_rate(NEUTRALS_GROUP).*(par.NR_FUND_IN_LOOP.*par.dt);
                % NCX_val(NEUTRALS_GROUP)=1-NCX_cum(NEUTRALS_GROUP);
               %                 
            end
            if par.CALCULATE_NDD
	            E_sigv(COLL_GROUP)=0.5*Ekin_val(COLL_GROUP)*1e-3;
                % cross section in center of mass for D
                s_E(COLL_GROUP)=(A1_NDD+E_sigv(COLL_GROUP).*(A2_NDD+E_sigv(COLL_GROUP).*(A3_NDD+E_sigv(COLL_GROUP).*(A4_NDD+E_sigv(COLL_GROUP).*A5_NDD))));
                sv_val(COLL_GROUP)=s_E(COLL_GROUP)./(E_sigv(COLL_GROUP).*exp(31.397./sqrt(E_sigv(COLL_GROUP))));  % in millibarns (10^-31 m^2)
                sv_val(COLL_GROUP)=sv_val(COLL_GROUP).*vnbiD(COLL_GROUP)*1e-31;
                sv_val(COLL_GROUP)=sv_val(COLL_GROUP).*(1+3*Ti(COLL_GROUP)./Ekin_val(COLL_GROUP)); % correction according to background Ti (approximation for Mikkelsen 1989)
                ndd_val(COLL_GROUP)=(par.MARKER_WEIGHT.*(1-par.Hisotope_frac).*ni(COLL_GROUP)).*sv_val(COLL_GROUP);       % neutron rate you would have over 1s
                ndd_val(~sn_ejected)=0;
                ndd_cum(COLL_GROUP)=ndd_val(COLL_GROUP).*(par.NR_FUND_IN_LOOP.*par.dt)+ndd_cum(COLL_GROUP);
            end
        end
% 		if par.CALCULATE_PDEP
% 			Eexch_e(~sn_ejected)=0;
% 			Eexch_i(~sn_ejected)=0;
% 		end
		
        %% Determine losses and additional output
        % in case of a Nan in coordinate, consider it ejected
		% ejected=ejected | isnan(x(:,1,:));
		if par.USE_VESSEL_LIMIT            
            RZout_check(sn_ejected)=ba_interp2(dim.vessel.wall_RZmap,Z_ind(sn_ejected),X_ind(sn_ejected),'linear'); % find limit value [0-1]
			ejected_wall(sn_ejected)    = logical(RZout_check(sn_ejected)>=EJECTED_THRESH | ...
			x(sn_ejected,2) <= dim.Zlow_lim | x(sn_ejected,2) >= dim.Zup_lim | ...
			x(sn_ejected,1) <= dim.Rlow_lim | x(sn_ejected,1) >= dim.Rup_lim);
			ejected(sn_ejected)    = ejected(sn_ejected) | ejected_wall(sn_ejected);
		else			
			psi_value(sn_ejected)=ba_interp2(maps(1).psi_XZ,Z_ind(sn_ejected),X_ind(sn_ejected),'linear'); % find psi-value
			if (dim.psi_scale(end)-dim.psi_scale(1))>0  		% psi_Scale increasing 
                ejected_wall(sn_ejected) = psi_value(sn_ejected) >= PSI_LIMIT_EJECTED;
				ejected(sn_ejected)=ejected(sn_ejected) | ejected_wall(sn_ejected);
			else 						% psi_scale decreasing 
                ejected_wall(sn_ejected) = psi_value(sn_ejected) <= PSI_LIMIT_EJECTED;
				ejected(sn_ejected)=ejected(sn_ejected) | ejected_wall(sn_ejected);
            end
        end
        
        %remove thermalized markers that have Ekin above 1.5*Ti
        if par.COULOMB_COLL
		    if par.remove_thermalized_markers
				THERMALIZED_MARKERS(sn_ejected)=(v_norm(sn_ejected)< 1.2247*vth_i(sn_ejected)) | (v_norm(sn_ejected)<=1.2247*MIN_COLL_V);
			else
			    % removing only very low velocities (to avoid negative energies)
				THERMALIZED_MARKERS(sn_ejected)=(v_norm(sn_ejected)<=1.2247*MIN_COLL_V);
			end
            THERMALIZED_MARKERS=logical(THERMALIZED_MARKERS);
            Ekin_val(THERMALIZED_MARKERS)=0.5.*(input.m/const.eV).*(v_norm(THERMALIZED_MARKERS)).^2;
            if sum(THERMALIZED_MARKERS)>=1 & TESTING_SD
%                 disp('removing slowing down ions with max Ekin....')
                disp(['removing ion with Ekin ' num2str(max(Ekin_val(THERMALIZED_MARKERS)))]);
            end            % updating output for time stamp right away
            output.Edep_th(THERMALIZED_MARKERS)=Ekin_val(THERMALIZED_MARKERS);
		    ejected_sd(sn_ejected) = THERMALIZED_MARKERS(sn_ejected);
            ejected(sn_ejected) = ejected(sn_ejected) | ejected_sd(sn_ejected);
            THERMALIZED_MARKERS=ejected*0;
            if ~par.CALCULATE_CX
                output.loss_wall(time_step)=sum(ejected & ~ ejected_sd );
            else
                output.loss_wall(time_step)=sum(ejected & ~ ejected_sd & ~ ejected_CX );
            end
        end
        
        if par.CALCULATE_CX
            % remove those hopeless markers who are really in CX loss region
            CX_NEUTRALS_BABIES(sn_ejected) = logical(CX_cum(sn_ejected)>CX_THRESH_EJECTED(sn_ejected));
            if sum(CX_NEUTRALS_BABIES)>=1
                output.x_CXn(CX_NEUTRALS_BABIES,:) = x(CX_NEUTRALS_BABIES,:); % record neutralization positions
                CX_THRESH_EJECTED(CX_NEUTRALS_BABIES) = rand(sum(CX_NEUTRALS_BABIES),1); % new life line for the potential reborn ions
                NCX_THRESH_EJECTED(CX_NEUTRALS_BABIES) = rand(sum(CX_NEUTRALS_BABIES),1); % new life line for the newly born neutrals
                CX_cum(CX_NEUTRALS_BABIES)=0;
                NCX_cum(CX_NEUTRALS_BABIES)=0;
                CX_val(CX_NEUTRALS_BABIES)=1;
                NCX_val(CX_NEUTRALS_BABIES)=1;
                CX_NEUTRALS(CX_NEUTRALS_BABIES)=1;
                CX_flag(sn_ejected) = CX_flag(sn_ejected) + CX_NEUTRALS_BABIES(sn_ejected); % counting number of neutralizations
                CX_NEUTRALS_BABIES(CX_NEUTRALS_BABIES)=0;
            end
            % definition of a neutral
            if par.ALLOW_REIONIZATIONS
                CX_IONS_BABIES(CX_NEUTRALS)=~logical(NCX_cum(CX_NEUTRALS)<NCX_THRESH_EJECTED(CX_NEUTRALS));
                CX_NEUTRALS(CX_NEUTRALS) = ~CX_IONS_BABIES(CX_NEUTRALS);
                % record the position of reionization for the show....
                if sum(CX_IONS_BABIES)>=1
                    output.x_CXi(CX_IONS_BABIES,:) = x(CX_IONS_BABIES,:); % record ionization positions
                end
                CX_IONS_BABIES(CX_IONS_BABIES)=0;
            end
                
                
            
            % the true ejection should wait for neutral collision with the
            % wall : there remains to do a special motion algorithm for
            % neutrals!
            ejected_CX(CX_NEUTRALS) = ejected_wall(CX_NEUTRALS);
%             ejected(sn_ejected) = ejected(sn_ejected) | ejected_CX(sn_ejected);
            if ~par.COULOMB_COLL
                output.loss_CX(time_step)=sum(ejected & CX_NEUTRALS); % some rough indicator of true CX losses
            else
                output.loss_CX(time_step)=sum(ejected & ~ ejected_sd & CX_NEUTRALS); % some rough indicator of true CX losses
            end
        end

        
        %% Determine new (normally larger) generation of markers
        if par.PROGRESSIVE_BIRTH
            time=time_step*par.dt*par.NR_FUND_IN_LOOP;
            ts_birth=round(time/par.birth_time_scale(2))+1;
            SIM_BORN=logical(squeeze(par.birth_matrix(ts_birth,:)))';
        end
        
        
        %%
        if (par.mode~=5) && all(ejected) 
            warning(['Every particle lost in process: ',num2str(par.PROCESS_NUMBER)])
%             break
        end
        
        part_just_lost=ejected & sn_ejected & SIM_BORN & ~part_lost_recently;		
		
        if length(find(part_just_lost))>0
		    part_just_lost_wall=ejected & ~ejected_CX & ~ ejected_sd & sn_ejected & SIM_BORN & ~part_lost_recently;
            output.flag_loss_counter(part_just_lost_wall)=1;
			% disp(['length(find(ejected & n_ejected)) = ' num2str(length(find(ejected & n_ejected)))]);
            % Store particles that have now been ejected (but not previously) in output
            output.x_ej_prev4(part_just_lost,:)=x_temp_prev3(part_just_lost,:);
            output.x_ej_prev3(part_just_lost,:)=x_temp_prev2(part_just_lost,:);
            output.x_ej_prev2(part_just_lost,:)=x_temp_prev(part_just_lost,:);
            output.x_ej_prev(part_just_lost,:)=x_temp(part_just_lost,:);
            output.x_ej(part_just_lost,:)=x(part_just_lost,:);
%             output.x_ej_next(part_just_lost,:)=x(part_just_lost,:);
            Bfield(part_just_lost)=sqrt(Bfield_sq(part_just_lost));
            % could use v_temp instead of v here if needed
            vpll_ej(part_just_lost)=dot(B(part_just_lost,:),v(part_just_lost,:),2)./Bfield(part_just_lost);
            output.Ekin_ej(part_just_lost)=0.5*(input.m/const.eV)*(dot(v(part_just_lost,:),v(part_just_lost,:),2));
            Eperp(part_just_lost)=output.Ekin_ej(part_just_lost)-0.5*(input.m/const.eV)*(vpll_ej(part_just_lost).^2);
            output.mm_ej(part_just_lost)=Eperp(part_just_lost)./Bfield(part_just_lost);
            output.vpll_ej(part_just_lost)=vpll_ej(part_just_lost);
            output.time_step_loss(part_just_lost)=time_step;
        end
        
        part_lost_recently=(output.flag_loss_counter>=1)&(output.flag_loss_counter<=par.RECORD_AFTER_LOSS);
        % removing escaped hopeless particles
        part_lost_recently(LOSS_POP_OUT)=0;
        part_lost_recently(NAN_OUTPUTS)=0;
        part_lost_recently(IMAG_OUTPUTS)=0;
        
        LOSS_POP=find(part_lost_recently);
        
        if length(LOSS_POP)>0
            for ii=1:length(LOSS_POP)
                output.x_ej_next(LOSS_POP(ii),:,output.flag_loss_counter(LOSS_POP(ii)))=x(LOSS_POP(ii),:);
            end
            output.flag_loss_counter(LOSS_POP)=output.flag_loss_counter(LOSS_POP)+1;
            
%         else
%             part_lost_recently=[];
        end
        
        output.loss(time_step)=sum(ejected);
       
        
        %% vpll crossing and midplane crossing
        
        % Increase the number of Z-crossings
        sign_Z_nu=sign(x_gc(:,2)-dim.func_Z_cross(x_gc(:,1)));
        
        expr=sign_Z_prev~=sign_Z_nu & ~ejected;
        output.nr_midplane_crossing(expr)=output.nr_midplane_crossing(expr)+1;
        
        % Increase the number of vpll-crossings
        sign_vpll_nu=sign(dot(v,B,2));
        
        expr=sign_vpll_prev~=sign_vpll_nu & ~ejected;
        output.nr_vpll_crossing(expr)=output.nr_vpll_crossing(expr)+1;
        
        %% Store position if POINCARE simulation
        if par.mode==5
            %evaluate distance travelled along the field line

			
			% Criterium for dot: change in tor. angle more than pi/2 (on domain [0,2*pi])
            
            % SWITCH BETWEEN BOTH THESE LINES TO EITHER PLOT WHEN CROSSING
            % PHI=0 OR PLOT EVERY LOOP (OR DECOMMENT THE LOOP OBTAIN)
            ind_part=find(n_ejected & abs(mod(x_temp(:,3),2*pi)-mod(x(:,3),2*pi))>pi/2); % Find indexes of particles crossing phi=0 boundary (note, cannot determine theta for ejected particles)
            %                 ind_part=find(~ejected); if mod(time_step,100)~=0; continue; end
            
            % Skip if no value is saved
            if ~isempty(ind_part)
                ind_time(ind_part)=ind_time(ind_part)+1;                         % Find next positions of output arrays
                if any(ind_time(ind_part)>size(output.x,3))                          % Check if still saveable with current array
                    warning('Simulation stopped since save array exceeded')
                    break
                end
                
                % Repmat all indexes by 3 and coordinates by # that crossed phi=0
                % to find indexes of particles
                ind_part_2D=repmat(ind_part,1,3); ind_coord_2D=repmat(1:3,length(ind_part),1); ind_out_2D=repmat(ind_time(ind_part),1,3);
                
                ind_xv          =sub2ind(size(x)       ,ind_part_2D(:),ind_coord_2D(:));                      %2D indexes to identify values to be saved
                ind_output_xv   =sub2ind(size(output.x),ind_part_2D(:),ind_coord_2D(:),ind_out_2D(:));        %Also take time as third dimension
                
                par.dt=0.5*par.dt;
                [~,v_n,B]=time_step_integration_GT_eq_struct(x,v,input.Fc_field);
                par.dt=par_dt_or;
                
                Bfield=dot(B,B,2);
                x_gc=x+bsxfun(@times,1./(qom*Bfield),cross(v_n,B,2));
                % Store the values (make vector copy)
                %output.x   (ind_output_xv)=x   (ind_xv);
                output.x_gc(ind_output_xv)=x_gc(ind_xv);
                %output.v   (ind_output_xv)=v   (ind_xv);
                
                % Store time
                ind_time_scale=sub2ind(size(par.time_scale),ind_part,ind_time(ind_part)); 	%2D indexes to identify values to be saved
                time=time_step*par.dt*par.NR_FUND_IN_LOOP;
                par.time_scale(ind_time_scale)=time;
                
                if par.Poincare_plot
                    % Delete previous Poincare plot
                    if exist('h','var')
                        delete(h);
                    end
                    switch par.poincare_type
                        case 'sfl'
                            % Indexes 2D
                            X_ind=((output.x_gc(:,1,:)-dim.R0)*dim.DX_inv)+dim.mid_Xzero;
                            Z_ind=( output.x_gc(:,2,:)        *dim.DZ_inv)+dim.mid_Z;
                            expr_nan=isnan(X_ind); expr_nan=expr_nan(:);
                            
                            % Find theta
                            theta_large(~expr_nan) = interpolate_theta_XZ(X_ind(~expr_nan),Z_ind(~expr_nan));
                            % Find psi (index)
                            psi_large(~expr_nan)=ba_interp2(maps(1).psi_XZ,Z_ind(~expr_nan),X_ind(~expr_nan),'linear'); % find psi-value
                            
                            % Fill with new Poincare plot
                            if prec_data
                                %                 h(:,1)=plot(ha,theta(prec.pop.ALL_TRAPPED,:)',psi_to_psi_overline(psi_large(prec.pop.ALL_TRAPPED,:))','+','displayname','trapped');
                                h(:,1)=plot(ha,squeeze(theta_large(~ejected&prec.pop.ALL_PASSING,:))',psi_to_psi_overline(squeeze(psi_large(n_ejected&prec.pop.ALL_PASSING,:)))','k.','displayname','passing');
                            else
                                h=plot(ha,squeeze(theta_large(~ejected,:,:))',squeeze(psi_to_psi_overline(psi_large(~ejected,:,:)))','.');
                            end
                        case 'tor'
                            % Fill with new Poincare plot
                            h=plot(ha,squeeze(output.x_gc(~ejected,1,:))',squeeze(output.x_gc(~ejected,2,:))','k');
                        case 'xy'
                            % Indexes 2D
                            X_ind=((output.x_gc(:,1,:)-dim.R0)*dim.DX_inv)+dim.mid_Xzero;
                            Z_ind=( output.x_gc(:,2,:)        *dim.DZ_inv)+dim.mid_Z;
                            expr_nan=isnan(X_ind); expr_nan=expr_nan(:);
                            
                            % Find theta
                            theta_large(~expr_nan) = interpolate_theta_XZ(X_ind(~expr_nan),Z_ind(~expr_nan));
                            chi=ba_interp2(maps(1).chi_XZ,Z_ind(~expr_nan),X_ind(~expr_nan),'linear');
                            r= sqrt(2*chi);
                            phi=output.x(:,3,:);
                            phi=phi(~expr_nan);
                            % Determine x and y-coordinates of particle in question
                            output.x_tilde(~expr_nan) = r .* cos(theta_large(~expr_nan)-phi);
                            output.y_tilde(~expr_nan) = r .* sin(theta_large(~expr_nan)-phi);
                            expr_3=~any((output.x_tilde.^2+output.y_tilde.^2)>dim.st.r2(301).^2,3);
                            
                            h=plot(ha,squeeze(output.x_tilde(expr_3,:,:))',squeeze(output.y_tilde(expr_3,:,:))','k');
                    end
                    drawnow
                end
            end
            % look at time evolution also for debugging and consistency
            if (mod(time_step,par.TIME_STAMP_PRECISION)==0)
                time_stamp=ceil((time_step)/par.TIME_STAMP_PRECISION);
                output.x(:,:,time_stamp)=x;
                output.v(:,:,time_stamp)=v;
            end

            %% Store x and v
        elseif (mod(time_step,par.TIME_STAMP_PRECISION)==0)
            % standard storage at each time stamp
            time_stamp=ceil((time_step)/par.TIME_STAMP_PRECISION);
            if par.CALCULATE_PDEP || par.CALCULATE_NDD
                output.psi_value_avg(:,time_stamp)=psi_value_avg./par.TIME_STAMP_PRECISION;
                psi_value_avg(COLL_GROUP)=0;
            end
            if par.CALCULATE_PDEP
                output.Pdep_e(:,time_stamp)=par.MARKER_WEIGHT.*const.eV.*Eexch_e./(par.time_scale(1));
                output.Pdep_i(:,time_stamp)=par.MARKER_WEIGHT.*const.eV.*Eexch_i./(par.time_scale(1));
                if par.IMPURITY_ZAVG>0
                    output.Pdep_imp(:,time_stamp)=par.MARKER_WEIGHT.*const.eV.*Eexch_imp./(par.time_scale(1));
                    Eexch_imp(:)=0;
					output.Etot_imp(:,time_stamp)=Etot_imp;
                end
                Eexch_e(:)=0;
                Eexch_i(:)=0;
                output.Etot_e(:,time_stamp)=Etot_e;
                output.Etot_i(:,time_stamp)=Etot_i;                
                output.delta_Ekin(:,time_stamp)=delta_Ekin;
                output.cum_tor_ang_mom(:,time_stamp)=mom_tor_tot;
            end
            if par.CALCULATE_CX
                output.CX_cum(:,time_stamp)=CX_cum;
                output.CX_flag(:,time_stamp)=CX_flag;
                % time stamp latest CX cross section
                output.CX_rate(sn_ejected&~NEUTRALS_GROUP,time_stamp)=CX_rate(sn_ejected&~NEUTRALS_GROUP);
%                 output.CX_flag(:,time_stamp)=CX_cum>=CX_THRESH_FLAG;
%                 output.ejected_CX(:,time_stamp)=ejected_CX;
            end
            if par.CALCULATE_NDD
				if par.RECORD_VSPACE_NDD
                    % using simple formula for ndd values for BP
                    % ndd_val_BP(COLL_GROUP)=ndd_val(COLL_GROUP);
					% generate a MB distribution of velocities
					sigma_MB(COLL_GROUP)=sqrt(Ti(COLL_GROUP).*(const.eV/const.mBulk));
					vXYZ_MB_Ti(COLL_GROUP,1)	= random('Normal',0,sigma_MB(COLL_GROUP));
					vXYZ_MB_Ti(COLL_GROUP,2)	= random('Normal',0,sigma_MB(COLL_GROUP));
					vXYZ_MB_Ti(COLL_GROUP,3)	= random('Normal',0,sigma_MB(COLL_GROUP));
					% calculating the velocities in cartesian coordinates 
					vXYZ_fast(COLL_GROUP,1) = v(COLL_GROUP,1).*cos(x(COLL_GROUP,3))-v(COLL_GROUP,3).*sin(x(COLL_GROUP,3));
					vXYZ_fast(COLL_GROUP,2) =-v(COLL_GROUP,1).*sin(x(COLL_GROUP,3))-v(COLL_GROUP,3).*cos(x(COLL_GROUP,3));
					vXYZ_fast(COLL_GROUP,3) = v(COLL_GROUP,2);
					% COM velocities
					vXYZ_COM_BP(COLL_GROUP,:) = (input.m.*vXYZ_fast(COLL_GROUP,:)+const.mBulk.*vXYZ_MB_Ti(COLL_GROUP,:))./(input.m+const.mBulk);
					v_rel_BP(COLL_GROUP) = sqrt(sum((vXYZ_fast(COLL_GROUP,:) - vXYZ_MB_Ti(COLL_GROUP,:)).^2,2));
					% neutron source in COM referential
					theta_iso(COLL_GROUP) = 2.*rand(sum(COLL_GROUP),1)-1;
					theta_iso(COLL_GROUP) = asin(theta_iso(COLL_GROUP));
					alpha_iso (COLL_GROUP)       =  (2*pi).*rand(sum(COLL_GROUP),1);
					vXYZ_neutrons(COLL_GROUP,1)    = cos(theta_iso(COLL_GROUP)).*cos(alpha_iso(COLL_GROUP));
					vXYZ_neutrons(COLL_GROUP,2)    = cos(theta_iso(COLL_GROUP)).*sin(alpha_iso(COLL_GROUP));
					vXYZ_neutrons(COLL_GROUP,3)    = sin(theta_iso(COLL_GROUP));
					vXYZ_neutrons(COLL_GROUP,:)    = const.v0_nDD.*vXYZ_neutrons(COLL_GROUP,:);
					% total velocity in laboratory referential [cartesians!]
                    v_n_BP(COLL_GROUP,:) = vXYZ_neutrons(COLL_GROUP,:) + vXYZ_COM_BP(COLL_GROUP,:);
					% recalculate the neutron rate for BP (true maxwellians) 
                    % formula : E = EA * mB / (mA+mB)
					E_sigv(COLL_GROUP)=(const.mBulk/(const.mBulk+input.m)).*(0.5*input.m/const.eV)*sum(v_rel_BP(COLL_GROUP,:).^2,2).*1e-3;
					% cross section in center of mass for D
					s_E(COLL_GROUP)=(A1_NDD+E_sigv(COLL_GROUP).*(A2_NDD+E_sigv(COLL_GROUP).*(A3_NDD+E_sigv(COLL_GROUP).*(A4_NDD+E_sigv(COLL_GROUP).*A5_NDD))));
					sv_val(COLL_GROUP)=s_E(COLL_GROUP)./(E_sigv(COLL_GROUP).*exp(31.397./sqrt(E_sigv(COLL_GROUP))));  % in millibarns (10^-31 m^2)
					sv_val(COLL_GROUP)=sv_val(COLL_GROUP).*v_rel_BP(COLL_GROUP)*1e-31;
					ndd_val_BP(COLL_GROUP)=(par.MARKER_WEIGHT.*(1-par.Hisotope_frac).*ni(COLL_GROUP)).*sv_val(COLL_GROUP);       % neutron rate you would have over 1s (absolute, not volumic)
					ndd_val_BP(~sn_ejected)=0;						

					
					% nDD values for Beam-Beam
					% only calculated if a previous step exits!
                    [~,~,BINVALS(COLL_GROUP)] = histcounts(psi_gc_value(COLL_GROUP),12);
                    % [~,BINVALS_PREV(COLL_GROUP)] = histc(psi_gc_value_prev(COLL_GROUP),dim.pn_nDD_bins);
					if sum(vXYZ_fast_prev(COLL_GROUP,3))~=0 && sum(BINVALS_PREV(COLL_GROUP))>0
					    % find the candidates for each radial bin!
						% find relevant bin for each ion						
						
						
                        for bb=1:12
                            BB_POP=find(bb==BINVALS.*COLL_GROUP);
                            if ~isempty(BB_POP)
                                RAND_POP=find(bb==BINVALS_PREV.*COLL_GROUP);
                                if ~isempty(RAND_POP)
                                    INDEXES_RAND_POP=randi(length(RAND_POP),length(BB_POP),1);
                                    SELECT_RAND_POP=RAND_POP(INDEXES_RAND_POP);
                                    try
                                        vXYZ_COM_BB(BB_POP,:) = 0.5.*(vXYZ_fast(BB_POP,:)+vXYZ_fast_prev(SELECT_RAND_POP,:));
                                    catch
                                        warning('BB_POP issue (again!)')
                                    end
                                end
                            end
                        end
						% neutron source in COM referential
                        theta_iso(COLL_GROUP) = 2.*rand(sum(COLL_GROUP),1)-1;
                        theta_iso(COLL_GROUP) = asin(theta_iso(COLL_GROUP));
						alpha_iso (COLL_GROUP)       =  (2*pi).*rand(sum(COLL_GROUP),1);
						vXYZ_neutrons(COLL_GROUP,1)    = cos(theta_iso(COLL_GROUP)).*cos(alpha_iso(COLL_GROUP));
						vXYZ_neutrons(COLL_GROUP,2)    = cos(theta_iso(COLL_GROUP)).*sin(alpha_iso(COLL_GROUP));
						vXYZ_neutrons(COLL_GROUP,3)    = sin(theta_iso(COLL_GROUP));
						vXYZ_neutrons(COLL_GROUP,:)    = const.v0_nDD.*vXYZ_neutrons(COLL_GROUP,:);
                        % total velocity in laboratory referential [cartesians!]
                        v_n_BB(COLL_GROUP,:) = vXYZ_neutrons(COLL_GROUP,:) + vXYZ_COM_BB(COLL_GROUP,:);
						v_rel_BB(COLL_GROUP) = sqrt(sum((vXYZ_fast(COLL_GROUP,:) - vXYZ_fast_prev(COLL_GROUP,:)).^2,2));
						
						% need also to calculate the neutron rate for BB
				        nfast(COLL_GROUP)=interp1qr(dim.pn,dim.nfast_prof,psi_gc_value(COLL_GROUP),psi_xi_pos(COLL_GROUP),psi_slope(COLL_GROUP));
                        % formula : E = EA * mB / (mA+mB)
                        E_COM_BB(COLL_GROUP)=0.5.*(0.5.*input.m/const.eV)*sum(v_rel_BB(COLL_GROUP,:).^2,2);
						E_sigv(COLL_GROUP)=E_COM_BB(COLL_GROUP)*1e-3;
					    % cross section in center of mass for D
						s_E(COLL_GROUP)=(A1_NDD+E_sigv(COLL_GROUP).*(A2_NDD+E_sigv(COLL_GROUP).*(A3_NDD+E_sigv(COLL_GROUP).*(A4_NDD+E_sigv(COLL_GROUP).*A5_NDD))));
						sv_val(COLL_GROUP)=s_E(COLL_GROUP)./(E_sigv(COLL_GROUP).*exp(31.397./sqrt(E_sigv(COLL_GROUP))));  % in millibarns (10^-31 m^2)
						sv_val(COLL_GROUP)=sv_val(COLL_GROUP).*v_rel_BB(COLL_GROUP)*1e-31;
						ndd_val_BB(COLL_GROUP)=(par.MARKER_WEIGHT.*nfast(COLL_GROUP)).*sv_val(COLL_GROUP);       % neutron rate you would have over 1s (absolute, not volumic)
						ndd_val_BB(~sn_ejected)=0;						
					end
					
					% recording the time stamp cartesian velocities for BB interations
					vXYZ_fast_prev(COLL_GROUP,:) = vXYZ_fast(COLL_GROUP,:);
					%psi_gc_value_prev(COLL_GROUP)=psi_gc_value(COLL_GROUP);
					BINVALS_PREV(COLL_GROUP)=BINVALS(COLL_GROUP);
					
					% recording fast neutrons data
					output.v_n_BP(COLL_GROUP,:,time_stamp)=v_n_BP(COLL_GROUP,:);
					output.v_n_BB(COLL_GROUP,:,time_stamp)=v_n_BB(COLL_GROUP,:);
					output.ndd_val_BP(:,time_stamp)=ndd_val_BP;
					output.ndd_val_BB(:,time_stamp)=ndd_val_BB;
					% ndd_val(COLL_GROUP)=ndd_val(COLL_GROUP)+ndd_val_BB(COLL_GROUP);
                end
                output.ndd_val(:,time_stamp)=ndd_val;
                output.ndd_cum(:,time_stamp)=ndd_cum;
            end
            if par.CALCULATE_NDD_THERMAL
                % update the output with extra ndd thermal info
                [output]=calculate_nDD_thermal(output,input,time_stamp);
                % fields required from initial particles: dim.x_ni dim.Ti_x_ni dim.ni_markers_weight               
            end
            if par.CALCULATE_CX
                disp(['length(find(ejected_CX )) = ' num2str(length(find(ejected_CX)))]);
            end
            if par.COULOMB_COLL
                disp(['length(find(ejected_sd )) = ' num2str(length(find(ejected_sd)))]);
            end
            disp(['length(find(ejected )) = ' num2str(length(find(ejected)))]);
            
            % further quantities are calculateed in evaluate_output
            output.x(:,:,time_stamp)=x;
            output.v(:,:,time_stamp)=v;
            output.x_gc(:,:,time_stamp)=x_gc;
           
            % Recorrect time-global for precision
            time=par.time_scale(time_stamp)+par.dt;
            
            if par.CALCULATE_PPHI_3D
                % Store analytic value
                output.pphi_an_temp(:,time_stamp)=Delta_pphi_an;
                Delta_pphi_an(:)=0;
            end
			if par.COULOMB_COLL & par.CALCULATE_DEVIATION
				output.deviation(:,time_stamp)=deviation;
			end
            disp(['storing x and v values at time stamp ' num2str(time_stamp) ' and time step ' num2str(time_step) '...' ])
        end
        
        %% SAVE DATA FILE intermediately
        if par.SAVE_DATA_FILE && mod(time_step,par.RECORD_PRECISION)==0 && par.SAVE_RAW_FILE
            disp('------------------------------------------------')
            toc(subtime)
            if par.mode~=5
                disp(strcat('time_step # ',num2str(time_step),' ; time_stamp # ',num2str(time_stamp),' ; time = ',num2str(par.time_scale(time_stamp))));
            else
                disp(['max phi crossings  ',num2str(max(ind_time)),' ; avg phi crossings  ',num2str(mean(ind_time(~ejected)))]);
            end
            disp(['Now at: ',num2str(100*time_step/par.NB_TIME_STEPS,'%2.0f'),' %'])
            disp(strcat('total # ejected particles =  ',num2str(sum(ejected))));
            save(par.SAVENAME_RAW,'-v7.3','input','output','par')
            
            subtime=tic;
        end
        
%         NAN_OUTPUTS(sn_ejected)=isnan(sum(v(sn_ejected,:),2)).*isnan(sum(x(sn_ejected,:),2));
%         if sum(NAN_OUTPUTS)>0
%             disp('error')
%         end

    end
    if exist('process_time_raw','var')
        process_time=process_time_raw;
    else
        process_time=toc(process_time);
    end
    if par.CALCULATE_PPHI_3D
        output.pphi_an_temp=par.dt*input.Z*cumsum(output.pphi_an_temp,2);
    end
	if par.COULOMB_COLL
	    output.ejected_sd=ejected_sd;
    end
    if par.CALCULATE_CX
	    output.ejected_CX=logical(ejected_CX);
        output.CX_NEUTRALS=CX_NEUTRALS;
    end
    output.ejected_wall=ejected_wall;
    
    disp('**************************************************************');
    disp(['Total execution time of process ',num2str(par.PROCESS_NUMBER),' is: ',num2str(process_time)]);


%% EVALUATE / SAVE DATA
if par.PROGRESSIVE_BIRTH
	par=remove_fields(par,{'birth_matrix'});              % Remove these fields
end
switch par.mode
    case 1 %TEST
        if par.SAVE_DATA_FILE && par.SAVE_RAW_FILE
            save(par.SAVENAME_RAW,'input','output','par','ejected','process_time','-v7.3') % Save in case evaluation fails
        end
            output=evaluate_output(input,output,ejected);
        if par.SAVE_DATA_FILE
            save(par.SAVENAME,'input','output','par','ejected','process_time','-v7.3')
        end
    case 2 % PREC
        if all(ejected)
            warning(['All particles ejected in precession simulation of process ',num2str(par.PROCESS_NUMBER)]);
        %    if par.SAVE_DATA_FILE
        %        save([par.SAVENAME_RAW,'prec'],'input','output','par','ejected','process_time','-v7.3')
        %    end
        %    exitcode=2;
        %    return
        end
        try
            % Get pop data
            output=evaluate_output(input,output,ejected);                    % Get more info, e.g. v_parallel
            prec = extract_precession_information_struct(input,output,ejected);    % Identify particles and orbit/precession frequency
%             output=remove_fields(output,[],{'x','v','vpll','nr_vpll_crossing','nr_midplane_crossing'});  % Remove any fields apart from...
            % the final gc position of a precession simulation is important
            % information and is kept as the new input position for the
            % next simulation
            input.x_gc_end=x_gc;
            input.x_end=x;   % squeeze(output.x_gc(:,:,end));
            input.v_end=v;		% squeeze(output.v(:,:,end));
			input=remove_fields(input,{'x_gc','v','Fc_field','x'});              % Remove these fields
			[~,B]=B_interpolation(x);
			Bfield_sq=dot(B,B,2);
			Bfield=sqrt(Bfield_sq);
			vpll=dot(B,v,2)./Bfield;
            input.vpll_end=vpll;
            input.Ekin_end=v_to_E(input.m,v);
            
            % a large amount of data is there but we keep the bare minimum
			% pphi_kin can be used to check time step precision
            output=remove_fields(output,[],{'x_gc','pphi_kin','v','time_step_loss','x_ej','loss','loss_wall','psi_outer','psi_gc_outer'});              % Remove any fields apart from...
%            	output=remove_fields(output,[],{'x','psi','v_plus_1','x_gc','pphi_kin','v','vpll','time_step_loss','x_ej','loss'});              % Remove any fields apart from...
%             	output=remove_fields(output,[],{'x_gc','pphi_kin','v','vpll','time_step_loss','x_ej','loss'});              % Remove any fields apart from...
            
            % Store the data of precession simulation in prec
            %prec.par=remove_fields(par,[],{'dt','NB_TIME_STAMPS'}); % Store NB_TIME_STAMPS and dt
            prec.par=par; % Store whole parameter set
            prec.par.end_time=par.time_scale(end);                  % Store the last tme
            prec.ejected=ejected;                                   % Store those ejected in 2D
 

			% resizing the file so that we have reasonable compromise between
			% information about what happened and size of data
			[par,output 	]=reduce_time_stamps(par,output); 
        catch err
            if ~exist([par.SAVENAME_RAW,'prec'],'file')
                warning('SAVING RAW DATA FROM PRECESSIONS IN MATPREC-FILE')
                if par.SAVE_DATA_FILE
                    save([par.SAVENAME_RAW,'prec'],'input','output','par','ejected','process_time','-v7.3')
                end
            else
                warning('RAW PRECESSION DATA FAILED TO BE ANALYZED')
            end
            rethrow(err)
        end
		
            % Save in case of only one file

            % if par.SAVE_DATA_FILE && par.NB_PROCESS==1 
            %    save(par.SAVENAME_STATS,'-v7.3','input','output','par','ejected','process_time','prec')
            % end
            % Delete a previous made .matprec-file if any
            if exist([par.SAVENAME_RAW,'prec'],'file')
                delete([par.SAVENAME_RAW,'prec'])
            end
			% split prec files
			if par.SAVE_DATA_FILE 
			    %save(par.SAVENAME_RAW,'-v7.3','input','output','par','ejected','process_time','prec')
                save(par.SAVENAME_STATS,'-v7.3','input','output','par','ejected','process_time','prec')
			    disp(['saved into ' par.SAVENAME_STATS]);
            end
			

        if all(ejected)
            warning(['All particles ejected in precession simulation of process ',num2str(par.PROCESS_NUMBER)]);
            if par.SAVE_DATA_FILE
                save([par.SAVENAME_RAW,'prec'],'input','output','par','ejected','process_time','-v7.3')
                save(par.SAVENAME_STATS,'-v7.3','input','output','par','ejected','process_time','prec')
           end
            %exitcode=2;
            %return
        end
        
        % Rerun for full simulation in mode 3
        disp('Precession information made')
        clearvars -EXCEPT par
        %exitcode=G_eq(num2str(par.PROCESS_NUMBER),num2str(par.NB_PROCESS));
        %disp('Continueing to FULL-simulation (mode 3)
        exitcode=2;
        return
    case {3,6} %FULL (after/without precession)
        % Save a RAW file before evaluation
        if par.SAVE_DATA_FILE && par.SAVE_RAW_FILE
            if exist('prec','var')
                save(par.SAVENAME_RAW,'input','output','par','ejected','process_time','prec','-v7.3')
            else
                save(par.SAVENAME_RAW,'input','output','par','ejected','process_time','-v7.3')
            end
        end
        
        % Evaluate the raw-data
        output_raw=output;
        output=evaluate_output(input,output,ejected); 
%       	output=remove_fields(output,[],{'ST_interaction','nr_vpll_crossing','nr_midplane_crossing','time_step_loss','loss','x_ej','v_ej','pphi_kin','mm','Delta_pphi','pphi_an','dpphi_dt','dmm_dt','nr_profile','wb','wd'}); % Remove everything but these fields.
%       	output=remove_fields(output,[],{'x','v','vpll','x_gc','ST_interaction','time_step_loss','loss','x_ej','mm_ej','vpll_ej','pphi_kin','mm','Delta_pphi','pphi_an','dpphi_dt','dmm_dt','nr_profile','wb','wd'}); % Remove everything but these fields.
        % next simulation
        input.x_gc_end=x_gc;
        input.x_end=x;   % squeeze(output.x_gc(:,:,end));
        input.v_end=v;		% squeeze(output.v(:,:,end));
		input=remove_fields(input,{'x_gc','v','Fc_field','x'});              % Remove these fields

		FIELDS_TO_SAVE={'x','v','vpll','x_gc','time_step_loss','loss','loss_wall','deviation',...
		 'x_ej','x_ej_next','Ekin_ej','vpll_ej','pphi_kin','mm',...
		 'psi_value_avg','ndd_val','ndd_cum','Edep_th','Pdep_e','Pdep_i','Etot_e','Etot_i','delta_Ekin','cum_tor_ang_mom',...
         'Delta_pphi','pphi_an','dpphi_dt','ejected_wall','ejected_sd'};
        if par.RECORD_LOSS_TRAJECTORY
             FIELDS_TO_SAVE{end+1}='x_ej_prev';
             FIELDS_TO_SAVE{end+1}='x_ej_prev2';
             FIELDS_TO_SAVE{end+1}='x_ej_prev3';
        end
        if par.CALCULATE_CX
             FIELDS_TO_SAVE{end+1}='x_CXn';
             FIELDS_TO_SAVE{end+1}='x_CXi';
             FIELDS_TO_SAVE{end+1}='CX_rate';
             FIELDS_TO_SAVE{end+1}='CX_flag';
             FIELDS_TO_SAVE{end+1}='ejected_CX';
             FIELDS_TO_SAVE{end+1}='CX_NEUTRALS';
        end
        if par.CALCULATE_NDD
             FIELDS_TO_SAVE{end+1}='ndd_val';
             FIELDS_TO_SAVE{end+1}='ndd_cum';
             if par.RECORD_VSPACE_NDD
                 FIELDS_TO_SAVE{end+1}='v_n_BP';
                 FIELDS_TO_SAVE{end+1}='v_n_BB';
                 FIELDS_TO_SAVE{end+1}='ndd_val_BP';
                 FIELDS_TO_SAVE{end+1}='ndd_val_BB';
             end
             if par.CALCULATE_NDD_THERMAL
                 FIELDS_TO_SAVE{end+1}='x_n_thermal';
                 FIELDS_TO_SAVE{end+1}='v_n_thermal';
                 FIELDS_TO_SAVE{end+1}='ndd_val_thermal';
             end
        end
        if par.IMPURITY_ZAVG>0
            FIELDS_TO_SAVE{end+1}='Pdep_imp';
            FIELDS_TO_SAVE{end+1}='Etot_imp';
        end
        output=remove_fields(output,[],FIELDS_TO_SAVE);

		input=remove_fields(input,{'x_gc','v','Fc_field'});              % Remove some fields ...
        
        % Reduce the number of time stamps
		[~ 	,output_raw	]=reduce_time_stamps(par,output_raw);
        [par,output 	]=reduce_time_stamps(par,output);  %#ok<ASGLU>
        
        % Save the data definitely
        if par.SAVE_DATA_FILE
            if exist('prec','var')
                save(par.SAVENAME,'input','output','par','ejected','process_time','prec','-v7.3')
            else
                save(par.SAVENAME,'input','output','par','ejected','process_time','-v7.3')
            end
        end
        
        output=output_raw;
        if par.SAVE_DATA_FILE && par.SAVE_RAW_FILE
            if exist('prec','var')
                save(par.SAVENAME_RAW,'input','output','par','ejected','process_time','prec','-v7.3')
            else
                save(par.SAVENAME_RAW,'input','output','par','ejected','process_time','-v7.3')
            end
        end
    case 5 %% POINCARE
        if par.SAVE_DATA_FILE
            save(par.SAVENAME_RAW,'input','output','par','ejected','process_time','-v7.3')
        end
        % Remove not used / obsolete fields in output array
        output.x_gc(:,:,max(ind_time)+1:end)=[];
        % keep trajectories for debugging and consistency
        %output.x(:,:,max(ind_time)+1:end)=[];
        %output.v(:,:,max(ind_time)+1:end)=[];
        par.time_scale(:,max(ind_time)+1:end)=[];
        
        % Indexes 2D
        X_ind=squeeze((output.x_gc(:,1,:)-dim.R0)*dim.DX_inv)+dim.mid_Xzero;
        Z_ind=squeeze( output.x_gc(:,2,:)        *dim.DZ_inv)+dim.mid_Z;
        
        %epxr_nan=expr_nan(:);
        
        % Store psi and theta in output-struct
        % Psi-surfaces
        output.psi=NaN(size(output.x_gc,1),size(output.x_gc,3));
        output.theta=output.psi;
        
		for part_index=1:size(X_ind,1)
			% Find theta
			expr_nan=isnan(X_ind(part_index,:));
			output.theta(part_index,~expr_nan) = interpolate_theta_XZ(X_ind(part_index,~expr_nan),Z_ind(part_index,~expr_nan));        
			% Find psi (index)
			output.psi(part_index,~expr_nan)=ba_interp2(maps(1).psi_XZ,Z_ind(part_index,~expr_nan),X_ind(part_index,~expr_nan),'cubic');
		end
        if par.SAVE_DATA_FILE
            save(par.SAVENAME,'input','output','par','ejected','process_time','-v7.3')
        end
    otherwise
        if par.SAVE_DATA_FILE
            save(par.SAVENAME_RAW,'input','output','par','ejected','process_time','-v7.3')
        end
        error('MODE NOT KNOWN OR OBSOLETE / SAVING RAW DATA-FILE')
end

%% EXIT
disp('**************************************************************');
disp('done');

exitcode=0; % Notation: exit prevents exitcode to work
disp('exitcode = ')
disp(exitcode)
end


%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_sat LOSS_POP_OUT]=apply_domain_limits(x,sn_ejected,dim)
        LOSS_POP_OUT=logical(sn_ejected*0);
        LOSS_POP_OUT(sn_ejected)=logical(x(sn_ejected,2) <= dim.Zlow_lim | x(sn_ejected,2) >= dim.Zup_lim |...
            x(sn_ejected,1) <= dim.Rlow_lim | x(sn_ejected,1) >= dim.Rup_lim );
        PHI_POP_OUT=logical(x(sn_ejected,3) <= 0 | x(sn_ejected,3) >= 2*pi);
		x(LOSS_POP_OUT,2)=max(x(LOSS_POP_OUT,2),dim.Zlow_lim);
		x(LOSS_POP_OUT,2)=min(x(LOSS_POP_OUT,2),dim.Zup_lim);
		x(LOSS_POP_OUT,1)=max(x(LOSS_POP_OUT,1),dim.Rlow_lim); 
		x(LOSS_POP_OUT,1)=min(x(LOSS_POP_OUT,1),dim.Rup_lim);	
        x(PHI_POP_OUT,3)=mod(x(PHI_POP_OUT,3),2*pi);
		x_sat=x;
end


