function define_basic_sim_parameters
global par

    % Please provide a unique ID 
    par.comment='COMPASS Upgrade NBI CX simulation';
    par.IMPURITY_ZAVG=0;        % set to 0 if impurity profiles not considered
    % par.mimp=  1.994307e-26;    % mass of Carbon ion

    %DEFINE_BASIC_SIM_PARAMETERS Summary of this function goes here
    par.PSI_LIMIT_EJECTED    = -0.35;
    par.TOL_PREC_DPPHI       = 0.02;   % sets up some tolerance to identify particles having "difficult" orbits (eg. close to edge)
    par.TOL_OUTER_HFS_THETA  = 0.22;   % sets up some tolerance for HFS crossing
    par.TOL_CORE_HFS_THETA   = 0.08;   % sets up some tolerancefor HFS crossing
    par.RELATIVISTIC_CORR    = 0;       % for kinetic energies beyond the speed of light

    par.COULOMB_COLL         = 1;       % for collisional effects using METIS profiles in 'pressure_profile.mat'
    par.SLOWING_DOWN         = 1;       % activate decrease of velocity for fast particles simulations
    par.PROGRESSIVE_BIRTH    = 1;       % used to split the input into smaller chunks inserted at equal time intervalls during the simulation
    par.NB_BIRTH_CHUNKS      = 100;
    par.CALCULATE_DEVIATION  = 0;
    par.SD_MARKERS_END       = 0;

    par.RECORD_AFTER_LOSS    = 1;      % number of time steps recorded after loss to the wall (if=0, only one value is recorded)

    par.INPUT_POWER          = 2e6;    % total injected power (for NBI simulation) - used in ndd and pdep calculation
    par.CALCULATE_NDD        = 0;       % apply calculation of DD neutron yield for the fast population against thermal D background
    par.CALCULATE_PDEP       = 0 &  (par.SLOWING_DOWN);       % yields Power deposited info (a bit demanding computationnaly)

    par.CALCULATE_CX         = 1;       % yields Power deposited info (a bit demanding computationnaly)
    par.N0_FAC_D             = 10.0;       % multiplicative factor on background D neutral density
    par.N0_FAC_D2            = 10.0;       % multiplicative factor on background D2 neutral density
    par.ALLOW_REIONIZATIONS  = 1;       % default = 1 [for testing of the multiples ionization functionality]
    par.USE_T0_TABLE         = 1;       % default = 1 [for testing of the neutral temperature table]

    par.USE_1DGC_POS         = 0;       % default = 0 [switch to GC value for 1D interpolation and psi avg calculations : found inaccurate.......]
    par.INPUT_TRUEPOS        = 1;       % default = 1 [used to force the code to use real positions, not gc, , as input]


    % Please provide the interpolation scheme   (advice: 3)
    % 0     -   self-written linear, with vector potential
    % 1     -   self-written linear
    % 2     -   interp linear               (ML standard)
    % 3     -   ba_interp linear            (File Exchange)
    % 4     -   griddedInterpolant linear   (makes use of lookup table)
    % 5     -   interp cubic                (ML standard)
    % 6     -   ba_interp cubic             (File Exchange)
    % 7     -   griddedInterpolant cubic    (makes use of lookup table)
    % 8     -   interp2 spline              (ML standard)
    % 9     -   griddedInterpolant spline   (makes use of lookup table - RAM intensive!)
    par.interp_scheme=3;

    switch par.interp_scheme
        case {2,5,8}
            error('The interpolation scheme isn''t allowed, since it is too slow');
        case {9}
            error('The interpolation scheme isn''t allowed, since the RAM required is very high and cannot be in shared RAM');
    end


    % Please enter scheme for the particle pusher (FabienB for optimized BORIS, BORIS with vector-operations, "old" Fabien approximation)
    % BORIS     - one-lined BORIS scheme without electric field  -- Better with interp schemes: 3 or 6
    % FabienB   - BORIS scheme based on older Fabien code which might work faster sometimes -- Better with interp schemes: 1 2 4 5 7 8 9
    switch par.interp_scheme
        case {3,6}
            par.scheme='BORIS';
        otherwise
            par.scheme='FabienB';
    end
    % par.scheme='FabienB';

    % Name to save files
    par.ID_NAME=strcat('G_eq_',par.ID);     % Name to save files
    par.coord_syst='toroidal';              % Switch for 3D coordinate system. Use 'flux': (theta,psi,phi) or 'toroidal': (R,Z,phi)

    % 3D field
    par.superimpose_2D_3D	=false; 					% Superimpose (if possible) the 2D and 3D fields, making a 2D interpolation obsolete and speeding up EBdyna_go. Disable to avoid numerical differences.

    par.APPLY_RMP           =false;                  % Resonant Magnetic Perturbations (RMP)
    if par.APPLY_RMP
        par.RMP_file        =['../data_tokamak/MARSF_15611_bres.mat']; 
    end
    par.APPLY_TFR           =false;                 % Toroidal Field Ripple
    if par.APPLY_TFR
        par.TFR_file        =['../data_tokamak/','TFR_',par.coord_syst,'_3D.mat']; % TFR_toroidal
        % difference between two files????
    %     par.TFR_file        =['../data_tokamak/','TF_',par.coord_syst,'_2019-04-26.mat'];
    end

    if (par.APPLY_RMP || par.APPLY_TFR) & par.mode~=2
        par.APPLY_3D=true;
    else
        par.APPLY_3D=false;
    end
    par.CALCULATE_PPHI_3D   =true...               % Determine pphi based on addition of (local) evolution
        & (par.APPLY_3D);

    % Sawtooth
    par.APPLY_SAWTOOTH      = false;
    if par.APPLY_SAWTOOTH
        par.st.t_reconnection = 1e-4;
        par.st.t_relaxation   = 1.8e-4;
        par.st.n_reconnection = 1e3;
        par.st.n_relaxation   = 1e2;
        par.st.DR=1e-9;
        par.st.DZ=1e-9;
        par.st.stable_time    = 10; % The first time point in the sawtooth that should be taken into account! Condition: (k-kr)*c > -1 and no numeric oscillations
        par.st.rec_end_time   = 3; % Parameter from which the reconnection should get kappa_c from interpolation rather than formula.
        if ~any(par.interp_scheme==[3,6]); error('Sawtooth only compatible with ba_interp-interpolation'); end;
        par.st.show_parameters = false & nargin==0;
    end

    % Simulation type / distribution
    par.GET_PREC_INFO       =false;                 % Find precession data in original 2D-field (trapped or not etc.)
    par.SAVE_DATA_FILE      =true;                  % Save output files (raw + full + possible prec)
    par.SAVE_RAW_FILE       =false;                 % Remove the infor from raw files (unnecessary in most cases)





    par.Poincare_plot=false;
    par.calculate_length_trajectory=false;

 
end

