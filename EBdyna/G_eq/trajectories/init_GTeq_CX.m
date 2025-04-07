    
    COLL_GROUP_D2  = ejected*0;

    % on neutral D
    A1_CX=3.2345;
    A2_CX=235.88;
    A3_CX=0.038371;
    A4_CX=3.8068e-6;
    A5_CX=1.1832e-10;
    A6_CX=2.3713;
    % on neutral D2
    b1_CX=17.3;
    b2_CX=105;
    b3_CX=2;
    b4_CX=10^4;
    b5_CX=-1.12;
    b6_CX=3.64*1e-4;
    b7_CX=0.9;
    b8_CX=5.03*1e-19;
    b9_CX=4;
    b10_CX=5.87*1e-28;
    b11_CX=5.5;
    
    %impact ionization on D+

    Eb_ionization=[0 ...
    5*10^2 ...
    1*10^03  ...
    2*10^03  ...
    5*10^03 ...
    1*10^04 ...
    2*10^04 ...
    5*10^04 ...
    1*10^05 ...
    2*10^05 ...
    5*10^05 ...
    1*10^06 ...
    2*10^06 ...
    5*10^06 ]'*2; % Deuterium : isotope mass = 2
    sigma_ionization_val=[0 ...
    1.46*10^-20 ...
    1.46*10^-19 ...
    1.02*10^-18 ...
    6.24*10^-18 ...
    1.94*10^-17 ...
    6.73*10^-17 ...
    1.43*10^-16 ...
    1.10*10^-16 ...
    6.99*10^-17 ...
    3.48*10^-17 ...
    1.94*10^-17 ...
    1.05*10^-17 ...
    4.62*10^-18]'*1e-4;  % (m^{2})    
%     Eb_pos   = ejected*0;
%     dEb_pos  = ejected*0;
%     Eb_slope = ejected*0;

    ni               = ejected*0;
    Ti                  = ejected*0;
    v_norm          = ejected*0;
    neutral_temp    = ejected*0;
    neutral_density = ejected*0;
    neutral_density_D2 = ejected*0;
    sigma_cx        = ejected*0;
    sigma_cx_D2     = ejected*0;
    CX_rate         = ejected*0;
    CX_cum          = ejected*0;
    CX_val          = 1.0+ejected*0;   
    NCX_rate         = ejected*0;     % for neutral CX rate
    NCX_cum          = ejected*0;     % for neutral life estiamtes
    NCX_val          = 1.0+ejected*0;     % for neutral life estiamtes
    NCX_THRESH_EJECTED         = rand(length(ejected),1);
    CX_THRESH_EJECTED= rand(length(ejected),1);   % we accumulate the CX collision until it reaches thresh then the marker is ejected
    CX_flag         = ejected*0; % flag to count # of neutralizations
    output.loss_CX    = zeros(par.NB_TIME_STEPS,1);
    CX_NEUTRALS      = logical(ejected*0);    % special motion for neutralized particles
    CX_NEUTRALS_BABIES = logical(ejected*0);  % temporary variable for newly born neutrals
    CX_IONS_BABIES = logical(ejected*0);  % temporary variable for newly born neutrals
    IONS_GROUP = logical(ejected*0+1);        % temporary variable for newly born neutrals
    NEUTRALS_GROUP = logical(ejected*0);      % temporary variable for newly born neutrals
    
    if par.USE_T0_TABLE
         E_sigv_log = ejected*0;  
         Elog_ind   = ejected*0;  
         T0_ind     = ejected*0;  
         Ti_ind     = ejected*0; 
    end