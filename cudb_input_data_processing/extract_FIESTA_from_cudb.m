time_stamp = TIME_EQUIL;                                   % Time stamp in seconds
SC_STRING = num2str(SHOT_NUMBER);                                 % Scenario identifier

if exist('fiesta_equilibrium') == 0
    run('/compass/Shared/Common/COMPASS-UPGRADE/RP1_Design/Scenarios/startup_FIESTA.m');
end

FIESTA = struct();
FIESTA.SC_STRING = SC_STRING;                       % Store scenario ID
FIESTA.flipud = 0;                                  % Flag to flip up/down (if needed)
FIESTA.time = time_stamp;                           % Store the selected time stamp

cdb = cdb_client();                                 % Connect to CDB
signal = cdb.get_signal(['a/Fiesta_OUT:' SC_STRING]);
FIESTA_TIME = signal.time_axis.data;
FIESTA.a = interp1(FIESTA_TIME, signal.data, time_stamp, 'nearest');

signal = cdb.get_signal(['Bt_vac_mag_axis/Fiesta_OUT:' SC_STRING]);
FIESTA.Bt_vac_mag_axis = interp1(FIESTA_TIME, signal.data, time_stamp, 'nearest');

signal = cdb.get_signal(['r_geom/Fiesta_OUT:' SC_STRING]);
FIESTA.r_geom = interp1(FIESTA_TIME, signal.data, time_stamp, 'nearest');

signal = cdb.get_signal(['kappa/Fiesta_OUT:' SC_STRING]);
FIESTA.kappa = interp1(FIESTA_TIME, signal.data, time_stamp, 'nearest');

signal = cdb.get_signal(['z_mag/Fiesta_OUT:' SC_STRING]);
FIESTA.z_mag = interp1(FIESTA_TIME, signal.data, time_stamp, 'nearest');

signal = cdb.get_signal(['r_mag/Fiesta_OUT:' SC_STRING]);
FIESTA.r_mag = interp1(FIESTA_TIME, signal.data, time_stamp, 'nearest');

signal = cdb.get_signal(['betan/Fiesta_OUT:' SC_STRING]);
FIESTA.betan = interp1(FIESTA_TIME, signal.data, time_stamp, 'nearest');

signal = cdb.get_signal(['betap/Fiesta_OUT:' SC_STRING]);
FIESTA.betap = interp1(FIESTA_TIME, signal.data, time_stamp, 'nearest');

signal = cdb.get_signal(['betat/Fiesta_OUT:' SC_STRING]);
FIESTA.betat = interp1(FIESTA_TIME, signal.data, time_stamp, 'nearest');

signal = cdb.get_signal(['Iplasma/Fiesta_OUT:' SC_STRING]);
FIESTA.ip = interp1(FIESTA_TIME, signal.data, time_stamp, 'nearest');

signal = cdb.get_signal(['psi_0/Fiesta_OUT:' SC_STRING]);
FIESTA.psi_axis = interp1(FIESTA_TIME, signal.data, time_stamp, 'nearest');

signal = cdb.get_signal(['psi_boundary/Fiesta_OUT:' SC_STRING]);
FIESTA.psi_boundary = interp1(FIESTA_TIME, signal.data, time_stamp, 'nearest');

signal = cdb.get_signal(['p/Fiesta_OUT:' SC_STRING]);
FIESTA.map2D.pressure = squeeze(interp1(FIESTA_TIME, signal.data, time_stamp, 'nearest'));
FIESTA.P_axis = interp2(signal.axis1.data, signal.axis2.data, FIESTA.map2D.pressure', FIESTA.r_mag, FIESTA.z_mag);
FIESTA.map2D.scale_R = signal.axis1.data;
FIESTA.map2D.scale_Z = signal.axis2.data;

signal = cdb.get_signal(['psi/Fiesta_OUT:' SC_STRING]);
FIESTA.map2D.psi = squeeze(interp1(FIESTA_TIME, signal.data, time_stamp, 'nearest'));

signal = cdb.get_signal(['psi_n/Fiesta_OUT:' SC_STRING]);
FIESTA.map2D.psi_n = squeeze(interp1(FIESTA_TIME, signal.data, time_stamp, 'nearest'));

signal = cdb.get_signal(['Bphi/Fiesta_OUT:' SC_STRING]);
FIESTA.map2D.Bphi = squeeze(interp1(FIESTA_TIME, signal.data, time_stamp, 'nearest'));

signal = cdb.get_signal(['Br/Fiesta_OUT:' SC_STRING]);
FIESTA.map2D.Br = squeeze(interp1(FIESTA_TIME, signal.data, time_stamp, 'nearest'));

signal = cdb.get_signal(['Bz/Fiesta_OUT:' SC_STRING]);
FIESTA.map2D.Bz = squeeze(interp1(FIESTA_TIME, signal.data, time_stamp, 'nearest'));

signal = cdb.get_signal(['f/Fiesta_OUT:' SC_STRING]);
FIESTA.map2D.F = squeeze(interp1(FIESTA_TIME, signal.data, time_stamp, 'nearest'));

signal = cdb.get_signal(['xp_lower_r/Fiesta_OUT:' SC_STRING]);
FIESTA.xp_lower_r = interp1(FIESTA_TIME, signal.data, time_stamp, 'nearest');

signal = cdb.get_signal(['xp_lower_z/Fiesta_OUT:' SC_STRING]);
FIESTA.xp_lower_z = interp1(FIESTA_TIME, signal.data, time_stamp, 'nearest');

signal = cdb.get_signal(['xp_upper_r/Fiesta_OUT:' SC_STRING]);
FIESTA.xp_upper_r = interp1(FIESTA_TIME, signal.data, time_stamp, 'nearest');

signal = cdb.get_signal(['xp_upper_z/Fiesta_OUT:' SC_STRING]);
FIESTA.xp_upper_z = interp1(FIESTA_TIME, signal.data, time_stamp, 'nearest');

if isnan(FIESTA.xp_upper_z)
    FIESTA.R_xpoint = FIESTA.xp_lower_r;            % Use lower x-point if upper is NaN
    FIESTA.Z_xpoint = FIESTA.xp_lower_z;
else
    FIESTA.R_xpoint = FIESTA.xp_upper_r;
    FIESTA.Z_xpoint = FIESTA.xp_upper_z;
end

signal = cdb.get_signal(['pprime/Fiesta_OUT:' SC_STRING]);
FIESTA.prof.pprime = squeeze(interp1(FIESTA_TIME, signal.data, time_stamp, 'nearest'));
FIESTA.prof.psiN = signal.axis1.data;

signal = cdb.get_signal(['ffprime/Fiesta_OUT:' SC_STRING]);
FIESTA.prof.ffprime = squeeze(interp1(FIESTA_TIME, signal.data, time_stamp, 'nearest'));

signal = cdb.get_signal(['f/Fiesta_OUT:' SC_STRING]);
FIESTA.prof.f = squeeze(interp1(FIESTA_TIME, signal.data, time_stamp, 'nearest'));

signal = cdb.get_signal(['q/Fiesta_OUT:' SC_STRING]);
FIESTA.prof.q = squeeze(interp1(FIESTA_TIME, signal.data, time_stamp, 'nearest'));

FIESTA.prof.psi = FIESTA.prof.psiN .* (FIESTA.psi_boundary - FIESTA.psi_axis) + FIESTA.psi_axis;  % Rescale psi

signal = cdb.get_signal(['boundary_closed_r/Fiesta_OUT:' SC_STRING]);
FIESTA.boundary_closed_r = squeeze(interp1(FIESTA_TIME, signal.data, time_stamp, 'nearest'));

signal = cdb.get_signal(['boundary_closed_z/Fiesta_OUT:' SC_STRING]);
FIESTA.boundary_closed_z = squeeze(interp1(FIESTA_TIME, signal.data, time_stamp, 'nearest'));

FIESTA.sign_Ip = sign(FIESTA.ip);                  % Current sign
FIESTA.sign_Bt = sign(FIESTA.Bt_vac_mag_axis);     % Field sign

save ([DATA_PLASMA_FOLDER,'FIESTA_equil.mat'],'FIESTA')                       % Save structure to MAT file

