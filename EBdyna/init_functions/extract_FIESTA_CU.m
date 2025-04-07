time_stamp=1550 % ms

if exist('fiesta_equilibrium')==0
    run('/compass/Shared/Common/IT/projects/fiesta/releases/V8.10/startup.m');
end
FIESTA_STRING='/compass/Shared/Common/COMPASS-UPGRADE/RP1_Design/Scenarios/32.4/Fiesta_32430/ITERlike_scenario_32430_v42_'

load([FIESTA_STRING num2str(time_stamp) 'ms.mat']);
disp('To remember to convert to tokamak coordinates for the magnetic field and poloidal flux !')
Sign_Ip=-1  % current should be ccw from top
Sign_Bt=-1  % I-mode with reversed (ccw) Bt 
FIESTA=struct();
FIESTA.flipud=0;

% get(equil)
grid=get(equil,'grid')
profiles=qprofile(equil);
profiles=rmfield(profiles,'data')

fn = fieldnames(profiles);
for k=1:numel(fn)
    prof=profiles.(fn{k});
    prof(1)=[];
    profiles.(fn{k})=prof;
end

ip=get(equil,'ip')
F2=get(equil,'F2')
jprofile=get(equil,'jprofile')
psi_axis=get(equil,'psi_axis')
psi_boundary=get(equil,'psi_boundary')
p_axis=get(equil,'p_axis')

Psi=get(equil,'Psi');
Psi_n=get(equil,'Psi_n');

Pressure=get(equil,'Pressure');
Bphi=get(equil,'Bphi');
Br=get(equil,'Br');
Bz=get(equil,'Bz');
jphi=get(equil,'J');

xpoint=get(equil,'xpoint');
boundary=get(equil,'boundary');

[rc, zc] = findclosed(equil, [psi_boundary]);

load([FIESTA_STRING num2str(time_stamp) 'ms_params.mat']);
%                irod: 2.2350e+07
%                itot: 1.9999e+06
%                psi0: -0.7339
%                psip: 4.1876
%                psia: -2.2576
%                 lif: 1.3515
%                 li1: 0.9099
%                 li2: 1.0267
%                 li3: 0.8214
%                  lp: 2.0939e-06
%                  r0: 0.8972
%                  z0: 0.0334
%                r0_j: 0.8754
%                z0_j: 0.0319
%             r0_geom: 0.8877
%             z0_geom: 0.0049
%               drsep: [NaN NaN]
%         aspectratio: 3.3864
%     aspectratio_fit: 3.3911
%              kappa0: 1.4210
%             kappa50: 1.4524
%             kappa95: 1.6853
%               kappa: 1.8124
%           kappa_fit: 1.8039
%                 f50: 4.5231
%           trapped50: 0.5736
%         delta_upper: 0.4583
%         delta_lower: 0.4583
%           delta_fit: 0.4150
%             delta95: 0.3438
%               decay: -1.1015
%                  j0: 1.1439e+07
%                  b0: 5.3308
%                  P0: 3.7475e+05
%                  q0: 0.8783
%                 q95: 2.4438
%                area: 0.3514
%              volume: 1.9048
%             surface: 12.2973
%                 rin: 0.6256
%                rout: 1.1499
%               betap: 0.3621
%               betat: 0.0168
%               betan: 0.0111
%              energy: 4.9764e+05

FIESTA.time_equil=time_stamp;
FIESTA.params=params;
FIESTA.prof=profiles;

FIESTA.ip=ip;
FIESTA.psi_axis=-Sign_Ip*psi_axis;
FIESTA.psi_boundary=-Sign_Ip*psi_boundary;
FIESTA.p_axis=p_axis;

FIESTA.R_xpoint=get(xpoint,'r');
FIESTA.Z_xpoint=get(xpoint,'z');

% FIESTA.R_boundary=get(boundary,'r');
% FIESTA.Z_boundary=get(boundary,'z');
FIESTA.R_boundary=rc;
FIESTA.Z_boundary=zc;



FIESTA.scale_R=get(grid,'r');
FIESTA.scale_Z=get(grid,'z');

FIESTA.F2=get(F2,'data','2D')';
FIESTA.Psi=-Sign_Ip*get(Psi,'data','2D')';
FIESTA.Psi_n=get(Psi_n,'data','2D')';
FIESTA.Bphi=Sign_Bt.*get(Bphi,'data','2D')';
FIESTA.Br=-Sign_Ip*get(Br,'data','2D')';
FIESTA.Bz=-Sign_Ip*get(Bz,'data','2D')';
FIESTA.jphi=get(jphi,'data','2D')';

FIESTA.params.psi0=-Sign_Ip*FIESTA.params.psi0;
FIESTA.params.psia=-Sign_Ip*FIESTA.params.psia;
FIESTA.params.psip=-Sign_Ip*FIESTA.params.psip;
FIESTA.prof.psi=-Sign_Ip.*FIESTA.prof.psi;
FIESTA.prof.jbar=-Sign_Ip.*FIESTA.prof.psi;
FIESTA.prof.jboot_int=-Sign_Ip.*FIESTA.prof.jboot_int;

FIESTA.params.b0=Sign_Bt.*FIESTA.params.b0;
FIESTA.prof.f=Sign_Bt.*FIESTA.prof.f;

FIESTA.sign_Ip=Sign_Ip;
FIESTA.sign_Bt=Sign_Bt;

save FIESTA_equil.mat FIESTA

