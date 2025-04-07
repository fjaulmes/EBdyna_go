

load('./NUR_data/NUR_out_CU_2NBI_phi45_Rt40.mat')


load('./data_tokamak/physics_constants.mat');
mass=mD
charge=1

SUBPOP=1
% INIPOP allows for new births during 3D simulation
INIPOP=1




Rm=NUR_out.R_i;
Zm=NUR_out.Z_i;
phim=NUR_out.phi_i;

vR=NUR_out.vR;
vZ=NUR_out.vZ;
vphi=NUR_out.vphi;

vnorm=sqrt(vR.^2+vZ.^2+vphi.^2);

E0=80*1e3;
v80keV=sqrt(2*eV*E0/(mD));


Rpos=squeeze(Rm(1,:));
Zpos=squeeze(Zm(1,:));
phipos=squeeze(phim(1,:));


% sub-sampling for short checks

vR=vR(INIPOP:SUBPOP:end);
vZ=vZ(INIPOP:SUBPOP:end);
vphi=vphi(INIPOP:SUBPOP:end);

Rpos=Rpos(INIPOP:SUBPOP:end);
Zpos=Zpos(INIPOP:SUBPOP:end);
phipos=phipos(INIPOP:SUBPOP:end);

% reshuffle everything since it is not systematically done in NUR
SHUFFLE_POS=(randperm(length(Rpos)));
Rpos=Rpos(SHUFFLE_POS);
Zpos=Zpos(SHUFFLE_POS);
phipos=phipos(SHUFFLE_POS);
vR=vR(SHUFFLE_POS);
vZ=vZ(SHUFFLE_POS);
vphi=vphi(SHUFFLE_POS);


input=struct();
input.x=[Rpos ; Zpos ; phipos]';
input.v=[vR ; vZ ; vphi]';
input.N_total=length(input.x);
input.m=mass;
input.Z=charge;

% outlier
Ekin_recalc=sum(input.v.^2,2);
Ekin_recalc=0.5*mD/eV*Ekin_recalc;

OUTLIER_BUGGY=find(Ekin_recalc<1e4)

if ~isempty(OUTLIER_BUGGY)
    GOOD_DATA=find(Ekin_recalc>1e4);
    input.x=input.x(GOOD_DATA,:);
    input.v=input.v(GOOD_DATA,:);
    input.N_total=length(input.x);
    warning('removed some buggy data from NUR output, careful!')
end

Ekin_recalc=sum(input.v.^2,2);
Ekin_recalc=0.5*mD/eV*Ekin_recalc;


save('./G_eq3D/input/NUR_markers_R40_600k.mat','input')
input

xgc=Rpos.*cos(phipos);
ygc=-Rpos.*sin(phipos);


%%
load('../../data_common/wallRZ_CU.mat');
    
figure;
hold on

plot(xgc,ygc,'r.')

% quiver(xgc(1:10:end),ygc(1:10:end),-vR(1:10:end).*cos(phipos(1:10:end))+vphi(1:10:end).*sin(phipos(1:10:end)),-vphi(1:10:end).*cos(phipos(1:10:end))-vR(1:10:end).*sin(phipos(1:10:end)),16)
quiver(xgc(1:10:end),ygc(1:10:end),vR(1:10:end).*cos(phipos(1:10:end))-vphi(1:10:end).*sin(phipos(1:10:end)),-vphi(1:10:end).*cos(phipos(1:10:end))-vR(1:10:end).*sin(phipos(1:10:end)),16)

plot(xgc(1:10:end),ygc(1:10:end),'r.')

plot_circle(0,0,Plasma.Rwall_min);
plot_circle(0,0,max(wall_CU.R));
axis equal

%%
figure
plot_CU_TOK_background
plot(Rpos,Zpos,'r.')
% plot(NUR_out.R_i,NUR_out.Z_i,'b.');hold on

%%
figure
hist(Ekin_recalc,81);



