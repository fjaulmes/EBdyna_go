
FILENAME
% load('./input_NUR/NUR_out_CU_2NBI_phi45_Rt40.mat')


load([DATA_PHYS_FOLDER,'physics_constants.mat']);
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



Ekin_recalc=sum(input.v.^2,2);
Ekin_recalc=0.5*mD/eV*Ekin_recalc;

save([DATA_INPUT_FOLDER,'NUR_markers_R0_300k.mat'],'input')
input

xgc=Rpos.*cos(phipos);
ygc=-Rpos.*sin(phipos);


%%
load('./data_common/compassu/wallRZ_CU.mat');

figure;
quiver(xgc,ygc,vR.*cos(phipos)-vphi.*sin(phipos),vphi.*cos(phipos)-vR.*sin(phipos),2)

hold on
plot(xgc,ygc,'r.')

plot_circle(0,0,min(wall_CU.R));
plot_circle(0,0,max(wall_CU.R));

%%
figure
plot_CU_TOK_background
plot(Rpos,Zpos,'r.')
% plot(NUR_out.R_i,NUR_out.Z_i,'b.');hold on

%%
figure
hist(Ekin_recalc,40);



