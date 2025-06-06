%test of suzuki crossection

fileID = fopen('/compass/home/zadvitskiy/graph_suzuki.dat');
paper=textscan(fileID,'%f%f','Delimiter',',');

Ea=logspace(1,3,20);
ne=[1e19 1e20 1e21];

BT=1;
Ap=1;

IMP=1;

for i=1:20
    
    E=Ea(i);
    if E<200
        Ti=ones(1,3)*E/10;
        Te=Ti;
    else
        Te=ones(1,3)*20;
        Te=Ti;
    end
    
    Zeff=ones(size(ne));

if (E<1e2)
    if (Ap==1)
        A=[-5.29e1 -1.36 7.19e-2 1.37e-2 4.54e-1 4.03e-1 -2.2e-1 6.66e-2 -6.77e-2 -1.48e-3];
    elseif (Ap==2)
        if (BT==5)
            A=[-6.79e1 -1.22 8.14e-2 1.39e-2 4.54e-1 4.65e-1 -2.73e-1 7.51e-2 -6.3e-2 -5.08e-4];
        elseif (BT==1)
            A=[-6.98e1 -1.21 8.33e-2 1.35e-2 4.45e-1 4.89e-1 -2.90e-1 7.86e-2 -6.3e-2 -5.12e-4];  
        else 
            error('B must be aproximated by 1 or 5 [1 5]')
        end
    else
    error('wrong Ab value')
    end
elseif (E>=1e2)
    if (Ap==1)
        A=[1.27e1 1.25 4.52e-1 1.05e-2 5.47e-1 -1.02e-1 3.60e-1 -2.98e-2 -9.59e-2 4.21e-3];
    elseif (Ap==2)
        if (BT==5)
            A=[1.41e1 1.11 4.08e-1 1.05e-2 5.47e-1 -4.03e-2 3.45e-1 -2.88e-2 -9.71e-2 4.74e-3];
        elseif (BT==1)
            A=[-2.41e1 -1.30 -1.54e-1 8.02e-3 4.81e-1 -1.49e-1 3.92e-1 -2.99e-2 -9.76e-2 4.79e-3];  
        else 
            error('B must be aproximated by 1 or 5 [1 5]')
        end
    else
    error('wrong Ab value')
    end
end

if (IMP==1)%for carbon and deuterium plasmas
B(1,1,1)=-1;
B(1,1,2)=-2.55e-2;
B(1,2,1)=-1.25e-1;
B(1,2,2)=-1.42e-2;
B(2,1,1)=3.88e-1;
B(2,1,2)=2.06e-2;
B(2,2,1)=2.97e-2;
B(2,2,2)=3.26e-3;
B(3,1,1)=-2.46e-2;
B(3,1,2)=-1.31e-3;
B(3,2,1)=-1.48e-3;
B(3,2,2)=-1.8e-4;
end


sigma(i,:)=sigma_Suzuki(E,ne,Te,Zeff,A,B);


end


figure;loglog(reshape(paper{1},[16 3]),reshape(paper{2},[16 3])*1e-16,'-*');hold;
loglog(Ea,sigma,'+-')

%% test of NUR vs BBNBI
%% read tokamak profiles
%load('/compass/Shared/Common/Theory-Modelling/NBI_FIESTA/EBdyna/data_tokamak/pressure_profile.mat')
%load('/compass/Shared/Common/Theory-Modelling/NBI_FIESTA/EBdyna/data_tokamak/flux_geometry.mat')
%load('/compass/Shared/Common/Theory-Modelling/NBI_FIESTA/EBdyna/data_tokamak/tokamak_map_dimensions.mat')

load('/compass/Shared/Common/COMPASS-UPGRADE/RP1_Design/Scenarios/4.4/EBdyna_minorradius_28cm/data_tokamak/pressure_profile.mat')
load('/compass/Shared/Common/COMPASS-UPGRADE/RP1_Design/Scenarios/4.4/EBdyna_minorradius_28cm/data_tokamak/flux_geometry.mat')
load('/compass/Shared/Common/COMPASS-UPGRADE/RP1_Design/Scenarios/4.4/EBdyna_minorradius_28cm/data_tokamak/tokamak_map_dimensions.mat')

%R0=0.55;%C
R0=0.894;%CU

psi_XZ_map=(psi_XZ_map-psi_scale(1))/range(psi_scale);
% figure;imagesc(Z_scale,R0+X_scale,psi_XZ_map)

nef=fit([psi_norm ]',[ne_prof ]','cubicinterp');
Tef=fit([psi_norm ]',[Te_prof ]','cubicinterp');
Zefff=fit([psi_norm]',[linspace(0,2,numel(ne_prof))]','cubicinterp');

ne2d=nef(psi_XZ_map);
ne2d=reshape(ne2d,size(psi_XZ_map));
ne2d(isnan(ne2d))=0;
ne2d(psi_XZ_map>1)=0;

Te2d=Tef(psi_XZ_map);
Te2d=reshape(Te2d,size(psi_XZ_map));
Te2d(isnan(Te2d))=0;
Te2d(psi_XZ_map>1)=0;

Zeff2d=Zefff(psi_XZ_map);
Zeff2d=reshape(Zeff2d,size(psi_XZ_map));
Zeff2d(isnan(Zeff2d))=0;
Zeff2d(psi_XZ_map>1)=0;

% make structures

Beam.bmwidra=0.1;  %source/(last grid) half width [m]
Beam.foclra=4; %distance from the source to optical focual point [m]
Beam.divra=2*0.00823688265; %source divergence [rad]
Beam.app(1)=1; %distance from the source to the apperture [m]
Beam.app(2)=1;% apperture half width [m]
Beam.R_t=0.65; %Tangency radius [m]
Beam.d_s=4.773;% distance from the source to the tangency point [m]
Beam.Z_s=0.41; %source elevation over midplane [m]
Beam.Z_a=0.322844; %aperture elevation over midplane [m]
Beam.Ab=2.01355; %beam speace atomic number
%Beam.Ab=2;
Beam.E=80e3;
%Beam.E=Beam.Emax; %beam Energy[eV]
Beam.frac_E1=0.58;
Beam.frac_E2=0.26;


param.nmark=1e5; %number of markers
param.nb=10000; %number of markers creted during one itaration (have toke it because of shinethrough)
param.nr=1000; %number of points along one ray to calculate probability function

nea=mean(ne_prof);
Plasma.ne=ne2d';%/nea*3.5e19; %density [m^-3]
Plasma.Te=Te2d'; %[eV]
Plasma.Zeff=ones(size(Te2d')) ;
Plasma.R=X_scale+R0; %[m] 1d array
Plasma.Z=Z_scale ;%[m] 1d array
Plasma.Ap=2; %plasma speace atomic number (int)
Plasma.IMP=1; %impurity index (for now IMP==1 CARBON impurity)
Plasma.BT=1; %Magnetic field

[R_i,Z_i,phi_i,vr,vz,vphi,E_i,Shth,X,Y,Z,sh]=NUR(Beam,Plasma,param);

%read BBNBI data
Rbb=ncread('/compass/Shared/Common/COMPASS-UPGRADE/RP1_Design/Scenarios/4.4/ASCOT/for_EBDYNA/for_EBDYNA/BBNBI_E80keV_P1MW_rtang0.65_torangle0.00_CW_incl5.00_z_0.41_rad_4.80.nc','/r-GC');
Zbb=ncread('/compass/Shared/Common/COMPASS-UPGRADE/RP1_Design/Scenarios/4.4/ASCOT/for_EBDYNA/for_EBDYNA/BBNBI_E80keV_P1MW_rtang0.65_torangle0.00_CW_incl5.00_z_0.41_rad_4.80.nc','/z-GC');
Ebb=ncread('/compass/Shared/Common/COMPASS-UPGRADE/RP1_Design/Scenarios/4.4/ASCOT/for_EBDYNA/for_EBDYNA/BBNBI_E80keV_P1MW_rtang0.65_torangle0.00_CW_incl5.00_z_0.41_rad_4.80.nc','/energy');

[Hbb, Rb_bb, Zb_bb]=hist2d(Rbb(find(Ebb==80e3)),Zbb(find(Ebb==80e3)),50);
[H, Rb, Zb]=hist2d(R_i(find(E_i==80e3)),Z_i(find(E_i==80e3)),50);

figure;
subplot(2,1,1);
surf(Rb_bb,Zb_bb,Hbb/sum(Hbb(:))/(range(Rb_bb)*range(Zb_bb)))
view(2);
shading interp
xlim([0.6 1.17]);ylim([-0.3 0.3]);
lim = caxis
title ('bbNBI')

subplot(2,1,2);
surf(Rb,Zb,H/sum(H(:))/(range(Rb)*range(Zb)))
view(2);
xlim([0.6 1.17]);ylim([-0.3 0.3]);
shading interp
caxis(lim)
title ('NUR')

[Hbb, Rb_bb, Zb_bb]=hist2d(Rbb(find(Ebb==80e3/2)),Zbb(find(Ebb==80e3/2)),50);
[H, Rb, Zb]=hist2d(R_i(find(E_i==80e3/2)),Z_i(find(E_i==80e3/2)),50);

figure;
subplot(2,1,1);
surf(Rb_bb,Zb_bb,Hbb/sum(Hbb(:))/(range(Rb_bb)*range(Zb_bb)))
view(2);
shading interp
xlim([0.6 1.17]);ylim([-0.3 0.3]);
lim = caxis
title ('bbNBI')

subplot(2,1,2);
surf(Rb,Zb,H/sum(H(:))/(range(Rb)*range(Zb)))
view(2);
xlim([0.6 1.17]);ylim([-0.3 0.3]);
shading interp
caxis(lim)
title ('NUR')


[Hbb, Rb_bb, Zb_bb]=hist2d(Rbb(find(Ebb==min(Ebb(:)))),Zbb(find(Ebb==min(Ebb(:)))),50);
[H, Rb, Zb]=hist2d(R_i(find(E_i==80e3/3)),Z_i(find(E_i==80e3/3)),50);

figure;
subplot(2,1,1);
surf(Rb_bb,Zb_bb,Hbb/sum(Hbb(:))/(range(Rb_bb)*range(Zb_bb)))
view(2);
shading interp
xlim([0.6 1.17]);ylim([-0.3 0.3]);
lim = caxis
title ('bbNBI')

subplot(2,1,2);
surf(Rb,Zb,H/sum(H(:))/(range(Rb)*range(Zb)))
view(2);
xlim([0.6 1.17]);ylim([-0.3 0.3]);
shading interp
caxis(lim)
title ('NUR')
