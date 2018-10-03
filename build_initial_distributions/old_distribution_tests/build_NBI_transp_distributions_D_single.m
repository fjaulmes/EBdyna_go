clear all;
initialize_folder_names;
filename=strcat(DATA_FOLDER,'q_profile.mat');
load(filename);
filename=strcat(DATA_FOLDER,'physics_constants.mat');
load(filename);
filename=strcat(DATA_FOLDER,'pressure_profile.mat');
load(filename);
filename=strcat(DATA_FOLDER,'motions_map_dimensions.mat');
load(filename);
filename=strcat(DATA_FOLDER,'flux_geometry.mat');
load(filename);
filename=strcat(DATA_FOLDER,'XZsmall_fields_tokamak_pre_collapse.mat');
load(filename);

% from eqdsk file
Raxis=R0+X_axis
rgeo=1.625;
zgeo=-0.006;
% from Benedikt
rmaxis=1.67;
zmaxis=0.0763;

rgeo_recalc=rmaxis-Raxis+R0;

% a compromise
rgeo=0.5*(rgeo+rgeo_recalc)

read_transp_data_file;

%%

%pause
close all

% Main important values for description of the distribution
mHe=mD
ZHe=1

NBI_RESCALE=1.0
DIM_RESCALE=1.015


% alphas_pos_x=zeros(Nalphas_simulated,1);
% alphas_pos_z=zeros(Nalphas_simulated,1);
% alphas_pos_phi=zeros(Nalphas_simulated,1);
% alphas_Ekin=zeros(Nalphas_simulated,1);
% alphas_mm=zeros(Nalphas_simulated,1);
% alphas_vpll=zeros(Nalphas_simulated,1);


% double data statistic for toroidla representation
NB_PROCESS=16
NB_SPLITS=round(1);

particles_weight=2.69992e+13
particles_weight=particles_weight/NB_SPLITS
PROCESS_NUMBER=1
% for PROCESS_NUMBER=1:NB_PROCESS
%     
    INIT_DATA_POS=mod(PROCESS_NUMBER-1,NB_SPLITS)+1
    SPLIT_RATIO=NB_SPLITS
    
    if mod(PROCESS_NUMBER-1,NB_SPLITS)==0
        disp('------- getting back to initial data -------')
        DIM_RESCALE=1.015;
    end
    
    d_pos_x=data(INIT_DATA_POS:SPLIT_RATIO:end,1)*0.01-rgeo;
    d_pos_z=data(INIT_DATA_POS:SPLIT_RATIO:end,2)*0.01-zmaxis;
    d_Ekin=NBI_RESCALE*data(INIT_DATA_POS:SPLIT_RATIO:end,4)*1e3;
    d_pitch=data(INIT_DATA_POS:SPLIT_RATIO:end,5);
    d_vtot=sqrt(2*(eV/mHe)*d_Ekin);
    % d_vpll=sign(d_pitch).*sqrt(2*(eV/mHe)*d_Ekin.*(d_pitch.^2)./(1+d_pitch.^2));
    d_vpll=-d_pitch.*d_vtot;
    
    alphas_pos_x=NBI_RESCALE*[d_pos_x ];
    alphas_pos_z=NBI_RESCALE*[d_pos_z ]+Z_axis;
    alphas_Ekin=[d_Ekin];
    alphas_vpll=[d_vpll];
    
    Epll=0.5*(mHe/eV)*alphas_vpll.^2;
    Eperp=max(alphas_Ekin-0.5*(mHe/eV)*alphas_vpll.^2,0);
    alphas_vperp=sqrt(2*(eV/mHe)*Eperp);
    alphas_vtot=sqrt(2*(eV/mHe)*alphas_Ekin);
    
    density_part_ratio=1
    Nalphas_simulated=length(alphas_Ekin)

    
    
    disp('number of particles generated');
    Nalphas_simulated
    
    alphas_pos_phi=zeros(Nalphas_simulated,1);
    
    
    % display particles
    figure(1);
    set(gca,'FontSize',22);
    
    hold on;
    axis xy square
    
    
    hold on; grid on
    contour(scale_X+R0,scale_Z,psi_XZsmall_map',psi_scale(2:22:end),'k')
    contour(scale_X+R0,scale_Z,psi_XZsmall_map',psi_scale(psi_rank_q1),'r','linewidth',4)
    contour(scale_X+R0,scale_Z,psi_XZsmall_map',[psi_scale(end) psi_scale(end)] ,'k','linewidth',4)
    
    plot(alphas_pos_x+R0,alphas_pos_z,'b.');
    
    pause(0.1);
    
    alphas_psi_value=interp2(scale_X,scale_Z,psi_XZsmall_map',alphas_pos_x,alphas_pos_z,'*linear');
    alphas_psi=interp1(psi_scale,1:257,alphas_psi_value);
    
    
    % filling up the uniformly distributes toroidal position values
    for(n=1:Nalphas_simulated)
        alphas_pos_phi(n)=my_rand(1)*2*pi;
    end
    
    alphas_mm=(alphas_Ekin-0.5*(mHe*alphas_vpll.^2)/eV)./interp2(scale_X,scale_Z,Btot_XZ_map',alphas_pos_x,alphas_pos_z,'*linear');
    
    NEGATIVE_MM=find(alphas_mm<0);
    
    disp('Beware of negative mm!!!!!!')
    disp('length(NEGATIVE_MM)');
    disp(length(NEGATIVE_MM));
    alphas_mm=max(alphas_mm,0);
    
    
%     FILENAME=strcat('initial_NBI60keV_transp_D_distribution',num2str(PROCESS_NUMBER),'.mat')
%     save (FILENAME,'alphas_pos_x','alphas_pos_z','alphas_pos_phi','alphas_Ekin','alphas_mm','alphas_vpll','Nalphas_simulated','particles_weight','mHe','ZHe');
%     DIM_RESCALE=DIM_RESCALE-0.01;
%     
% end

figure(2)
PART_POP=find((d_pitch>0.4).*(d_pitch<0.6));
hist(alphas_Ekin(PART_POP),20)


%%
PART_POP=find((alphas_pos_z>-0.1).*(alphas_pos_z<0.1).*(alphas_pos_x>X_axis-0.1).*(alphas_pos_x<X_axis+0.1));
EKIN_BIN_SIZE=8*1e3;
EKIN_BINS=(0:EKIN_BIN_SIZE:80*1e3);
Ekin_values=EKIN_BINS(1:end-1)+0.5*EKIN_BIN_SIZE;

PICH_BIN_SIZE=0.1;
PITCH_BINS=(-1.0:PICH_BIN_SIZE:1.0);
pitch_values=PITCH_BINS(1:end-1)+0.5*PICH_BIN_SIZE;
mHist = hist2d ([alphas_Ekin(PART_POP) d_pitch(PART_POP)], EKIN_BINS, PITCH_BINS);
contourf ( Ekin_values,pitch_values, mHist'/max(max(mHist)),(0 :0.05:1));
xlabel('Ekin')
ylabel('pitch')

%%
figure(3)
PART_POP=find((alphas_pos_z>-0.2).*(alphas_pos_z<0.2));
hist(alphas_pos_x(PART_POP)+R0,60)