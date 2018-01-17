function TF_to_TFR(TF,coord_syst)
%TF_to_TFR Determines from the total field (with arbitrary amplitude) the
%ripple
%   Detailed explanation goes here
global m glob

switch nargin
    case 0
        coord_syst='toroidal';
        warning(['Assuming a ',coord_syst,' coordinate system!'])
        TF=load(['./output/TF_',coord_syst,'_',datestr(now,'yyyy-mm-dd'),'.mat']);
    case 1
        coord_syst='toroidal';
        warning(['Assuming a ',coord_syst,' coordinate system!'])
    case 2
    otherwise
        error('Input arguments invalid')
end

%% LOAD FIELD / VARIABLES / FIELDNAMES
if isempty(TF)
    TF=load(['./output/TF_',coord_syst,'_',datestr(now,'yyyy-mm-dd'),'.mat']);
end
BA_fields=fieldnames(TF)';
BA_fields(strcmp(BA_fields,'symm'))=[]; % Remove the symm from the BA_fields

%% Load 2D field
par.paths=initialize_folder_names_struct;
filename=[par.paths.DATA_FOLDER,'tokamak_PR_map.mat'];
m=load(filename,'BX_PR_map','BZ_PR_map');
filename=[par.paths.DATA_FOLDER,'B_fields.mat'];
m2=load(filename,'Btor_PR_map'); m=combine_structs(m,m2); clear m2

filename=[par.paths.DATA_FOLDER,'flux_geometry.mat'];
glob=load(filename,'Z_PR_map','X_PR_map');
filename=[par.paths.DATA_FOLDER,'motions_map_dimensions.mat'];
d2=load(filename,'scale_X','scale_Z','Z_axis','R0','Raxis','a'); glob=combine_structs(glob,d2); clear d2

%% 1. Average |B| of 3D field (find 1/R from vacuum field)
for BA=BA_fields
	if any(strcmp(BA{1},{'dAR_dphi','dAZ_dphi','dAphi_dphi'})); continue; end;  % One doesn't need a 2D-variant of phi-derivatives (=0), but it needs to be scaled
	TF_2D.(BA{1})=mean(TF.(BA{1})(:,:,2:end),3);
end

%% 2. Scale TF with |B| on edge (no plasma response)
if strcmp(coord_syst,'flux')
	% The values for scaling
	mag_TF_on_edge=TF_2D.Bphi(:,end);   % Neglect flucuations of BR and BZ
	mag_B_on_edge=sqrt(m.BX_PR_map(:,end).^2+m.BZ_PR_map(:,end).^2+m.Btor_PR_map(:,end).^2);
	
	% The scaling factor
	sca_fac=mag_B_on_edge./mag_TF_on_edge;
	disp(['relative std. deviation scaling factor: ',num2str(std(sca_fac)/mean(sca_fac))])
	if ~exist('./output_TF','file'); mkdir('./output_TF'); end
	save(['output_TF/sca_fac_TF_',datestr(now,'yyyy-mm-dd')],'-v7.3','sca_fac'); 
    sca_fac=mean(sca_fac);
else
    warning('Choose to load scaling factor from a flux calculation of the same equilibrium!')
    sca_fac=load(['output_TF/sca_fac_TF_',datestr(now,'yyyy-mm-dd')],'sca_fac'); sca_fac=sca_fac.sca_fac; 
    disp(['relative std. deviation scaling factor: ',num2str(std(sca_fac)/mean(sca_fac))])
    sca_fac=mean(sca_fac);
end
% The scaled 2D and 3D field
for BA=BA_fields
    TF.   (BA{1})=sca_fac*TF.   (BA{1});
    if any(strcmp(BA{1},{'dAR_dphi','dAZ_dphi','dAphi_dphi'})); continue; end; % One doesn't need a 2D-variant of phi-derivatives (=0), but it needs to be scaled
	TF_2D.(BA{1})=sca_fac*TF_2D.(BA{1});
end

%% 3. Remove the 2D field from the 3D to have only the ripple
TFR=TF;
for BA=BA_fields
    if any(strcmp(BA{1},{'dAR_dphi','dAZ_dphi','dAphi_dphi'})); continue; end; % One doesn't need a 2D-variant of phi-derivatives (=0), but it needs to be scaled
	TFR.(BA{1})=bsxfun(@minus,TF.(BA{1}),TF_2D.(BA{1}));
end

%% Save output
save(['./output/TFR_',coord_syst,'_',datestr(now,'yyyy-mm-dd'),'.mat'],'-v7.3','-struct','TFR');
 
%% Plot the resulting ripple
if false
    if isfield(TF,'symm')
        % Increase TF size for plotting
        for BA=BA_fields
            TF.(BA{1})=cat(3,repmat(TF.(BA{1})(:,:,1:end-1),[1 1 TF.symm]),TF.(BA{1})(:,:,1));
        end
    end
    
    plot_TFR(TF,TF_2D,coord_syst)
end
clearvars -global
end

function plot_TFR(TF_3D,TF_2D,coord_syst)
global glob
 %% Determine ripple and other quantities
    % B (and expected B0 relation)
    B_2D=sqrt(TF_2D.BR.^2+TF_2D.BZ.^2+TF_2D.Bphi.^2);
    B_3D=sqrt(TF_3D.BR.^2+TF_3D.BZ.^2+TF_3D.Bphi.^2);
    Bmax=max(B_3D,[],3); Bmin=min(B_3D,[],3);
    ripple=bsxfun(@times,bsxfun(@minus,B_3D,0.5*(Bmax+Bmin)),(0.5*(Bmax+Bmin)).^-1);
switch coord_syst
    case 'flux'
        B0=mean(B_2D(:,1));
        R=linspace(glob.R0-glob.a,glob.R0+glob.a,1e3);
        B_mag_check= B0*glob.Raxis./R;
    
        % Find theta which on HFS mid-plane
        psi_norms=2:size(B_2D,2);   % Skip the center, i.e. skip psi=psi_0
        thetas=zeros(size(psi_norms));
        theta_range=[pi/2 3*pi/2];
        
        [~,b_ind]=min(abs(linspace(0,2*pi,size(TF_2D.BR,1))-theta_range(1)));
        [~,e_ind]=min(abs(linspace(0,2*pi,size(TF_2D.BR,1))-theta_range(2)));
        for psi_norm=psi_norms
            [~,thetas(psi_norm-1)]=min(abs(glob.Z_PR_map(b_ind:e_ind,psi_norm)));
        end
        % Flip, so go from HFS to mid. and correct for boundary of thetas
        thetas=fliplr(thetas)+b_ind-1;
        
        % Extend to LFS
        thetas=cat(2,thetas,ones(1,size(B_2D,2)));
        psi_norms=[size(B_2D,2):-1:2 1:size(B_2D,2)];
        
        % For the theta=pi side
        ind_2D=sub2ind(size(TF_2D.BR),thetas,psi_norms);
        B_midplane_2D=B_2D(ind_2D);
        R_midplane_2D=glob.X_PR_map(ind_2D)+glob.R0;
        
        % For all other theta's
        R_2D=glob.X_PR_map(1:4:end,1:4:end)+glob.R0;
        B_2D_vec=B_2D(1:4:end,1:4:end);
        B_and_R_2D=sortrows([R_2D(:) B_2D_vec(:)],1);
    %% Figures
    delete(findall(0,'type','figure','tag',mfilename))
    
    % Plot the 1/R relation
    figure('tag',mfilename)
    ha=axes; hold (ha,'on')
    set(ha,'FontSize',38); %,'ticklabelinterpreter','latex')
    xlabel(ha,'$R$ [m]')
    ylabel(ha,'$\left|B\right|/\left|B_0\right|$')
    
    plot(ha,R,B_mag_check/B0,'g--','displayname','$\left|B\right|=\left|B_0\right|R_0/R$');
    plot(ha,R_midplane_2D,B_midplane_2D/B0,'r','displayname','$\left|B\right|$ midplane');
    plot(ha,B_and_R_2D(:,1),B_and_R_2D(:,2)/B0,'color',[0.75 0.75 0.75],'displayname','$\left|B\right|$ full 2D');
    
    % Relative error according to the 1/R relation
    figure('tag',mfilename)
    ha=axes; hold (ha,'on')
    set(ha,'FontSize',38); %,'ticklabelinterpreter','latex')
    xlabel(ha,'$R$ [m]')
    ylabel(ha,'$\Delta \left|\vec{B}\right|/\left|\vec{B}_{\propto R{-1}}\right|$')
    plot(ha,B_and_R_2D(:,1),B_and_R_2D(:,2)./(B0*glob.Raxis./B_and_R_2D(:,1))-1,'b.','linestyle','none','displayname','relative error over whole 2D domain');
    
    % Top view with mid-plane
    figure('tag',mfilename)
    ha=axes; hold (ha,'on')
    set(ha,'FontSize',38); %,'ticklabelinterpreter','latex')
    xlabel(ha,'$x$ [m]')
    ylabel(ha,'$y$ [m]')
    zlabel(ha,'$z$ [m]')
    
    phis=1:size(ripple,3);    % For all phi-positions
    % The coordinates in 2x2 (each poloidal point);
    [thetas_3D,phis_3D]     =ndgrid(thetas,phis);
    [psi_norms_3D,~]        =ndgrid(psi_norms  ,phis);
    
    ind_midplane=sub2ind(size(ripple),thetas_3D,psi_norms_3D,phis_3D);
    
    ripple_midplane_3D=ripple(ind_midplane);
    ripple_midplane_3D=reshape(ripple_midplane_3D,[size(thetas_3D,1),size(thetas_3D,2)]);
    
    R=repmat(R_midplane_2D',1,length(phis));
	contourf(R,phis_3D,ripple_midplane_3D,100,'linestyle','none','parent',ha);
    h=get(ha,'children');
    
	phis_3D_actual=zeros(size(phis_3D));
	phis_3D_actual(:)=interp1(1:size(B_3D,3),linspace(0,2*pi,size(B_3D,3)),phis_3D(:));
	
    x_midplane= R.*cos(phis_3D_actual);
    y_midplane=-R.*sin(phis_3D_actual);
	xy_lim=max(abs([x_midplane(:);y_midplane(:)]));
    set(h,'XData',x_midplane,'YData',y_midplane);
    set(ha,'xlim',[-xy_lim xy_lim],'ylim',[-xy_lim xy_lim])
    axis(ha,'equal');
          
    % The ripple in 2D with check
    figure('tag',mfilename)
    ha=axes; hold(ha,'on')
    set(gca,'FontSize',38)
    xlabel('$R$ [m]')
    ylabel('$\eta_{max}$ [\%]')
    
    if ~exist('./input/ripple_1D_comparison_digitalized.csv','file')
        error('Check with vacuum ripple field not found')
    end
    
    check_data=csvread('./input/ripple_1D_comparison_digitalized.csv');
    ripple_2D=max(ripple,[],3);
    ripple_2D=ripple_2D(ind_2D)*100;
    plot(ha,check_data(:,1),check_data(:,2),'r+-')
    plot(ha,R_midplane_2D,ripple_2D,'bx')
    
    ylim(ha,[0 max(ripple_2D)*1.1])
   case 'toroidal'
        %% Figures
    delete(findall(0,'type','figure','tag',mfilename))
    
    % Plot the 1/R relation
    figure('tag',mfilename)
    ha=axes; hold (ha,'on')
    set(ha,'FontSize',38); %,'ticklabelinterpreter','latex')
    xlabel(ha,'$R$ [m]')
    ylabel(ha,'$\left|B\right|/\left|B_0\right|$')
    
    scale_R=glob.scale_X+glob.R0;
    ind_Xaxis=interp1(scale_R,1:length(scale_R),glob.Raxis);
    ind_Zaxis=interp1(glob.scale_Z,1:length(glob.scale_Z),glob.Z_axis);
    
    B0=interp2(B_2D,ind_Zaxis,ind_Xaxis);
    R=linspace(glob.R0-glob.a,glob.R0+glob.a,1e3);
    B_mag_check= B0*glob.Raxis./R;
    
    B_and_R_2D=sortrows([repmat(scale_R(:),[length(glob.scale_Z),1]) B_2D(:)],1);
    
	B_midplane_2D=interp2(B_2D,ind_Zaxis,1:length(scale_R));
    plot(ha,R,B_mag_check/B0,'g--','displayname','$\left|B\right|=\left|B_0\right|R_0/R$');
    plot(ha,scale_R,B_midplane_2D/B0,'r','displayname','$\left|B\right|$ midplane');
    plot(ha,B_and_R_2D(:,1),B_and_R_2D(:,2)/B0,'color',[0.75 0.75 0.75],'displayname','$\left|B\right|$ full 2D');
    
    % Relative error according to the 1/R relation
    figure('tag',mfilename)
    ha=axes; hold (ha,'on')
    set(ha,'FontSize',38); %,'ticklabelinterpreter','latex')
    xlabel(ha,'$R$ [m]')
    ylabel(ha,'$\Delta \left|\vec{B}\right|/\left|\vec{B}_{\propto R{-1}}\right|$')
    plot(ha,B_and_R_2D(:,1),B_and_R_2D(:,2)./(B0*glob.Raxis./B_and_R_2D(:,1))-1,'b.','linestyle','none','displayname','relative error over whole 2D domain');
    
    % Top view with mid-plane
    figure('tag',mfilename)
    ha=axes; hold (ha,'on')
    set(ha,'FontSize',38); %,'ticklabelinterpreter','latex')
    xlabel(ha,'$x$ [m]')
    ylabel(ha,'$y$ [m]')
    zlabel(ha,'$z$ [m]')
    
    phis=linspace(0,2*pi,size(ripple,3));    % For all phi-positions
    
    % The coordinates in 2x2 (each poloidal point);
    ripple_midplane_3D=squeeze(interp3(ripple,ind_Zaxis,1:length(scale_R),1:length(phis)));
    x_midplane=bsxfun(@times,scale_R',cos(phis));
    y_midplane=bsxfun(@times,scale_R',-sin(phis));
        
	[~,h]=contourf(x_midplane,y_midplane,ripple_midplane_3D,100,'linestyle','none','parent',ha);
    
    xy_lim=max(abs([x_midplane(:);y_midplane(:)]));
    set(h,'XData',x_midplane,'YData',y_midplane);
    set(ha,'xlim',[-xy_lim xy_lim],'ylim',[-xy_lim xy_lim])
    axis(ha,'equal');
    drawnow; saveas(ha,'TEMP.fig')      
    
    % The ripple in 2D with check
    figure('tag',mfilename)
    ha=axes; hold(ha,'on')
    set(gca,'FontSize',38)
    xlabel('$R$ [m]')
    ylabel('$\eta_{max}$ [\%]')
    
    if ~exist('./input/ripple_1D_comparison_digitalized.csv','file')
        error('Check with vacuum ripple field not found')
    end
     
    check_data=csvread('./input/ripple_1D_comparison_digitalized.csv');
    ripple_2D=max(ripple,[],3);
    ripple_2D=interp2(ripple_2D,ind_Zaxis,1:length(scale_R))*100;
    plot(ha,check_data(:,1),check_data(:,2),'r+-')
    plot(ha,scale_R,ripple_2D,'bx')
    
    ylim(ha,[0 max(ripple_2D)*1.1])
    drawnow;
    saveas(ha,'TEMP2.fig')  
end
end