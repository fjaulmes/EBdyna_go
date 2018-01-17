function Fourier_3D_field
%Fourier_3D_field See penetration of mode in psi-m plane with m poloidal wave number
%   Detailed explanation goes here

global par maps dim
%% CLEAN
delete(findall(0,'type','figure','tag',mfilename));
unit_test=false;
%% UNIT TEST
if unit_test
    warning('Executing a UNIT TEST!') %#ok<*UNRCH>
    
    L=2^11;
    t = linspace(0,1,L);
    
    % Modes
    x=@(m) cos(2*pi*m*t);         % Basic mode
    X = [x(0);x(1);x(2);x(3);x(4);x(5);x(6);x(7);x(8);x(9);x(10);x(11);x(12);x(13);x(14);x(15)];
    
    % Make 100 rows with each a superposition of a discrete # modes
    row_nr=3;
    rng(1);
    nr_modes=1+round(rand(row_nr,1)*14);
    which_modes=round(rand(row_nr,15)*15);
    
    B_rad=zeros(row_nr,size(X,2));
    for i=1:row_nr
        B_rad(i,:)=sum(X(1+which_modes(i,1:nr_modes(i)),:),1);
    end
    B_rad=B_rad';  % Put in order theta , psi
    psi_star=linspace(0,1,size(B_rad,2)); % make psi-vector23
else
    %% Load field
    % 3D field
    flip_theta=false;       % Change direction of poloidal coordinate
    if true || isempty(maps) || isempty(dim)
        par.APPLY_RMP     =true;                  % Resonant Magnetic Perturbations (RMP)
        par.RMP_file        =['../data_tokamak/','RMP_n=2_-90_flux_2017-01-05.mat'];
        par.APPLY_TFR           =false;                 % Toroidal Field Ripple
        par.APPLY_3D    = par.APPLY_RMP | par.APPLY_TFR;
        par.APPLY_SAWTOOTH=false; 
        par.superimpose_2D_3D=false;
        
        par.CALCULATE_PPHI_3D   =false...               % Determine pphi based on addition of (local) evolution
            & (par.APPLY_RMP | par.APPLY_TFR);
        par.mode=3;
        par.interp_scheme=3;
        par.paths=initialize_folder_names_struct;
        par.coord_syst='flux';
        
        par.step.R=1; par.step.Z=1;  par.step.phi=1; % Parameters for coarser grid
        
        [maps,dim]=initialize_maps(true,false);
        maps2=load([par.paths.DATA_FOLDER,'dl_values.mat'],'dl_X','dl_Z');
        maps=combine_structs(maps,maps2); clear maps2
        
        % Increase size in toroidal direction if symmetry is used in storage
        if dim.n3D.symm~=1
            maps.B_3D=cat(3,repmat(maps.B_3D(:,:,1:end-1,:),[1 1 dim.n3D.symm 1]),maps.B_3D(:,:,1,:));
            dim.n3D.size_3D(3)=size(maps.B_3D,3);
        end
        % Optionally try changing theta or phi definition        
        if flip_theta;
            maps.B_3D=maps.B_3D(end:-1:1,:,:,:); % theta direction
        end
       % maps.B_3D=maps.B_3D(:,:,end:-1:1,:); % phi direction
    else
        warning('Re-using already present 3D map')
    end
    
    %% Get radial B field
    psi_star=(1-dim.psi_scale/dim.psi_scale(1));
    psi_star=interp1(1:length(dim.psi_scale),psi_star,1:dim.n3D.size_3D(2));
    
    %% Configure B-maps
    size_2D=size(maps.dl_Z); %theta, psi
    
    if flip_theta
        theta_ind=fliplr(linspace(1,size_2D(1),dim.n3D.size_3D(1)));    % Points where 3D-field is known (theta)
    else
        theta_ind=linspace(1,size_2D(1),dim.n3D.size_3D(1));
    end
    psi_ind=linspace(1,size_2D(2),dim.n3D.size_3D(2));      % Points where 3D-field is known (psi)
    
    dr_rad_x=-maps.dl_Z(theta_ind,psi_ind);                 % x-component e_theta
    dr_rad_z= maps.dl_X(theta_ind,psi_ind);                 % z-component e_theta
    n_dr=sqrt(dr_rad_x.^2+dr_rad_z.^2);
    
    dr_rad_x=dr_rad_x./n_dr;
    dr_rad_z=dr_rad_z./n_dr;
    
    dr_rad_x(isnan(dr_rad_x))=0;
    dr_rad_z(isnan(dr_rad_z))=0;
    
    %% Determine radial component
    B_rad=bsxfun(@plus,bsxfun(@times,maps.B_3D(:,:,:,1),dr_rad_x),bsxfun(@times,maps.B_3D(:,:,:,2),dr_rad_z));
    
    B_rad=B_rad(1:1:end-1,:,1:end-1); % Shrink to remove overlap of last element (i.e. domain t=linspace(0,1,513); t(end)=[];)
    
end
%% Determine FFT
if size(B_rad,3)>1
    % As example say they're 9 original phi-positions, 8 unique.
    % The Fourier transform will return the modes n=[0 1 2 3 4 -3 -2 -1]
    % Since m-n-phase space is point symetric, we remove negative n, 
    % and compensate this by doubling the amplitude for [1 2 3]
    
    N=2^nextpow2(size(B_rad,3));      % Number of points in array [8]
    max_n=N/2;                        % (max mode -> Nyquist -> N/2) [4]
    B_ft_n  = (fft(B_rad ,N,3))/(N);  % FFT in toroidal direction (normalized)
    
    B_ft_n  = B_ft_n(:,:,1:max_n+1);                % Remove negative n-phases [-3 -2 -1]
    B_ft_n(:,:,2:end-1) = 2*B_ft_n(:,:,2:end-1);	% Double to compensate for removing negative mode numbers [1 2 3]
else
    B_ft_n  = B_rad;
    warning('Only one toroidal mode. Not FFT in phi!')
end

M = 2^nextpow2(size(B_ft_n,1));     % Number of points in array (poloidal)
[m,B_rad_fourier]=fft_mode_with_m(B_ft_n,M,1);

% Note B_rad_fourier is not a point-symmetric matrix in n and m (still depending
% on the radial coordinate). 

%% Radial coordinate
% Choose q-surface
q=4;

% Find correct radial index
if ~isfield(dim,'q_initial_profile')
    dim2=load([par.paths.DATA_FOLDER,'q_profile.mat'],'q_initial_profile'); dim=combine_structs(dim,dim2); clear dim2
end
psi_norm=interp1(dim.q_initial_profile,1:dim.NB_PSI,q);     % Radial index in FINESSE maps
psi_norm=round(dim.n3D.ind_2D_to_3D(psi_norm));             % Radial index in 3D maps

%% PLOTTING
if unit_test
    n_modes=0;
else
    n_modes=[2];
end

for n_mode=n_modes
    %% Plot FFT
    hf=figure('name',['n=',num2str(n_mode)],'tag','Fourier components radial perturbation');
    ha1=subplot(1,1,1,'parent',hf);
    xlabel(ha1,'$m$','FontSize',26)
    ylabel(ha1,'$\psi^*$','FontSize',26)
    set(ha1,'FontSize',22)
    hold(ha1,'on')
    set(ha1,'XLim',[-15 15]);
    set(ha1,'YLim',[0 1]);
    set(ha1,'YDir','normal')
    
    [~,h]=contourf(m,psi_star,B_rad_fourier(:,:,n_mode+1)',100,'parent',ha1); set(h,'LineStyle','none');
    %     h=imagesc(m,psi_star,B_rad_fourier(:,:,n_mode+1)','parent',ha1);
    title(['$B_r$ at $n=$',num2str(n_mode)],'FontSize',26)
    
    cl=colorbar('peer',ha1);
    title(cl,'$\left|\right.b_{nm}\left.\right|$','interpreter','latex','FontSize',26)
    %     colormap winter
    colormap jet
    try
        
        set(ha1,'ticklabelinterpreter','latex');
        set(cl,'ticklabelinterpreter','latex');
    catch
    end
    
    % ADD Q-profile dots
    q_a=max(dim.q_initial_profile);
    m_max=floor(n_mode*q_a);
    psi_norm_rational_q=interp1(dim.q_initial_profile,1:length(dim.q_initial_profile),(1:m_max)/n_mode);
    psi_star_rational_q=(psi_norm_rational_q-1)/(length(dim.q_initial_profile)-1);
    
    plot(ha1,(1:m_max)/n_mode,psi_star_rational_q,'w.','MarkerSize',50)
    
    %% Plot poloidal mode number spectra at q=4
    hf=figure('name',['n=',num2str(n_mode)],'tag',mfilename);
    ha1=subplot(1,1,1,'parent',hf);
    xlabel(ha1,'$m$','FontSize',26)
    ylabel(ha1,'$B_r$ [T]','FontSize',26)
    set(ha1,'FontSize',22)
    hold(ha1,'on')
    set(ha1,'XLim',[-25 25]);
    
    x_ax=m(abs(m)<=25);
    y_ax=B_rad_fourier(:,psi_norm,n_mode+1);
    y_ax=y_ax(abs(m)<=25);
    
    % Resonance need to be red
    res=q*n_mode;
    expr=   x_ax==res | x_ax==-res;
    
    bar(x_ax,y_ax,'parent',ha1,'FaceColor','k');
    y_ax(~expr)=0;
    bar(x_ax,y_ax,'parent',ha1,'FaceColor','r');
    title(['$B_r$ at $n=$',num2str(n_mode),' at $q$=',num2str(q)],'FontSize',26)
    try
        set(ha1,'ticklabelinterpreter','latex');
        set(cl,'ticklabelinterpreter','latex');
    catch
    end
end
%% Plot B_rad
% Make figure
hf=figure('name','B_rad','tag',mfilename);
ha1=subplot(1,1,1,'parent',hf);
xlabel(ha1,'$\varphi$ [rad]','FontSize',26)
ylabel(ha1,'$\theta$ [rad]','FontSize',26)
set(ha1,'FontSize',22)
hold(ha1,'on')
set(ha1,'XLim',[0 2*pi]); set(ha1,'YLim',[-pi pi]);

% Plot
phi=linspace(0,2*pi,dim.n3D.size_3D(3));
theta=linspace(-pi,pi,dim.n3D.size_3D(1));


B_image=squeeze(B_rad(:,psi_norm,:));
mirror_B_image=1+(size(B_image,1))*0.5;
B_image=cat(1,B_image(mirror_B_image:end,:),B_image(1:(mirror_B_image),:));
B_image(:,end+1)=B_image(:,1);

%imagesc(phi,theta,B_image,'parent',ha1);
contourf(phi,theta,B_image,100,'parent',ha1,'LineStyle','none');
title(['$B_r$ at $q=$',num2str(q)],'FontSize',26)
cl=colorbar('peer',ha1);
colormap(ha1,'jet')
title(cl,'[T]','interpreter','latex','FontSize',26)
set(ha1,'YDir','reverse')
try
    set(ha1,'ticklabelinterpreter','latex');
    set(cl,'ticklabelinterpreter','latex');
catch
end
drawnow;
end

function [n,P1]=fft_mode_with_m(X,N,dim)
% EXAMPLE:  We have N=8 -> N/2 =4 as the maximum mode.
%           The modes are in the order: [0 1 2 3 4 -3 -2 -1];
%           We sort them to [-4 -3 -2 -1 0 1 2 3 4], so 0 is in the center
%           Note that the amplitude of mode 4 is divided by 2 since it is
%           twice in the output

if mod(log2(N),1)~=0
    error('Requested Fourier transform with something different than a power of 2')
end
    
%% FFT
N_normal=2^nextpow2(size(X,dim));
Y = fft(X,N,dim);           % FFT, with frequency in rad/s
n= -N/2:N/2 ;               % Negative and positive mode numbers (e.g. -4 to 4)
n=n * N_normal/N;           % Normalize the mode-axis to the correct domain (e.g. if N=16 was requested, then n=-4:0.5:4 instead n=-4:4)

%% Plot quantitues
P2=abs(Y)/(N);

ind_N_2=N/2+1;              % Position of half of the domain (doubled)  (e.g. 4)
ind_pos=2:ind_N_2-1;        % Positions of positive frequencies         (e.g. 1 to 3)
ind_neg=ind_N_2+1:N;        % Positions of negative frequencies         (e.g. 5 to 7 -> -3 to -1 )

% e.g. [ -4                 -3:-1           0           1:3             4];
switch dim
    case 1
        P1=[0.5*	P2(ind_N_2  ,:,:) ;...      % mode -N/2
                    P2(ind_neg  ,:,:) ;...      % mode 0<-n<N/2
                    P2(1        ,:,:) ;...      % mode 0
                    P2(ind_pos  ,:,:) ;...      % mode 0<n<N/2
           0.5*     P2(ind_N_2  ,:,:)];         % mode N/2
    case 2
        P1=[0.5*	P2(:,ind_N_2    ,:) ;...
                    P2(:,ind_neg    ,:) ;...
                    P2(:,1          ,:) ;...
                    P2(:,ind_pos    ,:) ;...
           0.5*     P2(:,ind_N_2    ,:)];
    case 3
        P1=cat(3,...
            0.5*    P2(:,:,ind_N_2    ) ,...
                    P2(:,:,ind_neg 	) ,...
                    P2(:,:,1      	) ,...
                    P2(:,:,ind_pos  ) ,...
            0.5*    P2(:,:,ind_N_2    ) );
    otherwise
        error('Dimension not programmed');
end

end

%% Script to `play' with FFT in ML
% Definition of modes
% M = 64;     % 64 sample points
% t=linspace(0,1,M+1); t(end)=[];   % Time sampling (not the double points
% at t=0 and t=1 (phi=0 and phi=2*pi)
% modes = linspace(0, M, 8001); modes(end)=[];  % For plotting of ideal FFT
% 
% mode_in=[10.5 15];              % Mode 10.5 and 15
% % Position of mode in modes-vector
% [~,ind(1)]=min(abs(modes-mode_in(1)));  
% [~,ind(2)]=min(abs(modes-mode_in(2)));
% % Set the (ideal) fft to 1 at the correct frequencies
% X_dtft = zeros(size(modes));                  
% X_dtft(ind)=1;
% plot(modes, abs(X_dtft),'displayname','ideal fft / active modes')
% xlabel('mode number'); ylabel('A');
% xlim([0 M])
% 
% xf=@(mode) cos(2*pi*mode*t);        % Just the regular expression for harmonics
% x=xf(mode_in(1))+xf(mode_in(2));  % Superposition of the various modes
% 
% % Fourier transform with equal of sample points
% X = fft(x);     
% m=0:(M-1) ;     % The integer number of modes
% X=X./(M);       % FFT normalized to the number of modes
% 
% % Fourier transform with higher order than sampling (aliasing)
% P = 64*M;       
% Xp = fft(x, P); 
% mp=(0:(P-1)) *(M / P);  % the horizontal axis, 1:M-1
% Xp=Xp./(M);
% 
% hold on
% plot(m, abs(X), 'r','displayname','straight forward FFT')
% plot(mp, abs(Xp), 'k','displayname','FFT of higher order')
% hold off