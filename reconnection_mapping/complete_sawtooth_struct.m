function complete_sawtooth_struct
%#ok<*AGROW>
%#ok<*NASGU>
 
clearvars -global
old_intp=get(0,'defaulttextinterpreter');
set(0,'defaulttextinterpreter','latex');

%% LOAD PARAMETERS
par.paths=initialize_folder_names_struct;
 
% const=load(strcat(par.paths.DATA_FOLDER,'physics_constants.mat'));
maps=load('../data_tokamak/flux_geometry.mat','Z_PR_map', 'X_PR_map');
maps2=load('../data_tokamak/XZsmall_fields_tokamak_pre_collapse.mat','psi_norm_XZsmall_map'); maps=combine_structs(maps,maps2); clear maps2

dim=load(strcat(par.paths.DATA_FOLDER,'psi_profiles.mat'),'psi_pol_initial_profile'); 
dim2=load(strcat(par.paths.DATA_FOLDER,'q_profile.mat'),'q_initial_profile'); dim=combine_structs(dim,dim2); clear dim2
dim2=load('../data_tokamak/motions_map_dimensions.mat','scale_X','scale_Z','R0'); dim=combine_structs(dim,dim2); clear dim2
if length(dim.psi_pol_initial_profile)~=length(dim.q_initial_profile)
    error('Length of psi and q profiles not the same. Please correct input in DATA TOKAMAK folder')
end

 if verLessThan('matlab','8.4')
     intp_method='cubic';
 else
     intp_method='phchip';
 end

%% Time settings
time_steps_reconnection=1e2;    % time steps of reconnection sawtooth   (r1 linearly decreasing)
time_steps_relaxation=1e2;      % time steps during relaxation          (k (deformation) linearly decreasing)

time_reconnection=1e-4;         % Time of the reconnection process (s)
time_relaxation=2e-4;            % Time of relaxation process (s)

% (Complete) time array
time=get_time(time_steps_reconnection,time_steps_relaxation,time_reconnection,time_relaxation,'linear');
time(time_steps_reconnection+1)=time_reconnection-1e-10;

%% SFL-coordinates (initial)
psi=max(dim.psi_pol_initial_profile)-dim.psi_pol_initial_profile;  % Define psi as 0 zero in center of plasma and finite at edge
chi=cumtrapz(psi,dim.q_initial_profile);          % toroidal flux (psi) = int_0si  q d(psi) with trapezium integration

psi_min=psi-chi;                              % for m=n=1, psi_star = psi-phi_H, with helical poloidal flux p_H = n/m chi

%% RADIAL POSITIONS
r=sqrt(2*chi);  % Definition of r-coordinate

% r at q=1 (with top of psi_star)
ind_r0      = interp1(dim.q_initial_profile,1:length(r),1,intp_method);   % Index q=1
ind_r0_int  = floor(ind_r0);                                        % Integer index q=1
h_ind_r0_int= floor(ind_r0/2);                                      % Integer index between center and q=1, for finding rmix
r0          = interp1(dim.q_initial_profile,r,1,intp_method);       % r-value at q=1

% Mixing radius (psi_star = 0). Interpoldation from h_ind_r0_int since we do not want the center point 
ind_rmix = interp1(psi_min(h_ind_r0_int:end), (h_ind_r0_int:length(dim.psi_pol_initial_profile)	),0,intp_method);     % index mixing radius
rmix     = interp1(psi_min(h_ind_r0_int:end),r(h_ind_r0_int:end                      ),0,intp_method);     % r-valua mixing radius

%% RECONNECTION r-parameters (based on interpolation / finding equal flux)
% Assumption: r1 (position from original axis to center `old' core) moves with constant speed 
r1_t=@(t) r0*(1-t/time_reconnection);
r1=r1_t(time(1:time_steps_reconnection+1));
r1      =cat(1,r1   ,zeros(time_steps_relaxation+1,1));

psi_plus=interp1(r(1:ind_r0_int+1),psi_min(1:ind_r0_int+1),r1,intp_method);

r2=interp1(psi_min(ind_r0_int+1:end),r(ind_r0_int+1:end),psi_plus,intp_method);
r2(1)=r1(1);

if any(r2<r1); error('radial positions badly conditioned'); end;
if (any(r2>rmix) || r2(end)~=rmix); error('Check r2 vector for conditioning'); end;

% Derived r-values
r3=sqrt(r2.^2-r1.^2);
eta=r2-r1;

%% FLUX VALUES
% 1-values to determine dpsi / dr (which is dpsi/dchi * dchi/dr, with chi=0.5*r^2)
q1 =    interp1(r,dim.q_initial_profile,r1,intp_method);
q2 =    interp1(r,dim.q_initial_profile,r2,intp_method);

psi_1_acc= q2.*(1-q1) .*r1;  % dpsi/dr *q1*q2
psi_2_acc= q1.*(1-q2) .*r2;  % dpsi/dr *q1*q2

%% RECONNECTION PARAMETERS  (K and KC)
% KC: reconnection parameter equal to deformation with continuous E-field / no signular potential.
% Note this parameter is defined by the commented line, yet calculated by
% the next to remove singularity at r1=0 (end reconnection).
kc2 =2./(r2+r1) .* (psi_1_acc + psi_2_acc)./(psi_2_acc./r2 - psi_1_acc./r1);
kc =2./(r2+r1) .* (psi_1_acc + psi_2_acc)./(q1 - q2);
kc(1)=0;

% CURRENT DIFFUSION WITH K -> 0
% Let k decrease to 0
k_t=@(t) kc(time_steps_reconnection+1)*(1-(t-time_reconnection)./(time_relaxation));
k	=cat(1,kc(1:time_steps_reconnection+1)    ,k_t(time(time_steps_reconnection+2:end)));
k(end)=0;

% Determine kr, c, rho, r+ based on complete vectors
kr=(r2-r1)./(r2+r1);
kr_inv=1./kr;

%% Time derivatime values
r1_dot=zeros(size(r1)); 
r1_dot(1:time_steps_reconnection+1)=-r0/time_reconnection;

r2_dot=r1_dot .* psi_1_acc./psi_2_acc;
r2_dot(1)=-r1_dot(1);

r3_dot=(r2_dot.*r2-r1_dot.*r1)./r3;

kr_dot = 2*(r1.*r2_dot-r2.*r1_dot)./(r2+r1).^2;
kr_dot (time_steps_reconnection+2:end)=0;

k_dot=gradient(k,time);
k_dot(time_steps_reconnection+2:end)=-kc(end)/time_relaxation;

%% Determine r_plus in space
% Points of plot / calculation space
x_vec=linspace(-rmix,rmix,5e2+1);
y_vec=linspace(-rmix,rmix,5e2+1);
[~,x,y]=ndgrid(1,x_vec,y_vec);

% Pole coordinates from X-point
R =sqrt(bsxfun(@plus,bsxfun(@minus,r2,x).^2,y.^2));
theta=atan(bsxfun(@times,y,1./bsxfun(@minus,r2,x)));
sin_theta=sin(theta);
sin_theta_sq=sin_theta.^2;
cos_theta_sq=1-sin_theta_sq;

%% Determine area
rr=sqrt(x.^2+y.^2);

expr_1=bsxfun(@ge,rr,r2) ;
expr_2=bsxfun(@lt,bsxfun(@plus,bsxfun(@minus,x,eta).^2,y.^2),r1.^2);
expr_3=~expr_1 & ~expr_2;

%% R_PLUS
c=bsxfun(@minus,kr_inv + k ,bsxfun(@times,R.^2./cos_theta_sq,1./r3.^2));
rho=c.*sqrt(cos_theta_sq)./(1 + sqrt(1+bsxfun(@times,k-kr,c)));

r_plus=bsxfun(@times,r3,sqrt(rho.^2+sin_theta_sq));
r_plus(expr_1 | expr_2)=NaN;

%% ELECTRIC POTENTIAL
Phi= bsxfun(@times,r2_dot,R.*sin_theta)...
    +0.5*theta.*bsxfun(@times,bsxfun(@minus,r_plus.^2,r3.^2),k_dot-kr_dot)...
    +bsxfun(@times,theta,(k-kc).*r3.*r3_dot)...             % Term should always remain 0
    -bsxfun(@times,(r1.*r1_dot+r2.*r2_dot),sin_theta.*sqrt(cos_theta_sq))...
    +bsxfun(@times,r3_dot.*r3,sin_theta).*rho;

Phi(expr_1)=NaN; % Actually 0, but this plots nicer

a=bsxfun(@times,R,r2_dot-r1_dot).*sin(theta);
Phi(expr_2)=a(expr_2);

%% PSI 
psi_star = zeros(size(Phi));
% Psi_star area 3
psi_star(expr_3) = interp1(r3(1:time_steps_reconnection+1),psi_plus(1:time_steps_reconnection+1),r_plus(expr_3),intp_method);

% Psi_star area 2
zeta = sqrt(bsxfun(@plus,bsxfun(@minus,x,eta).^2,y.^2));
psi_star(expr_2) = interp1(r,psi_min,zeta(expr_2),intp_method);

r_original = sqrt(x.^2+y.^2);
Psi_original =  repmat(interp1(r,psi_min,r_original),[size(r_plus,1),1,1]);
psi_star(expr_1) = Psi_original(expr_1);

Chi= 0.5*r_original.^2;
Psi=bsxfun(@plus,psi_star,Chi);

psi_star(psi_star(:)<0)=NaN;
%% SFC
tor_pos=pi/2;

r_ind=interp1(chi,1:length(chi),Chi);
pol_angle            =atan2(y,x);
omega=pol_angle-tor_pos;
omega = mod(omega,2*pi);

pol_ind=linspace(0,2*pi,size(maps.X_PR_map,1)); 
pol_ind=interp1(pol_ind,1:length(pol_ind),omega);

R_pos = dim.R0+ba_interp2(maps.X_PR_map,r_ind,pol_ind);
Z_pos = ba_interp2(maps.Z_PR_map,r_ind,pol_ind);

%% PLOTTING
colors=linspecer(2);
% Make figures
delete(findall(0,'type','figure','tag',mfilename))
hf1=figure('name','r-space','tag',mfilename,'windowstyle','normal','units','normalized','position',[0 0 1 1],'Color','w');
hf2=figure('name','RZ-space','tag',mfilename,'windowstyle','normal','units','normalized','position',[0 0 1 1],'Color','w');

% Make axes figure 1
plot_quant={'r_plus','psi_star','Phi'};
plot_labels={'$r_+$','$\psi_*$','$\Phi$'};

% Draw a subplot on the first figure
ha1=subplot(2,2,1,'parent',hf1);
drawnow; pos_1{1}=get(ha1,'position'); delete(ha1)
[ha1,h1,h2]=plotyy(1:10,1:10,1:10,1:10,'parent',hf1);
set(ha1(1),'position',pos_1{1},'parent',hf1); set(ha1(2),'position',pos_1{1});
delete(h1); delete(h2); linkaxes(ha1,'x');
xlim(ha1(1),[-rmix rmix]); ylim(ha1(1),[-rmix rmix]); 
hold(ha1(1),'on'); hold(ha1(2),'on');
% Plot psi and q-profiles
% Parameter for x-axis
middle_y_ax= y_vec==0;

% psi - profiles
plot(ha1(1),[-fliplr(r), r] ,[fliplr(psi_min),psi_min],'displayname','initial','linestyle','- ','color',colors(1,:));
plot(ha1(1),[-flipud(r3);r3],[flipud(psi_plus);psi_plus]   ,'displayname','final'  ,'linestyle','--','color',colors(1,:));
ylim(ha1(1),[0 max(psi_min)])

% q-profiles
pos_x_axis=x_vec>=0;
hq_ini=plot(ha1(2),[-fliplr(r), r] ,[fliplr(dim.q_initial_profile),dim.q_initial_profile],'displayname','initial','linestyle','-' ,'color',colors(2,:));
chi_final = squeeze(Chi(end,pos_x_axis,middle_y_ax));
q_final = gradient(chi_final,squeeze(Psi(end,pos_x_axis,middle_y_ax)));
r_ori   =   [-fliplr(r_original(end,pos_x_axis,middle_y_ax)), r_original(end,pos_x_axis,middle_y_ax)];
hq_final = plot(ha1(2),r_ori,[fliplr(q_final),q_final],'displayname','final ','linestyle','--','color',colors(2,:));
ylim(ha1(2),[min(dim.q_initial_profile) dim.q_initial_profile(floor(ind_rmix)+1)])

hasbehavior(hq_ini,'legend',false)
hasbehavior(hq_final,'legend',false)

set(ha1(2),'ytick',1)

ht1(1)=title(ha1(1),['$t$ = '  ,num2str(0,'%1.2e'),' [s]'],'interpreter','latex'); 

xlabel(ha1(1),'$r$ [$r_{mix}$]','interpreter','latex')
ylabel(ha1(1),'$\psi$ [Wb]','interpreter','latex'); ylabel(ha1(2),'$q$','interpreter','latex')
set(ha1(1),'YColor',colors(1,:)); set(ha1(2),'YColor',colors(2,:))
hold(ha1(1),'on');   hold(ha1(2),'on'); 

%% Add axis for other plots
% Make axes figure 1
for i=1:3
    ha2(i)=subplot(2,2,i+1,'parent',hf1);
    axis(ha2(i),'equal');
    hold(ha2(i),'on');
    xlim(ha2(i),[-1 1]);
    ylim(ha2(i),[-1 1]);
    xlabel(ha2(i),'$\tilde{x}$','interpreter','latex'); ylabel(ha2(i),'$\tilde{y}$','interpreter','latex'); 
    ht2(i)=title(ha2(i),plot_labels{i},'interpreter','latex');  
    pos_2{i}=get(ha2(i),'position'); 
    set(ha2(i),'FontSize',28,'tickLabelInterpreter','latex')
end
linkaxes(ha2,'xy')

% Make axes figure 2
if exist('hf2','var')
    c=contourc(dim.R0+dim.scale_X,dim.scale_Z,maps.psi_norm_XZsmall_map',[513 513]); % LCFS contour
    for i=1:3
        ha3(i)=subplot(1,3,i,'parent',hf2);
        axis(ha3(i),'equal'); hold(ha3(i),'on');
        plot(ha3(i),c(1,2:end),(c(2,2:end)),'r','displayname','LCFS');
        xlabel(ha3(i),'$R$ [m]','interpreter','latex'); 
        ht3(i)=title(ha3(i),plot_labels{i},'interpreter','latex'); 
        pos_3{i}=get(ha3(i),'position'); 
        set(ha3(i),'FontSize',28,'tickLabelInterpreter','latex')
    end
    linkaxes(ha3,'xy');
    xlim(ha3(1),[min(dim.scale_X) max(dim.scale_X)]+dim.R0); 
    ylim(ha3(1),[min(dim.scale_Z) max(dim.scale_Z)]);
    ylabel(ha3(1),'$Z$ [m]','interpreter','latex'); 
end

%% PLOT THE MOVIE
a=0;
%F(size(r_plus,1)) = struct('cdata',[],'colormap',[]);
for i=1:size(r_plus,1)
% figure 1
if i~=1
    delete(h1);
    delete(h2);   
end
h1(1)=plot(ha1(1),x_vec,psi_star(i,:,middle_y_ax),'linestyle','-.','displayname','current time','color',colors(1,:));

for j=1:3
    quant=squeeze(eval([plot_quant{j},'(i,:,:)']));
    if all(isnan(quant(:)))
        continue
    end
    if strcmp(plot_quant{j},'psi_star')
        quant=quant/psi_plus(1);
        set(ha2(j),'CLim',[0 1]);
    elseif strcmp(plot_quant{j},'r_plus')
        quant=quant/rmix;
        set(ha2(j),'CLim',[0 1]);
    end
    [~,h2(j)]=contourf(squeeze(x)/rmix,squeeze(y)/rmix,quant,10,'parent',ha2(j),'linestyle','none');
    set(ha2(j),'position',pos_2{j});
    cl1(j)=colorbar('peer',ha2(j));
    colormap(ha2(j),'jet');
end

% Adjust figure / title ect,
set(ht1(1),'string',['t = '  ,num2str(time(i),'%1.2e'),' [s]']); 

try
    title(cl1(1),'[$\sqrt{2\pi Wb}$]','interpreter','latex')
    title(cl1(2),'[$\psi_{*0}$]','interpreter','latex')
    title(cl1(3),'[V]','interpreter','latex')
catch
end
    caxis(ha2(1),[0 1]);
    caxis(ha2(2),[0 1]);

% figure 2
if i~=1
    delete(h3);
end

for j=1:3
    quant=squeeze(eval([plot_quant{j},'(i,:,:)']));
    if all(isnan(quant(:)))
        continue
    end
    if strcmp(plot_quant{j},'psi_star')
        quant=quant/psi_plus(1);
        set(ha3(j),'CLim',[0 1]);
    elseif strcmp(plot_quant{j},'r_plus')
        quant=quant/rmix;
        set(ha3(j),'CLim',[0 1]);
    end
    [~,h3(j)]=contourf(squeeze(R_pos),squeeze(Z_pos),quant,10,'parent',ha3(j),'linestyle','none');
    set(ha3(j),'position',pos_3{j});
%     cl2(j)=colorbar('peer',ha3(j));
    colormap(ha3(j),'jet');
end

drawnow; 
if false
    savefig(hf2,[datestr(now,'yyyy-mm-dd'),'_ST_RZ_frame_',num2str(i)],'compact');
end
%F(i) = getframe(hf2);
end

%% Wrap up
set(0,'defaulttextinterpreter',old_intp)
delete(hf1)
delete(hf2)

%% Save video
if false && exist('F','var')
    writerObj = VideoWriter([datestr(now,'yyyy-mm-dd'),'_sawtooth_RZ.mp4'],'MPEG-4');
    writerObj.FrameRate = 10;
    open(writerObj);
    for i=1:length(F)
        writeVideo(writerObj,F(i));
    end
    close(writerObj);
end
end

%% Provide time array
function time=get_time(n_1,n_2,trec,trel,option)
switch option
    case 'linear'
        %% Linear time growth
        time=[linspace(0,trec,n_1+1),linspace(trec,trec+trel,n_2+1)]';
        %time(n_1+1)=[];     % Remove the duplicate value at the transition of reconnections and relaxation
    case 'normal'
        %% Normal time growth (i.e. more points at start and transition
        norm_dist=@(mu,sigma,x) 1/sqrt(2*pi) / sigma * exp(-(x-mu).^2/(2*sigma^2)); % Normal distribution function
        
        % Distribution function with 0.5 constant and 0.5 normally distributed functions
        f_1=@(t) 0.5/trec+0.5*norm_dist(0,0.1*trec,t)+0.5*norm_dist(trec,0.1*trec,t);
        f_2=@(t) 0.5/trel+norm_dist(trec,0.1*trel,t);
        
        % Time array reconnection
        t_1=linspace(0,trec,1e7);   % Start with linear time distribution
        for i=1:10
            c=cumtrapz(t_1,f_1(t_1));   % CDF (Cumulative Distribution function
            c(end)=1;                   % Hardcode such that all time points are included
            
            t_1=interp1(c,t_1,linspace(0,1,n_1+1));   % Linearly in CDF distributed time points
        end
        if any(c>1); error('Not enough original points for proper CDF approximation'); end
        
        % Time array relaxation
        t_2=linspace(trec,trec+trel,1e7);  % Start with linear time distribution
        for i=1:10
            c=cumtrapz(t_2,f_2(t_2));   % CDF (Cumulative Distribution function
            c(end)=1;                   % Hardcode such that all time points are included
            
            t_2=interp1(c,t_2,linspace(0,1,n_2+1));   % Linearly in CDF distributed time points
        end
        if any(c>1); error('Not enough original points for proper CDF approximation'); end
        
        time=[t_1, t_2]';    
        
end
end