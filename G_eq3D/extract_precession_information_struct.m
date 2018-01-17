function [prec] = extract_precession_information_struct(input,output,ejected )
%extract_precession_information_struct Identifies particles in categories
%   Detailed explanation goes here
global par dim
N_job=size(output.x,1);
N_not_ejected=sum(~ejected);



if ~isfield(output,'theta_gc')
    %% Indexes 2D
    X_ind=((output.x_gc(:,1,:)-dim.R0)*dim.DX_inv)+dim.mid_Xzero;
    Z_ind=( output.x_gc(:,2,:)        *dim.DZ_inv)+dim.mid_Z;
    
    %% Indexes 3D
    % Find theta
    if par.interp_scheme==1
        [X_ind,Z_ind] = interp_index_list_2D ([dim.size_X,dim.size_Z],X_ind,Z_ind);     % Return the indexes / slopes as reps. X_ind and Z_ind
    end
    output.theta_gc = squeeze(interpolate_theta_XZ(X_ind,Z_ind));
end
%% Parallel velocity
[max_vpll,index_max]=max(abs(output.vpll),[],2);
sign_vpll=max_vpll*0;
for j=1:length(index_max)
    sign_vpll(j)=sign(output.vpll(j,index_max(j)));
end

%% Radial position
if ~isfield(dim,'radial_r_value_flux')
    dim2=load(strcat(par.paths.DATA_FOLDER,'tokamak_PR_map.mat'),'radial_r_value_flux');
    r=interp1(1:length(dim2.radial_r_value_flux),dim2.radial_r_value_flux,output.psi_norm,'*cubic');
else
    r=interp1(1:length(dim.radial_r_value_flux),dim.radial_r_value_flux,output.psi_norm,'*cubic');
end

% Find r-value
r_avg=mean(r,2);



%% Cross HFS
[HFS_theta,HFS_theta_ind]=min(abs(output.theta_gc-pi),[],2);
theta_prev=zeros(N_job,1);
theta_next=zeros(N_job,1);
for i=1:N_job
    theta_prev(i)=output.theta_gc(i,min(HFS_theta_ind(i)  ,size(output.theta_gc,2)))-output.theta_gc(i,max(HFS_theta_ind(i)-5,1));
    theta_next(i)=output.theta_gc(i,min(HFS_theta_ind(i)+5,size(output.theta_gc,2)))-output.theta_gc(i,max(HFS_theta_ind(i)  ,1));
end
minX=min(output.x(:,1,:),[],3)-dim.R0;
maxX=max(output.x(:,1,:),[],3)-dim.R0;

% two conditions need to be verified
HFS_crossing=sign(theta_prev)~=sign(theta_next);
HFS_crossing=(minX<=dim.X_axis) & HFS_crossing;

HFS_crossing=(HFS_crossing & ((HFS_theta<=(0.0012-0.0011*sqrt(r_avg/max(r_avg))))));
sigma_LFS=maxX>=dim.X_axis;
sigma_HFS=minX<dim.X_axis;

%% Bounce
sign_change=any(bsxfun(@ne,sign(output.vpll),sign(output.vpll(:,1))),2);

%% Categorizing
% nb_STAGNATION_HFS=sum (~sign_change &   HFS_crossing & ~sigma_LFS);

prec.pop.ALL_TRAPPED=        sign_change &  ~HFS_crossing                        ;
prec.pop.ALL_PASSING=       ~sign_change &   HFS_crossing &  sigma_LFS           ...
                     |(~sign_change &  ~HFS_crossing &                  sigma_HFS);
                 
prec.pop.TRAPPED_MINUS=prec.pop.ALL_TRAPPED   &                                                   (sign_vpll==-1)      ;
prec.pop.TRAPPED_PLUS= prec.pop.ALL_TRAPPED   &                                                   (sign_vpll== 1)      ;

prec.pop.POTATOES=           sign_change &   HFS_crossing &  sigma_LFS           ;

prec.pop.STAGNATION=       (~sign_change &  ~HFS_crossing &  sigma_LFS &    ~sigma_HFS)...
                     |(~sign_change &   HFS_crossing & ~sigma_LFS)          ;

prec.pop.CO_PASSING=       (~sign_change &	HFS_crossing &  sigma_LFS &                     (sign_vpll== 1))...
                |     (~sign_change &  ~HFS_crossing &  sigma_LFS &     sigma_HFS &     (sign_vpll== 1));

prec.pop.COUNTER_PASSING=  (~sign_change &   HFS_crossing &  sigma_LFS &                     (sign_vpll==-1))...
                |     (~sign_change &  ~HFS_crossing &  sigma_LFS &     sigma_HFS &     (sign_vpll==-1));
            
%% Filter only ejected

fnames=fieldnames(prec.pop);
for i=1:length(fnames)
    prec.pop.(fnames{i})=prec.pop.(fnames{i})&~ejected;
end

%% Find indexes (obsolete, easily extracted in later stage)
% prec.ind.ALL_TRAPPED=        input.particle_nr(prec.pop.ALL_TRAPPED);
% prec.ind.ALL_PASSING=        input.particle_nr(prec.pop.ALL_PASSING);
% prec.ind.TRAPPED_PLUS=       input.particle_nr(prec.pop.TRAPPED_PLUS);
% prec.ind.TRAPPED_MINUS=      input.particle_nr(prec.pop.TRAPPED_MINUS);
% prec.ind.POTATOES=           input.particle_nr(prec.pop.POTATOES);
% prec.ind.STAGNATION=         input.particle_nr(prec.pop.STAGNATION);
% prec.ind.CO_PASSING=         input.particle_nr(prec.pop.CO_PASSING);
% prec.ind.COUNTER_PASSING=    input.particle_nr(prec.pop.COUNTER_PASSING);

%% Display
disp('---------------------------------------------');
disp('-- GROUPS OF PARTICLE ORBITS (PERCENTAGES) --');
disp('---------------------------------------------');
disp('POTATOES_LIST=');
disp(100*sum(prec.pop.POTATOES)       /N_not_ejected);
disp('STAGNATION_LIST=');
disp(100*sum(prec.pop.STAGNATION)     /N_not_ejected);

disp('TRAPPED_MINUS=');
disp(100*sum(prec.pop.TRAPPED_MINUS)  /N_not_ejected);
disp('TRAPPED_PLUS=');
disp(100*sum(prec.pop.TRAPPED_PLUS)   /N_not_ejected);

disp('CO_PASSING=');
disp(100*sum(prec.pop.CO_PASSING)     /N_not_ejected);
disp('COUNTER_PASSING=');
disp(100*sum(prec.pop.COUNTER_PASSING)/N_not_ejected);

disp('TOTAL=');
disp(100/N_not_ejected*sum(...
        prec.pop.COUNTER_PASSING...
    |   prec.pop.CO_PASSING ...
    |   prec.pop.TRAPPED_PLUS ...
    |   prec.pop.TRAPPED_MINUS ...
    |   prec.pop.STAGNATION ...
    |   prec.pop.POTATOES));

%% Average q/psi
prec.psi_avg=mean(output.psi_gc,2);
prec.q_avg=mean(output.q_gc,2);

%% Precession / bounce frequencies (for trapped particles)
% (moved to evaluate_output, so it is also determined in 3D fields)
prec.wb=output.wb;
prec.wd=output.wd;

end
