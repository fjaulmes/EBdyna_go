function build_homogenuous_dist(N_part)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if nargin==0
    N_part=6e6;
end



seed=1;
% seed='shuffle';
rng(seed);
%% Load default parameters
par.paths=initialize_folder_names_struct;   % Folder names
[const,maps,dim]=load_distr_maps(par);

E_to_v=@(m,E) sqrt(2*E.*(const.eV/m));
%% Positions
x=zeros(N_part+1,1);
i=2;
while N_part~=size(x,1)
    clear x
    R=dim.R0+dim.scale_X(1)+(dim.scale_X(end)-dim.scale_X(1))*rand(N_part*i,1);
    Z=dim.scale_Z(1)+(dim.scale_Z(end)-dim.scale_Z(1))*rand(N_part*i,1);

    x(:,1)=R;
    x(:,2)=Z;
    X_ind=((x(:,1,:)-dim.R0)*1./dim.DX)+dim.mid_Xzero;
    Z_ind=( x(:,2,:)        *1./dim.DX)+dim.mid_Z;
    psi=interp2(maps(1).psi_XZ,Z_ind,X_ind,'linear'); % find psi-value
    ejected=psi<=0 | isnan(psi);
    x(ejected,:)=[];
    if size(x,1)<N_part
        i=i+1;
        disp('Increasing random sampling since not enough particles initialized within LCFS')
        continue
    else
        x=x(1:N_part,:);
        break
    end
end

x(:,3)=rand(N_part,1)*2*pi;
%% 
input.m=const.mD;
input.Z=1;

input.Ekin=20+(rand(size(x,1),1)*40)*1e3;

%% 
vtot=E_to_v(input.m,input.Ekin);
vpll=(rand(N_part,1)*2-1).*vtot;


%% SAVE
save('../particles_equilibrium/input/G_eq_19xx_particles.mat','-v7.3','input','x','vpll');
end

