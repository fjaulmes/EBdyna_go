function [x,v,B,E]=time_step_integration_GT_eq_struct(x,v,Fc_field)
%%time_step_integration_GT_eq_struct
% Multifunctional timestep with BORIS or Fabien method.
global qom par

% for distance(Delta_l) calculation
x_prev=x;

%% Get B field
switch par.scheme
    case {'FabienB'}
        [E,BR,BZ,Bphi]=B_interpolation(x);
        Bfield_sq=BR.^2+BZ.^2+Bphi.^2;
    case {'BORIS'}
        [E,B]=B_interpolation(x);
        Bfield_sq=dot(B,B,2);    
end
Bfield=sqrt(Bfield_sq);

% in case the toroidal rotation is part of the input data
if par.APPLY_FC==1
    Ecx=x(:,1).*Fc_field;
    E(:,1)=E(:,1)+Ecx;
end
    
    
%% Multiplication factors
ft=0.5*qom*par.dt;
f=tan(ft*Bfield)./Bfield; % tan(a)/a-correction
pre_fac=2*f./(1+f.^2.*Bfield_sq);

switch par.scheme
    case 'FabienB'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        v=v+E*ft;
        fB2=f.*Bfield_sq;
        %% Component wise 2f/(1+fB^2) *[vxB-f|B|^2 v + fB(v.B)]
        M1=-fB2+f.*BR.^2;
        M2=Bphi+f.*BR.*BZ;
        M3=(-BZ)+f.*BR.*Bphi;
        v_X=M1.*v(:,1)+M2.*v(:,2)+M3.*v(:,3);
        
        M1=(-Bphi)+f.*BZ.*BR;
        M2=-fB2+f.*BZ.^2;
        M3=BR+f.*BZ.*Bphi;
        v_Z=M1.*v(:,1)+M2.*v(:,2)+M3.*v(:,3);
        
        M1=BZ+f.*Bphi.*BR;
        M2=(-BR)+f.*Bphi.*BZ;
        M3=-fB2+f.*Bphi.^2;
        v_phi=M1.*v(:,1)+M2.*v(:,2)+M3.*v(:,3);
        
        % next velocity half time step values
        v(:,1)=pre_fac.*v_X+v(:,1);
        v(:,2)=pre_fac.*v_Z+v(:,2);
        v(:,3)=pre_fac.*v_phi+v(:,3);
        
        v=v+E*ft;
        % Add B=output if nessecary
        if nargout>2
            B=cat(2,BR,BZ,Bphi);
        end
    case 'BORIS'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Velocity (n+1/2)
        v=v+E*ft;
        v=v+bsxfun(@times,pre_fac,...
            (fast_cross(v,B,2)-bsxfun(@times,f.*Bfield_sq,v)+bsxfun(@times,f.*dot(v,B,2),B)));
        v=v+E*ft;        
    otherwise
        error('Algorithm not known')
end

%% Position update (n+1)
R_update    =x(:,1)+v(:,1)*par.dt;  % new radius
Rphi_update =       v(:,3)*par.dt;  % excursion in phi direction

% new radius
x(:,1)=sqrt(R_update.^2+Rphi_update.^2);    % new x (n+1)
R_inv=1./x(:,1);                            % inverse radios

% new Z
x(:,2)=x(:,2)+par.dt*v(:,2);

% new phi
x(:,3)=x(:,3)+asin(Rphi_update.*R_inv);
% x(:,3)=mod(x(:,3),2*pi);  (Leave this out. Why would you anyway?)

%% Cylindrical correction (rotate velocities)
v_tmp=v;
v(:,1)=( R_update   .*v_tmp(:,1)+   Rphi_update.*v_tmp(:,3)).*R_inv;
v(:,3)=(-Rphi_update.*v_tmp(:,1)+	R_update   .*v_tmp(:,3)).*R_inv;



return
end

%% FAST CROSS PRODUCT
function cr=fast_cross(a,b,dim)
% No (inverse) permutation stuff, saves ~25% in cross product
if isequal(size(a),size(b)) && size(a,dim)==3
    switch dim
        case 1
            cr = [	a(2,:,:).*b(3,:,:)-a(3,:,:).*b(2,:,:)
                    a(3,:,:).*b(1,:,:)-a(1,:,:).*b(3,:,:)
                    a(1,:,:).*b(2,:,:)-a(2,:,:).*b(1,:,:)];
        case 2
            cr = [  a(:,2,:).*b(:,3,:)-a(:,3,:).*b(:,2,:),...
                    a(:,3,:).*b(:,1,:)-a(:,1,:).*b(:,3,:),...
                    a(:,1,:).*b(:,2,:)-a(:,2,:).*b(:,1,:)];
        otherwise
            cr=cross(a,b,dim);
    end
elseif size(a,dim)==3 && size(b,dim)==1 && isequal(size(a),size(permute(b,[2 1 3])))
    warning('Transposing field since only 1 particle is being pushed')
    cr=fast_cross(a,permute(b,[2 1 3]),dim);
else
    cr=cross(a,b,dim);
end
end