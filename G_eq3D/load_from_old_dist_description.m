    %% 1.1.1. Split
    input.N_total=dist.Nalphas_simulated;
    expr_job=get_expr_job([input.N_total 1],par.NB_PROCESS,par.PROCESS_NUMBER);  % Get expression which particles this job will simulate
    input.particle_nr=find(expr_job);
    input.N_job=sum(expr_job(:)); % Number of particles in simulation
    %% 1.1.2. Load general parameters
    input.Ekin=dist.alphas_Ekin(expr_job);%Kinetic energy
    input.m=dist.mHe;                   % Mass
    input.Z=dist.ZHe;                   % Charge
    qom=input.Z*const.eV/input.m;       %Kinectic constant
    %% 1.1.3. Location
    x=zeros(input.N_job,3);
    
    x(:,1)=dim.R0+dist.alphas_pos_x(expr_job);
    x(:,2)=dist.alphas_pos_z(expr_job);
    x(:,3)=dist.alphas_pos_phi(expr_job);
    input.x=x;
%     input.x_gc=x;
    
    
    %% 1.1.4. Velocity
    if isfield(dist,'v_X')
        % 1.1.4.1.1 Directly from file
        v=zeros(input.N_job,3);
        v(:,1)=dist.v_X(expr_job);
        v(:,2)=dist.v_Z(expr_job);
        v(:,3)=dist.v_phi(expr_job);
    elseif isfield(dist,'alphas_vpll')
        %1.1.4.2.1 Find perpendicular velocity
        vpll=dist.alphas_vpll(expr_job);
        vtot=E_to_v(input.m,input.Ekin);
        vperp=sqrt(vtot.^2-vpll.^2);
        % 1.1.4.2.2. Determine v by direction of b
        [N,b]=make_normal_vector(input,x,input.N_job); % Vector of perpendicular velocity
        v=bsxfun(@times,vpll,b)+bsxfun(@times,vperp,N);
    else
        error('Distribution file has improper velocity data')
    end
    input.v=v;

    % consider plasma rotation effects
    par.APPLY_FC=0;
    if isfield(dist,'alphas_momentum')
        input.ang_rot=dist.alphas_momentum(expr_job);
		input.Fc_field=input.ang_rot.^2/qom;
        par.APPLY_FC=1;
	else
		input.Fc_field=x(:,1)*0;
    end
    if isfield(dist,'alphas_weight')
        input.weight=dist.alphas_weight(expr_job);
	else
		input.weight=x(:,1)*0+1;
    end