function [ output ] = evaluate_output(input,output,ejected,fields_retained)
global const par dim maps
%evaluate_output Finds certain physical parameters after simulation
%   Detailed explanation goes here
narginchk(3,4)
if ~exist('ejected','var') || isempty(ejected)
    ejected=isnan(output.x(:,1,end));
end

ejected_straight_away=ejected*0;
X_ind=((output.x(:,1,1)-dim.R0)*dim.DX_inv)+dim.mid_Xzero;
Z_ind=( output.x(:,2,1)        *dim.DZ_inv)+dim.mid_Z;
psi_value=ejected*0;
psi_value(~ejected)=ba_interp2(maps(1).psi_XZ,Z_ind(~ejected),X_ind(~ejected),'linear'); % find psi-value
if dim.psi_scale(1)<0  		% psi_Scale increasing to 0
	ejected_straight_away=ejected_straight_away | psi_value >= 0;
else 						% psi_scale decreasing to 0
	ejected_straight_away=ejected_straight_away | psi_value <= 0;
end

%% TRUCTATION (optional)
truncate_output=false;
% One can choose to `truncate' the arrays (which had to be used in older versions of evaluate_output.m)
% In that case the particles which have been ejected before the first time
% step will not be `evaluated'. However, these do not pop up

% ejected_straight_away=ejected_straight_away | isnan(output.x(:,1,1));
% if any(ejected_straight_away)
%     warning('Truncating output during evaluation since some are ejected before first time stamp')
%     truncate_output=true;
%     output=truncate(output,ejected_straight_away);
%	 % for delta calculations, these particles need also to be removed from the input
%     input=truncate(input,ejected_straight_away);
%     ejected=ejected(~ejected_straight_away);
%     input.N_job=sum(~ejected_straight_away);
% else
%    truncate_output=false;
% end

%% Speed at n+1
par.dt=0.5*par.dt;
Fc_field=input.Fc_field;
[~,output.v_plus_1,B]=time_step_integration_GT_eq_struct(output.x,output.v,Fc_field);
par.dt=2*par.dt;

%% Find sfl-coordinates in case no par.coord_syst isn't 'flux' and hence not already in previous output
X_ind=((output.x(:,1,:)-dim.R0)*dim.DX_inv)+dim.mid_Xzero;
Z_ind=( output.x(:,2,:)        *dim.DZ_inv)+dim.mid_Z;
output.psi_norm=ba_interp2(maps(1).psi_norm_XZ,Z_ind,X_ind,'linear');

if par.interp_scheme==1
    [X_ind,Z_ind] = interp_index_list_2D ([dim.size_X,dim.size_Z],X_ind,Z_ind);     % Return the indexes / slopes as reps. X_ind and Z_ind
end
expr_calc=isfinite(output.x(:,1,:));
output.theta = interpolate_theta_XZ(X_ind,Z_ind,expr_calc);

%% B field and speed
Bfield=sqrt(dot(B,B,2));
b=bsxfun(@times,Bfield.^-1,B);

% Parallel speed
output.vpll=            dot(output.v_plus_1 ,b              ,2);
v_plus_1_sq=            dot(output.v_plus_1 ,output.v_plus_1,2);
v_perp_sq  =v_plus_1_sq-dot(output.vpll     ,output.vpll    ,2);

%% GC_of orbit 
output.x_gc=output.x+bsxfun(@times,(input.m/const.eV)*1./(input.Z*Bfield),cross(output.v_plus_1,b,2));

%% Radial index gyro center
X_ind_gc=((output.x_gc(:,1,:)-dim.R0)*dim.DX_inv)+dim.mid_Xzero;
Z_ind_gc=( output.x_gc(:,2,:)        *dim.DZ_inv)+dim.mid_Z;

output.psi_norm_gc=ba_interp2(maps(1).psi_norm_XZ,Z_ind_gc,X_ind_gc,'linear');
if par.interp_scheme==1
    [X_ind_gc,Z_ind_gc] = interp_index_list_2D ([dim.size_X,dim.size_Z],X_ind_gc,Z_ind_gc);     % Return the indexes / slopes as reps. X_ind and Z_ind
end
output.theta_gc = interpolate_theta_XZ(X_ind_gc,Z_ind_gc,expr_calc);

% Make upper limit to prevent invalid 1D interpolation
output.psi_norm   (output.psi_norm   >dim.NB_PSI)=dim.NB_PSI; 
output.psi_norm_gc(output.psi_norm_gc>dim.NB_PSI)=dim.NB_PSI;

%% Psi-surfaces
output.psi=interp1(1:dim.NB_PSI,dim.psi_scale,output.psi_norm,'pchip');
output.psi_gc=interp1(1:dim.NB_PSI,dim.psi_scale,output.psi_norm_gc,'pchip');

%% q-values
q_prf=load(strcat(par.paths.DATA_FOLDER,'q_profile.mat'),'q_initial_profile');
output.q   =interp1(1:dim.NB_PSI,q_prf.q_initial_profile,output.psi_norm,'pchip');
output.q_gc=interp1(1:dim.NB_PSI,q_prf.q_initial_profile,output.psi_norm_gc,'pchip');

%% Pphi-calculation
norm_fac=input.Z*...
    (maps(1).psi_global); % normalization factor for pphi = qRA_phi (center - edge)
if ~isfield(output,'pphi_kin')
    [ output.pphi_kin ] = get_pphi_kin(input,output.x,output.v_plus_1,output.psi);
    output.pphi_kin=output.pphi_kin/norm_fac;
end

% Analytical value (0 in 2D)
if par.CALCULATE_PPHI_3D && isfield(output,'pphi_an_temp')
    output.pphi_an =output.pphi_an_temp /norm_fac;
    output.pphi_an=bsxfun(@plus,squeeze(output.pphi_kin(:,1)),output.pphi_an);
    output=rmfield(output,'pphi_an_temp');
end

% Delta pphi
output.Delta_pphi=output.pphi_kin(:,:,end)-input.pphi_kin;

%% Larmor frequency
output.larmor_freq=input.Z*const.eV*Bfield/input.m;

%% Magnetic moment
Eperp=0.5*(input.m/const.eV)*v_perp_sq;
output.mm=Eperp./Bfield;

%% Kinetic energy
output.Ekin=0.5*(input.m/const.eV)*v_plus_1_sq;

output.Ekin_end=output.Ekin(:,:,end);

%% Profiles
if verLessThan('matlab','8.4')
    output.nr_profile=histc(output.psi_norm,0:dim.NB_PSI,1); 	% Histogram count for E-profile ()
    output.nr_profile(end,:)=[]; 								% Particles at psi_norm>=dim.psi do not count
else
    output.nr_profile=zeros(dim.NB_PSI,par.NB_TIME_STAMPS);
    make_profile_2D=false;
    if make_profile_2D && ~verLessThan('matlab','8.6')
        output.profile_2D=zeros(dim.size_X,dim.size_Z,par.NB_TIME_STAMPS);
    end
	% counting
    for i=1:par.NB_TIME_STAMPS
        output.nr_profile(:,i)=histcounts(output.psi_norm(:,:,i),0:dim.NB_PSI);
        if make_profile_2D && ~verLessThan('matlab','8.6')
            edges_X=[dim.scale_X-dim.DX dim.scale_X(end)+dim.DX]+dim.R0;
            edges_Z=[dim.scale_Z-dim.DZ dim.scale_Z(end)+dim.DZ];
            [output.profile_2D(:,:,i)] = histcounts2(output.x(:,1,i),output.x(:,2,i),edges_X,edges_Z);
        end
    end
end

%% Squeeze sizes
fnames_output=fieldnames(output);
for i=1:length(fnames_output)
    output.(fnames_output{i})=squeeze(output.(fnames_output{i}));
end

%% dpphi_dt and dmm_dt
waitbar_switch=false; % Switch to have a waitbar during fitting of pphi and mm

if waitbar_switch
    delete(findall(0,'tag',mfilename));
    set(0,'defaulttextinterpreter','tex');
    wb=waitbar(0,'Fitting pphi and mm','CreateCancelBtn','setappdata(gcbf,''cancelling'',1)','tag',mfilename);
    up_wb=round(input.N_job/1000);
end
output.dpphi_dt=zeros(input.N_job,2);
output.dmm_dt=zeros(input.N_job,2);

output.sim_exit=sum(isfinite(output.pphi_kin),2)+1;
output.sim_exit(all(isfinite(output.pphi_kin),2))=NaN;

for i=1:input.N_job
    if ~ejected(i)
        output.dmm_dt  (i,:)=polyfit(par.time_scale,output.mm      (i,:),1);
        output.dpphi_dt(i,:)=polyfit(par.time_scale,output.pphi_kin(i,:),1);
    elseif output.sim_exit(i)>2
        output.dmm_dt  (i,:)=polyfit(par.time_scale(1:(output.sim_exit(i)-1)),output.mm      (i,1:(output.sim_exit(i)-1)),1);
        output.dpphi_dt(i,:)=polyfit(par.time_scale(1:(output.sim_exit(i)-1)),output.pphi_kin(i,1:(output.sim_exit(i)-1)),1);
    end
    if waitbar_switch && mod(i,up_wb)==0
        waitbar(i/input.N_job,wb,...
            ['Fitting pphi and mm @ ',num2str(i/input.N_job*100,'%10.1f'),'%'])
        if getappdata(wb,'cancelling')
            delete(wb)
            set(0,'defaulttextinterpreter','latex')
            error('Cancel button pressed')
        end
    end
end
% Remove the constant term in the linear fit
output.dmm_dt(:,2)=[];
output.dpphi_dt(:,2)=[];

if waitbar_switch
    delete(wb); drawnow
    set(0,'defaulttextinterpreter','latex');
end

%% Identify if the lost particles were prompt or not
if isfield(output,'nr_midplane_crossing')
    output.promp_loss=output.nr_midplane_crossing<2 & ejected;
end

%% Identify which particles might have had a ST interaction
if par.APPLY_SAWTOOTH
    output.ST_interaction=any(output.psi_norm<=1.1*dim.st.ind_rmix,2);
end

%% New precession / bounce frequency
output.wd=NaN(input.N_job,2);
output.wb=NaN(input.N_job,1);

% Precession / bounce frequencies (for trapped particles)
% Particles which are trapped (based on 2 vpll crossings) 
sign_initial=sign(output.vpll(:,1)); % Sign of first vpll-crossing
pop_cross=bsxfun(@eq,-sign_initial*2,diff(sign(output.vpll),1,2));
full_orbit=sum(pop_cross,2)>1; %True if particle has at least one full orbit from vpll_crossings
ind_full_orbit=find(full_orbit);

for i=1:input.N_job
    % wd
    output.wd(i,:)=polyfit(par.time_scale',squeeze(output.x(i,3,:)),1);
    
    if any(i==ind_full_orbit)
        % wb
        ind_cross=find(pop_cross(i,:)); % Crossing same direction
        T=(par.time_scale(ind_cross(end))-par.time_scale(ind_cross(1)))/(length(ind_cross)-1);
        output.wb(i)=2*pi/T;
    end
end

output.wd(:,2)=[];

%% Identify locally trapped particles
output.ripple_trapped=output.nr_vpll_crossing>3 & output.nr_midplane_crossing==0;

%% Remove any fields which aren't requested
if exist('fields_retained','var')
    output=remove_fields(output,[],fields_retained);
end

%% Re-size the output if truncation has been done
if truncate_output
    output_names=fieldnames(output); %#ok<UNRCH>
    for i=1:length(output_names)
        if size(output.(output_names{i}),1)~=input.N_job;  continue; end;
        temp_output=output.(output_names{i});           % Temporary store truncated output
        truncated_size=size(output.(output_names{i}));  % The truncated size
        original_size=truncated_size;                   % Original size 2, 3 etc. set correctly
        original_size(1)=length(ejected_straight_away); % Original size 1 is the number of original particles
        
        % Reset output value
        if islogical(output.(output_names{i}))      
            output.(output_names{i})=true(original_size);   % True, either prompt loss / ejected
        else
            output.(output_names{i})=NaN(original_size);
        end
        % Fill output with temporary stored array
        output.(output_names{i})(~ejected_straight_away,:,:)=temp_output;
    end
%     
%     input.N_job=length(ejected_straight_away);
%     
%     ejected_temp = ejected;
%     ejected = ejected_straight_away;
%     ejected(~ejected_straight_away) = ejected_temp;
end
end

%% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Truncate the output since some are ejected before first time stamp
function struct=truncate(struct,expr)
N_job=length(expr);
names=fieldnames(struct);
for i=1:length(names)
    if size(struct.(names{i}),1)==N_job
        struct.(names{i})=struct.(names{i})(~expr,:,:);
    end
end
end