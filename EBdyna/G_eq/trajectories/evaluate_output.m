function [ output ] = evaluate_output(input,output,ejected,fields_retained)
global const par dim maps
if par.CALCULATE_CX
CX_NEUTRALS=output.CX_NEUTRALS;
else
CX_NEUTRALS=ejected*0;
end
%evaluate_output Finds certain physical parameters after simulation
%   Detailed explanation goes here
narginchk(3,4)
if ~exist('ejected','var') || isempty(ejected)
    ejected=isnan(output.x(:,1,end));
end

ejected_straight_away=ejected*0;
%X_ind=((output.x(:,1,1)-dim.R0)*dim.DX_inv)+dim.mid_Xzero;
%Z_ind=( output.x(:,2,1)        *dim.DZ_inv)+dim.mid_Z;
psi_value=ejected*0;

%security
output.x(:,2,:)=max(output.x(:,2,:),dim.Zlow_lim);
output.x(:,2,:)=min(output.x(:,2,:),dim.Zup_lim);
output.x(:,1,:)=max(output.x(:,1,:),dim.Rlow_lim); 
output.x(:,1,:)=min(output.x(:,1,:),dim.Rup_lim); 

			
%psi_value=ba_interp2(maps(1).psi_XZ,Z_ind(),X_ind(),'cubic'); % find psi-value

par.RECORD_PRECISION=round(par.NB_TIME_STEPS/par.NB_TIME_STAMPS);
ejected_straight_away=ejected_straight_away | output.time_step_loss<par.RECORD_PRECISION;


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
output.v_plus_1=output.v;
if ~par.CALCULATE_CX
[~,output.v_plus_1,B]=time_step_integration_GT_eq_struct(output.x,output.v,Fc_field);
else
[~,output.v_plus_1(~CX_NEUTRALS,:,:),B(~CX_NEUTRALS,:,:)]=time_step_integration_GT_eq_struct(output.x(~CX_NEUTRALS,:,:),output.v(~CX_NEUTRALS,:,:),Fc_field(~CX_NEUTRALS));
output.v_plus_1(CX_NEUTRALS,:,:)=output.v(CX_NEUTRALS,:,:);
end
par.dt=2*par.dt;

%% Find sfl-coordinates in case no par.coord_syst isn't 'flux' and hence not already in previous output
X_ind=((output.x(:,1,:)-dim.R0)*dim.DX_inv)+dim.mid_Xzero;
Z_ind=( output.x(:,2,:)        *dim.DZ_inv)+dim.mid_Z;
if isfield(maps(1),'psi_norm_XZ')  % for sawteeth: psi between 1 and dim.NB_PSI
    output.psi_norm=ba_interp2(maps(1).psi_norm_XZ,Z_ind,X_ind,'cubic');
end

if par.interp_scheme==1
    [X_ind,Z_ind] = interp_index_list_2D ([dim.size_X,dim.size_Z],X_ind,Z_ind);     % Return the indexes / slopes as reps. X_ind and Z_ind
end
expr_calc=isfinite(output.x(:,1,:));
if isfield (maps(1),'theta_XZ')
output.theta = interpolate_theta_XZ(X_ind,Z_ind,expr_calc);
end

%% B field and speed
% for t_ind = size(output.x,3)
%     B_t_ind=ba_interp2(maps(1).B_2D       ,squeeze(Z_ind(:,t_ind)),squeeze(X_ind(:,t_ind)),'cubic');
%     B(:,:,t_ind)=squeeze(B_t_ind);
% end
B=ba_interp2(maps(1).B_2D       ,Z_ind,X_ind,'cubic');
B=squeeze(B);
B=permute(B,[1 3 2]);

Bfield=sqrt(dot(B,B,2));
b=bsxfun(@times,Bfield.^-1,B);

% Parallel speed
output.vpll=            dot(output.v_plus_1 ,b              ,2);
v_plus_1_sq=            dot(output.v_plus_1 ,output.v_plus_1,2);
v_perp_sq  =v_plus_1_sq-dot(output.vpll     ,output.vpll    ,2);

%% GC_of orbit : causes trouble when neutrals are present .....
gc_x_corr=output.x*0;
if ~par.CALCULATE_CX
% supposed to be better than the real time calculation, since it uses cubic interpolation
gc_x_corr=bsxfun(@times,(input.m/const.eV)*1./(input.Z*Bfield),cross(output.v_plus_1,b,2));
output.x_gc=output.x+gc_x_corr;
else
gc_x_corr(~CX_NEUTRALS,:,:)=bsxfun(@times,(input.m/const.eV)*1./(input.Z*Bfield(~CX_NEUTRALS,:,:)),cross(output.v_plus_1(~CX_NEUTRALS,:,:),b(~CX_NEUTRALS,:,:),2));
output.x_gc(~CX_NEUTRALS,:,:)=output.x(~CX_NEUTRALS,:,:)+gc_x_corr(~CX_NEUTRALS,:,:);
output.x_gc(CX_NEUTRALS,:,:)=output.x(CX_NEUTRALS,:,:);
end
% still some outliers seem to cause problem!
output.x_gc=real(output.x_gc);
%security
output.x_gc(:,2,:)=max(output.x_gc(:,2,:),dim.Zlow_lim);
output.x_gc(:,2,:)=min(output.x_gc(:,2,:),dim.Zup_lim);
output.x_gc(:,1,:)=max(output.x_gc(:,1,:),dim.Rlow_lim); 
output.x_gc(:,1,:)=min(output.x_gc(:,1,:),dim.Rup_lim); 


%% Radial index gyro center
X_ind_gc=((output.x_gc(:,1,:)-dim.R0)*dim.DX_inv)+dim.mid_Xzero;
Z_ind_gc=( output.x_gc(:,2,:)        *dim.DZ_inv)+dim.mid_Z;

if isfield(maps(1),'psi_norm_XZ')
    output.psi_norm_gc=ba_interp2(maps(1).psi_norm_XZ,Z_ind_gc,X_ind_gc,'linear');
end
if par.interp_scheme==1
    [X_ind_gc,Z_ind_gc] = interp_index_list_2D ([dim.size_X,dim.size_Z],X_ind_gc,Z_ind_gc);     % Return the indexes / slopes as reps. X_ind and Z_ind
end
if isfield (maps(1),'theta_XZ')
output.theta_gc = interpolate_theta_XZ(X_ind_gc,Z_ind_gc,expr_calc);
end

B=ba_interp2(maps(1).B_2D       ,Z_ind_gc,X_ind_gc,'cubic');
B=squeeze(B);
B=permute(B,[1 3 2]);

Bfield=sqrt(dot(B,B,2));
b=bsxfun(@times,Bfield.^-1,B);

% Parallel speed at gc
output.vpll=            dot(output.v_plus_1 ,b              ,2);
v_plus_1_sq=            dot(output.v_plus_1 ,output.v_plus_1,2);
v_perp_sq  =v_plus_1_sq-dot(output.vpll     ,output.vpll    ,2);

% Make upper limit to prevent invalid 1D interpolation
% output.psi_norm   (output.psi_norm   >dim.NB_PSI)=dim.NB_PSI; 
% output.psi_norm_gc(output.psi_norm_gc>dim.NB_PSI)=dim.NB_PSI;

%% Psi-surfaces
% output.psi_norm=squeeze(output.psi_norm)
% output.psi_norm_gc=squeeze(output.psi_norm_gc)
% output.psi=interp1(1:dim.NB_PSI,dim.psi_scale,output.psi_norm,'pchip');
% output.psi_gc=interp1(1:dim.NB_PSI,dim.psi_scale,output.psi_norm_gc,'pchip');

% extended interpolation to the whole domain for psi!
% output.psi=squeeze(output.psi_norm)*0;
% output.psi_gc=squeeze(output.psi_norm)*0;
output.psi=output.vpll*0;
output.psi_gc=output.vpll*0;

for TS=1:size(output.vpll,3)
%     output.psi(:,TS)=interp2(dim.scale_X,dim.scale_Z,maps(1).psi_XZ',squeeze(output.x(:,1,TS))-dim.R0,squeeze(output.x(:,2,TS)),'cubic');
%     output.psi_gc(:,TS)=interp2(dim.scale_X,dim.scale_Z,maps(1).psi_XZ',squeeze(output.x_gc(:,1,TS))-dim.R0,squeeze(output.x_gc(:,2,TS)),'cubic');
    % unfortunate dimensions because pphi evaluation functin used further down
    output.psi(:,1,TS)=interp2(dim.scale_X,dim.scale_Z,maps(1).psi_XZ',output.x(:,1,TS)-dim.R0,output.x(:,2,TS),'cubic');
    output.psi_gc(:,1,TS)=interp2(dim.scale_X,dim.scale_Z,maps(1).psi_XZ',output.x_gc(:,1,TS)-dim.R0,output.x_gc(:,2,TS),'cubic');
end

%% q-values
% q_prf=load(strcat(par.folders.DATA_SHOT,'q_profile.mat'),'q_initial_profile');
% output.q   =interp1(1:dim.NB_PSI,q_prf.q_initial_profile,output.psi_norm,'pchip');
% output.q_gc=interp1(1:dim.NB_PSI,q_prf.q_initial_profile,output.psi_norm_gc,'pchip');

% find outermost flux surface position
psi_outer=max(squeeze(output.psi)')';
psi_gc_outer=max(squeeze(output.psi_gc)')';

output.psi_outer=psi_outer;
output.psi_gc_outer=psi_gc_outer;

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
output.Delta_pphi2=squeeze(output.pphi_kin(:,:,end)-output.pphi_kin(:,:,1));

%% Larmor frequency
output.larmor_freq=input.Z*const.eV*Bfield/input.m;

%% Magnetic moment
Eperp=0.5*(input.m/const.eV)*v_perp_sq;
output.mm=Eperp./Bfield;

%% Kinetic energy
output.Ekin=0.5*(input.m/const.eV)*v_plus_1_sq;

output.Ekin_end=output.Ekin(:,:,end);

%% Profiles
% using guiding center positions
if verLessThan('matlab','8.4')
    output.nr_profile=histc(output.psi_gc,dim.scale_psi_mp,1); 	    % Histogram count for mid-plane profile
    output.nr_profile(end,:)=[]; 								% Particles at psi_gc>=scale_psi_mp(end) do not count
else
    make_profile_2D=false; % don't do that
    if make_profile_2D && ~verLessThan('matlab','8.6')
        output.profile_2D=zeros(dim.size_X,dim.size_Z,par.NB_TIME_STAMPS);
    end
	% counting
    if make_profile_2D
        output.nr_profile=zeros(length(dim.scale_psi_mp)-1,par.NB_TIME_STAMPS);
        for i=1:par.NB_TIME_STAMPS
    %         output.nr_profile(:,i)=histcounts(output.psi_gc(:,:,i),dim.scale_psi_mp);
            if make_profile_2D && ~verLessThan('matlab','8.6')
                edges_X=[dim.scale_X-dim.DX dim.scale_X(end)+dim.DX]+dim.R0;
                edges_Z=[dim.scale_Z-dim.DZ dim.scale_Z(end)+dim.DZ];
                [output.profile_2D(:,:,i)] = histcounts2(output.x_gc(:,1,i),output.x_gc(:,2,i),edges_X,edges_Z);
            end
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