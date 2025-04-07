%% CONTINUE FROM RAW
% Check if a raw file exist
if exist(par.SAVENAME_RAW,'file') && ~(par.mode==2 && exist([par.SAVENAME_RAW,'prec'],'file'))
    load_raw=true;                  % Set as true, flip to false if parameters do not match
    raw=load(par.SAVENAME_RAW);
    
    % First check if the mode is consistent
    if par.mode~=raw.par.mode
        load_raw=false;
    else
        % Check other fields that need to coincide
        fields_to_check = {'PROCESS_NUMBER','NB_PROCESS','TEST_DIST','ID','interp_scheme','coord_syst','APPLY_RMP','APPLY_TFR','APPLY_SAWTOOTH','DISTNAME','dt'};
        for i=1:length(fields_to_check)
            if ~isequal(par.(fields_to_check{i}),raw.par.(fields_to_check{i}))
                warning(['Continueing from raw file failed since ',fields_to_check{i},' do not coincide'])
                load_raw=false;
            end
        end
        % Check if the input is equal
        if isfield(raw,'input') && ~isequal(input.particle_nr,raw.input.particle_nr)
            warning('Particles not equal in RAW file')
            load_raw=false;
        end
        % Replace the data with the data from the RAW file
        if load_raw
            warning('Continueing from old raw file')
            
            % Load output
            output=raw.output;
            if par.mode==5 % POINCARE
                par.time_scale=raw.par.time_scale;
                % Find last saved position of each particle
                for i=1:input.N_job
                    ind_time_i=find(isfinite(output.x(i,1,:)),1,'last');
                    if ~isempty(ind_time_i)
                        ind_time(i)=ind_time_i;
                        % Get x and v only if one position is stored
                        x(i,:)=output.x(i,:,ind_time(i));
                        v(i,:)=output.v(i,:,ind_time(i));
                    elseif isfield(raw,'input')
                        ind_time(i)=0;
                        x(i,:)=raw.input.x(i,:);
                        v(i,:)=raw.input.v(i,:);
                    end
                end
                if strcmp(par.poincare_type,'xy')
                    X_ind=((output.x_gc(:,1,:)-dim.R0)*dim.DX_inv)+dim.mid_Xzero;
                    Z_ind=( output.x_gc(:,2,:)        *dim.DZ_inv)+dim.mid_Z;
                    expr_nan=isnan(X_ind); expr_nan=expr_nan(:);
                    
                    % Find theta
                    theta_large(~expr_nan) = interpolate_theta_XZ(X_ind(~expr_nan),Z_ind(~expr_nan));
                    chi=ba_interp2(maps(1).chi_XZ,Z_ind(~expr_nan),X_ind(~expr_nan),'linear');
                    r= sqrt(2*chi);
                    phi=output.x(:,3,:);
                    phi=phi(~expr_nan);
                    % Determine x and y-coordinates of particle in question
                    output.x_tilde(~expr_nan) = r .* cos(theta_large(~expr_nan)-phi);
                    output.y_tilde(~expr_nan) = r .* sin(theta_large(~expr_nan)-phi);
                end
            else        % ALL OTHER
                if any(size(par.time_scale)~=size(raw.par.time_scale)) || any(par.time_scale~=raw.par.time_scale) || par.dt~=raw.par.dt
                    error('Times of RAW file don''t match required simulation')
                end
                if ~any(isfinite(output.x(:)))
                    error('Raw file has no output data (all lost?')
                end
                par.comment=raw.par.comment;
                if isfield(raw,'process_time')
                    disp(['Process ',num2str(par.PROCESS_NUMBER),' simulation recovered'])
                    time_stamp=length(par.time_scale);
                    time_step=par.NB_TIME_STEPS;
                    process_time_raw=raw.process_time;
                    ejected=raw.ejected;
                else
                    time_stamp=find(any(isfinite(output.x(:,1,:)),1),1,'last');
                    % Find those who have been ejected
                    ejected=isnan(output.x(:,1,time_stamp));
                end
                
                x=output.x(:,:,time_stamp);
                v=output.v(:,:,time_stamp);
                
                % Time is set to continue from this point on
                start_time_step=round(time_stamp*par.TIME_STAMP_PRECISION)+1;
                time=par.time_scale(time_stamp);
            end
        end
    end
    clear raw
end