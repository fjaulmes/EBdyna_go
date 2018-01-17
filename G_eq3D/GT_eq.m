function [exitcode]=GT_eq(x,v,input,output,ejected)
global maps dim par time qom
%GT_eq simulation code for integration of particle motion
%   Processes particles according to parameters, maps, dim etc.
%
%   Works with modes:
% [] - DEFAULT (first 2 then 3)
% 1  - TEST
% 2  - PRECESSION
%       Will use a (shorter) simulation in 2D field and categorize particles.
%       Determines bounce / drift frequencies
%       Calls G_eq with same 'box'-nr for full simulation
%       Produces:
%           RAW     - file:     In which par and output are stored
%           STATS   - file:     With precession data.
%               .mat        - succesvol
%               .matprec    - failed. Can be used for new precession or is used if all particles in equilibrium were ejected
% 3/4  - FULL (after precession/without precession)
%       Longer simulation in required field (e.g. RMP-field, with ST)
%       Produces:
%           RAW     - file:     In which par and output are stored
%           FULL    - file:     After evaluation (maybe with precession data)
% 5  - POINCARE
%       Plots the poincare map at phi==0 in a figure if GUI is used. Stores
%       position based on phi=0 crossing
%       Produces:
%           RAW     - file:     In which par and output are stored.
% NOTE: THE RAW FILE COULD BE USED TO CONTINUE FROM THIS POINT. DELETE RAW
% FILES WITH INCORRECT PARAMETERS OR EBdyna_go MIGH CRASH / PRODUCE FAULTY
% OUTPUT

%#ok<*NASGU>  % Removes warning in Matlab editor that some values are not used (they are often saved in the output file(s)).

%% Default start parameters 
% (Set first time step to 1, but this could be changed when continueing from a RAW-file
start_time_step=1;
time=0;

%% Load PREC-data if requested
if par.GET_PREC_INFO && isfield(input,'prec')
    prec=input.prec; input=rmfield(input,'prec');
	prec_data=true; % For Poincare plotting in particle categories
else
	prec_data=false; % For Poincare plotting in particle categories
end

%% Add data for specific mode
if par.mode==5
    % Alter output for variable saving of parameters
    output.x=NaN(input.N_job,3,1e4);    output.x_gc=output.x;	output.v=output.x;
    ind_time=zeros(input.N_job,1);   % Vector counting number of indexes already stored
    par.time_scale=NaN(input.N_job,size(output.x,3)); % Time of each save
    
    switch par.poincare_type
        case 'sfl'
            % Pre-allocate and determine size of psi's to be plotted
            psi_large  =NaN(size(output.x,1),1,size(output.x,3));
            theta_large=NaN(size(output.x,1),1,size(output.x,3));
        case 'xy'
            theta_large=NaN(size(output.x,1),1,size(output.x,3));
            output.x_tilde  =NaN(size(output.x,1),1,size(output.x,3));
            output.y_tilde=NaN(size(output.x,1),1,size(output.x,3));
    end
    if par.Poincare_plot
        [ha,~,psi_to_psi_overline]=make_poincare;
    end
    par_dt_or=par.dt;
end

%% Add additional output
%Parameters to determine analytic value of pphi
if par.CALCULATE_PPHI_3D
    output.pphi_an_temp=zeros(input.N_job,par.NB_TIME_STAMPS);
    Delta_pphi_an=zeros(input.N_job,1);
end

% call this script (continue_from_raw) if you have to continue from raw files
% it is quite annoying and useless most of the time so now
% I have cpmmented this part out and pushed it to a separate script

% continue_from_raw;


%% Display input / output
disp('**************************************************************');
disp('PARAMETERS:')
disp(par)
disp('**************************************************************');
disp('INPUT PARTICLES')
disp(input)
disp('**************************************************************');
disp('RAW CONSISTS OF:')
disp(output)
% benchmarking time between two data recordings
disp('**************************************************************');
disp('**************************************************************');

%% SET START PARAMETERS (Temporarily stored parameters etc.)
[~,B]=B_interpolation(x);
Bfield_sq=dot(B,B,2);
x_gc=x+bsxfun(@times,1./(qom*Bfield_sq),cross(v,B,2));



sign_Z_nu=sign(x_gc(:,2)-dim.func_Z_cross(x_gc(:,1))); % Add shift for current equilibrium found empirically! Change when using other equilibrium
sign_vpll_nu=sign(dot(v,B,2));

process_time=tic;
subtime=tic;

%% MAIN SIMULATION LOOP


% Use matprec-file if code has produced a matprec file previouslly, but
% extract_precession failed.
if par.mode==2 && exist([par.SAVENAME_RAW,'prec'],'file')
    disp('prec file already exists in output folder !')
    TEMP_NAME=[par.SAVENAME_STATS(1:end-4),'_matprec.mat'];
    copyfile([par.SAVENAME_RAW,'prec'],TEMP_NAME)
    %load(TEMP_NAME);
    delete([par.SAVENAME_RAW,'prec']);
end
    for time_step=start_time_step:par.NB_TIME_STEPS
        %% Update in advance of next iteration loop
        %Store not ejected variable (to compare which particle have been
        %lost in this loop)
        n_ejected=~ejected;
        
        % Temporarily store position for recording where particles are lost
        x_temp=x;
        v_temp=v;
        
        % Temporarily store signs of current vpll and Z
        sign_vpll_prev=sign_vpll_nu;
        sign_Z_prev=sign_Z_nu;
        		
        %% Time step
        if ~par.CALCULATE_PPHI_3D
            for i=1:par.NR_FUND_IN_LOOP
                [x(n_ejected,:),v(n_ejected,:)]=time_step_integration_GT_eq_struct(x(n_ejected,:),v(n_ejected,:),input.Fc_field(n_ejected,:));
                time=time+par.dt;
            end
        else
            for i=1:par.NR_FUND_IN_LOOP
                v_temp_2=v;
                field_3D=find_3D_Afield(x(n_ejected,:),{'dAphi_dphi','dAR_dphi','dAZ_dphi'});
                
                % The time step
                [x(n_ejected,:),v(n_ejected,:)]=time_step_integration_GT_eq_struct(x(n_ejected,:),v(n_ejected,:),input.Fc_field(n_ejected,:));
                time=time+par.dt;
                % Estimate velocity at n
                v_n=0.5*(v(n_ejected,:)+v_temp_2(n_ejected,:));
                
                % Evaluate pphi
                Delta_pphi_an(n_ejected)=Delta_pphi_an(n_ejected)...
                    +v_n(:,1).*field_3D.dAR_dphi...
                    +v_n(:,2).*field_3D.dAZ_dphi...
                    +v_n(:,3).*field_3D.dAphi_dphi;
            end
        end
        
        %% Determine losses and additional output
        % in case of a Nan in coordinate, consider it ejected
		% ejected=ejected | isnan(x(:,1,:)) ;
		
		X_ind=((x(:,1)-dim.R0)*dim.DX_inv)+dim.mid_Xzero;
        Z_ind=( x(:,2)        *dim.DZ_inv)+dim.mid_Z;
		psi_value=ejected*0;
		psi_value(~ejected)=ba_interp2(maps(1).psi_XZ,Z_ind(~ejected),X_ind(~ejected),'linear'); % find psi-value
		if dim.psi_scale(1)<0  		% psi_Scale increasing to 0
			ejected=ejected | psi_value >= 0;
		else 						% psi_scale decreasing to 0
			ejected=ejected | psi_value <= 0;
		end
        
        if all(ejected)
            warning(['Every particle lost in process: ',num2str(par.PROCESS_NUMBER)])
            break
        end
		
		% debugging
        if length(find(ejected & n_ejected))>0
			disp(['length(find(ejected & n_ejected)) = ' num2str(length(find(ejected & n_ejected)))]);
% 			x_temp(ejected & n_ejected,:)
		end
        % Store particles that have now been ejected (but not previously) in output
        output.x_ej(ejected & n_ejected,:)=x_temp(ejected & n_ejected,:);
        output.v_ej(ejected & n_ejected,:)=x_temp(ejected & n_ejected,:);
        output.time_step_loss(ejected & n_ejected)=time_step;
        output.loss(time_step)=sum(ejected);
        
        %% vpll crossing and midplane crossing
        [~,B]=B_interpolation(x);
        Bfield_sq=dot(B,B,2);
        x_gc=x+bsxfun(@times,1./(qom*Bfield_sq),cross(v,B,2));
        
        % Increase the number of Z-crossings
        sign_Z_nu=sign(x_gc(:,2)-dim.func_Z_cross(x_gc(:,1)));
        
        expr=sign_Z_prev~=sign_Z_nu & ~ejected;
        output.nr_midplane_crossing(expr)=output.nr_midplane_crossing(expr)+1;
        
        % Increase the number of vpll-crossings
        sign_vpll_nu=sign(dot(v,B,2));
        
        expr=sign_vpll_prev~=sign_vpll_nu & ~ejected;
        output.nr_vpll_crossing(expr)=output.nr_vpll_crossing(expr)+1;
        
        %% Store position if POINCARE simulation
        if par.mode==5
            % Criterium for dot: change in tor. angle more than pi/2 (on domain [0,2*pi])
            
            % SWITCH BETWEEN BOTH THESE LINES TO EITHER PLOT WHEN CROSSING
            % PHI=0 OR PLOT EVERY LOOP (OR DECOMMENT THE LOOP OBTAIN)
            ind_part=find(n_ejected & abs(mod(x_temp(:,3),2*pi)-mod(x(:,3),2*pi))>pi/2); % Find indexes of particles crossing phi=0 boundary (note, cannot determine theta for ejected particles)
            %                 ind_part=find(~ejected); if mod(time_step,100)~=0; continue; end
            
            % Skip if no value is saved
            if ~isempty(ind_part)
                ind_time(ind_part)=ind_time(ind_part)+1;                         % Find next positions of output arrays
                if any(ind_time(ind_part)>size(output.x,3))                          % Check if still saveable with current array
                    warning('Simulation stopped since save array exceeded')
                    break
                end
                
                % Repmat all indexes by 3 and coordinates by # that crossed phi=0
                % to find indexes of particles
                ind_part_2D=repmat(ind_part,1,3); ind_coord_2D=repmat(1:3,length(ind_part),1); ind_out_2D=repmat(ind_time(ind_part),1,3);
                
                ind_xv          =sub2ind(size(x)       ,ind_part_2D(:),ind_coord_2D(:));                      %2D indexes to identify values to be saved
                ind_output_xv   =sub2ind(size(output.x),ind_part_2D(:),ind_coord_2D(:),ind_out_2D(:));        %Also take time as third dimension
                
                par.dt=0.5*par.dt;
                [~,v_n,B]=time_step_integration_GT_eq_struct(x,v,input.Fc_field);
                par.dt=par_dt_or;
                
                Bfield=dot(B,B,2);
                x_gc=x+bsxfun(@times,1./(qom*Bfield),cross(v_n,B,2));
                % Store the values (make vector copy)
                %output.x   (ind_output_xv)=x   (ind_xv);
                output.x_gc(ind_output_xv)=x_gc(ind_xv);
                %output.v   (ind_output_xv)=v   (ind_xv);
                
                % Store time
                ind_time_scale=sub2ind(size(par.time_scale),ind_part,ind_time(ind_part)); 	%2D indexes to identify values to be saved
                time=time_step*par.dt*par.NR_FUND_IN_LOOP;
                par.time_scale(ind_time_scale)=time;
                
                if par.Poincare_plot
                    % Delete previous Poincare plot
                    if exist('h','var')
                        delete(h);
                    end
                    switch par.poincare_type
                        case 'sfl'
                            % Indexes 2D
                            X_ind=((output.x_gc(:,1,:)-dim.R0)*dim.DX_inv)+dim.mid_Xzero;
                            Z_ind=( output.x_gc(:,2,:)        *dim.DZ_inv)+dim.mid_Z;
                            expr_nan=isnan(X_ind); expr_nan=expr_nan(:);
                            
                            % Find theta
                            theta_large(~expr_nan) = interpolate_theta_XZ(X_ind(~expr_nan),Z_ind(~expr_nan));
                            % Find psi (index)
                            psi_large(~expr_nan)=ba_interp2(maps(1).psi_XZ,Z_ind(~expr_nan),X_ind(~expr_nan),'linear'); % find psi-value
                            
                            % Fill with new Poincare plot
                            if prec_data
                                %                 h(:,1)=plot(ha,theta(prec.pop.ALL_TRAPPED,:)',psi_to_psi_overline(psi_large(prec.pop.ALL_TRAPPED,:))','+','displayname','trapped');
                                h(:,1)=plot(ha,squeeze(theta_large(~ejected&prec.pop.ALL_PASSING,:))',psi_to_psi_overline(squeeze(psi_large(n_ejected&prec.pop.ALL_PASSING,:)))','k.','displayname','passing');
                            else
                                h=plot(ha,squeeze(theta_large(~ejected,:,:))',squeeze(psi_to_psi_overline(psi_large(~ejected,:,:)))','.');
                            end
                        case 'tor'
                            % Fill with new Poincare plot
                            h=plot(ha,squeeze(output.x_gc(~ejected,1,:))',squeeze(output.x_gc(~ejected,2,:))','k');
                        case 'xy'
                            % Indexes 2D
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
                            expr_3=~any((output.x_tilde.^2+output.y_tilde.^2)>dim.st.r2(301).^2,3);
                            
                            h=plot(ha,squeeze(output.x_tilde(expr_3,:,:))',squeeze(output.y_tilde(expr_3,:,:))','k');
                    end
                    drawnow
                end
            end
            % look at time evolution also for debugging and consistency
            if (mod(time_step,par.TIME_STAMP_PRECISION)==0)
                time_stamp=ceil((time_step)/par.TIME_STAMP_PRECISION);
                output.x(:,:,time_stamp)=x;
                output.v(:,:,time_stamp)=v;
            end

            %% Store x and v
        elseif (mod(time_step,par.TIME_STAMP_PRECISION)==0)
            time_stamp=ceil((time_step)/par.TIME_STAMP_PRECISION);
            output.x(:,:,time_stamp)=x;
            output.v(:,:,time_stamp)=v;
            
            % Recorrect time-global for precision
            time=par.time_scale(time_stamp)+par.dt;
            
            if par.CALCULATE_PPHI_3D
                % Store analytic value
                output.pphi_an_temp(:,time_stamp)=Delta_pphi_an;
                Delta_pphi_an(:)=0;
            end
            disp(['storing x and v values at time stamp ' num2str(time_stamp) ' and time step ' num2str(time_step) '...' ])
        end
        
        %% SAVE DATA FILE intermediately
        if par.SAVE_DATA_FILE && mod(time_step,par.RECORD_PRECISION)==0
            disp('------------------------------------------------')
            toc(subtime)
            if par.mode~=5
                disp(strcat('time_step # ',num2str(time_step),' ; time_stamp # ',num2str(time_stamp),' ; time = ',num2str(par.time_scale(time_stamp))));
            else
                disp(['max phi crossings  ',num2str(max(ind_time)),' ; avg phi crossings  ',num2str(mean(ind_time(~ejected)))]);
            end
            disp(['Now at: ',num2str(100*time_step/par.NB_TIME_STEPS,'%2.0f'),' %'])
            disp(strcat('total # ejected particles =  ',num2str(sum(ejected))));
            save(par.SAVENAME_RAW,'-v7.3','input','output','par')
            
            subtime=tic;
        end
    end
    if exist('process_time_raw','var')
        process_time=process_time_raw;
    else
        process_time=toc(process_time);
    end
    if par.CALCULATE_PPHI_3D
        output.pphi_an_temp=par.dt*input.Z*cumsum(output.pphi_an_temp,2);
    end
    disp('**************************************************************');
    disp(['Total execution time of process ',num2str(par.PROCESS_NUMBER),' is: ',num2str(process_time)]);


%% EVALUATE / SAVE DATA
switch par.mode
    case 1 %TEST
        if par.SAVE_DATA_FILE
            save(par.SAVENAME_RAW,'input','output','par','ejected','process_time','-v7.3') % Save in case evaluation fails
        end
            output=evaluate_output(input,output,ejected);
        if par.SAVE_DATA_FILE
            save(par.SAVENAME,'input','output','par','ejected','process_time','-v7.3')
        end
    case 2 % PREC
        if all(ejected)
            warning(['All particles ejected in precession simulation of process ',num2str(par.PROCESS_NUMBER)]);
            if par.SAVE_DATA_FILE
                save([par.SAVENAME_RAW,'prec'],'input','output','par','ejected','process_time','-v7.3')
            end
            exitcode=2;
            return
        end
        try
            % Get pop data
            output=evaluate_output(input,output,ejected);                    % Get more info, e.g. v_parallel
            prec = extract_precession_information_struct(input,output,ejected);    % Identify particles and orbit/precession frequency
%             output=remove_fields(output,[],{'x','v','vpll','nr_vpll_crossing','nr_midplane_crossing'});  % Remove any fields apart from...
            % the final gc position of a precession simulation is important
            % information and is kept as the new input position for the
            % next simulation
            input.x_ini=x;
            input.v_ini=v;
            input.x=x_gc;   % squeeze(output.x_gc(:,:,end));
            input.v=v;		%squeeze(output.v(:,:,end));
            input.vpll=squeeze(output.vpll(:,end));
                
            % a large amount of data is there but we keep the bare minimum
			% pphi_kin can be used to check time step precision
            output=remove_fields(output,[],{'x_gc','pphi_kin','vpll','time_step_loss','x_ej','loss'});              % Remove any fields apart from...
            
            % Store the data of precession simulation in prec
            prec.par=remove_fields(par,[],{'dt','NB_TIME_STAMPS'}); % Store NB_TIME_STAMPS and dt
            prec.par.end_time=par.time_scale(end);                  % Store the last tme
            prec.ejected=ejected;                                   % Store those ejected in 2D
 

			% resizing the file so that we have reasonable compromise between
			% information about what happened and size of data
			[par,output 	]=reduce_time_stamps(par,output); 
 
            % Save in case of only one file
            if par.SAVE_DATA_FILE && par.NB_PROCESS==1
                save(par.SAVENAME_STATS,'-v7.3','input','output','par','ejected','process_time','prec')
            end
            % Delete a previous made .matprec-file if any
            if exist([par.SAVENAME_RAW,'prec'],'file')
                delete([par.SAVENAME_RAW,'prec'])
            end
			% split prec files
			if par.SAVE_DATA_FILE && par.NB_PROCESS>1
			    save(par.SAVENAME_RAW,'-v7.3','input','output','par','ejected','process_time','prec')
            end
			
        catch err
            if ~exist([par.SAVENAME_RAW,'prec'],'file')
                warning('SAVING RAW DATA FROM PRECESSIONS IN MATPREC-FILE')
                if par.SAVE_DATA_FILE
                    save([par.SAVENAME_RAW,'prec'],'input','output','par','ejected','process_time','-v7.3')
                end
            else
                warning('RAW PRECESSION DATA FAILED TO BE ANALYZED')
            end
            rethrow(err)
        end
        if all(ejected)
            warning(['All particles ejected in precession simulation of process ',num2str(par.PROCESS_NUMBER)]);
            if par.SAVE_DATA_FILE
                save([par.SAVENAME_RAW,'prec'],'input','output','par','ejected','process_time','-v7.3')
            end
            exitcode=2;
            return
        end
        
        % Rerun for full simulation in mode 3
        disp('Precession information made')
        clearvars -EXCEPT par
        %exitcode=G_eq(num2str(par.PROCESS_NUMBER),num2str(par.NB_PROCESS));
        %disp('Continueing to FULL-simulation (mode 3)
        exitcode=2;
        return
    case {3,6} %FULL (after/without precession)
        % Save a RAW file before evaluation
        if par.SAVE_DATA_FILE
            if exist('prec','var')
                save(par.SAVENAME_RAW,'input','output','par','ejected','process_time','prec','-v7.3')
            else
                save(par.SAVENAME_RAW,'input','output','par','ejected','process_time','-v7.3')
            end
        end
        
        % Evaluate the raw-data
        output_raw=output;
        output=evaluate_output(input,output,ejected); 
%         output=remove_fields(output,[],{'ST_interaction','nr_vpll_crossing','nr_midplane_crossing','time_step_loss','loss','x_ej','v_ej','pphi_kin','mm','Delta_pphi','pphi_an','dpphi_dt','dmm_dt','nr_profile','wb','wd'}); % Remove everything but these fields.
        output=remove_fields(output,[],{'x','v','x_gc','ST_interaction','time_step_loss','loss','x_ej','v_ej','pphi_kin','mm','Delta_pphi','pphi_an','dpphi_dt','dmm_dt','nr_profile','wb','wd'}); % Remove everything but these fields.
        
        % Reduce the number of time stamps
		[~ 	,output_raw	]=reduce_time_stamps(par,output_raw);
        [par,output 	]=reduce_time_stamps(par,output);  %#ok<ASGLU>
        
        % Save the data definitely
        if par.SAVE_DATA_FILE
            if exist('prec','var')
                save(par.SAVENAME,'input','output','par','ejected','process_time','prec','-v7.3')
            else
                save(par.SAVENAME,'input','output','par','ejected','process_time','-v7.3')
            end
        end
        
        output=output_raw;
        if par.SAVE_DATA_FILE
            if exist('prec','var')
                save(par.SAVENAME_RAW,'input','output','par','ejected','process_time','prec','-v7.3')
            else
                save(par.SAVENAME_RAW,'input','output','par','ejected','process_time','-v7.3')
            end
        end
    case 5 %% POINCARE
        if par.SAVE_DATA_FILE
            save(par.SAVENAME_RAW,'input','output','par','ejected','process_time','-v7.3')
        end
        % Remove not used / obsolete fields in output array
        output.x_gc(:,:,max(ind_time)+1:end)=[];
        % keep trajectories for debugging and consistency
        %output.x(:,:,max(ind_time)+1:end)=[];
        %output.v(:,:,max(ind_time)+1:end)=[];
        par.time_scale(:,max(ind_time)+1:end)=[];
        
        % Indexes 2D
        X_ind=squeeze((output.x_gc(:,1,:)-dim.R0)*dim.DX_inv)+dim.mid_Xzero;
        Z_ind=squeeze( output.x_gc(:,2,:)        *dim.DZ_inv)+dim.mid_Z;
        
        expr_nan=isnan(X_ind);
        %epxr_nan=expr_nan(:);
        
        % Store psi and theta in output-struct
        % Psi-surfaces
        output.psi=NaN(size(output.x_gc,1),size(output.x_gc,3));
        output.theta=output.psi;
        
        % Find theta
        theta(~expr_nan) = interpolate_theta_XZ(X_ind(~expr_nan),Z_ind(~expr_nan));
        
        % Find psi (index)
        psi=ba_interp2(maps(1).psi_XZ,Z_ind(~expr_nan),X_ind(~expr_nan),'cubic');
        if par.SAVE_DATA_FILE
            save(par.SAVENAME,'input','output','par','ejected','process_time','-v7.3')
        end
    otherwise
        if par.SAVE_DATA_FILE
            save(par.SAVENAME_RAW,'input','output','par','ejected','process_time','-v7.3')
        end
        error('MODE NOT KNOWN OR OBSOLETE / SAVING RAW DATA-FILE')
end

%% EXIT
disp('**************************************************************');
disp('done');

exitcode=0; % Notation: exit prevents exitcode to work
disp('exitcode = ')
disp(exitcode)
end

%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAKE A FIGURE FOR POINCARE-PLOT
function [ha,hq,psi_to_psi_overline]=make_poincare
global par dim maps
q_to_psi = @(q) interp1(1:length(dim.psi_scale),dim.psi_scale,interp1(dim.q_initial_profile,1:length(dim.q_initial_profile),q),'cubic'); % psi from q
psi_to_q = @(psi) interp1(dim.psi_scale,dim.q_initial_profile,psi,'cubic');
psi_to_psi_overline=[];
% Plot lines of rational values with RMP_mode
RMP_mode=2;
rational_q_surfaces=4:(-1/RMP_mode):1; % For q=1 to q=4

%% Make a figure
try set(groot,'defaulttextinterpreter','latex');    catch err; end;
delete(findall(0,'type','figure','tag','Poincare_map'))
hf=figure('Name','Poincare map','tag','Poincare_map');

open_fig_ST_xy='./2016-12-13_ST_xy_frame_4.fig';
open_fig_ST_RZ='./2016-12-13_ST_RZ_frame_4.fig';

switch par.poincare_type
    case 'sfl'
        %% psi-theta=space
        psi_to_psi_overline = @(psi) 1-psi/maps(1).psi_global;   % normalized psi from psi
        q_to_psi_overline = @(q) psi_to_psi_overline(q_to_psi(q));                  % normalized psi from q
        
        [AX,h1,h2]=plotyy([0 1],[0 1],[0 1],[0 1],'parent',hf);                     % Make an double y-axis figure
        delete(h1); delete(h2); linkaxes(AX,'xy'); clear h1 h2
        ha=AX(1);   hold(ha,'on'); set(ha,'FontSize',20);                           % Have an psi-theta plot figure
        hq=AX(2);   hold(hq,'on'); set(hq,'FontSize',20);                           % Have an q-theta plot figure
        set(ha,'YColor','k'); set(hq,'YColor','r');
        xlabel(ha,'$\theta$ [rad]','FontSize',24,'interpreter','latex');
        ylabel(ha,'$\overline\psi$','FontSize',24,'interpreter','latex');
        ylabel(hq,'$q$','FontSize',24,'interpreter','latex');
        set(ha,'XLim',[0 2*pi],'Ylim',[0 psi_to_psi_overline(0)]);
        
        % Determine psi on these surfaces
        psi_overline_q_surfaces=q_to_psi_overline(rational_q_surfaces);
        
        % Make lines on q-values and label them
        q_surf_name{length(rational_q_surfaces)}=[];
        for i=1:length(rational_q_surfaces)
            if mod(rational_q_surfaces(i),1)==0
                line([0 2*pi],[psi_overline_q_surfaces(i) psi_overline_q_surfaces(i)],'color','r','linestyle','-','LineWidth',1.5,'parent',hq);
                q_surf_name{i}=[num2str(rational_q_surfaces(i)*RMP_mode),'/',num2str(RMP_mode)];
            else
                line([0 2*pi],[psi_overline_q_surfaces(i) psi_overline_q_surfaces(i)],'color','r','linestyle',':','LineWidth',1.5,'parent',hq);
            end
        end
        set(hq,'YTick',fliplr(psi_overline_q_surfaces),'YtickLabel',fliplr(q_surf_name))
        set(ha,'YTick',fliplr(psi_overline_q_surfaces));
        
    case 'tor'
        hq=[];
        q_XZ=psi_to_q(maps(1).psi_XZ);
        
        ha=axes('parent',hf);
        hold(ha,'on'); set(ha,'FontSize',20);
        axis(ha,'equal','xy')
        xlabel(ha,'$R$','FontSize',24,'interpreter','latex');
        ylabel(ha,'$Z$','FontSize',24,'interpreter','latex');
        if par.APPLY_SAWTOOTH && exist(open_fig_ST_RZ,'file')
            hf_temp=openfig(open_fig_ST_RZ);
            ha_temp=findall(hf_temp,'type','axes');
            ha_temp=ha_temp(2);
            h_temp=findall(ha_temp,'type','contour');
            set(flipud(h_temp),'parent',ha);
            close(hf_temp)
            colormap(ha,'jet')
            caxis(ha,[0 max(h_temp.ZData(:))]);
        end
        % Plot LCFS
        c=contourc(dim.R0+dim.scale_X,dim.scale_Z,maps(1).psi_norm_XZ',[513 513]);
        plot(ha,c(1,2:end),(c(2,2:end)),'r','displayname','LCFS','LineWidth',2)
        contour(dim.R0+dim.scale_X,dim.scale_Z,q_XZ',rational_q_surfaces,'ShowText','on','parent',ha,'color','k','linewidth',1);
        %         for i=1:length(rational_q_surfaces)
        %             c=contourc(dim.R0+dim.scale_X,dim.scale_Z,maps(1).psi_XZ',[psi_q_surfaces(i) psi_q_surfaces(i)]);
        %             if mod(rational_q_surfaces(i),1)==0
        %                 arg_line_color_style='r';
        %             else
        %                 arg_line_color_style='r:';
        %             end
        %             plot(ha,c(1,2:end),(c(2,2:end)),arg_line_color_style,'displayname',['resonant surface, q= ',[num2str(rational_q_surfaces(i)*RMP_mode),'/',num2str(RMP_mode)]],'LineWidth',1)
        %         end
    case 'xy'
        hq=[];
        q_XZ=psi_to_q(maps(1).psi_XZ);
        
        ha=axes('parent',hf);
        hold(ha,'on'); set(ha,'FontSize',20);
        axis(ha,'equal','xy')
        xlabel(ha,'$x$','FontSize',24,'interpreter','latex');
        ylabel(ha,'$y$','FontSize',24,'interpreter','latex');
        if par.APPLY_SAWTOOTH && exist(open_fig_ST_RZ,'file')
            hf_temp=openfig(open_fig_ST_xy);
            ha_temp=findall(hf_temp,'type','axes');
            ha_temp=ha_temp(2);
            h_temp=findall(ha_temp,'type','contour');
            %             set(flipud(h_temp),'parent',ha);
            [~,h_new]=contourf(h_temp.XData,h_temp.YData,h_temp.ZData,200,'linestyle',':','linecolor','r','linewidth',1,'parent',ha);
            close(hf_temp)
            colormap(ha,'summer')
            caxis(ha,[0 max(h_new.ZData(:))]);
            xlim([-1 1]*max(abs(h_new.XData(:))));
            ylim([-1 1]*max(abs(h_new.YData(:))));
        end
        % Plot q-values
        r=sqrt(2*maps(1).chi_XZ);
        theta= maps(1).theta_normal_XZ;
        x=r.*cos(theta);
        y=r.*sin(theta);
        contour(x,y,q_XZ,rational_q_surfaces,'ShowText','on','parent',ha,'color','b','linewidth',1);
    otherwise
        error('plotting of poincare not specified')
end

end
