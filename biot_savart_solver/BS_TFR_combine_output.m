function [  ] = BS_TFR_combine_output(  )
%BS_TFR_combine_output Combines files ending with p1 to a total file
%   After parallel processing of BS3D we would like to make maps of the
%   total fields. Therefore we combine the files

%% Get into the output folder
save_folder='./output/';


files=dir(save_folder);
fnames=[];
for i=3:length(files)
    if ~isempty(strfind(files(i).name,'individual')) && ~isempty(strfind(files(i).name,'TF'))
        fnames{end+1}=strcat(save_folder,files(i).name); %#ok<AGROW>
    end
end

%% Load first file
base=load(fnames{1}); % Skip 1 and 2 (. and ..)
par=base.par;
% par.NB_TF_COILS=14

if length(fnames)~=par.N_PROCESS
     error('Failed to match number of processes with files')
end

% if ~strcmp(par.coord_sys,'flux')
% 	warning('Please redo the calculation with a flux-option, with 2 radial points, to get the B field at the LCFS for scaling and ripple calculation')
% end

% Save name final file
save_name=strcat('./output_TF/TF_',par.coord_sys,'_',datestr(now,'yyyy-mm-dd'),'.mat')

switch par.nr_coils
    case 1
        disp('Only 1 coil calculated. Will by periodically extended by par.NB_TF_COILS')
        coil_nr='single_coil';
        
        if mod(par.size_total(3)-1,par.NB_TF_COILS)~=0;    error('Positions TF coils / phi-positions don''t match');   end;
    case 16
        disp('Full 16 coils calculated')
        coil_nr='16_coils';
    otherwise
        error('Number of coils not recognized')
end

%% Loop over all the rest and add them in the same field
BA_fields=fieldnames(base.field)';
sca_fac=5.0/3.8

for i=1:length(fnames)
    new=load(fnames{i}); % Load a coil file
        
    for BA=BA_fields
        if any(strcmp(BA{1},{'BR','BZ','Bphi','AR','AZ','Aphi','dAphi_dphi','dAR_dphi','dAZ_dphi'}))
            if i==1
                comb_field.(BA{1})=zeros(base.par.size_total); % Pre-allocate
            end
            comb_field.(BA{1})(new.par.indexes)=new.field.(BA{1})(:);
        end
    end
    disp(['Succesfully added ',num2str(i),' of ',num2str(base.par.N_PROCESS),' at: ',num2str(i*100/length(fnames),'%2.0f'),'%'])
end
par.paths=initialize_folder_names_struct;
filename=[par.paths.DATA_FOLDER,'motions_map_dimensions.mat'];
dim_map_RZ=load(filename,'scale_X','scale_Z','Z_axis','R0','Raxis','a');

BA_fields=fieldnames(comb_field)';
%% Determine complete field (of all 16 coils)
switch coil_nr
    case '16_coils'
        % Reduce size to 1/par.th since periodic
        if mod(size(comb_field,3),16)==0
            ind_end=(size(comb_field,3)-1)/16+1;
            comb_field=comb_field(:,:,ind_end);
        end
        
        % Single save the field
        save(save_name,'-v7.3','-struct','comb_field');
        total_field=comb_field; clear comb_field;
        total_field.symm=16;
    case 'single_coil'
        % Ooh ppff..  lets try to periodically combine them, but this takes
        % some RAM
        try
            % Okay, now we need to translate the field over some phi-angle to
            % get the BA-fields for the a coil
                   
            total_phi_elements=par.size_total(3);                  % # total elements                      [9]
            translate_phi=(total_phi_elements-1)/par.NB_TF_COILS               % # phi-positions between TF coils     [8/4=2]
            
            % Let's make even more use of symmetry, by storing only the
            % mode
            symm=1;
            for mode_trial=[par.NB_TF_COILS 8 4 2]
                if mod(total_phi_elements-1,mode_trial)==0
                    symm=mode_trial;
                    break
                end
            end
			size_total=par.size_total;
			size_total(3)=(total_phi_elements-1)/symm+1;     % Cut the total combined field to this value [5 with n=2] Other half is symmetric
            
            total_field.symm=symm
            % Pre-allocate the total field
            for BA=BA_fields
                total_field.(BA{1})=zeros(size_total);
            end
            
            for c=1:par.NB_TF_COILS
                phi_translation= (c-1) * translate_phi;        % Phi-positions that are translated in respect to the first
                indexes_translated = mod( (1:size_total(3)) - phi_translation , total_phi_elements-1);  % Start with [1:5], reverse by 2, with correct indexes
                indexes_translated(indexes_translated==0)=total_phi_elements-1;
                % Superposition of each field
                for j=1:length(BA_fields)	
                    total_field.(BA_fields{j})=total_field.(BA_fields{j})+comb_field.(BA_fields{j})(:,:,indexes_translated);
                end
            end

			%'BR','BZ','Bphi','AR','AZ','Aphi','dAphi_dphi','dAR_dphi','dAZ_dphi'

			disp('Original output Bphi0 value : ')
			mid_pos=round(0.5*size(total_field.Bphi()));
			[ ~ , Rpos]=min(abs(dim_map_RZ.scale_X+dim_map_RZ.R0-0.894));
			Bphi_mid=total_field.Bphi(Rpos,mid_pos(2),1)
			% if rescaling needed
			if sca_fac~=1
				warning(['Rescaling the calculatied TF fields by a factor ' num2str(sca_fac)])
				for j=1:length(BA_fields)
					total_field.(BA_fields{j})=sca_fac*(total_field.(BA_fields{j}));           
				end
			else
				disp('using simulation TF field amplitude : check the input TF coils currents to be sure!')	
			end

			disp('adjusting final value of field to have Bphi0 at about : ')
			Bphi_mid=total_field.Bphi(Rpos,mid_pos(2),1)
				
            save(save_name,'-v7.3','-struct','total_field');
            disp('Great succes!')
        catch err
            % save('./output_TF/TF_single_coil_only.mat','-v7.3','-struct','comb_field');
            rethrow(err)
        end
    otherwise
		warning('Something went wrong with the definition of the number of coils')
        save(save_name,'-v7.3','-struct','comb_field');
		save(savename,'-append','symm');
end
disp('Combined output file created')

clearvars -except total_field par
TF_to_TFR(total_field,par.coord_sys);
end